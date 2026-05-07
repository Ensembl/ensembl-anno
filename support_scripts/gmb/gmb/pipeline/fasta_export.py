#!/usr/bin/env python3
"""FASTA export for Gene Model Builder.

Python replacement for gtf_to_seq.pl.  Produces:
  - cdna.fa   : spliced cDNA for all selected transcripts
  - prot.fa   : protein for transcripts with CDS
  - cds.fa    : CDS nucleotide sequence (optional)

All sequences use stable IDs matching the consensus GFF3.
"""

from __future__ import annotations

import os
from typing import IO, TYPE_CHECKING, Any

if TYPE_CHECKING:
    from gmb.pipeline.config import PipelineConfig

import pandas as pd

from gmb.pipeline.annotate_cds_utrs import (
    build_spliced_seq,
    find_best_orf,
    map_cds_to_genomic,
    reverse_complement,
    translate,
)


def _build_cds_seq_from_intervals(
    cds_intervals: list[tuple[int, int]],
    strand: str,
    chrom_seq: str,
) -> str:
    """Build CDS nucleotide sequence from genomic CDS intervals.

    Exons must be concatenated in ascending genomic order before
    reverse-complementing.  RC(X+Y) == RC(Y)+RC(X), so concatenating in
    ascending order and then RC'ing the whole string naturally places the
    highest-coordinate exon (5' for minus-strand) first in the result.
    Concatenating in descending order before RC'ing inverts the exon order.
    """
    parts = [chrom_seq[cs:ce] for cs, ce in sorted(cds_intervals)]
    cds_nuc = "".join(parts)
    if strand == "-":
        cds_nuc = reverse_complement(cds_nuc)
    return cds_nuc


def export_fasta(
    loci: list[tuple[Any, list[dict]]],
    genome: dict[str, str],
    helixer_cds: pd.DataFrame | None,
    config: PipelineConfig,
    output_dir: str,
) -> dict[str, int]:
    """Export cDNA, protein, and CDS FASTA files.

    Parameters
    ----------
    loci : list of (cluster_id, list of model dicts)
        Selected isoforms from the scoring stage.
        Each model dict has keys: id, chrom, strand, df (exon DataFrame),
        combined_evidence, and a stable output ID assigned during GFF3 writing.
    genome : dict
        ``{chrom: sequence}`` from ``load_genome``.
    helixer_cds : pd.DataFrame or None
        CDS rows for Helixer models.
    config : PipelineConfig
    output_dir : str
        Output directory path.

    Returns
    -------
    dict
        Export statistics with keys ``total_transcripts``, ``with_cds``,
        ``with_protein``, ``partial_5``, ``partial_3``.
    """
    ecfg = config.export
    ocfg = config.orf

    # Build Helixer CDS lookup
    cds_by_tid = {}
    if helixer_cds is not None and not helixer_cds.empty:
        for tid, grp in helixer_cds.groupby("transcript_id"):
            cds_by_tid[tid] = sorted(zip(grp["Start"].values, grp["End"].values))

    os.makedirs(output_dir, exist_ok=True)

    cdna_path = os.path.join(output_dir, "cdna.fa")
    prot_path = os.path.join(output_dir, "prot.fa")
    cds_path = os.path.join(output_dir, "cds.fa")

    stats = {
        "total_transcripts": 0,
        "with_cds": 0,
        "with_protein": 0,
        "partial_5": 0,
        "partial_3": 0,
    }

    cdna_fh = open(cdna_path, "w") if ecfg.write_cdna else None
    prot_fh = open(prot_path, "w") if ecfg.write_protein else None
    cds_fh = open(cds_path, "w") if ecfg.write_cds else None

    try:
        for cluster_id, models in loci:
            gene_id = models[0].get("gene_id", f"GENE_{cluster_id}")

            for i, m in enumerate(models):
                tid = m.get("output_tid", f"{gene_id}.{i + 1}")
                chrom = m["chrom"]
                strand = m["strand"]
                evidence = m.get("combined_evidence", m.get("source", ""))

                if chrom not in genome:
                    continue

                chrom_seq = genome[chrom]
                exons = sorted(zip(m["df"]["Start"].values, m["df"]["End"].values))

                stats["total_transcripts"] += 1

                # --- cDNA ---
                cdna = build_spliced_seq(exons, strand, chrom_seq)
                if cdna_fh:
                    header = (
                        f">{tid} gene={gene_id} "
                        f"evidence={evidence} "
                        f"loc={chrom}:{exons[0][0]+1}-{exons[-1][1]}"
                        f"({strand})"
                    )
                    cdna_fh.write(f"{header}\n")
                    _write_seq(cdna_fh, cdna)

                # --- CDS + protein ---
                cds_intervals = None
                is_partial_5 = False
                is_partial_3 = False

                # Check for existing CDS (Helixer)
                existing_cds = cds_by_tid.get(m["id"])
                if existing_cds:
                    cds_intervals = existing_cds
                else:
                    # Predict ORF
                    orf = find_best_orf(cdna, min_codons=ocfg.min_codons)
                    if orf is not None:
                        cds_start, cds_end, is_partial_5, is_partial_3 = orf
                        cds_intervals = map_cds_to_genomic(cds_start, cds_end, exons, strand)

                if cds_intervals:
                    stats["with_cds"] += 1
                    if is_partial_5:
                        stats["partial_5"] += 1
                    if is_partial_3:
                        stats["partial_3"] += 1

                    cds_nuc = _build_cds_seq_from_intervals(cds_intervals, strand, chrom_seq)

                    protein = translate(cds_nuc)
                    if protein.endswith("*"):
                        protein = protein[:-1]

                    if protein:
                        stats["with_protein"] += 1

                    # Partial markers
                    orf_status = ""
                    if is_partial_5:
                        orf_status += " partial_5prime"
                    if is_partial_3:
                        orf_status += " partial_3prime"
                    if not is_partial_5 and not is_partial_3:
                        orf_status = " complete"

                    # Write CDS FASTA
                    if cds_fh:
                        cds_header = (
                            f">{tid} gene={gene_id} " f"cds_len={len(cds_nuc)}" f"{orf_status}"
                        )
                        cds_fh.write(f"{cds_header}\n")
                        _write_seq(cds_fh, cds_nuc)

                    # Write protein FASTA
                    if prot_fh and protein:
                        if not ecfg.include_partial and (is_partial_5 or is_partial_3):
                            continue
                        prot_header = (
                            f">{tid} gene={gene_id} " f"aa_len={len(protein)}" f"{orf_status}"
                        )
                        prot_fh.write(f"{prot_header}\n")
                        _write_seq(prot_fh, protein)

    finally:
        if cdna_fh:
            cdna_fh.close()
        if prot_fh:
            prot_fh.close()
        if cds_fh:
            cds_fh.close()

    # Remove empty optional files
    if ecfg.write_cds and stats["with_cds"] == 0 and os.path.exists(cds_path):
        os.remove(cds_path)

    print(
        f"  FASTA export: {stats['total_transcripts']} transcripts, "
        f"{stats['with_cds']} with CDS, {stats['with_protein']} with protein"
    )

    return stats


def _write_seq(fh: IO[str], seq: str, line_width: int = 80) -> None:
    """Write a sequence in wrapped FASTA format.

    Delegates to :func:`gmb.utils.fasta.write_seq`.
    """
    from gmb.utils.fasta import write_seq
    write_seq(fh, seq, line_width)
