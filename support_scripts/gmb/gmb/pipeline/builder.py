#!/usr/bin/env python3
"""Gene Model Builder orchestrator.

Loads evidence, applies filters, performs optional protein validation,
and generates consensus models.

DataFrame Schema Contract
-------------------------
Evidence DataFrames passed through the pipeline share these columns:

* ``Chromosome`` : str -- sequence/contig name
* ``Start`` : int -- 0-based start (half-open, pyranges convention)
* ``End`` : int -- 1-based end (half-open, pyranges convention)
* ``Strand`` : str -- ``"+"`` or ``"-"``
* ``transcript_id`` : str -- unique transcript identifier (source-prefixed)
* ``gene_id`` : str -- gene-level grouping identifier (source-prefixed)
* ``Source`` : str -- evidence origin, e.g. ``"Scallop"``, ``"Helixer"``
* ``Feature`` : str -- GFF3/GTF feature type (``"exon"``, ``"CDS"``, etc.)

Additional columns may be present depending on the source format
(``Score``, ``Coverage``, ``Identity``, ``combined_evidence``, etc.).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict
from typing import TYPE_CHECKING

import yaml

if TYPE_CHECKING:
    from gmb.pipeline.config import PipelineConfig

import numpy as np
import pandas as pd
import pyranges as pr

from gmb.pipeline.annotate_cds_utrs import (
    annotate_all_transcripts,
    build_spliced_seq,
    load_genome,
    reverse_complement,
    translate,
)
from gmb.pipeline.config import dump_config, list_build_presets, load_config
from gmb.pipeline.dedup_genes import dedup_genes
from gmb.pipeline.duplicate_transcript_collapse import collapse_exact_duplicate_transcripts
from gmb.pipeline.evidence_filter import (
    filter_backbone_models,
    filter_chimeras,
    filter_protein_evidence,
    split_mega_transcripts,
)
from gmb.pipeline.gff3_validate import validate_and_fix_gff3
from gmb.pipeline.protein_validation import batch_score_proteins, check_dependencies
from gmb.pipeline.scoring import select_isoforms
from gmb.pipeline.subset_utils import (
    add_subset_args,
    build_mapping,
    remap_df_seqnames,
    remap_genome_seqnames,
    resolve_subset_regions,
    subset_df_by_regions,
    write_subset_manifest,
)
from gmb.utils.intervals import same_strand_overlap_ids
from gmb.utils.logging import resolve_log_file, setup_logging


def regenerate_final_fasta(
    gff_rows: list[dict],
    genome_dict: dict[str, str],
    output_dir: str,
) -> dict[str, int]:
    """Derive cdna.fa, cds.fa, and prot.fa from the final GFF3 rows.

    Must be called AFTER all post-processing (validation, dedup, collapse)
    so that the FASTA sequences match the final GFF3 exon/CDS coordinates.
    Coordinates in *gff_rows* must be 0-based half-open (internal convention).
    Terminal stop codons are excluded from prot.fa (standard Ensembl convention).
    """
    from collections import defaultdict as _dd

    from gmb.utils.fasta import write_seq

    children: dict[str, list[dict]] = _dd(list)
    mrnas: dict[str, dict] = {}
    for r in gff_rows:
        feat = r.get("Feature", "")
        if feat == "mRNA":
            mrnas[r["ID"]] = r
        elif r.get("Parent"):
            children[r["Parent"]].append(r)

    stats = {"cdna": 0, "cds": 0, "prot": 0}
    cdna_path = os.path.join(output_dir, "cdna.fa")
    cds_path = os.path.join(output_dir, "cds.fa")
    prot_path = os.path.join(output_dir, "prot.fa")

    with (
        open(cdna_path, "w") as cdna_fh,
        open(cds_path, "w") as cds_fh,
        open(prot_path, "w") as prot_fh,
    ):
        for tid in sorted(mrnas):
            m = mrnas[tid]
            chrom = m["Chromosome"]
            strand = m["Strand"]
            if chrom not in genome_dict:
                continue
            chrom_seq = genome_dict[chrom]
            kids = children.get(tid, [])

            exon_ivs = sorted([(r["Start"], r["End"]) for r in kids if r["Feature"] == "exon"])
            cds_ivs = sorted([(r["Start"], r["End"]) for r in kids if r["Feature"] == "CDS"])

            if exon_ivs:
                cdna = build_spliced_seq(exon_ivs, strand, chrom_seq)
                cdna_fh.write(f">{tid}\n")
                write_seq(cdna_fh, cdna)
                stats["cdna"] += 1

            if cds_ivs:
                cds_parts = [chrom_seq[s:e] for s, e in cds_ivs]
                cds_nuc = "".join(cds_parts)
                if strand == "-":
                    cds_nuc = reverse_complement(cds_nuc)

                cds_fh.write(f">{tid}\n")
                write_seq(cds_fh, cds_nuc)
                stats["cds"] += 1

                protein = translate(cds_nuc)
                if protein.endswith("*"):
                    protein = protein[:-1]

                prot_fh.write(f">{tid}\n")
                write_seq(prot_fh, protein)
                stats["prot"] += 1

    print(
        f"  Final FASTA: {stats['cdna']} cDNA, {stats['cds']} CDS, "
        f"{stats['prot']} protein records"
    )
    return stats


def compute_cds_phases(cds_intervals: list, strand: str) -> list:
    """Compute GFF3 CDS phase values for a list of CDS intervals.

    GFF3 phase is the number of bases at the biological 5' start of each CDS
    segment that must be skipped to reach the first base of the next complete
    codon.  Phase 0 at the transcript's biological 5' CDS start.  Values 0,
    1, or 2 only; never ".".

    Parameters
    ----------
    cds_intervals : list of (start, end)
        CDS intervals as 0-based half-open genomic coordinates, in any order.
    strand : str
        "+" or "-".  Any other value returns an empty list (caller must guard
        against unresolved strand before calling).

    Returns
    -------
    list of int
        Phase per CDS interval in ascending genomic coordinate order
        (as required for GFF3 output).  Same length as cds_intervals.
    """
    if not cds_intervals or strand not in ("+", "-"):
        return []

    sorted_asc = sorted(cds_intervals)

    # Biological translation order: ascending for "+", descending for "-"
    bio_order = sorted_asc if strand == "+" else list(reversed(sorted_asc))

    cumulative_bp = 0
    bio_phases: list[int] = []
    for s, e in bio_order:
        bio_phases.append((3 - (cumulative_bp % 3)) % 3)
        cumulative_bp += e - s

    # Map back to ascending genomic (GFF3) order
    if strand == "+":
        return bio_phases
    else:
        return list(reversed(bio_phases))


def compute_percentile_guardrails(
    locus_df: pd.DataFrame,
    config: PipelineConfig,
    protein_supported_tids: set[str],
) -> dict[str, int]:
    """Compute effective guardrails based on high-confidence candidates.

    Parameters
    ----------
    locus_df : pd.DataFrame
        Exon-level DataFrame for clustered loci, containing at least
        ``transcript_id``, ``Source``, ``Start``, ``End``, and optionally
        ``combined_evidence``.
    config : PipelineConfig
        Pipeline configuration with ``validation`` section.
    protein_supported_tids : set of str
        Transcript IDs that overlap protein evidence.

    Returns
    -------
    dict
        Keys ``"effective_max_exon_len_bp"`` and
        ``"effective_max_transcript_span_bp"`` with integer limits.
    """
    val_cfg = config.validation

    # Start with configured base limits
    runtime_params = {
        "effective_max_exon_len_bp": val_cfg.max_exon_len_bp,
        "effective_max_transcript_span_bp": val_cfg.max_transcript_span_bp,
    }

    if (
        val_cfg.max_exon_len_mode != "percentile"
        and val_cfg.max_transcript_span_mode != "percentile"
    ):
        return runtime_params

    if val_cfg.max_exon_len_reference == "reference":
        print(
            "Warning: percentile reference mode 'reference' is not fully supported yet. Falling back to candidates_supported."
        )

    # Identify high-confidence transcripts: protein supported OR multi-source
    high_conf_tids = set()
    for (source, tid), grp in locus_df.groupby(["Source", "transcript_id"]):
        if tid in protein_supported_tids:
            high_conf_tids.add(tid)
            continue

        combined_ev = (
            grp["combined_evidence"].iloc[0] if "combined_evidence" in grp.columns else source
        )
        if len(combined_ev.split(",")) >= 2:
            high_conf_tids.add(tid)

    if not high_conf_tids:
        print(
            "Note: No high-confidence candidates found for percentile calculation. Using fixed limits."
        )
        return runtime_params

    hc_exons = locus_df[locus_df["transcript_id"].isin(high_conf_tids)]

    # Calculate exon lengths
    if val_cfg.max_exon_len_mode == "percentile":
        exon_lengths = hc_exons["End"] - hc_exons["Start"]
        if len(exon_lengths) > 0:
            pct_val = np.percentile(exon_lengths, val_cfg.max_exon_len_percentile)
            calc_limit = int(pct_val * val_cfg.max_exon_len_factor)
            # Use the more permissive of configured vs calculated
            runtime_params["effective_max_exon_len_bp"] = max(val_cfg.max_exon_len_bp, calc_limit)

    # Calculate transcript spans
    if val_cfg.max_transcript_span_mode == "percentile":
        spans = []
        for _tid, grp in hc_exons.groupby("transcript_id"):
            s = grp["End"].max() - grp["Start"].min()
            spans.append(s)
        if spans:
            pct_val = np.percentile(spans, val_cfg.max_transcript_span_percentile)
            calc_limit = int(pct_val * val_cfg.max_transcript_span_factor)
            runtime_params["effective_max_transcript_span_bp"] = max(
                val_cfg.max_transcript_span_bp, calc_limit
            )

    return runtime_params


def compute_utr_end_support(
    model: dict,
    locus_df: pd.DataFrame,
    config: PipelineConfig,
) -> dict[str, object]:
    """Determine if 5' and 3' ends are supported by other sources.

    Parameters
    ----------
    model : dict
        Candidate gene model dict with keys ``id``, ``strand``, ``start``,
        ``end``, and optionally ``protein_coding_score``.
    locus_df : pd.DataFrame
        Exon-level DataFrame for the enclosing locus.
    config : PipelineConfig
        Pipeline configuration with ``utr`` and ``protein_validation``
        sections.

    Returns
    -------
    dict
        Support results with keys ``supported_5p``, ``supported_3p``,
        ``reason_5p``, ``reason_3p``, ``action_5p``, ``action_3p``.
    """
    utr_cfg = config.utr
    res = {
        "supported_5p": False,
        "supported_3p": False,
        "reason_5p": "untested",
        "reason_3p": "untested",
        "action_5p": "kept",
        "action_3p": "kept",
    }

    if not utr_cfg.require_end_support:
        res["supported_5p"] = True
        res["supported_3p"] = True
        res["reason_5p"] = "support_not_required"
        res["reason_3p"] = "support_not_required"
        return res

    mode = utr_cfg.end_support_mode
    tid = model["id"]
    target_sources = set(utr_cfg.end_support_sources)
    tol = utr_cfg.end_tolerance_bp

    # Extract end coordinate candidates from locus (excluding self)
    other_models = locus_df[locus_df["transcript_id"] != tid]
    other_ends_5p = []
    other_ends_3p = []

    for _, grp in other_models.groupby(["Source", "transcript_id"]):
        src = grp["Source"].iloc[0]
        if src not in target_sources:
            continue
        strand = grp["Strand"].iloc[0]
        if strand != model["strand"]:
            continue

        # Determine 5p/3p based on strand
        min_c = grp["Start"].min()
        max_c = grp["End"].max()
        if strand == "+":
            other_ends_5p.append(min_c)
            other_ends_3p.append(max_c)
        else:
            other_ends_5p.append(max_c)
            other_ends_3p.append(min_c)

    # Model's own ends
    strand = model["strand"]
    min_c = model["start"]
    max_c = model["end"]
    m_5p = min_c if strand == "+" else max_c
    m_3p = max_c if strand == "+" else min_c

    # multisource_end_agreement
    multi_5p = any(abs(e - m_5p) <= tol for e in other_ends_5p)
    multi_3p = any(abs(e - m_3p) <= tol for e in other_ends_3p)

    # protein_validated
    prot_supported = False
    if config.protein_validation.enabled and "protein_coding_score" in model:
        if model["protein_coding_score"] >= config.protein_validation.min_score:
            prot_supported = True

    for end_type, is_multi, req_multi in [
        ("5p", multi_5p, utr_cfg.require_multisource_for_utr_5p),
        ("3p", multi_3p, utr_cfg.require_multisource_for_utr_3p),
    ]:
        supported = False
        reason = "unsupported"

        if mode == "multisource_end_agreement":
            if not req_multi or is_multi:
                supported = True
                reason = "multi_source_agreement" if req_multi else "multisource_not_required"
            else:
                reason = "no_end_agreement"
        elif mode == "protein_validated":
            if prot_supported:
                supported = True
                reason = "protein_validated"
            else:
                reason = "protein_score_low"
        elif mode == "either":
            if (not req_multi or is_multi) or prot_supported:
                supported = True
                reason = "either_rule_met"
            else:
                reason = "neither_rule_met"
        elif mode == "off":
            supported = True
            reason = "support_off"

        res[f"supported_{end_type}"] = supported
        res[f"reason_{end_type}"] = reason

        if not supported:
            if utr_cfg.fallback_policy_when_unsupported == "drop_utr":
                res[f"action_{end_type}"] = "dropped"
            elif utr_cfg.fallback_policy_when_unsupported == "hard_cap":
                res[f"action_{end_type}"] = "capped"
            elif utr_cfg.fallback_policy_when_unsupported == "drop_transcript":
                res[f"action_{end_type}"] = "drop_transcript"

    return res


def _write_resolved_config(cfg, output_dir: str) -> None:
    """Write the fully-resolved configuration to the output directory.

    Produces two files:
    - ``resolved_config.yaml`` — human-readable YAML of every setting
    - ``resolved_config_sha256`` — SHA-256 hex digest of that file
    """
    import hashlib

    resolved_path = os.path.join(output_dir, "resolved_config.yaml")
    sha_path = os.path.join(output_dir, "resolved_config_sha256")
    data = dump_config(cfg)
    text = yaml.dump(data, default_flow_style=False, sort_keys=True)
    with open(resolved_path, "w") as fh:
        fh.write(text)
    digest = hashlib.sha256(text.encode()).hexdigest()
    with open(sha_path, "w") as fh:
        fh.write(digest + "\n")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_setup_args = parser.add_argument_group("Setup")
    parser.add_setup_args.add_argument(
        "--config",
        action="append",
        help="Path to a YAML config override. May be repeated to layer several "
        "files in order (each later --config overrides earlier ones on the "
        "same keys); a single --config continues to work exactly as before.",
    )
    parser.add_setup_args.add_argument(
        "--preset",
        default="fungi",
        help=(
            "Clade preset to load before any --config overrides "
            "(default: 'fungi').  Use --list-presets to see what is "
            "installed.  Pass 'none' to use the neutral standard base only."
        ),
    )
    parser.add_setup_args.add_argument(
        "--list-presets",
        action="store_true",
        help="Print available clade presets and exit.",
    )
    parser.add_setup_args.add_argument(
        "--check-deps", action="store_true", help="Check external tool dependencies and exit"
    )

    parser.add_input_args = parser.add_argument_group("Inputs")
    parser.add_input_args.add_argument("--genome", help="Genome FASTA")
    parser.add_input_args.add_argument("--scallop", help="Scallop GTF")
    parser.add_input_args.add_argument("--stringtie", help="StringTie GTF")
    parser.add_input_args.add_argument(
        "--minimap2", help="Minimap2 long-read transcript alignments (GTF), optional"
    )
    parser.add_input_args.add_argument("--helixer", help="Helixer GFF3 (ab initio backbone)")
    parser.add_input_args.add_argument(
        "--tiberius",
        help="Tiberius GTF (ab initio backbone, alternative to --helixer; "
        "mutually exclusive with --helixer)",
    )
    parser.add_input_args.add_argument("--orthodb", help="OrthoDB GTF")
    parser.add_input_args.add_argument("--uniprot", help="UniProt GTF")
    parser.add_input_args.add_argument("--genblast", help="GenBlast protein alignment GTF")

    parser.add_output_args = parser.add_argument_group("Outputs")
    parser.add_output_args.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_output_args.add_argument("--gene-prefix", default="GENE", help="Prefix for new IDs")
    parser.add_output_args.add_argument(
        "--log-file",
        help="Pipeline log path. Defaults to <output-dir>/gmb.log",
    )
    parser.add_output_args.add_argument(
        "--no-log-file",
        action="store_true",
        help="Disable default file logging",
    )
    parser.add_output_args.add_argument(
        "--validate-fasta",
        action="store_true",
        help="Run FASTA QC checks at end of pipeline",
    )

    add_subset_args(parser)

    return parser.parse_args()


def load_evidence(
    path: str | None,
    source_label: str,
) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    """Load exon and CDS rows from a GTF or GFF3 file.

    Parameters
    ----------
    path : str or None
        Filesystem path to the annotation file. ``None`` or non-existent
        paths return ``(None, None)``.
    source_label : str
        Label to assign as the ``Source`` column (e.g. ``"Scallop"``).

    Returns
    -------
    tuple of (DataFrame or None, DataFrame or None)
        ``(exons_df, cds_df)``. Both are ``None`` when *path* is missing.
    """
    if not path or not os.path.exists(path):
        return None, None
    print(f"Loading {source_label} from {path}...")
    if path.endswith(".gtf"):
        df = pr.read_gtf(path).df
    else:
        df = pr.read_gff3(path).df

    df["Source"] = source_label

    if "transcript_id" not in df.columns:
        df["transcript_id"] = pd.NA

    m_mask = df["Feature"].isin(["mRNA", "transcript"])
    e_mask = df["Feature"].isin(["exon", "CDS"])

    if "ID" in df.columns:
        df.loc[m_mask, "transcript_id"] = df.loc[m_mask, "transcript_id"].fillna(
            df.loc[m_mask, "ID"]
        )
    if "Parent" in df.columns:
        df.loc[e_mask, "transcript_id"] = df.loc[e_mask, "transcript_id"].fillna(
            df.loc[e_mask, "Parent"]
        )

    # Cascade any remaining
    if "Parent" in df.columns:
        df["transcript_id"] = df["transcript_id"].fillna(df["Parent"])
    if "ID" in df.columns:
        df["transcript_id"] = df["transcript_id"].fillna(df["ID"])
    if "gene_id" in df.columns:
        df["transcript_id"] = df["transcript_id"].fillna(df["gene_id"])

    df["transcript_id"] = df["transcript_id"].fillna("unknown")

    if "gene_id" not in df.columns:
        if "Parent" in df.columns:
            df["gene_id"] = df["Parent"]
        else:
            df["gene_id"] = df["transcript_id"]

    # Prefix with source label to prevent ID collisions across tools
    df["transcript_id"] = f"{source_label}_" + df["transcript_id"].astype(str)
    df["gene_id"] = f"{source_label}_" + df["gene_id"].astype(str)

    exons = df[df["Feature"] == "exon"].copy()
    cds = df[df["Feature"] == "CDS"].copy()
    return exons, cds


def main() -> None:
    """Entry point for the Gene Model Builder pipeline.

    Parses CLI arguments, loads evidence, applies filtering, scoring,
    validation, deduplication, and writes final GFF3 / FASTA outputs.
    """
    args = parse_args()

    if args.list_presets:
        presets = list_build_presets()
        if presets:
            for p in presets:
                print(p)
        else:
            print("No presets installed.  Run 'pip install -e .' to install package data.")
        return

    os.makedirs(args.output_dir, exist_ok=True)
    log_file = resolve_log_file(args.output_dir, args.log_file, args.no_log_file)
    setup_logging(log_file=log_file, capture_stdio=log_file is not None)
    if log_file:
        print(f"Logging to {log_file}")

    if args.helixer and args.tiberius:
        sys.exit("ERROR: pass only one ab initio backbone: --helixer or --tiberius, not both.")
    backbone_path, backbone_label = (
        (args.tiberius, "Tiberius") if args.tiberius else (args.helixer, "Helixer")
    )

    config = load_config(args.config, args.preset)
    # Sync the scoring gate's backbone string to whichever ab initio track was
    # actually loaded, so weighting/single-exon-without-support logic (which
    # matches on this label) applies correctly regardless of --helixer vs
    # --tiberius.
    config.scoring.backbone_label = backbone_label

    # Write the fully-resolved configuration so runs are reproducible.
    _write_resolved_config(config, args.output_dir)

    if args.check_deps:
        check_dependencies(config.protein_validation)
        print("Dependency check passed.")
        sys.exit(0)

    print("Loading genome...")
    genome_dict = load_genome(args.genome)

    print("Loading transcriptomic evidence...")
    scallop_exons, scallop_cds = load_evidence(args.scallop, "Scallop")
    stringtie_exons, stringtie_cds = load_evidence(args.stringtie, "StringTie")
    minimap2_exons, minimap2_cds = load_evidence(args.minimap2, "Minimap2")

    print("Loading ab initio evidence...")
    helixer_exons, helixer_cds = load_evidence(backbone_path, backbone_label)

    print("Loading protein evidence...")
    orthodb_exons, orthodb_cds = load_evidence(args.orthodb, "OrthoDB")
    uniprot_exons, uniprot_cds = load_evidence(args.uniprot, "UniProt")
    genblast_exons, genblast_cds = load_evidence(args.genblast, "GenBlast")

    stats = {}

    # --- Seqname mapping ---
    mapping = build_mapping(
        assembly_report=getattr(args, "assembly_report", None),
        seqname_map=getattr(args, "seqname_map", None),
    )
    if mapping:
        print("Applying seqname mapping to evidence...")
        try:
            genome_dict = remap_genome_seqnames(genome_dict, mapping)
        except ValueError as exc:
            sys.exit(
                f"ERROR: {exc}\n\n"
                "This usually means --assembly-report was used for an assembly where "
                "multiple sequence records are assigned to the same chromosome label. "
                "Do not use that assembly report as a simple seqname map for GMB build. "
                "Use accession seqnames end-to-end, or provide a true pseudomolecule "
                "FASTA/GFF set with transformed coordinates."
            )
        scallop_exons = remap_df_seqnames(scallop_exons, mapping, "Scallop")
        scallop_cds = remap_df_seqnames(scallop_cds, mapping)
        stringtie_exons = remap_df_seqnames(stringtie_exons, mapping, "StringTie")
        stringtie_cds = remap_df_seqnames(stringtie_cds, mapping)
        minimap2_exons = remap_df_seqnames(minimap2_exons, mapping, "Minimap2")
        minimap2_cds = remap_df_seqnames(minimap2_cds, mapping)
        helixer_exons = remap_df_seqnames(helixer_exons, mapping, backbone_label)
        helixer_cds = remap_df_seqnames(helixer_cds, mapping)
        orthodb_exons = remap_df_seqnames(orthodb_exons, mapping, "OrthoDB")
        orthodb_cds = remap_df_seqnames(orthodb_cds, mapping)
        uniprot_exons = remap_df_seqnames(uniprot_exons, mapping, "UniProt")
        uniprot_cds = remap_df_seqnames(uniprot_cds, mapping)
        genblast_exons = remap_df_seqnames(genblast_exons, mapping, "GenBlast")
        genblast_cds = remap_df_seqnames(genblast_cds, mapping)

    print("Filtering Evidence...")
    tx_frames = [
        df
        for df in [scallop_exons, stringtie_exons, minimap2_exons]
        if df is not None and not df.empty
    ]
    tx_exons = pd.concat(tx_frames, ignore_index=True) if tx_frames else pd.DataFrame()
    tx_exons_filtered = filter_chimeras(tx_exons, config, stats)

    if config.transcript_splitting.split_enabled:
        tx_exons_filtered = split_mega_transcripts(tx_exons_filtered, config, stats)

    h_exons_filt, h_cds_filt = filter_backbone_models(
        helixer_exons if helixer_exons is not None else pd.DataFrame(),
        helixer_cds if helixer_cds is not None else pd.DataFrame(),
        config,
        stats,
    )

    prot_frames = [
        df
        for df in [orthodb_exons, uniprot_exons, genblast_exons]
        if df is not None and not df.empty
    ]
    prot_exons = pd.concat(prot_frames, ignore_index=True) if prot_frames else pd.DataFrame()
    prot_exons_filt = filter_protein_evidence(prot_exons, config, stats, tx_exons_filtered)

    # --- Region subsetting (fast test mode) ---
    subset_regions = None
    if getattr(args, "sample_loci", None) is not None:
        # Build preliminary loci from candidate exons for sampling

        _all_exons_for_loci = []
        if not tx_exons_filtered.empty:
            _all_exons_for_loci.append(tx_exons_filtered)
        if not h_exons_filt.empty:
            _all_exons_for_loci.append(h_exons_filt)
        if _all_exons_for_loci:
            from gmb.pipeline.subset_utils import _build_loci_from_exons

            _combined_ex = pd.concat(_all_exons_for_loci, ignore_index=True)
            _loci_df = _build_loci_from_exons(_combined_ex)
            print(f"  Built {len(_loci_df)} preliminary loci for sampling")
            subset_regions = resolve_subset_regions(args, loci_df=_loci_df)
    else:
        subset_regions = resolve_subset_regions(args)

    if subset_regions:
        _n_before = {
            "transcriptomic": (
                tx_exons_filtered["transcript_id"].nunique() if not tx_exons_filtered.empty else 0
            ),
            "helixer": h_exons_filt["transcript_id"].nunique() if not h_exons_filt.empty else 0,
            "protein": (
                prot_exons_filt["transcript_id"].nunique() if not prot_exons_filt.empty else 0
            ),
        }
        print(f"  Subsetting to {len(subset_regions)} region(s)...")
        tx_exons_filtered = subset_df_by_regions(tx_exons_filtered, subset_regions)
        h_exons_filt = subset_df_by_regions(h_exons_filt, subset_regions)
        h_cds_filt = subset_df_by_regions(h_cds_filt, subset_regions)
        prot_exons_filt = subset_df_by_regions(prot_exons_filt, subset_regions)
        # Also subset CDS tracks
        if scallop_cds is not None and not scallop_cds.empty:
            scallop_cds = subset_df_by_regions(scallop_cds, subset_regions)
        if stringtie_cds is not None and not stringtie_cds.empty:
            stringtie_cds = subset_df_by_regions(stringtie_cds, subset_regions)
        if minimap2_cds is not None and not minimap2_cds.empty:
            minimap2_cds = subset_df_by_regions(minimap2_cds, subset_regions)
        _n_after = {
            "transcriptomic": (
                tx_exons_filtered["transcript_id"].nunique() if not tx_exons_filtered.empty else 0
            ),
            "helixer": h_exons_filt["transcript_id"].nunique() if not h_exons_filt.empty else 0,
            "protein": (
                prot_exons_filt["transcript_id"].nunique() if not prot_exons_filt.empty else 0
            ),
        }
        for track, n_b in _n_before.items():
            n_a = _n_after[track]
            print(f"    {track}: {n_b} → {n_a} transcripts")
        manifest_path = os.path.join(args.output_dir, "subset_regions.tsv")
        write_subset_manifest(subset_regions, getattr(args, "seed", 1), manifest_path)

    all_dfs = []
    if not tx_exons_filtered.empty:
        all_dfs.append(tx_exons_filtered)
    if not h_exons_filt.empty:
        all_dfs.append(h_exons_filt)

    if not all_dfs:
        print("No evidence left after filtering. Exiting.")
        sys.exit(0)

    candidate_exons = pd.concat(all_dfs, ignore_index=True)

    cds_dfs = []
    for d in [scallop_cds, stringtie_cds, minimap2_cds, h_cds_filt]:
        if d is not None and not d.empty:
            cds_dfs.append(d)
    candidate_cds = pd.concat(cds_dfs, ignore_index=True) if cds_dfs else None

    print("Translating candidates for scoring/validation...")
    annotations = annotate_all_transcripts(
        candidate_exons, genome_dict, candidate_cds, min_codons=config.orf.min_codons
    )

    validation_scores = {}
    validation_details = {}
    if config.protein_validation.enabled:
        print("Running Protein Validation Stage...")
        check_dependencies(config.protein_validation)
        protein_dict = {tid: ann["protein"] for tid, ann in annotations.items() if ann["protein"]}
        validation_scores, validation_details = batch_score_proteins(protein_dict, config)

    if validation_scores:
        candidate_exons["protein_coding_score"] = (
            candidate_exons["transcript_id"].map(validation_scores).fillna(0.0)
        )

    # Protein-alignment support. `protein_supported_tids` is the boolean signal
    # that drives scoring (protein_overlap_bonus) and the retention gates in
    # select_isoforms. `protein_support_sources` records WHICH track supplied
    # that support, using the identical overlap criterion, purely so the support
    # can be reported in evidence_attribution.tsv -- previously this information
    # was computed and used but never surfaced anywhere.
    # Support requires the alignment to be on the SAME STRAND as the candidate.
    # A plain `pr_tx.overlap(pr_prot)` silently matches antisense alignments
    # because pyranges reports every read_gtf-loaded track as unstranded (see
    # gmb.utils.intervals.same_strand_overlap_ids for the mechanism).
    protein_supported_tids = set()
    protein_support_sources: dict[str, set[str]] = defaultdict(set)
    if not prot_exons_filt.empty and not candidate_exons.empty:
        protein_supported_tids = same_strand_overlap_ids(candidate_exons, prot_exons_filt)
        # Same criterion, per protein track, so attribution can name the source.
        for prot_label in prot_exons_filt["Source"].astype(str).unique():
            track_df = prot_exons_filt[prot_exons_filt["Source"].astype(str) == prot_label]
            if track_df.empty:
                continue
            for supported_tid in same_strand_overlap_ids(candidate_exons, track_df):
                protein_support_sources[supported_tid].add(prot_label)
    stats["protein_supported_candidates"] = len(protein_supported_tids)

    print("Clustering loci...")
    pr_candidates = pr.PyRanges(candidate_exons)
    clustered = pr_candidates.cluster(slack=0, count=True)
    cluster_df = clustered.df

    print("Computing guardrails...")
    runtime_params = compute_percentile_guardrails(cluster_df, config, protein_supported_tids)
    stats["effective_max_exon_len_bp"] = runtime_params["effective_max_exon_len_bp"]
    stats["effective_max_transcript_span_bp"] = runtime_params["effective_max_transcript_span_bp"]

    print("Scoring and Selecting Isoforms...")
    selected_gff_rows = []
    selected_cdna_fa = []
    selected_prot_fa = []
    unstranded_exclusion_rows: list[dict] = []

    gene_counter = 1

    for _cid, locus_df in cluster_df.groupby("Cluster"):
        genes = select_isoforms(
            locus_df,
            config,
            protein_supported_tids,
            genome_dict,
            protein_support_sources=protein_support_sources,
        )
        if not genes:
            continue

        for gene_models in genes:
            gene_id = f"{args.gene_prefix}_{gene_counter:05d}"
            gene_chrom = gene_models[0]["chrom"]
            gene_strand = gene_models[0]["strand"]

            # Exclude loci with unresolved strand: CDS coordinates and protein
            # translations cannot be computed reliably without a known strand.
            # Record every excluded transcript for attribution before skipping.
            if gene_strand not in ("+", "-"):
                for model in gene_models:
                    unstranded_exclusion_rows.append(
                        {
                            "gene_id": gene_id,
                            "candidate_transcript_id": model["id"],
                            "chrom": gene_chrom,
                            "strand": gene_strand,
                            "evidence_sources": model.get("combined_evidence", ""),
                            "rejection_reason": "UNRESOLVED_STRAND",
                        }
                    )
                gene_counter += 1
                continue

            # Collect all mRNA rows first so we can recompute gene span from children
            gene_mrna_rows = []
            gene_insert_idx = len(selected_gff_rows)  # position for gene row

            for i, model in enumerate(gene_models):
                tid = model["id"]
                new_tid = f"{gene_id}.t{i + 1}"
                ann = annotations.get(tid)

                # Use annotation exons (same ones used for ORF/CDS mapping) as source of truth
                if ann and ann.get("exons"):
                    exons_sorted = sorted(ann["exons"])
                else:
                    ex_df = model["df"]
                    exons_sorted = sorted(zip(ex_df["Start"].values, ex_df["End"].values))

                # Run UTR Support validation if UTRs exist
                drop_whole_transcript = False

                # Assign default UTR support values if compute_utr_end_support fails or hasn't run
                utr_support = {
                    "supported_5p": True,
                    "supported_3p": True,
                    "action_5p": "kept",
                    "action_3p": "kept",
                    "reason_5p": "default",
                    "reason_3p": "default",
                }

                try:
                    utr_support = compute_utr_end_support(model, cluster_df, config)
                except Exception as e:
                    print(f"Warning: UTR support computation failed for {tid}: {e}")

                # Check 5p
                if not utr_support["supported_5p"] and utr_support["action_5p"] == "dropped":
                    if ann:
                        ann["five_prime_utr"] = []
                elif (
                    not utr_support["supported_5p"]
                    and utr_support["action_5p"] == "drop_transcript"
                ):
                    drop_whole_transcript = True

                # Check 3p
                if not utr_support["supported_3p"] and utr_support["action_3p"] == "dropped":
                    if ann:
                        ann["three_prime_utr"] = []
                elif (
                    not utr_support["supported_3p"]
                    and utr_support["action_3p"] == "drop_transcript"
                ):
                    drop_whole_transcript = True

                if drop_whole_transcript:
                    stats.setdefault("utr_drop_transcripts", 0)
                    stats["utr_drop_transcripts"] += 1
                    continue

                # Save attribution support for QC output later
                model["utr_support"] = utr_support

                # Collect all child feature coordinates to recompute mRNA span
                all_child_coords = list(exons_sorted)
                if ann:
                    all_child_coords.extend(ann.get("cds", []))
                    all_child_coords.extend(ann.get("five_prime_utr", []))
                    all_child_coords.extend(ann.get("three_prime_utr", []))

                # mRNA span = union of all children
                if all_child_coords:
                    mrna_start = min(s for s, e in all_child_coords)
                    mrna_end = max(e for s, e in all_child_coords)
                else:
                    mrna_start = model["start"]
                    mrna_end = model["end"]

                mrna_row = {
                    "Chromosome": gene_chrom,
                    "Source": "GMB",
                    "Feature": "mRNA",
                    "Start": mrna_start,
                    "End": mrna_end,
                    "Score": ".",
                    "Strand": gene_strand,
                    "Frame": ".",
                    "ID": new_tid,
                    "Parent": gene_id,
                    "Evidence": model.get("combined_evidence", ""),
                    "ProteinEvidence": model.get("protein_evidence", ""),
                    "gmb_score": model.get("score"),
                }

                # Attach here (not re-looked-up later) because `annotations`
                # and `validation_details` are keyed by the original
                # candidate `tid`, not the renamed `new_tid` written to the
                # final GFF3 -- a later re-lookup by `new_tid` would silently
                # miss every row.
                if ann:
                    mrna_row["orf_label"] = ann.get("orf_label")
                    mrna_row["is_partial_5"] = ann.get("is_partial_5")
                    mrna_row["is_partial_3"] = ann.get("is_partial_3")
                    mrna_row["internal_stop_count"] = (ann.get("protein") or "").count("*")
                    mrna_row["protein_length"] = len(ann.get("protein") or "")
                mrna_row["protein_coding_score"] = validation_scores.get(tid)
                # `tid` is the ORIGINAL candidate transcript ID; `new_tid` is
                # the final renamed one. Both are recorded so the sidecar can
                # be joined either way (see protein_validation.tsv below).
                mrna_row["candidate_transcript_id"] = tid
                if tid in validation_details:
                    mrna_row["protein_validation_detail"] = validation_details[tid]

                if "utr_support" in model:
                    mrna_row["utr_support"] = model["utr_support"]

                gene_mrna_rows.append(mrna_row)
                selected_gff_rows.append(mrna_row)

                for j, (ex_s, ex_e) in enumerate(exons_sorted):
                    selected_gff_rows.append(
                        {
                            "Chromosome": gene_chrom,
                            "Source": "GMB",
                            "Feature": "exon",
                            "Start": ex_s,
                            "End": ex_e,
                            "Score": ".",
                            "Strand": gene_strand,
                            "Frame": ".",
                            "ID": f"{new_tid}.exon{j + 1}",
                            "Parent": new_tid,
                        }
                    )

                if ann:
                    cds_phases = compute_cds_phases(ann["cds"], gene_strand)
                    for j, (s, e) in enumerate(ann["cds"]):
                        selected_gff_rows.append(
                            {
                                "Chromosome": gene_chrom,
                                "Source": "GMB",
                                "Feature": "CDS",
                                "Start": s,
                                "End": e,
                                "Score": ".",
                                "Strand": gene_strand,
                                "Frame": str(cds_phases[j]) if j < len(cds_phases) else ".",
                                "ID": f"{new_tid}.cds{j + 1}",
                                "Parent": new_tid,
                            }
                        )
                    for j, (s, e) in enumerate(ann["five_prime_utr"]):
                        selected_gff_rows.append(
                            {
                                "Chromosome": gene_chrom,
                                "Source": "GMB",
                                "Feature": "five_prime_UTR",
                                "Start": s,
                                "End": e,
                                "Score": ".",
                                "Strand": gene_strand,
                                "Frame": ".",
                                "ID": f"{new_tid}.5utr{j + 1}",
                                "Parent": new_tid,
                            }
                        )
                    for j, (s, e) in enumerate(ann["three_prime_utr"]):
                        selected_gff_rows.append(
                            {
                                "Chromosome": gene_chrom,
                                "Source": "GMB",
                                "Feature": "three_prime_UTR",
                                "Start": s,
                                "End": e,
                                "Score": ".",
                                "Strand": gene_strand,
                                "Frame": ".",
                                "ID": f"{new_tid}.3utr{j + 1}",
                                "Parent": new_tid,
                            }
                        )

                    if ann["cdna"]:
                        selected_cdna_fa.append(f">{new_tid}\n{ann['cdna']}")
                    if ann["protein"]:
                        selected_prot_fa.append(f">{new_tid}\n{ann['protein']}")

            # Gene span = union of all mRNA children
            gene_start = min(r["Start"] for r in gene_mrna_rows)
            gene_end = max(r["End"] for r in gene_mrna_rows)

            # Insert gene row before its children
            selected_gff_rows.insert(
                gene_insert_idx,
                {
                    "Chromosome": gene_chrom,
                    "Source": "GMB",
                    "Feature": "gene",
                    "Start": gene_start,
                    "End": gene_end,
                    "Score": ".",
                    "Strand": gene_strand,
                    "Frame": ".",
                    "ID": gene_id,
                    "Parent": "",
                },
            )

            gene_counter += 1

    print("Writing Outputs...")

    # --- Post-processing: Validation + Dedup ---
    print("  Validating GFF3 structural integrity...")
    selected_gff_rows, val_stats = validate_and_fix_gff3(selected_gff_rows, config, runtime_params)
    stats.update({f"validation_{k}": v for k, v in val_stats.items()})
    if val_stats["violations_found"] > 0:
        print(
            f"    {val_stats['violations_found']} violations found, "
            f"{val_stats['transcripts_fixed']} fixed, "
            f"{val_stats['transcripts_dropped']} dropped, "
            f"{val_stats['exons_synthesized']} exons synthesized, "
            f"{val_stats['utrs_trimmed']} UTRs trimmed"
        )
    else:
        print("    All transcripts passed validation")

    print("  Deduplicating overlapping genes...")
    selected_gff_rows, dedup_stats = dedup_genes(selected_gff_rows, config)
    stats.update({f"dedup_{k}": v for k, v in dedup_stats.items()})
    if dedup_stats.get("genes_merged", 0) + dedup_stats.get("genes_dropped", 0) > 0:
        print(
            f"    {dedup_stats.get('genes_merged', 0)} merged, "
            f"{dedup_stats.get('genes_dropped', 0)} dropped, "
            f"{dedup_stats.get('genes_output', 0)} genes remaining"
        )
    else:
        print(f"    No duplicates found ({dedup_stats.get('genes_output', 0)} genes)")

    # --- Exact-duplicate transcript collapse ---
    # Runs strictly after dedup_genes() (which brings duplicate fragments
    # together as isoforms of one gene via gene-level overlap merging) and
    # before FASTA/evidence-attribution output, on the final isoform set --
    # see gmb.pipeline.duplicate_transcript_collapse's module docstring for the
    # full root-cause account and insertion-point rationale.
    print("  Collapsing exact-duplicate transcripts...")
    protein_by_tid = {
        rec.split("\n", 1)[0].lstrip(">"): rec.split("\n", 1)[1]
        for rec in selected_prot_fa
        if "\n" in rec
    }
    selected_gff_rows, collapse_log_rows, collapse_stats = collapse_exact_duplicate_transcripts(
        selected_gff_rows,
        config,
        protein_by_tid=protein_by_tid,
        genome_dict=genome_dict,
        protein_supported_tids=protein_supported_tids,
    )
    stats.update({f"duplicate_collapse_{k}": v for k, v in collapse_stats.items()})
    if collapse_stats["transcripts_removed"] > 0:
        print(
            f"    {collapse_stats['duplicate_groups_found']} exact-duplicate group(s), "
            f"{collapse_stats['transcripts_removed']} transcript(s) collapsed across "
            f"{collapse_stats['genes_with_collapse']} gene(s)"
        )
        collapse_path = os.path.join(args.output_dir, "collapsed_duplicate_transcripts.tsv")
        pd.DataFrame(collapse_log_rows).to_csv(collapse_path, sep="\t", index=False)
        print(f"    Collapse log: {collapse_path}")
    else:
        print("    No exact-duplicate transcripts found")

    # Final gene count
    final_gene_count = sum(1 for r in selected_gff_rows if r.get("Feature") == "gene")
    stats["total_loci"] = final_gene_count

    surviving_tids = {r["ID"] for r in selected_gff_rows if r.get("Feature") == "mRNA"}

    gff3_path = os.path.join(args.output_dir, "consensus.gff3")
    out_df = pd.DataFrame(selected_gff_rows)

    if not out_df.empty:
        # Convert to 1-based start for GFF3. pyranges is 0-based half-open (start is 0-based, end is 1-based)
        out_df["Start"] = out_df["Start"] + 1
        out_df = out_df.sort_values(["Chromosome", "Start"])
        with open(gff3_path, "w") as fh:
            fh.write("##gff-version 3\n")
            for _, r in out_df.iterrows():
                attr = f"ID={r['ID']}"
                if r["Parent"]:
                    attr += f";Parent={r['Parent']}"
                if "Evidence" in r and pd.notna(r["Evidence"]) and r["Evidence"] != "":
                    attr += f";Evidence={r['Evidence']}"
                if (
                    "ProteinEvidence" in r
                    and pd.notna(r["ProteinEvidence"])
                    and r["ProteinEvidence"] != ""
                ):
                    attr += f";ProteinEvidence={r['ProteinEvidence']}"
                if (
                    "CollapsedFrom" in r
                    and pd.notna(r["CollapsedFrom"])
                    and r["CollapsedFrom"] != ""
                ):
                    attr += f";CollapsedFrom={r['CollapsedFrom']}"
                fh.write(
                    f"{r['Chromosome']}\t{r['Source']}\t{r['Feature']}\t{r['Start']}\t{r['End']}\t{r['Score']}\t{r['Strand']}\t{r['Frame']}\t{attr}\n"
                )
    else:
        with open(gff3_path, "w") as fh:
            fh.write("##gff-version 3\n")

    print("  Regenerating final FASTA from post-processed GFF3...")
    regenerate_final_fasta(selected_gff_rows, genome_dict, args.output_dir)

    # --- Evidence attribution TSV ---
    evidence_rows = []

    # Pre-calculate children for fast access
    by_parent = defaultdict(list)
    for r in selected_gff_rows:
        if r.get("Parent"):
            by_parent[r["Parent"]].append(r)

    # Iterate ALL output mRNAs directly
    output_mrnas = [r for r in selected_gff_rows if r.get("Feature") == "mRNA"]

    # Extract original models to map UTR support metadata
    for _cid, locus_df in cluster_df.groupby("Cluster"):
        genes = select_isoforms(
            locus_df,
            config,
            protein_supported_tids,
            genome_dict,
            protein_support_sources=protein_support_sources,
        )
        if genes:
            for g in genes:
                for _idx, _m in enumerate(g):
                    # Emitted ID matches `gene_id.t{idx+1}` logic
                    # To accurately find it here we need to reconstruct what it was called or just use order, but the logic in select isoforms creates it.
                    pass  # We will instead pass it via a side-channel or Evidence field since tid changes

    # Better approach: parse Evidence field or just store it in mRNA row temporarily
    {m["ID"]: m for m in output_mrnas}

    for m in output_mrnas:
        tid = m["ID"]
        gene_id = m["Parent"]

        children = by_parent[tid]
        cds_rows = [c for c in children if c["Feature"] == "CDS"]
        exon_rows = [c for c in children if c["Feature"] == "exon"]
        utr5_rows = [c for c in children if c["Feature"] == "five_prime_UTR"]
        utr3_rows = [c for c in children if c["Feature"] == "three_prime_UTR"]

        cds_bp = sum(c["End"] - c["Start"] for c in cds_rows)
        utr5_bp = sum(c["End"] - c["Start"] for c in utr5_rows)
        utr3_bp = sum(c["End"] - c["Start"] for c in utr3_rows)

        # Max exon and Intron len
        exon_lens = [c["End"] - c["Start"] for c in exon_rows]
        max_exon_len_bp = max(exon_lens) if exon_lens else 0

        sorted_exons = sorted(exon_rows, key=lambda x: x["Start"])
        intron_lens = []
        for i in range(len(sorted_exons) - 1):
            intron_lens.append(sorted_exons[i + 1]["Start"] - sorted_exons[i]["End"])
        max_intron_len_bp = max(intron_lens) if intron_lens else 0

        transcript_span_bp = (
            sorted_exons[-1]["End"] - sorted_exons[0]["Start"] if sorted_exons else 0
        )

        evidence_sources = m.get("Evidence", "")
        # Protein-alignment tracks that support this transcript. These are
        # supporting evidence, never candidate models, so they are reported in
        # their own column rather than folded into evidence_sources (which
        # names the tracks the structure was BUILT from).
        protein_alignment_sources = m.get("ProteinEvidence", "")

        row_dict = {
            "gene_id": gene_id,
            "transcript_id": tid,
            "evidence_sources": evidence_sources,
            "protein_alignment_sources": protein_alignment_sources,
            "has_protein_alignment_support": bool(protein_alignment_sources),
            "exon_count": len(exon_rows),
            "cds_bp": cds_bp,
            "utr_5p_bp": utr5_bp,
            "utr_3p_bp": utr3_bp,
            "max_exon_len_bp": max_exon_len_bp,
            "max_intron_len_bp": max_intron_len_bp,
            "transcript_span_bp": transcript_span_bp,
            "gmb_score": m.get("gmb_score"),
        }

        # Extract support from the mRNA dict (injected during gene_model_builder formulation)
        if "utr_support" in m:
            s_dict = m["utr_support"]
            row_dict.update(
                {
                    "utr_5p_supported": s_dict.get("supported_5p", True),
                    "utr_3p_supported": s_dict.get("supported_3p", True),
                    "utr_5p_action": s_dict.get("action_5p", "kept"),
                    "utr_3p_action": s_dict.get("action_3p", "kept"),
                    "utr_5p_reason": s_dict.get("reason_5p", "default"),
                    "utr_3p_reason": s_dict.get("reason_3p", "default"),
                }
            )

        evidence_rows.append(row_dict)

    if evidence_rows:
        ev_df = pd.DataFrame(evidence_rows)
        ev_path = os.path.join(args.output_dir, "evidence_attribution.tsv")
        ev_df.to_csv(ev_path, sep="\t", index=False)
        print(f"  Evidence attribution: {ev_path}")

    if unstranded_exclusion_rows:
        excl_df = pd.DataFrame(unstranded_exclusion_rows)
        excl_path = os.path.join(args.output_dir, "unstranded_exclusions.tsv")
        excl_df.to_csv(excl_path, sep="\t", index=False)
        n_excl = len(unstranded_exclusion_rows)
        stats["unstranded_exclusions"] = n_excl
        print(
            f"  Unstranded exclusions: {n_excl} transcript(s) excluded (UNRESOLVED_STRAND). "
            f"See {excl_path}"
        )

    # --- Protein validation TSV (per-transcript DIAMOND/Psauron detail) ---
    # Written whenever the stage ran, independent of whether any row_dict
    # above needed it -- this is the compact table canonical-transcript
    # selection (gmb.pipeline.canonical_selection) and manual QC read from.
    #
    # Percentage-scale note: diamond_pident/diamond_qcov/diamond_scov are
    # 0-100 (DIAMOND's own native percentage convention); psauron_score and
    # protein_coding_score are 0-1 (Psauron's own convention, and the
    # diamond/psauron weighted-combination result -- see
    # ProteinValidationConfig.diamond_weight/psauron_weight in config.py).
    # gmb_score is open-ended (can be negative) and is not a percentage at
    # all -- it's the additive scoring.py isoform score.
    if config.protein_validation.enabled:
        protein_val_rows = []
        for m in output_mrnas:
            tid = m["ID"]
            gene_id = m["Parent"]
            detail = m.get("protein_validation_detail", {})
            protein_val_rows.append(
                {
                    "gene_id": gene_id,
                    "transcript_id": tid,
                    # Original pre-rename candidate ID: lets this sidecar be
                    # joined to anything keyed on candidate IDs, and lets a
                    # consumer verify which candidate produced the result.
                    "candidate_transcript_id": m.get("candidate_transcript_id"),
                    # Sequence-identity reuse key -- stable across the
                    # candidate->final transcript rename, and the key used to
                    # deduplicate identical proteins for InterPro review.
                    "protein_sha256": detail.get("protein_sha256"),
                    "diamond_hit": detail.get("diamond_hit"),
                    "diamond_pident": detail.get("diamond_pident"),
                    "diamond_qcov": detail.get("diamond_qcov"),
                    "diamond_scov": detail.get("diamond_scov"),
                    "diamond_bitscore": detail.get("diamond_bitscore"),
                    "diamond_evalue": detail.get("diamond_evalue"),
                    "psauron_score": detail.get("psauron_score"),
                    "protein_length": m.get("protein_length"),
                    "orf_label": m.get("orf_label"),
                    "is_partial_5": m.get("is_partial_5"),
                    "is_partial_3": m.get("is_partial_3"),
                    "internal_stop_count": m.get("internal_stop_count"),
                    "protein_coding_score": m.get("protein_coding_score"),
                    "gmb_score": m.get("gmb_score"),
                    "evidence_sources": m.get("Evidence", ""),
                    # Provenance: what produced this row. A consumer can
                    # decide whether the result is still reusable without
                    # rerunning DIAMOND/psauron.
                    "validation_status": detail.get("validation_status"),
                    "validation_reason": detail.get("validation_reason"),
                    # Whether DIAMOND/psauron actually ran for THIS row, or
                    # whether the row is sharing another transcript's result
                    # because they translate to an identical protein (see
                    # protein_sha256 above, and protein_validation.py's
                    # batch_score_proteins doc for the "compute once, consume
                    # many" contract this makes auditable per-row).
                    "protein_validation_source": detail.get("protein_validation_source"),
                    "protein_validation_reused_from": detail.get("protein_validation_reused_from"),
                    "diamond_version": detail.get("diamond_version"),
                    "psauron_version": detail.get("psauron_version"),
                    "diamond_db": detail.get("diamond_db"),
                }
            )
        if protein_val_rows:
            pv_df = pd.DataFrame(protein_val_rows)
            pv_path = os.path.join(args.output_dir, "protein_validation.tsv")
            pv_df.to_csv(pv_path, sep="\t", index=False)
            print(f"  Protein validation detail: {pv_path}")

    def _default_encoder(obj):
        if hasattr(obj, "item"):
            return obj.item()
        raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

    summary_out = {
        "summary": stats,
        "filtering": stats,
        "utr": {"transcripts_dropped": stats.get("utr_drop_transcripts", 0)},
    }

    with open(os.path.join(args.output_dir, "summary.json"), "w") as fh:
        json.dump(summary_out, fh, indent=2, default=_default_encoder)

    with open(os.path.join(args.output_dir, "summary.tsv"), "w") as fh:
        fh.write("Metric\tValue\n")
        for k, v in stats.items():
            fh.write(f"{k}\t{v}\n")

    if getattr(args, "validate_fasta", False):
        from gmb.pipeline.fasta_qc import print_report, validate_fasta

        genome_path = getattr(args, "genome", None)
        qc_report = validate_fasta(args.output_dir, genome_path)
        print_report(qc_report)
        report_path = os.path.join(args.output_dir, "fasta_qc_report.json")
        with open(report_path, "w") as fh:
            json.dump(qc_report, fh, indent=2)
        if not qc_report.get("pass", False):
            failed = qc_report.get("failed_checks", [])
            print(f"ERROR: FASTA QC failed: {', '.join(failed)}")
            sys.exit(1)

    print("Done!")


if __name__ == "__main__":
    main()
