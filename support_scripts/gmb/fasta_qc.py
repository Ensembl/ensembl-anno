#!/usr/bin/env python3
"""FASTA QC validation for Gene Model Builder outputs.

Checks that prot.fa and cdna.fa are consistent with consensus.gff3:
  - Coverage: every CDS-bearing mRNA has a protein; every mRNA has a cDNA
  - ID mapping: FASTA headers match GFF transcript IDs exactly
  - No ghost records: no FASTA entries for transcripts absent from GFF
  - No duplicates: no repeated FASTA headers
  - Sequence correctness: reconstructed translations match prot.fa

Can be run standalone or called from gene_model_builder with --validate-fasta.
"""

from __future__ import annotations

import json
import os
import re
import sys
from collections import Counter

# Allow imports from the gmb package directory
sys.path.insert(0, os.path.dirname(__file__))

from annotate_cds_utrs import build_spliced_seq, reverse_complement, translate


def parse_gff3(gff3_path: str) -> dict:
    """Parse consensus GFF3 into structured transcript data.

    Returns
    -------
    dict with keys:
        genes: dict of gene_id -> {chrom, strand, start, end}
        transcripts: dict of tid -> {gene_id, chrom, strand, exons, cds}
        n_genes, n_transcripts, n_cds_transcripts
    """
    genes = {}
    transcripts = {}
    exons_by_parent = {}
    cds_by_parent = {}

    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, _src, feat, start, end, _score, strand, frame, attrs_str = parts
            start = int(start)
            end = int(end)

            attrs = {}
            for kv in attrs_str.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k] = v

            if feat == "gene":
                genes[attrs["ID"]] = {
                    "chrom": chrom,
                    "strand": strand,
                    "start": start,
                    "end": end,
                }
            elif feat == "mRNA":
                tid = attrs["ID"]
                parent = attrs.get("Parent", "")
                transcripts[tid] = {
                    "gene_id": parent,
                    "chrom": chrom,
                    "strand": strand,
                    "exons": [],
                    "cds": [],
                }
            elif feat == "exon":
                parent = attrs.get("Parent", "")
                exons_by_parent.setdefault(parent, []).append(
                    (start - 1, end)  # Convert to 0-based half-open
                )
            elif feat == "CDS":
                parent = attrs.get("Parent", "")
                frame_val = frame if frame != "." else "0"
                cds_by_parent.setdefault(parent, []).append(
                    (start - 1, end, int(frame_val))  # 0-based half-open + phase
                )

    # Attach exons and CDS to transcripts
    for tid, tx in transcripts.items():
        tx["exons"] = sorted(exons_by_parent.get(tid, []))
        raw_cds = sorted(cds_by_parent.get(tid, []), key=lambda x: x[0])
        tx["cds"] = [(s, e) for s, e, _p in raw_cds]
        tx["cds_phases"] = [p for _s, _e, p in raw_cds]

    cds_tids = {tid for tid, tx in transcripts.items() if tx["cds"]}

    return {
        "genes": genes,
        "transcripts": transcripts,
        "n_genes": len(genes),
        "n_transcripts": len(transcripts),
        "n_cds_transcripts": len(cds_tids),
        "cds_transcript_ids": cds_tids,
    }


def parse_fasta_ids(fasta_path: str) -> list[str]:
    """Extract record IDs (first whitespace-delimited token after '>') from FASTA."""
    ids = []
    if not os.path.exists(fasta_path):
        return ids
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0].strip())
    return ids


def parse_fasta_records(fasta_path: str) -> dict[str, str]:
    """Parse FASTA into {id: sequence} dict."""
    records = {}
    if not os.path.exists(fasta_path):
        return records
    current_id = None
    parts = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = "".join(parts)
                current_id = line[1:].split()[0].strip()
                parts = []
            else:
                parts.append(line.strip())
    if current_id is not None:
        records[current_id] = "".join(parts)
    return records


def load_genome(genome_path: str) -> dict[str, str]:
    """Load genome FASTA into {chrom: sequence} dict."""
    genome = {}
    current_id = None
    parts = []
    with open(genome_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current_id is not None:
                    genome[current_id] = "".join(parts)
                current_id = line[1:].split()[0].strip()
                parts = []
            else:
                parts.append(line.strip().upper())
    if current_id is not None:
        genome[current_id] = "".join(parts)
    return genome


def run_coverage_checks(
    gff_data: dict,
    prot_ids: list[str],
    cdna_ids: list[str],
) -> dict:
    """Part A1: coverage and ID consistency checks."""
    gff_mrna_ids = set(gff_data["transcripts"].keys())
    cds_tids = gff_data["cds_transcript_ids"]
    prot_id_set = set(prot_ids)
    cdna_id_set = set(cdna_ids)

    # Missing proteins: CDS-bearing mRNAs not in prot.fa
    missing_proteins = sorted(cds_tids - prot_id_set)

    # Extra proteins: in prot.fa but not an mRNA in GFF
    extra_proteins = sorted(prot_id_set - gff_mrna_ids)

    # Missing cDNA: mRNAs not in cdna.fa
    missing_cdna = sorted(gff_mrna_ids - cdna_id_set)

    # Extra cDNA
    extra_cdna = sorted(cdna_id_set - gff_mrna_ids)

    # Duplicate headers
    prot_dupes = {k: v for k, v in Counter(prot_ids).items() if v > 1}
    cdna_dupes = {k: v for k, v in Counter(cdna_ids).items() if v > 1}

    return {
        "n_genes": gff_data["n_genes"],
        "n_transcripts": gff_data["n_transcripts"],
        "n_cds_transcripts": gff_data["n_cds_transcripts"],
        "n_prot_records": len(prot_ids),
        "n_cdna_records": len(cdna_ids),
        "missing_proteins": missing_proteins[:50],
        "missing_proteins_total": len(missing_proteins),
        "extra_proteins": extra_proteins[:50],
        "extra_proteins_total": len(extra_proteins),
        "missing_cdna": missing_cdna[:50],
        "missing_cdna_total": len(missing_cdna),
        "extra_cdna": extra_cdna[:50],
        "extra_cdna_total": len(extra_cdna),
        "duplicate_prot_headers": dict(list(prot_dupes.items())[:50]),
        "duplicate_cdna_headers": dict(list(cdna_dupes.items())[:50]),
    }


def run_sequence_checks(
    gff_data: dict,
    prot_records: dict[str, str],
    cdna_records: dict[str, str],
    genome: dict[str, str],
    max_checks: int = 0,
) -> dict:
    """Part A2: sequence correctness checks.

    Parameters
    ----------
    max_checks : int
        Max transcripts to check. 0 = check all.
    """
    transcripts = gff_data["transcripts"]
    cds_tids = gff_data["cds_transcript_ids"]

    protein_mismatches = []
    cdna_mismatches = []
    internal_stops = []
    frame_issues = []

    tids_to_check = sorted(cds_tids)
    if max_checks > 0 and len(tids_to_check) > max_checks:
        import random

        random.seed(42)
        tids_to_check = sorted(random.sample(tids_to_check, max_checks))

    for tid in tids_to_check:
        tx = transcripts[tid]
        chrom = tx["chrom"]
        strand = tx["strand"]
        exons = tx["exons"]
        cds_intervals = tx["cds"]

        if chrom not in genome:
            continue

        chrom_seq = genome[chrom]

        # --- cDNA check ---
        if tid in cdna_records and exons:
            expected_cdna = build_spliced_seq(exons, strand, chrom_seq)
            observed_cdna = cdna_records[tid]
            if expected_cdna != observed_cdna:
                mismatch_pos = next(
                    (i for i, (a, b) in enumerate(zip(expected_cdna, observed_cdna)) if a != b),
                    min(len(expected_cdna), len(observed_cdna)),
                )
                cdna_mismatches.append(
                    {
                        "transcript_id": tid,
                        "expected_length": len(expected_cdna),
                        "observed_length": len(observed_cdna),
                        "first_mismatch_pos": mismatch_pos,
                    }
                )

        # --- Protein check ---
        if tid in prot_records and cds_intervals:
            # Reconstruct CDS sequence from genome
            # Concatenate in ascending genomic order, then RC for minus strand.
            # RC(X+Y)==RC(Y)+RC(X): ascending order correctly places the
            # highest-coordinate (5') exon first after RC for minus-strand CDS.
            cds_parts = [chrom_seq[s:e] for s, e in sorted(cds_intervals)]
            cds_nuc = "".join(cds_parts)
            if strand == "-":
                cds_nuc = reverse_complement(cds_nuc)

            expected_protein = translate(cds_nuc)
            if expected_protein.endswith("*"):
                expected_protein = expected_protein[:-1]

            observed_protein = prot_records[tid]

            # Check for internal stops
            if "*" in observed_protein:
                internal_stops.append(
                    {
                        "transcript_id": tid,
                        "stop_positions": [
                            i for i, c in enumerate(observed_protein) if c == "*"
                        ],
                    }
                )

            if expected_protein != observed_protein:
                mismatch_pos = next(
                    (
                        i
                        for i, (a, b) in enumerate(zip(expected_protein, observed_protein))
                        if a != b
                    ),
                    min(len(expected_protein), len(observed_protein)),
                )
                protein_mismatches.append(
                    {
                        "transcript_id": tid,
                        "expected_length": len(expected_protein),
                        "observed_length": len(observed_protein),
                        "first_mismatch_pos": mismatch_pos,
                        "strand": strand,
                    }
                )

    return {
        "transcripts_checked": len(tids_to_check),
        "protein_mismatches": protein_mismatches[:50],
        "protein_mismatches_total": len(protein_mismatches),
        "cdna_mismatches": cdna_mismatches[:50],
        "cdna_mismatches_total": len(cdna_mismatches),
        "internal_stops": internal_stops[:50],
        "internal_stops_total": len(internal_stops),
        "frame_issues": frame_issues[:50],
    }


def validate_fasta(
    output_dir: str,
    genome_path: str | None = None,
    max_seq_checks: int = 0,
) -> dict:
    """Run all FASTA QC checks and return report dict.

    Parameters
    ----------
    output_dir : str
        Directory containing consensus.gff3, prot.fa, cdna.fa.
    genome_path : str or None
        Path to genome FASTA. Required for sequence correctness checks.
    max_seq_checks : int
        Max transcripts for sequence comparison. 0 = all.

    Returns
    -------
    dict : QC report.
    """
    gff3_path = os.path.join(output_dir, "consensus.gff3")
    prot_path = os.path.join(output_dir, "prot.fa")
    cdna_path = os.path.join(output_dir, "cdna.fa")

    if not os.path.exists(gff3_path):
        return {"error": f"consensus.gff3 not found in {output_dir}"}

    # Parse inputs
    gff_data = parse_gff3(gff3_path)
    prot_ids = parse_fasta_ids(prot_path)
    cdna_ids = parse_fasta_ids(cdna_path)

    # A1: Coverage checks
    report = run_coverage_checks(gff_data, prot_ids, cdna_ids)

    # A2: Sequence checks (if genome provided)
    if genome_path and os.path.exists(genome_path):
        genome = load_genome(genome_path)
        prot_records = parse_fasta_records(prot_path)
        cdna_records = parse_fasta_records(cdna_path)
        seq_report = run_sequence_checks(
            gff_data, prot_records, cdna_records, genome, max_seq_checks
        )
        report["sequence_checks"] = seq_report

    # Determine pass/fail
    report["pass"] = (
        report["missing_proteins_total"] == 0
        and report["extra_proteins_total"] == 0
        and report["missing_cdna_total"] == 0
        and report["extra_cdna_total"] == 0
        and len(report["duplicate_prot_headers"]) == 0
        and len(report["duplicate_cdna_headers"]) == 0
    )

    return report


def print_report(report: dict) -> None:
    """Print a human-readable summary to stdout."""
    if "error" in report:
        print(f"ERROR: {report['error']}")
        return

    status = "PASS" if report["pass"] else "FAIL"
    print(f"\n=== FASTA QC Report [{status}] ===")
    print(f"  Genes in GFF:            {report['n_genes']}")
    print(f"  Transcripts in GFF:      {report['n_transcripts']}")
    print(f"  CDS-bearing transcripts: {report['n_cds_transcripts']}")
    print(f"  Protein FASTA records:   {report['n_prot_records']}")
    print(f"  cDNA FASTA records:      {report['n_cdna_records']}")

    if report["missing_proteins_total"] > 0:
        print(f"\n  MISSING proteins ({report['missing_proteins_total']} total):")
        for tid in report["missing_proteins"][:10]:
            print(f"    - {tid}")
        if report["missing_proteins_total"] > 10:
            print(f"    ... and {report['missing_proteins_total'] - 10} more")

    if report["extra_proteins_total"] > 0:
        print(f"\n  EXTRA proteins ({report['extra_proteins_total']} total):")
        for tid in report["extra_proteins"][:10]:
            print(f"    - {tid}")
        if report["extra_proteins_total"] > 10:
            print(f"    ... and {report['extra_proteins_total'] - 10} more")

    if report["missing_cdna_total"] > 0:
        print(f"\n  MISSING cDNA ({report['missing_cdna_total']} total):")
        for tid in report["missing_cdna"][:10]:
            print(f"    - {tid}")

    if report["extra_cdna_total"] > 0:
        print(f"\n  EXTRA cDNA ({report['extra_cdna_total']} total):")
        for tid in report["extra_cdna"][:10]:
            print(f"    - {tid}")

    if report["duplicate_prot_headers"]:
        print(f"\n  DUPLICATE prot headers: {len(report['duplicate_prot_headers'])}")

    if report["duplicate_cdna_headers"]:
        print(f"\n  DUPLICATE cdna headers: {len(report['duplicate_cdna_headers'])}")

    if "sequence_checks" in report:
        sc = report["sequence_checks"]
        print(f"\n  Sequence checks ({sc['transcripts_checked']} transcripts):")
        print(f"    Protein mismatches:  {sc['protein_mismatches_total']}")
        print(f"    cDNA mismatches:     {sc['cdna_mismatches_total']}")
        print(f"    Internal stops:      {sc['internal_stops_total']}")
        if sc["protein_mismatches"]:
            print("    Sample protein mismatches:")
            for m in sc["protein_mismatches"][:5]:
                print(
                    f"      {m['transcript_id']}: expected {m['expected_length']}aa, "
                    f"got {m['observed_length']}aa, first diff at pos {m['first_mismatch_pos']}"
                )

    print()


def main():
    import argparse

    parser = argparse.ArgumentParser(description="FASTA QC for Gene Model Builder outputs")
    parser.add_argument("output_dir", help="Directory with consensus.gff3, prot.fa, cdna.fa")
    parser.add_argument("--genome", help="Genome FASTA for sequence correctness checks")
    parser.add_argument(
        "--max-seq-checks",
        type=int,
        default=0,
        help="Max transcripts for sequence checks (0=all)",
    )
    args = parser.parse_args()

    report = validate_fasta(args.output_dir, args.genome, args.max_seq_checks)
    print_report(report)

    # Write JSON report
    report_path = os.path.join(args.output_dir, "fasta_qc_report.json")
    with open(report_path, "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"Report written to {report_path}")

    sys.exit(0 if report.get("pass", False) else 1)


if __name__ == "__main__":
    main()
