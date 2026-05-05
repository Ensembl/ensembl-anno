#!/usr/bin/env python3
"""Re-translate CDS from consensus GFF3 using the fixed ascending-order code.

Usage:
    python retranslate_from_gff3.py \
        --gff3  output_full/consensus.gff3 \
        --genome z_tritici/zymoseptoria_tritici.fa \
        --out   output_fixed/prot.fa
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from collections import defaultdict

# Use the fixed annotate_cds_utrs functions from the gmb directory
sys.path.insert(0, os.path.dirname(__file__))
from annotate_cds_utrs import reverse_complement, translate


def load_genome(fasta_path: str) -> dict[str, str]:
    seq: dict[str, str] = {}
    name = None
    parts: list[str] = []
    opener = gzip.open if fasta_path.endswith(".gz") else open
    with opener(fasta_path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seq[name] = "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        seq[name] = "".join(parts).upper()
    print(f"  Loaded {len(seq)} sequences from genome", file=sys.stderr)
    return seq


def parse_gff3(gff3_path: str):
    """Yield (chrom, strand, tid, gene_id, evidence, cds_intervals_0based)."""
    opener = gzip.open if gff3_path.endswith(".gz") else open

    # Pass 1: collect mRNA metadata
    mrna_meta: dict[str, dict] = {}  # tid -> {chrom, strand, gene_id, evidence}
    # Pass 1 & 2 in one scan: collect mRNA meta then CDS
    cds_by_tid: dict[str, list[tuple[int, int]]] = defaultdict(list)

    with opener(gff3_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feat, start_s, end_s, score, strand, phase, attrs = parts
            start_1 = int(start_s)
            end_1 = int(end_s)
            # Convert to 0-based half-open
            start_0 = start_1 - 1
            end_0 = end_1  # end is already 1-based inclusive → 0-based exclusive = same value

            # Parse attrs
            attr_dict: dict[str, str] = {}
            for kv in attrs.split(";"):
                kv = kv.strip()
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attr_dict[k] = v

            if feat == "mRNA":
                tid = attr_dict.get("ID", "")
                gene_id = attr_dict.get("Parent", "")
                evidence = attr_dict.get("Evidence", "")
                mrna_meta[tid] = {
                    "chrom": chrom,
                    "strand": strand,
                    "gene_id": gene_id,
                    "evidence": evidence,
                    "start": start_0,
                    "end": end_0,
                }
            elif feat == "CDS":
                parent = attr_dict.get("Parent", "")
                cds_by_tid[parent].append((start_0, end_0))

    return mrna_meta, cds_by_tid


def build_cds_seq(cds_intervals: list[tuple[int, int]], strand: str, chrom_seq: str) -> str:
    """Concatenate CDS intervals in ascending genomic order, then RC for minus strand."""
    parts = [chrom_seq[cs:ce] for cs, ce in sorted(cds_intervals)]
    cds_nuc = "".join(parts)
    if strand == "-":
        cds_nuc = reverse_complement(cds_nuc)
    return cds_nuc


def write_seq(fh, seq: str, line_width: int = 80) -> None:
    for i in range(0, len(seq), line_width):
        fh.write(seq[i : i + line_width] + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff3", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-aa", type=int, default=10)
    args = ap.parse_args()

    print("Loading genome…", file=sys.stderr)
    genome = load_genome(args.genome)

    print("Parsing GFF3…", file=sys.stderr)
    mrna_meta, cds_by_tid = parse_gff3(args.gff3)
    print(f"  {len(mrna_meta)} mRNAs, {len(cds_by_tid)} with CDS", file=sys.stderr)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    stats = {"total": 0, "with_cds": 0, "with_protein": 0, "stops_removed": 0}

    with open(args.out, "w") as fh:
        for tid, meta in mrna_meta.items():
            intervals = cds_by_tid.get(tid)
            if not intervals:
                continue

            chrom = meta["chrom"]
            strand = meta["strand"]
            gene_id = meta["gene_id"]
            evidence = meta["evidence"]

            if chrom not in genome:
                continue

            stats["total"] += 1
            stats["with_cds"] += 1

            chrom_seq = genome[chrom]
            cds_nuc = build_cds_seq(intervals, strand, chrom_seq)
            protein = translate(cds_nuc)

            # Strip trailing stop
            if protein.endswith("*"):
                protein = protein[:-1]
                stats["stops_removed"] += 1

            if len(protein) < args.min_aa:
                continue

            stats["with_protein"] += 1
            aa_len = len(protein)
            header = (
                f">{tid} gene={gene_id} aa_len={aa_len} "
                f"evidence={evidence} strand={strand} "
                f"cds_exons={len(intervals)}"
            )
            fh.write(header + "\n")
            write_seq(fh, protein)

    print(
        f"Done: {stats['total']} CDS transcripts → {stats['with_protein']} proteins written to {args.out}",
        file=sys.stderr,
    )
    print(f"  ({stats['stops_removed']} trailing stop codons stripped)", file=sys.stderr)


if __name__ == "__main__":
    main()
