#!/usr/bin/env python3
"""
Produce a fixed prot.fa by re-translating Helixer-sourced CDS with the corrected
ascending-sort code, while keeping the original proteins for non-Helixer transcripts.

The bug (descending sort before RC on minus strand) only affects transcripts whose
CDS intervals come from Helixer.  For ORF-predicted (transcriptome-only) models the
CDS is derived from the correctly-oriented cdna, so the sort order had no visible
effect on those proteins.

Usage:
    python make_fixed_prot.py \
        --gff3     output_full/consensus.gff3 \
        --genome   z_tritici/zymoseptoria_tritici.fa \
        --old-prot output_full/prot.fa \
        --out      output_fixed/prot_fixed.fa
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(__file__))
from annotate_cds_utrs import reverse_complement, translate


# ── helpers ──────────────────────────────────────────────────────────────────

def load_genome(path: str) -> dict[str, str]:
    seq: dict[str, str] = {}
    name = None
    parts: list[str] = []
    with open(path) as fh:
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
    return seq


def load_old_prot(path: str) -> dict[str, str]:
    """Return {tid: sequence} from a FASTA file (single or multi-line)."""
    prots: dict[str, str] = {}
    tid = None
    parts: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if tid is not None:
                    prots[tid] = "".join(parts)
                tid = line.split()[0][1:]
                parts = []
            else:
                parts.append(line)
    if tid is not None:
        prots[tid] = "".join(parts)
    return prots


def parse_gff3(path: str):
    """Return (mrna_meta, cds_by_tid).

    mrna_meta  : {tid: {chrom, strand, gene_id, evidence}}
    cds_by_tid : {tid: [(start_0, end_0), ...]}  (0-based half-open, ascending)
    """
    mrna_meta: dict[str, dict] = {}
    cds_by_tid: dict[str, list[tuple[int, int]]] = defaultdict(list)

    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, _src, feat, s1, e1, _sc, strand, _ph, attrs = cols
            # GFF3 1-based inclusive → 0-based half-open
            s0, e0 = int(s1) - 1, int(e1)

            attr_dict = {}
            for kv in attrs.split(";"):
                kv = kv.strip()
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attr_dict[k] = v

            if feat == "mRNA":
                tid = attr_dict.get("ID", "")
                mrna_meta[tid] = {
                    "chrom": chrom,
                    "strand": strand,
                    "gene_id": attr_dict.get("Parent", ""),
                    "evidence": attr_dict.get("Evidence", ""),
                }
            elif feat == "CDS":
                parent = attr_dict.get("Parent", "")
                cds_by_tid[parent].append((s0, e0))

    return mrna_meta, cds_by_tid


def build_cds_seq(intervals: list[tuple[int, int]], strand: str, chrom_seq: str) -> str:
    """Concatenate in ascending genomic order, then RC for minus strand (fixed)."""
    parts = [chrom_seq[s:e] for s, e in sorted(intervals)]
    nuc = "".join(parts)
    if strand == "-":
        nuc = reverse_complement(nuc)
    return nuc


def write_seq(fh, seq: str, width: int = 80) -> None:
    for i in range(0, len(seq), width):
        fh.write(seq[i : i + width] + "\n")


def is_helixer(evidence: str) -> bool:
    return "Helixer" in evidence


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff3",     required=True)
    ap.add_argument("--genome",   required=True)
    ap.add_argument("--old-prot", required=True)
    ap.add_argument("--out",      required=True)
    args = ap.parse_args()

    print("Loading genome…", file=sys.stderr)
    genome = load_genome(args.genome)

    print("Loading original proteins…", file=sys.stderr)
    old_prots = load_old_prot(args.old_prot)
    print(f"  {len(old_prots)} proteins in original prot.fa", file=sys.stderr)

    print("Parsing GFF3…", file=sys.stderr)
    mrna_meta, cds_by_tid = parse_gff3(args.gff3)
    print(f"  {len(mrna_meta)} mRNAs, {len(cds_by_tid)} with CDS", file=sys.stderr)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    stats = {
        "helixer_retranslated": 0,
        "non_helixer_kept": 0,
        "helixer_fixed_stops": 0,
        "no_old_protein": 0,
    }

    with open(args.out, "w") as fh:
        for tid, meta in mrna_meta.items():
            intervals = cds_by_tid.get(tid)
            if not intervals:
                continue

            evidence = meta["evidence"]
            chrom    = meta["chrom"]
            strand   = meta["strand"]
            gene_id  = meta["gene_id"]

            if not is_helixer(evidence):
                # Keep original protein unchanged
                old_prot = old_prots.get(tid)
                if old_prot is None:
                    stats["no_old_protein"] += 1
                    continue
                stats["non_helixer_kept"] += 1
                fh.write(f">{tid}\n")
                write_seq(fh, old_prot)
                continue

            # Helixer gene → re-translate with fixed code
            if chrom not in genome:
                continue
            chrom_seq = genome[chrom]
            cds_nuc   = build_cds_seq(intervals, strand, chrom_seq)
            protein   = translate(cds_nuc)
            if protein.endswith("*"):
                protein = protein[:-1]

            # Track improvement vs old
            old_prot = old_prots.get(tid, "")
            if old_prot and "*" in old_prot and "*" not in protein:
                stats["helixer_fixed_stops"] += 1

            stats["helixer_retranslated"] += 1
            aa_len = len(protein)
            fh.write(
                f">{tid} gene={gene_id} aa_len={aa_len} evidence={evidence} "
                f"strand={strand} cds_exons={len(intervals)}\n"
            )
            write_seq(fh, protein)

    print("Done.", file=sys.stderr)
    print(f"  Helixer re-translated : {stats['helixer_retranslated']}", file=sys.stderr)
    print(f"  Non-Helixer kept      : {stats['non_helixer_kept']}", file=sys.stderr)
    print(
        f"  Internal-stop fixes   : {stats['helixer_fixed_stops']} "
        "(Helixer genes where old had stop, new does not)",
        file=sys.stderr,
    )
    if stats["no_old_protein"]:
        print(f"  Non-Helixer missing   : {stats['no_old_protein']}", file=sys.stderr)


if __name__ == "__main__":
    main()
