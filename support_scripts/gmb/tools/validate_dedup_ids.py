"""Validate the cross-chromosome gene_id/transcript_id namespacing fix.

Standalone (non-pytest) demonstration for a GTF file: shows the raw
duplicate-ID collision, then loads the file through
`gmb.compare.annotation_loader.load_annotation` and confirms:
  - original duplicated IDs are detected
  - internal IDs become unique
  - gene/transcript/exon/CDS counts are preserved
  - features remain linked to the correct gene/transcript

Usage:
    export PYTHONPATH=/path/to/support_scripts/gmb:$PYTHONPATH
    python3 tools/validate_dedup_ids.py tiberius_annotations/GCA_009914755.4_tiberius.gtf
"""

import re
import sys
from collections import defaultdict

from gmb.compare.annotation_loader import load_annotation

_GTF_ATTR = re.compile(r'(\w+)\s+"([^"]*)"')


def raw_duplicate_summary(path):
    """Count raw gene_id/transcript_id values that occur on more than one seqid."""
    gene_id_seqs = defaultdict(set)
    tx_id_seqs = defaultdict(set)
    feat_counts = defaultdict(int)

    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, feat, attrs_str = parts[0], parts[2], parts[8]
            feat_counts[feat] += 1
            attrs = dict(_GTF_ATTR.findall(attrs_str))
            if feat == "gene" and attrs.get("gene_id"):
                gene_id_seqs[attrs["gene_id"]].add(chrom)
            if feat == "transcript" and attrs.get("transcript_id"):
                tx_id_seqs[attrs["transcript_id"]].add(chrom)

    dup_genes = {g: s for g, s in gene_id_seqs.items() if len(s) > 1}
    dup_tx = {t: s for t, s in tx_id_seqs.items() if len(s) > 1}
    return feat_counts, gene_id_seqs, tx_id_seqs, dup_genes, dup_tx


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <gtf_path>")
        sys.exit(1)

    path = sys.argv[1]

    print("=" * 70)
    print(f"Raw file inspection: {path}")
    print("=" * 70)
    feat_counts, gene_id_seqs, tx_id_seqs, dup_genes, dup_tx = raw_duplicate_summary(path)
    print(f"gene features:       {feat_counts.get('gene', 0)}")
    print(f"transcript features: {feat_counts.get('transcript', 0) + feat_counts.get('mRNA', 0)}")
    print(f"exon features:       {feat_counts.get('exon', 0)}")
    print(f"CDS features:        {feat_counts.get('CDS', 0)}")
    print(f"unique raw gene_id values:       {len(gene_id_seqs)}")
    print(f"raw gene_id values on >1 chrom:  {len(dup_genes)}")
    print(f"unique raw transcript_id values: {len(tx_id_seqs)}")
    print(f"raw transcript_id values on >1 chrom: {len(dup_tx)}")

    print()
    print("=" * 70)
    print("Loading through gmb.compare.annotation_loader.load_annotation ...")
    print("=" * 70)
    genes, mrna, exons, cds = load_annotation(path, "Validation")

    print()
    print("=" * 70)
    print("Post-load checks")
    print("=" * 70)
    ok = True

    def check(label, cond):
        nonlocal ok
        status = "OK" if cond else "FAIL"
        if not cond:
            ok = False
        print(f"[{status}] {label}")

    # GFF3 inputs use "mRNA" rather than "transcript" for the transcript
    # feature line; count either so this check works for both formats.
    raw_tx_count = feat_counts.get("transcript", 0) + feat_counts.get("mRNA", 0)

    check(
        f"gene feature count preserved ({len(genes)} == {feat_counts.get('gene', 0)})",
        len(genes) == feat_counts.get("gene", 0),
    )
    check(
        f"transcript feature count preserved ({len(mrna)} == {raw_tx_count})",
        len(mrna) == raw_tx_count,
    )
    check(
        f"exon feature count preserved ({len(exons)} == {feat_counts.get('exon', 0)})",
        len(exons) == feat_counts.get("exon", 0),
    )
    check(
        f"CDS feature count preserved ({len(cds)} == {feat_counts.get('CDS', 0)})",
        len(cds) == feat_counts.get("CDS", 0),
    )

    if dup_genes:
        check(
            f"internal gene_id is fully unique ({genes['gene_id'].nunique()} == {len(genes)} rows, "
            f"vs {len(gene_id_seqs)} raw unique values)",
            genes["gene_id"].nunique() == len(genes),
        )
        check(
            f"internal transcript_id is fully unique ({mrna['transcript_id'].nunique()} == {len(mrna)} rows)",
            mrna["transcript_id"].nunique() == len(mrna),
        )
        check(
            "original_gene_id column preserves the raw (colliding) IDs",
            genes["original_gene_id"].nunique() == len(gene_id_seqs),
        )

        # Spot-check: every exon/CDS row's gene_id must match its own gene's
        # internal (namespaced) ID, not some other chromosome's gene sharing
        # the same original ID.
        gene_lookup = dict(zip(genes["gene_id"], genes["Chromosome"]))
        mismatches = 0
        for _, row in exons.iterrows():
            gid = row["gene_id"]
            if gid in gene_lookup and gene_lookup[gid] != row["Chromosome"]:
                mismatches += 1
        check(f"every exon's gene_id resolves to a gene on the same chromosome ({mismatches} mismatches)", mismatches == 0)
    else:
        print("[INFO] No cross-chromosome ID collisions detected in this file — "
              "namespacing was a no-op, as expected for genome-wide-unique IDs.")
        check(
            "gene_id left untouched (original_gene_id == gene_id)",
            (genes["gene_id"] == genes["original_gene_id"]).all() if not genes.empty else True,
        )

    print()
    print("RESULT:", "ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
