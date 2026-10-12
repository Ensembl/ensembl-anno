#!/usr/bin/env python3
"""Gene Model Builder — Python integration example.

Shows how an orchestration pipeline calls GMB through its stable API,
``run_gene_model_builder`` (the same function `gmb-run` uses). It runs automatic
long-read consensus (when per-read alignments are given), preflight, build and
finalise, and writes <output-dir>/run_summary.json.

    python examples/run_gene_model_builder.py \\
        --genome /data/genome.fa --helixer /data/helixer.gff3 \\
        --short-read /data/scallop.gtf /data/stringtie.gtf \\
        --protein /data/orthodb.gtf /data/genblast.gtf \\
        --long-read-alignments /data/minimap2_annotation.gtf \\
        --output-dir /work/gmb --gene-prefix XXGMB
"""

from __future__ import annotations

import argparse
import sys

from gmb import run_gene_model_builder


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--genome", required=True)
    p.add_argument("--helixer", required=True, help="Helixer GFF3 (genome sequence names)")
    p.add_argument("--short-read", nargs="*", default=[], help="short-read assembly GTFs (max 2)")
    p.add_argument("--protein", nargs="*", default=[], help="protein alignment GTFs (max 3)")
    p.add_argument("--long-read-alignments", help="per-read long-read GTF (collapsed by GMB)")
    p.add_argument("--preset", default="fungi_latest_settings",
                   help="default: the recommended fungal settings")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--gene-prefix", default="GMB")
    args = p.parse_args(argv)

    result = run_gene_model_builder(
        genome=args.genome,
        backbone=args.helixer,
        backbone_kind="helixer",
        short_read=args.short_read,
        protein_alignment=args.protein,
        long_read_alignments=args.long_read_alignments,
        preset=args.preset,
        output_dir=args.output_dir,
        gene_prefix=args.gene_prefix,
    )
    if not result.ok:
        print(f"GMB failed at {result.failed_stage}: {result.failure_reason}", file=sys.stderr)
        print(f"See {result.output_dir}/run_summary.json and {result.output_dir}/logs/",
              file=sys.stderr)
        return 1
    print("Handover:", result.handover_dir)
    print("  annotation      ", result.annotation)
    print("  all isoforms    ", result.annotation_all_isoforms)
    print("  proteins        ", result.proteins)
    print("  long-read models", result.long_read_models)
    return 0


if __name__ == "__main__":
    sys.exit(main())
