#!/usr/bin/env python3
"""Gene Model Builder — Python integration example.

Shows how an orchestration pipeline calls GMB through its stable API. Nothing
here imports an internal module: ``run_gene_model_builder`` is the whole
interface.

    python examples/run_gene_model_builder.py \
        --genome /data/genome.fa \
        --backbone /data/backbone.gff3 --backbone-kind helixer \
        --short-read /data/scallop.gtf /data/stringtie.gtf \
        --protein /data/orthodb.gtf \
        --preset fungi --output-dir /work/gmb --gene-prefix XXGMB
"""

from __future__ import annotations

import argparse
import json
import sys

from gmb import run_gene_model_builder


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--genome", required=True)
    p.add_argument("--backbone", required=True)
    p.add_argument(
        "--backbone-kind", default="helixer", choices=["helixer", "tiberius", "backbone"]
    )
    p.add_argument(
        "--short-read", nargs="*", default=[], help="assembled short-read transcript GTFs"
    )
    p.add_argument(
        "--long-read", nargs="*", default=[], help="collapsed long-read transcript GTFs (optional)"
    )
    p.add_argument("--protein", nargs="*", default=[], help="protein-to-genome alignment GTFs")
    p.add_argument("--preset", default="standard", choices=["standard", "apicomplexa", "fungi"])
    p.add_argument("--config", nargs="*", default=[], help="YAML overlays, last wins")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--gene-prefix", default="GMB")
    args = p.parse_args(argv)

    result = run_gene_model_builder(
        genome=args.genome,
        backbone=args.backbone,
        backbone_kind=args.backbone_kind,
        short_read=args.short_read,
        long_read=args.long_read or None,
        protein_alignment=args.protein,
        preset=args.preset,
        config=args.config or None,
        output_dir=args.output_dir,
        gene_prefix=args.gene_prefix,
    )

    print(f"stages run        : {', '.join(result.stages_run)}")
    print(f"preflight verdict : {result.preflight_verdict}")
    print(f"QC passed         : {result.qc_passed}")
    print(f"overall ok        : {result.ok}")

    # Preflight FAIL stops the run before the expensive build.
    if result.preflight_verdict == "FAIL":
        print(
            "\nPreflight FAILED — the evidence bundle is unfit to annotate from.", file=sys.stderr
        )
        for check in (result.preflight_report or {}).get("checks", []):
            if check.get("verdict") == "FAIL":
                target = f" [{check['target']}]" if check.get("target") else ""
                print(f"  - {check['name']}{target}: {check['message']}", file=sys.stderr)
        return 1

    # A QC failure means the annotation exists but must not be handed over.
    if result.qc_passed is False:
        print("\nQC FAILED — do not hand over this build.", file=sys.stderr)
        print(
            json.dumps((result.qc_report or {}).get("failed_checks", []), indent=1),
            file=sys.stderr,
        )
        print(f"Provenance for the failing run: {result.run_manifest}", file=sys.stderr)
        return 1

    if not result.ok:
        print("\nA stage failed; see the logs under " f"{args.output_dir}/logs/", file=sys.stderr)
        return 1

    print("\nHandover files")
    for label, path in (
        ("annotation (canonical)", result.annotation),
        ("annotation (all isoforms)", result.annotation_all_isoforms),
        ("proteins", result.proteins),
        ("CDS", result.cds),
        ("cDNA", result.cdna),
        ("evidence attribution", result.evidence_attribution),
        ("resolved config", result.resolved_config),
        ("run manifest", result.run_manifest),
        ("handover manifest", result.handover_manifest),
    ):
        if path:
            print(f"  {label:26s} {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
