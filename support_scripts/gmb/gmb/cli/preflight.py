#!/usr/bin/env python3
"""``gmb-preflight`` — validate an evidence bundle before building.

Takes the same input arguments as ``gmb-build`` so a production pipeline can run
it with the identical argument list, then run the build only if it passes::

    gmb-preflight --preset fungi --genome "$GENOME" --helixer "$BACKBONE" \\
        --scallop "$SCALLOP" --stringtie "$STRINGTIE" --orthodb "$ORTHODB" \\
        --output-dir "$OUT/preflight" \\
      && gmb-build --preset fungi --genome "$GENOME" ... --output-dir "$OUT/build"

Exit codes
----------
``0``  PASS or WARN — safe to build.
``1``  FAIL — at least one check says the build should not proceed as configured.
       Override with ``--allow-fail`` when you know better; the verdict is still
       written to the report either way.
``2``  preflight could not run at all (unreadable genome, bad arguments).
"""

from __future__ import annotations

import argparse
import json
import os
import sys

from gmb.preflight import FAIL, run_preflight
from gmb.pipeline.config import load_config, validate_selection_policy
from gmb.utils.logging import setup_logging


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="gmb-preflight",
        description=(
            "Validate a GMB evidence bundle before an expensive build. Reports "
            "seqid compatibility, ID collisions, strand completeness, structural "
            "distributions, per-track splice quality, resolved evidence roles and "
            "weights, and which optional selection policies will actually fire. "
            "Never reads a reference annotation."
        ),
    )
    setup = p.add_argument_group("Setup")
    setup.add_argument("--config", action="append", default=None,
                       help="YAML config override; repeatable, last wins.")
    setup.add_argument("--preset", default="standard",
                       help="Clade preset (default: standard, the neutral base).")

    inputs = p.add_argument_group("Inputs (same names as gmb-build)")
    inputs.add_argument("--genome", required=True, help="Genome FASTA")
    inputs.add_argument("--scallop", help="Short-read transcript assembly GTF")
    inputs.add_argument("--stringtie", help="Short-read transcript assembly GTF")
    inputs.add_argument("--minimap2", help="Long-read transcript models GTF (optional)")
    inputs.add_argument("--helixer", help="Helixer GFF3 (ab initio backbone)")
    inputs.add_argument("--tiberius", help="Tiberius GTF (ab initio backbone)")
    inputs.add_argument("--orthodb", help="OrthoDB protein-to-genome GTF")
    inputs.add_argument("--uniprot", help="UniProt protein-to-genome GTF")
    inputs.add_argument("--genblast", help="GenBlast protein-to-genome GTF")

    out = p.add_argument_group("Output")
    out.add_argument("--output-dir", help="Directory for preflight_report.{json,txt}")
    out.add_argument("--json", action="store_true",
                     help="Print the report as JSON on stdout instead of text.")
    out.add_argument("--allow-fail", action="store_true",
                     help="Exit 0 even when a check FAILs (verdict still reported).")
    out.add_argument("--log-file", help="Log path (default: <output-dir>/gmb-preflight.log)")
    out.add_argument("--no-log-file", action="store_true", help="Disable file logging.")
    return p


def _tracks_from_args(args) -> list:
    """Map CLI arguments onto (label, path) pairs.

    The label is what the role resolver sees, so it must match the labels the
    builder uses for the same flags -- this is the single place the two are
    kept in step.
    """
    pairs = [
        ("Scallop", args.scallop),
        ("StringTie", args.stringtie),
        ("Minimap2", args.minimap2),
        ("Helixer", args.helixer),
        ("Tiberius", args.tiberius),
        ("OrthoDB", args.orthodb),
        ("UniProt", args.uniprot),
        ("GenBlast", args.genblast),
    ]
    return [{"label": label, "path": path} for label, path in pairs if path]


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    if args.helixer and args.tiberius:
        print("error: --helixer and --tiberius are mutually exclusive "
              "(exactly one backbone).", file=sys.stderr)
        return 2

    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
    log_file = None if args.no_log_file else (
        args.log_file or (os.path.join(args.output_dir, "gmb-preflight.log")
                          if args.output_dir else None))
    if log_file:
        setup_logging(log_file=log_file)

    config = load_config(args.config, args.preset)

    # The backbone flag decides which label carries the backbone role, exactly
    # as in gmb-build.
    if args.tiberius:
        config.scoring.backbone_label = "Tiberius"
    elif args.helixer:
        config.scoring.backbone_label = "Helixer"

    for warning in validate_selection_policy(config):
        print(f"config warning: {warning}", file=sys.stderr)

    report = run_preflight(
        {"genome": args.genome, "tracks": _tracks_from_args(args)}, config)

    if args.json:
        print(json.dumps(report.to_dict(), indent=1))
    else:
        print(report.to_text())

    if args.output_dir:
        with open(os.path.join(args.output_dir, "preflight_report.json"), "w") as fh:
            json.dump(report.to_dict(), fh, indent=1)
        with open(os.path.join(args.output_dir, "preflight_report.txt"), "w") as fh:
            fh.write(report.to_text() + "\n")

    if report.verdict == FAIL and not args.allow_fail:
        print("\npreflight FAILED — fix the inputs or the configuration, or "
              "re-run with --allow-fail to proceed anyway.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
