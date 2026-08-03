#!/usr/bin/env python3
"""Standalone long-read consensus tool.

Converts a raw Minimap2 per-read alignment GTF (one record per aligned read)
into a small, provenance-preserving consensus transcript track suitable for:

  * inspection as a standalone consensus track;
  * use as ``--minimap2`` evidence in ``gmb-build``;
  * downstream use by any tool that consumes transcript GTF.

This tool does NOT require the GMB builder, protein validation, canonical
selection, or InterPro.

Usage — pure long-read mode (no short-read support required):

    gmb-longread-consensus \\
        --input raw_minimap2.gtf \\
        --output-dir out/ \\
        --preset pfalciparum_pure

Usage — short-read-assisted mode (structure-aware rescue):

    gmb-longread-consensus \\
        --input raw_minimap2.gtf \\
        --output-dir out/ \\
        --preset pfalciparum_assisted \\
        --scallop tissue1.scallop.gtf \\
        --scallop tissue2.scallop.gtf \\
        --stringtie tissue1.stringtie.gtf

Split-only (reuse the split later):

    gmb-longread-consensus \\
        --input raw_minimap2.gtf \\
        --output-dir out/ \\
        --split-only

Reuse a previous split:

    gmb-longread-consensus \\
        --reuse-split-dir out/_split_by_seqname \\
        --output-dir out2/ \\
        --preset pfalciparum_pure

Process one chromosome (Slurm array / local testing):

    gmb-longread-consensus \\
        --input-split-gtf split/1.gtf \\
        --seqname 1 \\
        --output-dir chr1_out/ \\
        --preset pfalciparum_pure
"""

from __future__ import annotations

import argparse
import os
import sys
import time

from gmb.pipeline.longread import (
    build_consensus_for_seqname,
    load_longread_consensus_config,
    split_by_seqname,
    validate_config,
    validate_split_dir,
    write_consensus_gtf,
    write_dropped_records_tsv,
    write_rescue_attribution_tsv,
    write_run_manifest,
    write_summary,
)
from gmb.utils.io import ensure_dir
from gmb.utils.logging import resolve_log_file, setup_logging


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="gmb-longread-consensus",
        description=(
            "Collapse a raw Minimap2 long-read alignment GTF into a small "
            "consensus transcript track.  Never loads the full raw file into "
            "memory — always splits by seqname first.\n\n"
            "See docs/longread_consensus.md for the full user guide."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # --- Input ---
    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument(
        "--input",
        metavar="RAW_GTF",
        help="Raw Minimap2 alignment GTF (full genome or subset).",
    )
    input_group.add_argument(
        "--reuse-split-dir",
        metavar="DIR",
        help=(
            "Reuse an existing split directory (must contain split_manifest.tsv). "
            "Skips Stage 1; validates checksums before Stage 2."
        ),
    )
    input_group.add_argument(
        "--input-split-gtf",
        metavar="GTF",
        help=(
            "Process exactly one pre-split seqname GTF (Stage 2 only). "
            "Must be paired with --seqname."
        ),
    )

    # --- Stage control ---
    parser.add_argument(
        "--split-only",
        action="store_true",
        help="Run Stage 1 only (split by seqname) and exit. Writes split_manifest.tsv.",
    )
    parser.add_argument(
        "--split-dir",
        default=None,
        metavar="DIR",
        help=(
            "Directory for split files (default: <output-dir>/_split_by_seqname). "
            "Only used with --input."
        ),
    )
    parser.add_argument(
        "--seqname",
        default=None,
        help=(
            "Restrict to a single seqname.  "
            "Required with --input-split-gtf; optional with --input to subset."
        ),
    )

    # --- Output ---
    parser.add_argument(
        "--output-dir",
        required=True,
        metavar="DIR",
        help="Output directory.  Created if it does not exist.",
    )

    # --- Config / preset ---
    parser.add_argument(
        "--config",
        default=None,
        metavar="YAML",
        help="YAML config file with threshold overrides.",
    )
    parser.add_argument(
        "--preset",
        default=None,
        metavar="NAME",
        help=(
            "Named preset from configs/longread_consensus/ "
            "(e.g. pfalciparum_pure, pfalciparum_assisted).  "
            "Merged before --config when both are given."
        ),
    )

    # --- Short-read rescue inputs (repeatable) ---
    parser.add_argument(
        "--scallop",
        dest="scallop",
        action="append",
        default=[],
        metavar="GTF",
        help="Scallop GTF for short-read rescue.  May be repeated for multiple samples.",
    )
    parser.add_argument(
        "--stringtie",
        dest="stringtie",
        action="append",
        default=[],
        metavar="GTF",
        help="StringTie GTF for short-read rescue.  May be repeated for multiple samples.",
    )

    # --- Logging ---
    parser.add_argument("--log-file", default=None, help="Log path.")
    parser.add_argument("--no-log-file", action="store_true")

    return parser.parse_args(argv)


def main(argv=None) -> None:
    start_time = time.time()
    args = parse_args(argv)

    ensure_dir(args.output_dir)
    log_file = resolve_log_file(args.output_dir, args.log_file, args.no_log_file)
    setup_logging(log_file=log_file, capture_stdio=log_file is not None)
    if log_file:
        print(f"Logging to {log_file}")

    # --- Validate input mode ---
    if args.input_split_gtf and not args.seqname:
        sys.exit("ERROR: --input-split-gtf requires --seqname.")
    if args.split_only and not args.input:
        sys.exit("ERROR: --split-only requires --input.")

    # --- Load and validate config ---
    try:
        cfg = load_longread_consensus_config(path=args.config, preset=args.preset)
    except (ValueError, FileNotFoundError) as exc:
        sys.exit(f"ERROR: Config error: {exc}")

    scallop_paths = [p for p in args.scallop if p]
    stringtie_paths = [p for p in args.stringtie if p]

    try:
        validate_config(cfg, scallop_paths, stringtie_paths)
    except ValueError as exc:
        sys.exit(f"ERROR: Config validation failed: {exc}")

    print(f"Config: rescue={'enabled' if cfg.shortread_rescue.enabled else 'disabled'}")
    if cfg.shortread_rescue.enabled:
        me = cfg.shortread_rescue.multi_exon
        se = cfg.shortread_rescue.single_exon
        if me.enabled:
            print(
                f"  multi-exon rescue: policy={me.policy}, junction_tol={me.junction_tolerance_bp}bp"
            )
        if se.enabled:
            print(
                f"  single-exon rescue: policy={se.policy}, "
                f"min_overlap={se.reciprocal_overlap_min}, "
                f"boundary_tol={se.boundary_tolerance_bp}bp"
            )
        if scallop_paths:
            print(f"  Scallop: {len(scallop_paths)} file(s)")
        if stringtie_paths:
            print(f"  StringTie: {len(stringtie_paths)} file(s)")

    # --- Stage 1: split (or reuse) ---
    seqname_paths: dict[str, str] = {}

    if args.input_split_gtf:
        # Single pre-split file mode (Slurm array / local test)
        if not os.path.exists(args.input_split_gtf):
            sys.exit(f"ERROR: --input-split-gtf not found: {args.input_split_gtf}")
        seqname_paths[args.seqname] = args.input_split_gtf
        print(f"Using pre-split file: {args.input_split_gtf} (seqname={args.seqname})")

    elif args.reuse_split_dir:
        # Reuse validated split directory
        print(f"Validating split directory: {args.reuse_split_dir}")
        try:
            seqname_paths = validate_split_dir(args.reuse_split_dir)
        except (FileNotFoundError, ValueError) as exc:
            sys.exit(f"ERROR: {exc}")
        if args.seqname:
            if args.seqname not in seqname_paths:
                sys.exit(
                    f"ERROR: seqname '{args.seqname}' not found in split dir "
                    f"(available: {sorted(seqname_paths)})"
                )
            seqname_paths = {args.seqname: seqname_paths[args.seqname]}
        print(f"Reusing {len(seqname_paths)} validated split file(s).")

    else:
        # Normal: split from raw input
        if not args.input:
            sys.exit("ERROR: one of --input, --reuse-split-dir, or --input-split-gtf is required.")
        if not os.path.exists(args.input):
            sys.exit(f"ERROR: --input file not found: {args.input}")

        split_dir = args.split_dir or os.path.join(args.output_dir, "_split_by_seqname")
        print(f"Input: {args.input} ({os.path.getsize(args.input):,} bytes)")
        print("Stage 1: splitting raw file by seqname…")
        seqname_paths = split_by_seqname(
            args.input,
            split_dir,
            only_seqname=args.seqname,
            write_manifest=True,
        )
        if not seqname_paths:
            sys.exit("ERROR: no records found (check --seqname matches the input file).")

        if args.split_only:
            print(f"Split-only mode: {len(seqname_paths)} seqname file(s) written to {split_dir}")
            print("Done.")
            return

    # --- Stage 2: consensus building ---
    print(f"Stage 2: building consensus for {len(seqname_paths)} seqname(s)…")
    all_records: list[dict] = []
    all_stats: list[dict] = []
    all_dropped: list[dict] = []

    for seqname in sorted(seqname_paths):
        t0 = time.time()
        records, stats, dropped = build_consensus_for_seqname(
            seqname,
            seqname_paths[seqname],
            cfg,
            scallop_paths=scallop_paths,
            stringtie_paths=stringtie_paths,
        )
        stats["elapsed_s"] = round(time.time() - t0, 1)
        all_records.extend(records)
        all_stats.append(stats)
        all_dropped.extend(dropped)
        print(
            f"  {seqname}: {stats['input_reads']:,} reads → "
            f"{stats['consensus_models_output']:,} consensus models "
            f"({stats['elapsed_s']}s)"
        )

    # --- Write outputs ---
    output_gtf = os.path.join(args.output_dir, "minimap2_consensus.gtf")
    write_consensus_gtf(all_records, output_gtf)
    print(f"Wrote {len(all_records):,} consensus transcripts → {output_gtf}")

    attribution_path = os.path.join(args.output_dir, "longread_consensus_rescue_attribution.tsv")
    n_rescued_attr = write_rescue_attribution_tsv(all_records, attribution_path)
    if n_rescued_attr:
        print(f"Rescue attribution: {n_rescued_attr} rows → {attribution_path}")
    else:
        write_rescue_attribution_tsv(all_records, attribution_path)

    dropped_path = os.path.join(args.output_dir, "longread_consensus_dropped.tsv")
    n_dropped = write_dropped_records_tsv(all_dropped, dropped_path)
    print(f"Dropped record log: {n_dropped} rows → {dropped_path}")

    args_dict = {
        "input": args.input,
        "reuse_split_dir": args.reuse_split_dir,
        "input_split_gtf": args.input_split_gtf,
        "output_dir": args.output_dir,
        "seqname": args.seqname,
        "config": args.config,
        "preset": args.preset,
        "scallop": scallop_paths,
        "stringtie": stringtie_paths,
    }
    summary = write_summary(all_stats, all_records, args_dict, cfg, args.output_dir, start_time)

    manifest_path = write_run_manifest(
        args.output_dir,
        input_path=args.input,
        scallop_paths=scallop_paths,
        stringtie_paths=stringtie_paths,
        output_files={
            "minimap2_consensus_gtf": output_gtf,
            "rescue_attribution_tsv": attribution_path,
            "dropped_records_tsv": dropped_path,
            "summary_json": os.path.join(args.output_dir, "longread_consensus_summary.json"),
            "summary_tsv": os.path.join(args.output_dir, "longread_consensus_summary.tsv"),
        },
        summary=summary,
        cfg=cfg,
        preset=args.preset,
        start_time=start_time,
    )
    print(f"Run manifest → {manifest_path}")
    print("Done.")


if __name__ == "__main__":
    main()
