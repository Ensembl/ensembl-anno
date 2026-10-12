"""gmb-run: the whole GMB workflow in one command.

    [long-read consensus] -> preflight -> build -> finalise

A thin wrapper around :func:`gmb.api.run_gene_model_builder`, so the command and
the Python API apply the same preparation policy. Non-interactive; suitable for
Slurm or workflow managers. Exit codes:

    0  success: finalised, QC-checked annotation in <output-dir>/finalise/
    1  unexpected error (see the traceback)
    2  invalid arguments
    3  long-read preparation failed (input rejected or consensus failed)
    4  preflight failed (inputs unusable)
    5  gene build failed
    6  sequence QC failed (annotation written to build/, not finalised)
    7  finalisation failed
"""

from __future__ import annotations

import argparse
import os
import sys
import time

EXIT_CODES = {"longread_consensus": 3, "preflight": 4, "build": 5, "sequence_qc": 6, "finalise": 7}
RECOMMENDED_FUNGAL_PRESET = "fungi_latest_settings"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="gmb-run",
        description="Run the full GMB workflow: automatic long-read consensus (when per-read "
        "alignments are given), preflight, build and finalise. For fungi use "
        f"--preset {RECOMMENDED_FUNGAL_PRESET} (RUN_FUNGI.md).",
        epilog="Exit codes: 0 ok, 2 bad arguments, 3 long-read preparation, 4 preflight, "
        "5 build, 6 sequence QC, 7 finalise, 1 unexpected error.",
    )
    s = p.add_argument_group("settings")
    s.add_argument(
        "--preset",
        required=True,
        help=f"configuration preset; for fungi: {RECOMMENDED_FUNGAL_PRESET}",
    )
    s.add_argument(
        "--config",
        action="append",
        default=[],
        help="extra YAML overlay (advanced; repeatable, last wins)",
    )
    i = p.add_argument_group("inputs (every file must use the genome's sequence names)")
    i.add_argument("--genome", required=True, help="genome FASTA")
    bb = i.add_mutually_exclusive_group(required=True)
    bb.add_argument("--helixer", help="Helixer GFF3 (ab initio backbone, with CDS)")
    bb.add_argument("--tiberius", help="Tiberius GTF (ab initio backbone)")
    bb.add_argument("--backbone", help="any other ab initio backbone (GTF/GFF3 with CDS)")
    i.add_argument("--backbone-label", help="source label for --backbone")
    i.add_argument("--scallop", help="short-read transcript assembly GTF (stranded)")
    i.add_argument("--stringtie", help="short-read transcript assembly GTF (stranded)")
    i.add_argument("--orthodb", help="protein-to-genome alignment GTF (e.g. OrthoDB genBlastG)")
    i.add_argument("--uniprot", help="protein-to-genome alignment GTF (e.g. UniProt genBlastG)")
    i.add_argument("--genblast", help="protein-to-genome alignment GTF")
    lr = i.add_mutually_exclusive_group()
    lr.add_argument(
        "--long-read-alignments",
        help="per-read long-read alignments (GTF, one transcript per read, strand "
        "set upstream); collapsed automatically with gmb-longread-consensus",
    )
    lr.add_argument(
        "--long-read-models",
        help="already-collapsed long-read transcript models (GTF); used as is",
    )
    i.add_argument(
        "--longread-consensus-preset",
        help="consensus preset (advanced; default follows --preset: fungi)",
    )
    i.add_argument("--longread-consensus-config", help="consensus YAML (advanced)")
    o = p.add_argument_group("output")
    o.add_argument(
        "--output-dir", required=True, help="run directory; the handover is <output-dir>/finalise/"
    )
    o.add_argument("--gene-prefix", default="GMB", help="prefix of new gene IDs (default GMB)")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    from gmb import __version__
    from gmb.api import run_gene_model_builder
    from gmb.pipeline.config import list_build_presets

    if args.preset not in [*list_build_presets(), "standard", "none"]:
        print(
            f"ERROR: unknown preset {args.preset!r}; available: {list_build_presets()}",
            file=sys.stderr,
        )
        return 2
    for label in (
        "genome",
        "helixer",
        "tiberius",
        "backbone",
        "scallop",
        "stringtie",
        "orthodb",
        "uniprot",
        "genblast",
        "long_read_alignments",
        "long_read_models",
        "longread_consensus_config",
    ):
        path = getattr(args, label)
        if path and not os.path.isfile(path):
            print(f"ERROR: --{label.replace('_', '-')} file not found: {path}", file=sys.stderr)
            return 2
    for path in args.config:
        if not os.path.isfile(path):
            print(f"ERROR: --config file not found: {path}", file=sys.stderr)
            return 2

    out = os.path.abspath(args.output_dir)
    os.makedirs(os.path.join(out, "logs"), exist_ok=True)
    log_fh = open(os.path.join(out, "logs", "gmb-run.log"), "a")

    def log(msg: str) -> None:
        line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
        print(line, flush=True)
        log_fh.write(line + "\n")
        log_fh.flush()

    log(f"gmb-run {__version__}: {' '.join(sys.argv[1:] if argv is None else argv)}")
    if args.preset == "fungi":
        log(
            f"NOTE: 'fungi' is the historical fungal baseline. The recommended fungal "
            f"settings are --preset {RECOMMENDED_FUNGAL_PRESET} (RUN_FUNGI.md)."
        )
    kind, backbone = next(
        (k, getattr(args, k)) for k in ("helixer", "tiberius", "backbone") if getattr(args, k)
    )
    try:
        result = run_gene_model_builder(
            genome=os.path.abspath(args.genome),
            backbone=os.path.abspath(backbone),
            backbone_kind=kind,
            backbone_label=args.backbone_label,
            short_read=[os.path.abspath(x) for x in (args.scallop, args.stringtie) if x],
            protein_alignment=[
                os.path.abspath(x) for x in (args.orthodb, args.uniprot, args.genblast) if x
            ],
            long_read=[os.path.abspath(args.long_read_models)] if args.long_read_models else None,
            long_read_alignments=(
                os.path.abspath(args.long_read_alignments) if args.long_read_alignments else None
            ),
            longread_consensus_preset=args.longread_consensus_preset,
            longread_consensus_config=args.longread_consensus_config,
            preset=args.preset,
            config=[os.path.abspath(c) for c in args.config],
            output_dir=out,
            gene_prefix=args.gene_prefix,
            log=log,
        )
    except ValueError as exc:
        log(f"ERROR: {exc}")
        return 2

    if result.ok:
        log("GMB run complete. Share these finalised files:")
        for label, path in (
            ("annotation (canonical tagged)", result.annotation),
            ("all isoforms", result.annotation_all_isoforms),
            ("proteins", result.proteins),
            ("CDS", result.cds),
            ("cDNA", result.cdna),
            ("handover manifest", result.handover_manifest),
        ):
            log(f"  {label:30s} {path}")
        log(f"  run summary                    {out}/run_summary.json")
        return 0
    stage = result.failed_stage or "build"
    log(
        f"GMB run FAILED at {stage} (exit {EXIT_CODES.get(stage, 1)}). "
        f"Logs: {out}/logs/. Summary: {out}/run_summary.json"
    )
    return EXIT_CODES.get(stage, 1)


if __name__ == "__main__":
    raise SystemExit(main())
