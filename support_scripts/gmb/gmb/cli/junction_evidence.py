#!/usr/bin/env python3
"""``gmb-junction-evidence`` -- splice-junction diagnostics for any annotation.

Loads a STAR junction manifest once and reports, for every transcript of one or
more GFF3/GTF files, what the junctions say about each CDS (or exon) intron. Use
it to compare Helixer, transcript assemblies, GMB output and a reference with the
same evidence representation as ``gmb-build --star-junctions``::

    gmb-junction-evidence --star-junctions manifest.tsv --genome genome.fa \\
        --annotation helixer=helixer.gff3 --annotation gmb=consensus.gff3 \\
        --output-dir junction_report/

Reads junction files for the selected roles only (default ``construction``).
Writes ``<label>_introns.tsv``, ``<label>_models.tsv``, ``libraries.tsv`` and
``junction_evidence_summary.json``. Diagnostics only; nothing is selected.
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import sys
import time

from gmb.pipeline.junction_diagnostics import (
    JunctionAnnotator,
    read_models,
    write_diagnostics,
    write_library_table,
)
from gmb.pipeline.star_junctions import JunctionIndex, JunctionManifestError


def add_junction_args(parser: argparse.ArgumentParser) -> None:
    """The junction-input options shared by gmb-build, gmb-preflight and this tool."""
    g = parser.add_argument_group("STAR junction evidence (optional, diagnostics only)")
    g.add_argument(
        "--star-junctions",
        metavar="MANIFEST",
        help="TSV manifest of per-library STAR SJ.out.tab files (docs/star_junctions.md)",
    )
    g.add_argument(
        "--star-junction-roles",
        default="construction",
        help="comma-separated manifest roles to load (default: construction; "
        "'all' loads every library)",
    )


def roles_of(value: str) -> tuple[str, ...]:
    return tuple(r.strip() for r in (value or "").split(",") if r.strip())


def load_genome(path: str | None) -> dict:
    if not path:
        return {}
    from gmb.utils.fasta import load_genome as _load

    return _load(path)


def peak_rss_mb() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return round(r / (1 << 20) if sys.platform == "darwin" else r / 1024, 1)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="gmb-junction-evidence", description=__doc__.split("\n\n")[0])
    add_junction_args(p)
    p.add_argument("--genome", help="genome FASTA (enables motif and non-canonical checks)")
    p.add_argument(
        "--annotation",
        action="append",
        default=[],
        metavar="LABEL=PATH",
        help="annotation to report (GFF3 or GTF); repeatable",
    )
    p.add_argument(
        "--feature",
        default="CDS",
        choices=("CDS", "exon"),
        help="which segments define a model's introns (default CDS)",
    )
    p.add_argument("--output-dir", required=True)
    args = p.parse_args(argv)
    if not args.star_junctions:
        p.error("--star-junctions is required")
    t0 = time.perf_counter()
    genome = load_genome(args.genome)
    t_genome = time.perf_counter() - t0
    try:
        index = JunctionIndex.from_manifest(
            args.star_junctions,
            roles_of(args.star_junction_roles),
            {k: len(v) for k, v in genome.items()} if genome else None,
        )
    except JunctionManifestError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    os.makedirs(args.output_dir, exist_ok=True)
    write_library_table(os.path.join(args.output_dir, "libraries.tsv"), index)
    annot = JunctionAnnotator(index, genome)
    reports = {}
    for spec in args.annotation:
        label, _, path = spec.partition("=")
        if not path:
            label, path = os.path.splitext(os.path.basename(spec))[0], spec
        t = time.perf_counter()
        models = read_models(path, source=label, feature=args.feature)
        info = write_diagnostics(args.output_dir, label, annot, models)
        info.update(path=os.path.abspath(path), seconds=round(time.perf_counter() - t, 3))
        reports[label] = info
    summary = {
        "index": index.stats,
        "feature": args.feature,
        "genome_seconds": round(t_genome, 3),
        "annotations": reports,
        "peak_rss_mb": peak_rss_mb(),
        "notes": junction_caveats(index),
    }
    with open(os.path.join(args.output_dir, "junction_evidence_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
    print(
        f"{index.stats['libraries_loaded']} libraries, {index.stats['junctions']:,} junctions; "
        f"wrote {args.output_dir}"
    )
    return 0


def junction_caveats(index: JunctionIndex) -> list[str]:
    units = index.stats.get("count_units", [])
    return [
        "Diagnostics only: no model was selected or changed using these junctions.",
        f"Read counts are in library count units {units}; 'mate' units can count one fragment twice.",
        "STAR column 6 (annotated) is ignored: in two-pass mode it is not independent evidence.",
        "n_libraries counts replicate libraries of the same condition; use n_groups for breadth.",
        "Absence of a junction is reported as unknown or unobservable, never as contradiction.",
        "If these junctions came from the same alignments as the transcript assemblies, agreement "
        "between them is not independent validation.",
    ]


if __name__ == "__main__":
    sys.exit(main())
