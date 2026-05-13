#!/usr/bin/env python3
"""Generate (or refresh) golden fixture files for the z_tritici region regression test.

Usage::

    # From support_scripts/gmb/
    python tests/generate_golden_fixtures.py

    # Optionally specify a pre-existing output dir (skip re-running pipeline)
    python tests/generate_golden_fixtures.py --from-dir /path/to/output

The script writes into ``tests/fixtures/expected/``:

    consensus_region1.gff3        -- copy of the consensus GFF3
    region1_gene_count.txt        -- integer gene count (one line)

These files are used by :class:`TestRegion1GoldenRegression` in
``test_z_tritici_subset.py``.

The pipeline is run on the pre-subsetted fixtures in
``tests/fixtures/z_tritici_region1/`` (first 500 kb of chr 1, ~8 seconds).
If those fixtures do not exist, run ``tests/create_fixtures.sh`` first.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

_HERE = Path(__file__).parent
_GMB_DIR = _HERE.parent
_REGION_FIXTURES = _HERE / "fixtures" / "z_tritici_region1"
_EXPECTED = _HERE / "fixtures" / "expected"


def _count_gff3_genes(gff3: Path) -> int:
    count = 0
    with open(gff3) as fh:
        for line in fh:
            parts = line.split("\t")
            if len(parts) >= 3 and parts[2] == "gene":
                count += 1
    return count


def run(from_dir: Path | None = None) -> None:
    if not _REGION_FIXTURES.exists():
        print(
            f"ERROR: Region fixtures not found at {_REGION_FIXTURES}.\n"
            "Run tests/create_fixtures.sh to generate them.",
            file=sys.stderr,
        )
        sys.exit(1)

    if from_dir:
        output_dir = from_dir
    else:
        tmp = tempfile.mkdtemp(prefix="golden_gen_")
        output_dir = Path(tmp)
        print(f"Running pipeline on pre-subsetted region fixtures into {output_dir} …")
        cmd = [
            sys.executable,
            "-m",
            "gmb.cli.build",
            "--scallop",
            str(_REGION_FIXTURES / "scallop_geneset.gtf"),
            "--stringtie",
            str(_REGION_FIXTURES / "stringtie_geneset.gtf"),
            "--helixer",
            str(_REGION_FIXTURES / "helixer_remapped.gff3"),
            "--orthodb",
            str(_REGION_FIXTURES / "orthodb_geneset.gtf"),
            "--uniprot",
            str(_REGION_FIXTURES / "uniprot_geneset.gtf"),
            "--genome",
            str(_REGION_FIXTURES / "genome.fa"),
            "--output-dir",
            str(output_dir),
            "--gene-prefix",
            "ZTGOLD",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(_GMB_DIR))
        if result.returncode != 0:
            print("PIPELINE FAILED:", result.stderr[-2000:], file=sys.stderr)
            sys.exit(1)

    gff3_src = output_dir / "consensus.gff3"
    if not gff3_src.exists():
        print(f"ERROR: {gff3_src} not found", file=sys.stderr)
        sys.exit(1)

    _EXPECTED.mkdir(parents=True, exist_ok=True)

    gff3_dst = _EXPECTED / "consensus_region1.gff3"
    shutil.copy(gff3_src, gff3_dst)
    print(f"Wrote {gff3_dst}")

    gene_count = _count_gff3_genes(gff3_src)
    count_dst = _EXPECTED / "region1_gene_count.txt"
    count_dst.write_text(str(gene_count) + "\n")
    print(f"Wrote {count_dst}  ({gene_count} genes)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--from-dir",
        type=Path,
        default=None,
        help="Use an existing pipeline output directory instead of re-running.",
    )
    args = parser.parse_args()
    run(from_dir=args.from_dir)
