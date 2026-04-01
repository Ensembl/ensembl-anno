#!/usr/bin/env python3
"""Generate (or refresh) golden fixture files for the z_tritici seqname-1 regression test.

Usage::

    # From support_scripts/gmb/
    python tests/generate_golden_fixtures.py

    # Optionally specify a pre-existing output dir (skip re-running pipeline)
    python tests/generate_golden_fixtures.py --from-dir /path/to/output

The script writes into ``tests/fixtures/expected/``:

    consensus_seqname1.gff3       -- copy of the consensus GFF3
    seqname1_gene_count.txt       -- integer gene count (one line)

These files are used by :class:`TestZTriticiGoldenRegression` in
``test_z_tritici_subset.py``.
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
_Z_TRITICI = _GMB_DIR / "z_tritici"
_FIXTURES = _HERE / "fixtures" / "expected"
_SCRIPT = _GMB_DIR / "gene_model_builder.py"


def _count_gff3_genes(gff3: Path) -> int:
    count = 0
    with open(gff3) as fh:
        for line in fh:
            parts = line.split("\t")
            if len(parts) >= 3 and parts[2] == "gene":
                count += 1
    return count


def run(from_dir: Path | None = None) -> None:
    if from_dir:
        output_dir = from_dir
    else:
        tmp = tempfile.mkdtemp(prefix="golden_gen_")
        output_dir = Path(tmp)
        print(f"Running pipeline into {output_dir} …")
        cmd = [
            sys.executable,
            str(_SCRIPT),
            "--scallop",
            str(_Z_TRITICI / "scallop_geneset.gtf"),
            "--stringtie",
            str(_Z_TRITICI / "stringtie_geneset.gtf"),
            "--helixer",
            str(_Z_TRITICI / "helixer_remapped.gff3"),
            "--orthodb",
            str(_Z_TRITICI / "orthodb_geneset.gtf"),
            "--uniprot",
            str(_Z_TRITICI / "uniprot_geneset.gtf"),
            "--genome",
            str(_Z_TRITICI / "zymoseptoria_tritici.fa"),
            "--assembly-report",
            str(_Z_TRITICI / "GCF_000219625.1_MYCGR_v2.0_assembly_report.txt"),
            "--output-dir",
            str(output_dir),
            "--gene-prefix",
            "ZTGOLD",
            "--seqname",
            "1",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print("PIPELINE FAILED:", result.stderr[-2000:], file=sys.stderr)
            sys.exit(1)

    gff3_src = output_dir / "consensus.gff3"
    if not gff3_src.exists():
        print(f"ERROR: {gff3_src} not found", file=sys.stderr)
        sys.exit(1)

    _FIXTURES.mkdir(parents=True, exist_ok=True)

    gff3_dst = _FIXTURES / "consensus_seqname1.gff3"
    shutil.copy(gff3_src, gff3_dst)
    print(f"Wrote {gff3_dst}")

    gene_count = _count_gff3_genes(gff3_src)
    count_dst = _FIXTURES / "seqname1_gene_count.txt"
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
