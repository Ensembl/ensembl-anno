#!/usr/bin/env python3
"""Subset integration tests using bundled Z. tritici example data.

Two test modes
--------------

**Fast mode (default, CI-ready, ~8 seconds):**
    Uses pre-subsetted fixture files in ``tests/fixtures/z_tritici_region1/``
    covering the first 500 kb of chromosome 1.  These were created by the
    ``create_fixtures.sh`` script (or reproduced with the awk/python commands
    in that script).  The pipeline runs without ``--seqname`` or ``--region``
    flags because the fixtures are already chr-1-only.

    All tests in ``TestZTriticiRegion1`` and ``TestRegion1GoldenRegression``
    use this mode.  They are marked ``@pytest.mark.integration``.

**Slow mode (local only, ~25 min):**
    Runs the full pipeline on the entire bundled z_tritici data using
    ``--seqname 1``.  This exercises the seqname-mapping and subsetting
    code paths with the large OrthoDB file (6.7 M lines).  Only runs when
    the environment variable ``RUN_SLOW_INTEGRATION=1`` is set.

Why pre-subsetted fixtures?
---------------------------
The bundled ``z_tritici/orthodb_geneset.gtf`` has 6.7 M lines.  Even when
the pipeline is given ``--seqname 1``, it must load the entire file before
filtering — because pyranges reads the file into memory first.  Chromosome 1
alone has 1.4 M lines, making the seqname-1 run take ~25 minutes on a laptop.

The 500 kb region subset contains ~47 K OrthoDB lines, ~700 Scallop/StringTie
lines, and ~1500 Helixer lines.  It produces 136 genes and finishes in <10 s.

Golden-file regression
----------------------
``tests/fixtures/expected/consensus_region1.gff3`` and
``tests/fixtures/expected/region1_gene_count.txt`` store expected outputs.
The structural comparison is sorted (gene loci by seqname/start/end/strand)
so minor ordering changes do not cause spurious failures.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_HERE = Path(__file__).parent
_GMB_DIR = _HERE.parent
_Z_TRITICI = _GMB_DIR / "z_tritici"
_FIXTURES = _HERE / "fixtures" / "expected"
_REGION_FIXTURES = _HERE / "fixtures" / "z_tritici_region1"

# Check for pyranges — required by the pipeline
try:
    import pyranges  # noqa: F401
except ImportError:
    pytest.skip("pyranges required", allow_module_level=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _run_pipeline_region(output_dir: Path) -> subprocess.CompletedProcess:
    """Run gene_model_builder on the pre-subsetted 500 kb region fixture."""
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
        "ZTTEST",
    ]
    return subprocess.run(cmd, capture_output=True, text=True, cwd=str(_GMB_DIR))


def _run_pipeline_seqname1(output_dir: Path) -> subprocess.CompletedProcess:
    """Run gene_model_builder on the full z_tritici dataset with --seqname 1.

    This exercises the seqname-mapping and subsetting code paths with the
    large OrthoDB file.  Takes ~25 minutes on a laptop — only for local runs.
    """
    cmd = [
        sys.executable,
        "-m",
        "gmb.cli.build",
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
        "ZTTEST",
        "--seqname",
        "1",
    ]
    return subprocess.run(cmd, capture_output=True, text=True, cwd=str(_GMB_DIR))


def _parse_gff3_genes(gff3_path: Path) -> list[dict]:
    """Return list of gene feature dicts from a GFF3 file."""
    genes = []
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] == "gene":
                genes.append(
                    {
                        "seqname": parts[0],
                        "start": int(parts[3]),
                        "end": int(parts[4]),
                        "strand": parts[6],
                        "attrs": parts[8],
                    }
                )
    return genes


def _parse_gff3_mrnas(gff3_path: Path) -> list[str]:
    """Return list of mRNA IDs from a GFF3 file."""
    mrnas = []
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "mRNA":
                continue
            for attr in parts[8].split(";"):
                if attr.startswith("ID="):
                    mrnas.append(attr[3:])
    return mrnas


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def region1_output(tmp_path_factory):
    """Run the pipeline once on the pre-subsetted 500 kb region; share output."""
    if not _REGION_FIXTURES.exists():
        pytest.skip(
            "Pre-subsetted region fixtures not found. "
            "Run tests/create_fixtures.sh to generate them."
        )
    output_dir = tmp_path_factory.mktemp("ztritici_region1")
    result = _run_pipeline_region(output_dir)
    return output_dir, result


# ---------------------------------------------------------------------------
# Fast integration tests (CI-ready, ~8 seconds)
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestZTriticiRegion1:
    """Smoke tests using the pre-subsetted 500 kb z_tritici region fixture."""

    def test_pipeline_exits_zero(self, region1_output):
        """Pipeline should complete without errors."""
        _, result = region1_output
        assert result.returncode == 0, (
            f"Pipeline failed (exit {result.returncode}):\n"
            f"STDOUT:\n{result.stdout[-3000:]}\n"
            f"STDERR:\n{result.stderr[-3000:]}"
        )

    def test_required_output_files_exist(self, region1_output):
        """All mandatory output files must be present."""
        output_dir, _ = region1_output
        for fname in ("consensus.gff3", "cdna.fa", "prot.fa", "summary.json", "summary.tsv"):
            assert (output_dir / fname).exists(), f"Missing required output: {fname}"

    def test_nonzero_gene_count(self, region1_output):
        """At least one gene must be called in the region."""
        output_dir, _ = region1_output
        gff3 = output_dir / "consensus.gff3"
        genes = _parse_gff3_genes(gff3)
        assert len(genes) > 0, "No genes in consensus.gff3 — pipeline produced empty output"

    def test_no_protein_evidence_as_mrna(self, region1_output):
        """Protein evidence (OrthoDB_/UniProt_) must never appear as mRNA features.

        Protein alignments are support-only evidence — they must not be
        promoted to mRNA features in the output GFF3.
        """
        output_dir, _ = region1_output
        gff3 = output_dir / "consensus.gff3"
        mrna_ids = _parse_gff3_mrnas(gff3)
        bad = [m for m in mrna_ids if m.startswith(("OrthoDB_", "UniProt_"))]
        assert bad == [], (
            f"Protein evidence IDs appeared as mRNA in consensus.gff3: {bad[:5]}"
        )

    def test_summary_json_structure(self, region1_output):
        """summary.json must have the expected top-level structure."""
        output_dir, _ = region1_output
        with open(output_dir / "summary.json") as fh:
            report = json.load(fh)

        assert "summary" in report, "summary.json missing 'summary' key"
        assert "total_loci" in report["summary"], "'summary.total_loci' key missing"
        assert report["summary"]["total_loci"] >= 1, "total_loci is 0"
        assert "filtering" in report, "summary.json missing 'filtering' key"

    def test_evidence_attribution_written(self, region1_output):
        """evidence_attribution.tsv must be present and non-empty."""
        output_dir, _ = region1_output
        tsv = output_dir / "evidence_attribution.tsv"
        if not tsv.exists():
            pytest.skip("evidence_attribution.tsv not written by this config")
        with open(tsv) as fh:
            lines = [line for line in fh if not line.startswith("#")]
        assert len(lines) > 1, "evidence_attribution.tsv has no data rows"

    def test_cdna_fasta_has_sequences(self, region1_output):
        """cdna.fa should contain at least one transcript sequence."""
        output_dir, _ = region1_output
        content = (output_dir / "cdna.fa").read_text()
        assert content.count(">") >= 1, "cdna.fa has no FASTA headers"
        seqs = [line for line in content.splitlines() if line and not line.startswith(">")]
        assert len("".join(seqs)) > 0, "cdna.fa has headers but no sequence data"

    def test_protein_fasta_has_sequences(self, region1_output):
        """prot.fa should contain at least one protein sequence."""
        output_dir, _ = region1_output
        content = (output_dir / "prot.fa").read_text()
        assert content.count(">") >= 1, "prot.fa has no FASTA headers"


# ---------------------------------------------------------------------------
# Golden regression tests (fast, uses pre-subsetted fixtures + stored expected)
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestRegion1GoldenRegression:
    """Golden-file regression: compare output against stored expected files.

    Tests are skipped if the fixture files do not exist.  To (re-)generate
    them, run::

        python tests/generate_golden_fixtures.py

    or copy outputs from a known-good run into tests/fixtures/expected/.
    """

    def test_gene_count_matches_golden(self, region1_output):
        """Gene count should match the stored golden count exactly."""
        golden = _FIXTURES / "region1_gene_count.txt"
        if not golden.exists():
            pytest.skip("Golden gene count not found — run generate_golden_fixtures.py")

        output_dir, _ = region1_output
        genes = _parse_gff3_genes(output_dir / "consensus.gff3")
        expected = int(golden.read_text().strip())
        assert len(genes) == expected, (
            f"Gene count changed: got {len(genes)}, expected {expected}. "
            "If intentional, re-run generate_golden_fixtures.py."
        )

    def test_gff3_loci_match_golden(self, region1_output):
        """Canonical (sorted) gene loci must match the stored golden GFF3."""
        golden = _FIXTURES / "consensus_region1.gff3"
        if not golden.exists():
            pytest.skip("Golden GFF3 not found — run generate_golden_fixtures.py")

        output_dir, _ = region1_output
        actual_genes = sorted(
            (g["seqname"], g["start"], g["end"], g["strand"])
            for g in _parse_gff3_genes(output_dir / "consensus.gff3")
        )
        expected_genes = sorted(
            (g["seqname"], g["start"], g["end"], g["strand"])
            for g in _parse_gff3_genes(golden)
        )
        assert actual_genes == expected_genes, (
            f"Gene loci changed vs golden. "
            f"Added: {set(actual_genes) - set(expected_genes)}, "
            f"Removed: {set(expected_genes) - set(actual_genes)}"
        )


# ---------------------------------------------------------------------------
# Slow seqname-1 test (local only, ~25 min, gated by env var)
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestZTriticiSeqname1Slow:
    """Full seqname-1 pipeline test — exercises seqname-mapping + subsetting code path.

    Gated by ``RUN_SLOW_INTEGRATION=1`` because the OrthoDB file has 6.7 M
    lines (chr-1 subset alone has 1.4 M), making this take ~25 min.

    Run locally with::

        RUN_SLOW_INTEGRATION=1 pytest tests/test_z_tritici_subset.py::TestZTriticiSeqname1Slow -v
    """

    @pytest.fixture(scope="class")
    def seqname1_output(self, tmp_path_factory):
        if not os.environ.get("RUN_SLOW_INTEGRATION"):
            pytest.skip("Set RUN_SLOW_INTEGRATION=1 to run this slow test (~25 min)")
        output_dir = tmp_path_factory.mktemp("ztritici_seqname1_slow")
        result = _run_pipeline_seqname1(output_dir)
        return output_dir, result

    def test_pipeline_exits_zero(self, seqname1_output):
        _, result = seqname1_output
        assert result.returncode == 0, (
            f"Pipeline failed:\nSTDOUT:{result.stdout[-2000:]}\nSTDERR:{result.stderr[-2000:]}"
        )

    def test_nonzero_gene_count(self, seqname1_output):
        output_dir, _ = seqname1_output
        genes = _parse_gff3_genes(output_dir / "consensus.gff3")
        assert len(genes) > 0

    def test_all_genes_on_seqname_1(self, seqname1_output):
        """Subsetting must produce only chr-1 features."""
        output_dir, _ = seqname1_output
        genes = _parse_gff3_genes(output_dir / "consensus.gff3")
        non_chr1 = [g for g in genes if g["seqname"] != "1"]
        assert non_chr1 == [], (
            f"Found {len(non_chr1)} genes not on chromosome '1': "
            f"{[g['seqname'] for g in non_chr1[:5]]}"
        )
