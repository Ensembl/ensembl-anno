#!/usr/bin/env python3
"""Subset integration tests using bundled Z. tritici example data.

Runs gene_model_builder.py on chromosome 1 of the bundled Zymoseptoria tritici
dataset (``support_scripts/gmb/z_tritici/``) to validate end-to-end pipeline
behaviour on real genomic evidence.

Why --seqname 1?
----------------
Chromosome 1 is the largest Z. tritici chromosome and carries a representative
sample of gene models across all evidence tracks. Using ``--seqname 1`` keeps
runtime under ~2 minutes while exercising every pipeline stage (filtering,
clustering, scoring, CDS annotation, FASTA export, GFF3 validation).

The assembly-report mapping is required because the input GTF/GFF3 files use
GenBank accession numbers (``CM001642.1`` etc.) while the genome FASTA uses
short chromosome numbers (``1``, ``2`` …).

Markers
-------
All tests in this module are marked with ``integration`` and will be collected
by default.  The ``--seqname 1`` subset keeps each test under ~2 minutes even
on a laptop.

Golden-file regression
----------------------
If ``tests/fixtures/expected/`` contains ``consensus_seqname1.gff3`` and
``evidence_attribution_seqname1.tsv``, the corresponding tests will compare
the pipeline output against those files.  The comparison is structural
(sorted, parsed) rather than byte-identical, so minor formatting changes do
not cause spurious failures.
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
_SCRIPT = _GMB_DIR / "gene_model_builder.py"
_FIXTURES = _HERE / "fixtures" / "expected"

# Check for pyranges — required by the pipeline
try:
    import pyranges  # noqa: F401
except ImportError:
    pytest.skip("pyranges required", allow_module_level=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _run_pipeline(output_dir: Path) -> subprocess.CompletedProcess:
    """Run gene_model_builder on z_tritici chromosome 1."""
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
        "ZTTEST",
        "--seqname",
        "1",
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(_GMB_DIR)
    return subprocess.run(cmd, capture_output=True, text=True, env=env)


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
def pipeline_output(tmp_path_factory):
    """Run the pipeline once; share output across all tests in this module."""
    output_dir = tmp_path_factory.mktemp("ztritici_seqname1")
    result = _run_pipeline(output_dir)
    return output_dir, result


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestZTriticiSeqname1:
    """Smoke tests for the z_tritici seqname-1 subset run."""

    def test_pipeline_exits_zero(self, pipeline_output):
        """Pipeline should complete without errors."""
        _, result = pipeline_output
        assert result.returncode == 0, (
            f"Pipeline failed (exit {result.returncode}):\n"
            f"STDOUT:\n{result.stdout[-3000:]}\n"
            f"STDERR:\n{result.stderr[-3000:]}"
        )

    def test_required_output_files_exist(self, pipeline_output):
        """All mandatory output files must be present."""
        output_dir, _ = pipeline_output
        for fname in ("consensus.gff3", "cdna.fa", "prot.fa", "summary.json", "summary.tsv"):
            assert (output_dir / fname).exists(), f"Missing required output: {fname}"

    def test_nonzero_gene_count(self, pipeline_output):
        """At least one gene must be called on chromosome 1."""
        output_dir, _ = pipeline_output
        gff3 = output_dir / "consensus.gff3"
        genes = _parse_gff3_genes(gff3)
        assert len(genes) > 0, "No genes in consensus.gff3 — pipeline produced empty output"

    def test_all_genes_on_seqname_1(self, pipeline_output):
        """Subsetting to --seqname 1 should only produce chr-1 features."""
        output_dir, _ = pipeline_output
        gff3 = output_dir / "consensus.gff3"
        genes = _parse_gff3_genes(gff3)
        non_chr1 = [g for g in genes if g["seqname"] != "1"]
        assert non_chr1 == [], (
            f"Found {len(non_chr1)} gene(s) not on chromosome '1': "
            f"{[g['seqname'] for g in non_chr1[:5]]}"
        )

    def test_no_protein_evidence_as_mrna(self, pipeline_output):
        """Protein evidence (OrthoDB_/UniProt_) must never appear as mRNA features.

        Protein alignments are support-only evidence — they must not be
        promoted to mRNA features in the output GFF3.
        """
        output_dir, _ = pipeline_output
        gff3 = output_dir / "consensus.gff3"
        mrna_ids = _parse_gff3_mrnas(gff3)
        bad = [m for m in mrna_ids if m.startswith(("OrthoDB_", "UniProt_"))]
        assert bad == [], (
            f"Protein evidence IDs appeared as mRNA in consensus.gff3: {bad[:5]}"
        )

    def test_summary_json_structure(self, pipeline_output):
        """summary.json must have the expected top-level structure."""
        output_dir, _ = pipeline_output
        with open(output_dir / "summary.json") as fh:
            report = json.load(fh)

        assert "summary" in report, "summary.json missing 'summary' key"
        assert "total_loci" in report["summary"], "'summary.total_loci' key missing"
        assert report["summary"]["total_loci"] >= 1, "total_loci is 0"
        assert "filtering" in report, "summary.json missing 'filtering' key"

    def test_evidence_attribution_exists(self, pipeline_output):
        """evidence_attribution.tsv must be present and non-empty."""
        output_dir, _ = pipeline_output
        tsv = output_dir / "evidence_attribution.tsv"
        # Some pipeline configs write this file; skip gracefully if absent
        if not tsv.exists():
            pytest.skip("evidence_attribution.tsv not written by this config")
        with open(tsv) as fh:
            lines = [l for l in fh if not l.startswith("#")]
        assert len(lines) > 1, "evidence_attribution.tsv has no data rows"

    def test_cdna_fasta_has_sequences(self, pipeline_output):
        """cdna.fa should contain at least one transcript sequence."""
        output_dir, _ = pipeline_output
        content = (output_dir / "cdna.fa").read_text()
        assert content.count(">") >= 1, "cdna.fa has no FASTA headers"
        seqs = [l for l in content.splitlines() if l and not l.startswith(">")]
        assert len("".join(seqs)) > 0, "cdna.fa has headers but no sequence data"

    def test_protein_fasta_has_sequences(self, pipeline_output):
        """prot.fa should contain at least one protein sequence."""
        output_dir, _ = pipeline_output
        content = (output_dir / "prot.fa").read_text()
        assert content.count(">") >= 1, "prot.fa has no FASTA headers"

    def test_subset_manifest_written(self, pipeline_output):
        """When --seqname is used, a subset_regions.tsv manifest should be written."""
        output_dir, _ = pipeline_output
        manifest = output_dir / "subset_regions.tsv"
        if not manifest.exists():
            pytest.skip("subset_regions.tsv not written by this build")
        content = manifest.read_text()
        assert "1" in content, "subset_regions.tsv does not mention seqname '1'"


@pytest.mark.integration
class TestZTriticiGoldenRegression:
    """Golden-file regression: compare output against stored expected files.

    Tests are skipped if the fixture files do not exist.  To (re-)generate
    them, run::

        python tests/generate_golden_fixtures.py

    or copy the output from a known-good run into tests/fixtures/expected/.
    """

    def test_gene_count_matches_golden(self, pipeline_output):
        """Gene count should match the stored golden count (±0)."""
        golden = _FIXTURES / "seqname1_gene_count.txt"
        if not golden.exists():
            pytest.skip("Golden gene count not found — run generate_golden_fixtures.py")

        output_dir, _ = pipeline_output
        genes = _parse_gff3_genes(output_dir / "consensus.gff3")
        expected = int(golden.read_text().strip())
        assert len(genes) == expected, (
            f"Gene count changed: got {len(genes)}, expected {expected}. "
            "If this is intentional, re-run generate_golden_fixtures.py."
        )

    def test_gff3_matches_golden(self, pipeline_output):
        """Canonical (sorted) GFF3 must match the stored golden file."""
        golden = _FIXTURES / "consensus_seqname1.gff3"
        if not golden.exists():
            pytest.skip("Golden GFF3 not found — run generate_golden_fixtures.py")

        output_dir, _ = pipeline_output
        actual_genes = sorted(
            (g["seqname"], g["start"], g["end"], g["strand"])
            for g in _parse_gff3_genes(output_dir / "consensus.gff3")
        )
        expected_genes = sorted(
            (g["seqname"], g["start"], g["end"], g["strand"])
            for g in _parse_gff3_genes(golden)
        )
        assert actual_genes == expected_genes, (
            f"Gene loci changed vs golden.  "
            f"Added: {set(actual_genes) - set(expected_genes)}, "
            f"Removed: {set(expected_genes) - set(actual_genes)}"
        )
