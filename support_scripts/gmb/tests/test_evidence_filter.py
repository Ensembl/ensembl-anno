#!/usr/bin/env python3
"""Tests for evidence_filter module."""

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from config import load_config
from evidence_filter import filter_chimeras, filter_helixer_models, filter_protein_evidence


@pytest.fixture
def config():
    return load_config()


# ---------------------------------------------------------------------------
# Protein evidence filtering
# ---------------------------------------------------------------------------


class TestFilterProteinEvidence:
    def _make_protein_df(self, transcripts):
        """Build exon-level DataFrame from transcript specs.

        transcripts: list of (tid, chrom, strand, exons)
        exons: list of (start, end)
        """
        rows = []
        for tid, chrom, strand, exons in transcripts:
            for s, e in exons:
                rows.append(
                    {
                        "Chromosome": chrom,
                        "Start": s,
                        "End": e,
                        "Strand": strand,
                        "transcript_id": tid,
                        "Source": "OrthoDB",
                        "Feature": "exon",
                    }
                )
        return pd.DataFrame(rows)

    def test_remove_short_fragments(self, config):
        """Fragments shorter than min_protein_aa * 3 bp should be removed."""
        # Fragment: 80bp span (< 90bp threshold for 30aa)
        # Adequate: 500bp span
        df = self._make_protein_df(
            [
                ("frag1", "1", "+", [(100, 180)]),  # 80bp
                ("good1", "1", "+", [(1000, 1500)]),  # 500bp
            ]
        )
        stats = {}
        result = filter_protein_evidence(df, config, stats=stats)
        tids = set(result["transcript_id"].unique())
        assert "frag1" not in tids
        assert "good1" in tids
        assert stats["protein_short_fragments_removed"] >= 1

    def test_keep_adequate_protein(self, config):
        """Protein models >= min_protein_aa should be retained."""
        df = self._make_protein_df(
            [
                ("prot1", "1", "+", [(100, 300)]),  # 200bp (66aa)
            ]
        )
        result = filter_protein_evidence(df, config)
        assert "prot1" in result["transcript_id"].values

    def test_collapse_redundant(self, config):
        """Two models at same locus with high overlap → collapse to one."""
        # Two nearly identical models
        df = self._make_protein_df(
            [
                ("prot1", "1", "+", [(1000, 1500)]),
                ("prot2", "1", "+", [(1010, 1510)]),
            ]
        )
        result = filter_protein_evidence(df, config)
        n = result["transcript_id"].nunique()
        assert n == 1, f"Expected 1 after collapse, got {n}"

    def test_locus_competition_top_k(self, config):
        """More than top_n models at one locus → keep top N."""
        # Create 10 overlapping models
        transcripts = []
        for i in range(10):
            offset = i * 5
            transcripts.append((f"prot{i}", "1", "+", [(1000 + offset, 1500 + offset)]))
        df = self._make_protein_df(transcripts)
        stats = {}
        result = filter_protein_evidence(df, config, stats=stats)
        n = result["transcript_id"].nunique()
        # After redundancy collapse + top_n, should be limited
        assert n <= config.protein_filter.top_n_per_locus + 2

    def test_remove_long_artifacts(self, config):
        """Protein spans > max_span_bp should be removed."""
        df = self._make_protein_df(
            [
                ("long1", "1", "+", [(100, 60000)]),  # > 50kb
                ("good1", "1", "+", [(200, 500)]),
            ]
        )
        stats = {}
        result = filter_protein_evidence(df, config, stats=stats)
        assert "long1" not in result["transcript_id"].values
        assert stats["protein_long_artifacts_removed"] >= 1

    def test_empty_input(self, config):
        """Empty DataFrame → empty result."""
        df = pd.DataFrame(
            columns=["Chromosome", "Start", "End", "Strand", "transcript_id", "Source", "Feature"]
        )
        result = filter_protein_evidence(df, config)
        assert result.empty


# ---------------------------------------------------------------------------
# Chimera detection
# ---------------------------------------------------------------------------


class TestFilterChimeras:
    def test_large_intron_removed(self, config):
        """Multi-exon transcript with intron > max_intron_length is removed."""
        df = pd.DataFrame(
            {
                "Chromosome": ["1", "1"],
                "Start": [100, 5000],
                "End": [200, 5100],
                "Strand": ["+", "+"],
                "transcript_id": ["chimera1", "chimera1"],
                "Source": ["Scallop", "Scallop"],
            }
        )
        stats = {}
        result = filter_chimeras(df, config, stats=stats)
        assert result.empty or "chimera1" not in result["transcript_id"].values
        assert stats["chimeras_large_intron"] >= 1

    def test_single_exon_not_penalized(self, config):
        """Single-exon transcripts should always be kept."""
        df = pd.DataFrame(
            {
                "Chromosome": ["1"],
                "Start": [100],
                "End": [200],
                "Strand": ["+"],
                "transcript_id": ["single1"],
                "Source": ["Scallop"],
            }
        )
        result = filter_chimeras(df, config)
        assert "single1" in result["transcript_id"].values

    def test_normal_intron_kept(self, config):
        """Multi-exon transcript with normal intron stays."""
        df = pd.DataFrame(
            {
                "Chromosome": ["1", "1"],
                "Start": [100, 300],
                "End": [200, 400],
                "Strand": ["+", "+"],
                "transcript_id": ["normal1", "normal1"],
                "Source": ["Scallop", "Scallop"],
            }
        )
        result = filter_chimeras(df, config)
        assert "normal1" in result["transcript_id"].values


# ---------------------------------------------------------------------------
# Helixer model filtering
# ---------------------------------------------------------------------------


class TestFilterHelixerModels:
    def test_short_cds_removed(self, config):
        """Helixer model with CDS < min_cds_bp is removed."""
        exons = pd.DataFrame(
            {
                "Chromosome": ["1"],
                "Start": [100],
                "End": [200],
                "Strand": ["+"],
                "transcript_id": ["hx1"],
                "Source": ["Helixer"],
            }
        )
        cds = pd.DataFrame(
            {
                "Chromosome": ["1"],
                "Start": [100],
                "End": [150],  # 50bp < 90bp
                "Strand": ["+"],
                "transcript_id": ["hx1"],
                "Source": ["Helixer"],
            }
        )
        stats = {}
        f_ex, f_cds = filter_helixer_models(exons, cds, config, stats=stats)
        assert f_ex.empty or "hx1" not in f_ex["transcript_id"].values

    def test_adequate_cds_kept(self, config):
        """Helixer model with CDS >= min_cds_bp is kept."""
        exons = pd.DataFrame(
            {
                "Chromosome": ["1"],
                "Start": [100],
                "End": [400],
                "Strand": ["+"],
                "transcript_id": ["hx2"],
                "Source": ["Helixer"],
            }
        )
        cds = pd.DataFrame(
            {
                "Chromosome": ["1"],
                "Start": [100],
                "End": [400],  # 300bp
                "Strand": ["+"],
                "transcript_id": ["hx2"],
                "Source": ["Helixer"],
            }
        )
        f_ex, f_cds = filter_helixer_models(exons, cds, config)
        assert "hx2" in f_ex["transcript_id"].values


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
