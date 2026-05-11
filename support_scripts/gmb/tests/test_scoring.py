#!/usr/bin/env python3
"""Tests for scoring module."""

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from gmb.pipeline.config import PipelineConfig, load_config
from gmb.pipeline.scoring import score_model, select_isoforms


@pytest.fixture
def config():
    return load_config()


def _make_model(tid, source, chrom="1", strand="+", exons=None, combined_evidence=None):
    """Create a model dict for testing."""
    if exons is None:
        exons = [(100, 500)]
    df = pd.DataFrame(
        {
            "Start": [s for s, e in exons],
            "End": [e for s, e in exons],
            "Chromosome": [chrom] * len(exons),
            "Strand": [strand] * len(exons),
        }
    )
    return {
        "id": tid,
        "source": source,
        "chrom": chrom,
        "strand": strand,
        "intron_chain": "single-exon" if len(exons) == 1 else "multi",
        "protein_support": False,
        "df": df,
        "start": exons[0][0],
        "end": exons[-1][1],
        "exon_count": len(exons),
        "combined_evidence": combined_evidence or source,
    }


def _make_locus_df(models):
    """Combine model dicts into a locus DataFrame."""
    dfs = []
    for m in models:
        df = m["df"].copy()
        df["Source"] = m["source"]
        df["transcript_id"] = m["id"]
        df["combined_evidence"] = m["combined_evidence"]
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


class TestScoreModel:
    def test_helixer_scores_higher(self, config):
        """Helixer model should score higher than single-source Scallop."""
        helixer = _make_model("hx1", "Helixer", combined_evidence="Helixer")
        scallop = _make_model("sc1", "Scallop", combined_evidence="Scallop")
        h_score = score_model(helixer, config, set())
        s_score = score_model(scallop, config, set())
        assert h_score > s_score

    def test_protein_support_bonus(self, config):
        """Model with protein support scores higher."""
        m1 = _make_model("tx1", "Scallop")
        s_no_prot = score_model(m1, config, set())
        s_with_prot = score_model(m1, config, {"tx1"})
        assert s_with_prot > s_no_prot
        assert s_with_prot - s_no_prot == config.scoring.protein_overlap_bonus

    def test_multi_source_bonus(self, config):
        """Multi-source model scores higher than single-source."""
        single = _make_model("tx1", "Scallop", combined_evidence="Scallop")
        multi = _make_model("tx2", "Scallop", combined_evidence="Scallop,StringTie")
        s1 = score_model(single, config, set())
        s2 = score_model(multi, config, set())
        assert s2 > s1


class TestSelectIsoforms:
    def test_max_isoforms_limit(self, config):
        """Only top-N isoforms should be selected per locus."""
        # Create models with different structures
        models = []
        for i in range(5):
            models.append(
                _make_model(
                    f"tx{i}",
                    "Helixer",
                    exons=[(100 + i * 10, 500 + i * 10)],
                    combined_evidence="Helixer",
                )
            )
        locus_df = _make_locus_df(models)
        selected = select_isoforms(locus_df, config, set())
        flat_selected = [m for g in selected for m in g]
        assert len(flat_selected) <= config.scoring.max_isoforms_per_locus

    def test_protein_supported_preferred(self, config):
        """Protein-supported model should be selected as primary."""
        m1 = _make_model("tx1", "Scallop", combined_evidence="Scallop")
        m2 = _make_model("tx2", "Scallop", exons=[(100, 500)], combined_evidence="Scallop")
        locus_df = _make_locus_df([m1, m2])
        selected = select_isoforms(locus_df, config, {"tx1"})
        flat_selected = [m for g in selected for m in g]
        # tx1 should be selected (protein support)
        selected_ids = {m["id"] for m in flat_selected}
        assert "tx1" in selected_ids

    def test_empty_locus(self, config):
        """Empty locus → empty result."""
        locus_df = pd.DataFrame(
            columns=[
                "Chromosome",
                "Start",
                "End",
                "Strand",
                "Source",
                "transcript_id",
                "combined_evidence",
            ]
        )
        selected = select_isoforms(locus_df, config, set())
        assert selected == []


class TestLowMinCodons:
    def test_short_orf_accepted_with_fungal_config(self):
        """ORFs of 40 codons should be accepted with fungal min_codons=33."""
        from gmb.pipeline.annotate_cds_utrs import find_best_orf

        # 40-codon ORF: ATG + 38 coding + stop = 40 codons
        seq = "ATG" + "GCT" * 38 + "TAA"  # 40 codons
        result = find_best_orf(seq, min_codons=33)
        assert result is not None
        assert (result[1] - result[0]) // 3 >= 33

    def test_very_short_orf_still_found(self):
        """Even 20-codon ORFs should be found (as fallback)."""
        from gmb.pipeline.annotate_cds_utrs import find_best_orf

        seq = "ATG" + "GCT" * 18 + "TAA"  # 20 codons
        result = find_best_orf(seq, min_codons=33)
        assert result is not None  # fallback should still return it


class TestSameGeneOverlapThreshold:
    """Verify that same_gene_overlap_threshold controls absorption of a small
    Helixer gene into a large transcriptome-spanning model.

    Geometry
    --------
    Large TX (StringTie, multi-exon, protein-supported): 1000–5000 (span 4000)
    Small Helixer (single-exon):                         4800–5600 (span  800)

    Overlap = 5000 − 4800 = 200 bp
    Overlap ratio relative to smaller model = 200 / 800 = 0.25

    At threshold 0.15  → 0.25 > 0.15 → same gene → Helixer absorbed (dropped)
    At threshold 0.30  → 0.25 < 0.30 → separate  → Helixer survives
    """

    def _make_locus_df(self):
        """Build a locus DataFrame with the two-model geometry described above."""
        tx_rows = pd.DataFrame(
            {
                "Start": [1000, 3000],
                "End": [2500, 5000],
                "Chromosome": ["1", "1"],
                "Strand": ["+", "+"],
                "Source": ["StringTie", "StringTie"],
                "transcript_id": ["tx1", "tx1"],
                # combined_evidence carries the multi-source tag used in scoring
                "combined_evidence": ["StringTie,Scallop", "StringTie,Scallop"],
            }
        )
        hx_rows = pd.DataFrame(
            {
                "Start": [4800],
                "End": [5600],
                "Chromosome": ["1"],
                "Strand": ["+"],
                "Source": ["Helixer"],
                "transcript_id": ["hx1"],
                "combined_evidence": ["Helixer"],
            }
        )
        return pd.concat([tx_rows, hx_rows], ignore_index=True)

    def _cfg(self, threshold: float) -> PipelineConfig:
        """Return a PipelineConfig with the given overlap threshold.

        max_isoforms_per_locus=1 ensures that once a primary is claimed at a
        locus there is no room for a secondary — any absorbed model is dropped.
        """
        cfg = PipelineConfig()
        cfg.scoring.same_gene_overlap_threshold = threshold
        cfg.scoring.max_isoforms_per_locus = 1
        return cfg

    def test_low_threshold_absorbs_helixer(self):
        """At threshold=0.15 the Helixer gene (overlap ratio 0.25) is merged
        into the TX locus and dropped because the locus is already full."""
        genes = select_isoforms(self._make_locus_df(), self._cfg(0.15), {"tx1"})
        selected_ids = {m["id"] for g in genes for m in g}

        assert "tx1" in selected_ids, "TX model should be selected"
        assert "hx1" not in selected_ids, (
            "Helixer gene should be absorbed into the TX locus at threshold=0.15"
        )

    def test_high_threshold_preserves_helixer(self):
        """At threshold=0.30 the same Helixer gene (overlap ratio 0.25) is
        treated as a separate locus and survives selection."""
        genes = select_isoforms(self._make_locus_df(), self._cfg(0.30), {"tx1"})
        selected_ids = {m["id"] for g in genes for m in g}

        assert "tx1" in selected_ids, "TX model should be selected"
        assert "hx1" in selected_ids, (
            "Helixer gene should survive as a separate locus at threshold=0.30"
        )
        assert len(genes) == 2, (
            f"Expected 2 separate gene loci, got {len(genes)}: "
            f"{[m['id'] for g in genes for m in g]}"
        )

    def test_threshold_is_a_strict_boundary(self):
        """A threshold set exactly equal to the overlap ratio (0.25) should
        NOT trigger a merge (> comparison, not >=)."""
        genes = select_isoforms(self._make_locus_df(), self._cfg(0.25), {"tx1"})
        selected_ids = {m["id"] for g in genes for m in g}

        assert "hx1" in selected_ids, (
            "Overlap ratio equal to threshold should not trigger merge "
            "(the comparison is strictly >, not >=)"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
