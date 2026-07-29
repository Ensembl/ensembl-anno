"""Tests for gmb.compare.tool_performance_analysis module."""

import json
import os
import tempfile

import pandas as pd
import pytest

from gmb.compare.tool_performance_analysis import (
    assign_strata,
    compute_boundary_analysis,
    compute_cross_ranking,
    compute_error_taxonomy,
    compute_intron_chain_analysis,
    compute_overall_profile,
    compute_reference_features,
    compute_stratified_metrics,
    load_comparison_details,
    load_comparison_summary,
    run_analysis,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def sample_ref_details():
    """Sample reference-perspective comparison details."""
    return pd.DataFrame(
        {
            "source": ["reference"] * 10,
            "gene_id": [f"GENE{i:03d}" for i in range(10)],
            "chrom": ["1"] * 10,
            "start": [1000 * i for i in range(10)],
            "end": [1000 * i + 900 for i in range(10)],
            "strand": ["+"] * 10,
            "classification": [
                "Exact_Match",
                "Exact_Match",
                "Structural_Mismatch",
                "Structural_Mismatch",
                "Partial_Match",
                "Partial_Match",
                "Partial_Match",
                "Missed",
                "Strand_Mismatch",
                "Exact_Match",
            ],
            "classification_cds": [
                "Exact_Match",
                "Exact_Match",
                "Exact_Match",
                "Structural_Mismatch",
                "Partial_Match",
                "Missed",
                "Partial_Match",
                "Missed",
                "Strand_Mismatch",
                "No_CDS",
            ],
            "matched_id": [
                "Q001",
                "Q002",
                "Q003",
                "Q003",
                "Q005",
                "Q006",
                "Q007",
                "",
                "Q009",
                "Q010",
            ],
            "best_match_transcript_id": [
                "Q001.1",
                "Q002.1",
                "Q003.1",
                "Q003.1",
                "Q005.1",
                "Q006.1",
                "Q007.1",
                "",
                "Q009.1",
                "Q010.1",
            ],
            "exon_overlap": [1.0, 0.95, 0.9, 0.85, 0.7, 0.5, 0.6, 0.0, 0.0, 0.98],
            "intron_chain_match": [
                True,
                True,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                True,
            ],
            "cds_overlap": [1.0, 0.95, 0.9, 0.8, 0.6, 0.0, 0.5, 0.0, 0.0, 0.0],
            "cds_intron_chain_match": [
                True,
                True,
                True,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
            ],
            "match_basis": ["exon"] * 10,
            "tool": ["TestTool"] * 10,
        }
    )


@pytest.fixture
def sample_cons_details():
    """Sample consensus-perspective comparison details."""
    return pd.DataFrame(
        {
            "source": ["consensus"] * 8,
            "gene_id": [f"Q{i:03d}" for i in range(1, 9)],
            "chrom": ["1"] * 8,
            "start": [1000 * i for i in range(8)],
            "end": [1000 * i + 900 for i in range(8)],
            "strand": ["+"] * 8,
            "classification": [
                "Exact_Match",
                "Exact_Match",
                "Structural_Mismatch",
                "Novel",
                "Partial_Match",
                "Partial_Match",
                "Partial_Match",
                "Novel",
            ],
            "classification_cds": [
                "Exact_Match",
                "Exact_Match",
                "Structural_Mismatch",
                "Novel",
                "Partial_Match",
                "Partial_Match",
                "Partial_Match",
                "Novel",
            ],
            "matched_id": [
                "GENE000",
                "GENE001",
                "GENE002",
                "",
                "GENE004",
                "GENE005",
                "GENE006",
                "",
            ],
            "best_match_transcript_id": [""] * 8,
            "exon_overlap": [1.0, 0.95, 0.9, 0.0, 0.7, 0.5, 0.6, 0.0],
            "intron_chain_match": [True, True, False, False, False, False, False, False],
            "cds_overlap": [1.0, 0.95, 0.9, 0.0, 0.6, 0.0, 0.5, 0.0],
            "cds_intron_chain_match": [True, True, True, False, False, False, False, False],
            "match_basis": ["exon"] * 8,
            "tool": ["TestTool"] * 8,
        }
    )


@pytest.fixture
def sample_features():
    """Sample reference feature DataFrame."""
    return pd.DataFrame(
        {
            "gene_length": [500, 2000, 5000, 10000, 50000, 1000, 3000, 80000, 15000, 900],
            "exon_count": [1, 2, 5, 12, 20, 1, 3, 15, 8, 1],
            "cds_length": [300, 1500, 3000, 8000, 20000, 0, 2000, 40000, 6000, 500],
            "transcript_count": [1, 3, 5, 10, 25, 1, 2, 30, 8, 1],
            "intron_count": [0, 1, 4, 11, 19, 0, 2, 14, 7, 0],
            "is_single_exon": [True, False, False, False, False, True, False, False, False, True],
        },
        index=[f"GENE{i:03d}" for i in range(10)],
    )
    # index name = gene_id


# ---------------------------------------------------------------------------
# Tests: compute_overall_profile
# ---------------------------------------------------------------------------


class TestOverallProfile:
    def test_basic_counts(self, sample_ref_details):
        profile = compute_overall_profile(sample_ref_details)
        assert profile["sampled_ref_genes"] == 10
        assert profile["exact_transcript_match_rate"] == 3 / 10
        assert profile["missed_rate"] == 1 / 10
        assert profile["strand_mismatch_rate"] == 1 / 10

    def test_locus_detection(self, sample_ref_details):
        profile = compute_overall_profile(sample_ref_details)
        # 10 - 1 missed = 9 detected
        assert profile["locus_detection_rate"] == 9 / 10

    def test_intron_chain(self, sample_ref_details):
        profile = compute_overall_profile(sample_ref_details)
        # 3 IC matches out of 9 detected
        assert profile["intron_chain_recovery_rate"] == 3 / 9

    def test_cds_exact(self, sample_ref_details):
        profile = compute_overall_profile(sample_ref_details)
        # 3 CDS exact matches (GENE000, GENE001, GENE002)
        assert profile["cds_exact_match_rate"] == 3 / 10

    def test_empty_df(self):
        empty = pd.DataFrame(
            columns=["classification", "classification_cds", "intron_chain_match", "tool"]
        )
        profile = compute_overall_profile(empty)
        assert profile == {}


# ---------------------------------------------------------------------------
# Tests: compute_error_taxonomy
# ---------------------------------------------------------------------------


class TestErrorTaxonomy:
    def test_basic_taxonomy(self, sample_ref_details, sample_cons_details):
        tax = compute_error_taxonomy(sample_ref_details, sample_cons_details)
        assert tax["tool"] == "TestTool"
        assert tax["exact_transcript"] == 3
        assert tax["missed_gene"] == 1
        assert tax["wrong_strand"] == 1

    def test_fragmentation_detected(self, sample_ref_details, sample_cons_details):
        # GENE002 and GENE003 both match Q003 -> fragmentation
        tax = compute_error_taxonomy(sample_ref_details, sample_cons_details)
        assert tax["fragmented_predictions"] >= 1

    def test_cds_exact_utr_diff(self, sample_ref_details, sample_cons_details):
        # GENE002 has CDS=Exact_Match but classification=Structural_Mismatch
        tax = compute_error_taxonomy(sample_ref_details, sample_cons_details)
        assert tax["exact_cds_utr_differs"] >= 1


# ---------------------------------------------------------------------------
# Tests: assign_strata
# ---------------------------------------------------------------------------


class TestAssignStrata:
    def test_exon_strata(self, sample_features):
        strata = assign_strata(sample_features)
        # GENE000 has 1 exon -> "1_exon"
        assert strata.loc["GENE000", "exon_stratum"] == "1_exon"
        # GENE001 has 2 exons -> "2-3_exon"
        assert strata.loc["GENE001", "exon_stratum"] == "2-3_exon"
        # GENE002 has 5 exons -> "4-10_exon"
        assert strata.loc["GENE002", "exon_stratum"] == "4-10_exon"
        # GENE003 has 12 exons -> ">10_exon"
        assert strata.loc["GENE003", "exon_stratum"] == ">10_exon"

    def test_intron_strata(self, sample_features):
        strata = assign_strata(sample_features)
        assert strata.loc["GENE000", "intron_stratum"] == "0_introns"
        assert strata.loc["GENE001", "intron_stratum"] == "1-2_introns"
        assert strata.loc["GENE002", "intron_stratum"] == "3-9_introns"
        assert strata.loc["GENE003", "intron_stratum"] == ">=10_introns"

    def test_transcript_strata(self, sample_features):
        strata = assign_strata(sample_features)
        assert strata.loc["GENE000", "transcript_stratum"] == "mono_tx"
        assert strata.loc["GENE001", "transcript_stratum"] == "2-5_tx"
        assert strata.loc["GENE003", "transcript_stratum"] == "6-20_tx"
        assert strata.loc["GENE004", "transcript_stratum"] == ">20_tx"

    def test_all_strata_assigned(self, sample_features):
        strata = assign_strata(sample_features)
        assert "exon_stratum" in strata.columns
        assert "intron_stratum" in strata.columns
        assert "gene_length_stratum" in strata.columns
        assert "cds_length_stratum" in strata.columns
        assert "transcript_stratum" in strata.columns
        # No NaN in exon_stratum
        assert strata["exon_stratum"].isna().sum() == 0


# ---------------------------------------------------------------------------
# Tests: compute_stratified_metrics
# ---------------------------------------------------------------------------


class TestStratifiedMetrics:
    def test_produces_rows(self, sample_ref_details, sample_features):
        strata = compute_stratified_metrics(sample_ref_details, sample_features)
        assert len(strata) > 0
        assert "tool" in strata.columns
        assert "stratum_type" in strata.columns
        assert "metric" in strata.columns
        assert "value" in strata.columns

    def test_contains_exon_strata(self, sample_ref_details, sample_features):
        strata = compute_stratified_metrics(sample_ref_details, sample_features)
        exon_rows = strata[strata["stratum_type"] == "exon"]
        assert len(exon_rows) > 0

    def test_metrics_in_valid_range(self, sample_ref_details, sample_features):
        strata = compute_stratified_metrics(sample_ref_details, sample_features)
        assert (strata["value"] >= 0).all()
        assert (strata["value"] <= 1).all()


# ---------------------------------------------------------------------------
# Tests: compute_intron_chain_analysis
# ---------------------------------------------------------------------------


class TestIntronChainAnalysis:
    def test_multi_exon_only(self, sample_ref_details, sample_features):
        ic = compute_intron_chain_analysis(sample_ref_details, sample_features)
        # Only genes with exon_count > 1 are multi-exon (7 of 10)
        assert ic["n_multi_exon_ref"] == 7

    def test_ic_rates(self, sample_ref_details, sample_features):
        ic = compute_intron_chain_analysis(sample_ref_details, sample_features)
        assert 0 <= ic["exact_intron_chain_rate"] <= 1


# ---------------------------------------------------------------------------
# Tests: compute_boundary_analysis
# ---------------------------------------------------------------------------


class TestBoundaryAnalysis:
    def test_basic_boundary(self, sample_ref_details):
        ref_genes = pd.DataFrame(
            {
                "ID": [f"GENE{i:03d}" for i in range(10)],
                "Chromosome": ["1"] * 10,
                "Start": [1000 * i for i in range(10)],
                "End": [1000 * i + 900 for i in range(10)],
                "Strand": ["+"] * 10,
            }
        )
        query_genes = pd.DataFrame(
            {
                "ID": [f"Q{i:03d}" for i in range(1, 11)],
                "Chromosome": ["1"] * 10,
                "Start": [1000 * i + 5 for i in range(10)],  # 5bp offset
                "End": [1000 * i + 905 for i in range(10)],
                "Strand": ["+"] * 10,
            }
        )
        bnd = compute_boundary_analysis(sample_ref_details, ref_genes, {"TestTool": query_genes})
        assert not bnd.empty
        # Should have rows for start, end, both at various thresholds
        assert "start" in bnd["boundary"].values
        assert "end" in bnd["boundary"].values
        assert "both" in bnd["boundary"].values

    def test_exact_boundaries(self):
        """When ref and query have identical boundaries, exact rate should be 1.0."""
        ref_details = pd.DataFrame(
            {
                "source": ["reference"] * 3,
                "gene_id": ["G1", "G2", "G3"],
                "classification": ["Exact_Match", "Structural_Mismatch", "Partial_Match"],
                "matched_id": ["Q1", "Q2", "Q3"],
                "tool": ["MyTool"] * 3,
            }
        )
        ref_genes = pd.DataFrame(
            {
                "ID": ["G1", "G2", "G3"],
                "Chromosome": ["1"] * 3,
                "Start": [100, 500, 900],
                "End": [200, 600, 1000],
                "Strand": ["+"] * 3,
            }
        )
        query_genes = pd.DataFrame(
            {
                "ID": ["Q1", "Q2", "Q3"],
                "Chromosome": ["1"] * 3,
                "Start": [100, 500, 900],
                "End": [200, 600, 1000],
                "Strand": ["+"] * 3,
            }
        )
        bnd = compute_boundary_analysis(ref_details, ref_genes, {"MyTool": query_genes})
        exact_start = bnd[(bnd["boundary"] == "start") & (bnd["threshold_bp"] == 0)]
        assert exact_start["rate_within"].iloc[0] == 1.0


# ---------------------------------------------------------------------------
# Tests: compute_cross_ranking
# ---------------------------------------------------------------------------


class TestCrossRanking:
    def test_ranking_order(self):
        profiles = {
            "ToolA": {
                "locus_detection_rate": 0.9,
                "exact_transcript_match_rate": 0.3,
                "cds_exact_match_rate": 0.5,
                "cds_any_match_rate": 0.8,
                "intron_chain_recovery_rate": 0.4,
                "missed_rate": 0.1,
            },
            "ToolB": {
                "locus_detection_rate": 0.85,
                "exact_transcript_match_rate": 0.4,
                "cds_exact_match_rate": 0.6,
                "cds_any_match_rate": 0.85,
                "intron_chain_recovery_rate": 0.5,
                "missed_rate": 0.15,
            },
        }
        strata = pd.DataFrame()
        ranking = compute_cross_ranking(strata, profiles)
        assert len(ranking) > 0
        # For exact_transcript_match_rate, ToolB (0.4) should rank first
        exact_row = ranking[ranking["metric"] == "exact_transcript_match_rate"].iloc[0]
        assert exact_row["rank_1"] == "ToolB"
        # For missed_rate (lower is better), ToolA (0.1) should rank first
        missed_row = ranking[ranking["metric"] == "missed_rate"].iloc[0]
        assert missed_row["rank_1"] == "ToolA"


# ---------------------------------------------------------------------------
# Tests: File I/O
# ---------------------------------------------------------------------------


class TestFileIO:
    def test_load_comparison_details(self, tmp_path):
        # Create a minimal comparison_details.tsv
        df = pd.DataFrame(
            {
                "source": ["reference", "consensus"],
                "gene_id": ["G1", "Q1"],
                "chrom": ["1", "1"],
                "start": [100, 100],
                "end": [200, 200],
                "strand": ["+", "+"],
                "classification": ["Exact_Match", "Exact_Match"],
                "classification_cds": ["Exact_Match", "Exact_Match"],
                "matched_id": ["Q1", "G1"],
                "best_match_transcript_id": ["Q1.1", "G1.1"],
                "exon_overlap": [1.0, 1.0],
                "intron_chain_match": [True, True],
                "cds_overlap": [1.0, 1.0],
                "cds_intron_chain_match": [True, True],
                "match_basis": ["exon", "exon"],
            }
        )
        df.to_csv(tmp_path / "comparison_details.tsv", sep="\t", index=False)

        loaded = load_comparison_details(str(tmp_path), "test_tool")
        assert len(loaded) == 2
        assert loaded["tool"].iloc[0] == "test_tool"

    def test_load_comparison_summary(self, tmp_path):
        summary = {"total_reference_genes": 100, "sensitivity": {"exact_match_rate": 0.5}}
        with open(tmp_path / "comparison_summary.json", "w") as f:
            json.dump(summary, f)

        loaded = load_comparison_summary(str(tmp_path))
        assert loaded["total_reference_genes"] == 100

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_comparison_details(str(tmp_path / "nonexistent"), "tool")


# ---------------------------------------------------------------------------
# Tests: Integration (run_analysis without reference)
# ---------------------------------------------------------------------------


class TestIntegration:
    def test_run_without_reference(self, tmp_path):
        """Test run_analysis with just comparison outputs, no reference GFF."""
        # Create minimal comparison outputs for two tools
        for tool_name in ["tool_a", "tool_b"]:
            tool_dir = tmp_path / f"compare_{tool_name}"
            tool_dir.mkdir()

            details = pd.DataFrame(
                {
                    "source": ["reference", "reference", "consensus", "consensus"],
                    "gene_id": ["G1", "G2", "Q1", "Q2"],
                    "chrom": ["1", "1", "1", "1"],
                    "start": [100, 500, 100, 500],
                    "end": [200, 600, 200, 600],
                    "strand": ["+", "+", "+", "+"],
                    "classification": ["Exact_Match", "Missed", "Exact_Match", "Novel"],
                    "classification_cds": ["Exact_Match", "Missed", "Exact_Match", "Novel"],
                    "matched_id": ["Q1", "", "G1", ""],
                    "best_match_transcript_id": ["Q1.1", "", "G1.1", ""],
                    "exon_overlap": [1.0, 0.0, 1.0, 0.0],
                    "intron_chain_match": [True, False, True, False],
                    "cds_overlap": [1.0, 0.0, 1.0, 0.0],
                    "cds_intron_chain_match": [True, False, True, False],
                    "match_basis": ["exon", "exon", "exon", "exon"],
                }
            )
            details.to_csv(tool_dir / "comparison_details.tsv", sep="\t", index=False)

            summary = {
                "total_reference_genes": 2,
                "total_consensus_genes": 2,
                "sensitivity": {"exact_match_rate": 0.5, "missed_rate": 0.5},
                "query_path": "fake.gff",
            }
            with open(tool_dir / "comparison_summary.json", "w") as f:
                json.dump(summary, f)

        out_dir = str(tmp_path / "output")
        result = run_analysis(
            comparisons={
                "tool_a": str(tmp_path / "compare_tool_a"),
                "tool_b": str(tmp_path / "compare_tool_b"),
            },
            reference_gff=None,
            output_dir=out_dir,
        )

        assert os.path.exists(os.path.join(out_dir, "tool_performance_summary.tsv"))
        assert os.path.exists(os.path.join(out_dir, "tool_performance_summary.json"))
        assert os.path.exists(os.path.join(out_dir, "tool_error_taxonomy.tsv"))
        assert os.path.exists(os.path.join(out_dir, "tool_cross_ranking.tsv"))

        # Check profiles
        assert "tool_a" in result["overall_profiles"]
        assert result["overall_profiles"]["tool_a"]["exact_transcript_match_rate"] == 0.5
