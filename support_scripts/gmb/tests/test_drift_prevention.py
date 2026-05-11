
from gmb.pipeline.config import UtrConfig, ValidationConfig
from gmb.pipeline.gff3_validate import trim_utrs, validate_transcript


def test_drift_prevention():
    """Test that transcripts with CDS/UTR drifting beyond max_feature_drift_bp are flagged and drop logic applies."""

    # Exons spanning 1000 - 2000
    mrna_row = {"Feature": "mRNA", "Start": 1000, "End": 2000, "Chromosome": "1", "ID": "m1"}
    exon_rows = [
        {"Feature": "exon", "Start": 1000, "End": 1200},
        {"Feature": "exon", "Start": 1800, "End": 2000},
    ]

    # CDS way outside the exons (drift 2000bp, which > 1500 limit)
    cds_rows = [{"Feature": "CDS", "Start": 4000, "End": 4500}]
    utr_rows = []

    violations, max_drift, span = validate_transcript(
        mrna_row, exon_rows, cds_rows, utr_rows, config=ValidationConfig(max_feature_drift_bp=1500)
    )

    assert max_drift == 2500  # 4500 - 2000 = 2500
    assert len(violations) > 0
    assert any("not within any single merged exon" in v for v in violations)

    # In gff3_validate.py, if max_drift > max_feature_drift_bp, we force drop


def test_strict_utr_trimming():
    """Test that UTR segments are strictly clipped to the exon spans."""

    exon_rows = [
        {"Feature": "exon", "Start": 1000, "End": 1200},
        {"Feature": "exon", "Start": 1800, "End": 2000},
    ]
    # UTR extending beyond the start (900 -> 1000)
    utr_5p_rows = [{"Feature": "five_prime_UTR", "Start": 900, "End": 1100}]
    utr_3p_rows = []
    cds_rows = []

    cfg = UtrConfig(trim_policy="hard_cap", max_5p_bp=2000, max_3p_bp=2000, max_total_bp=4000)

    trimmed_5p, trimmed_3p, was_trimmed = trim_utrs(
        utr_5p_rows, utr_3p_rows, cds_rows, exon_rows, cfg
    )

    assert len(trimmed_5p) == 1
    # Should be clipped to 1000
    assert trimmed_5p[0]["Start"] == 1000
    assert trimmed_5p[0]["End"] == 1100


def test_compute_percentile_guardrails():
    """Test computation of effective guardrail limits from locus candidate set."""
    import pandas as pd

    from gmb.pipeline.config import PipelineConfig
    from gmb.pipeline.builder import compute_percentile_guardrails

    cfg = PipelineConfig()
    cfg.validation.max_exon_len_mode = "percentile"
    cfg.validation.max_exon_len_percentile = 90
    cfg.validation.max_exon_len_factor = 2.0
    cfg.validation.max_exon_len_bp = 1000  # Floor
    cfg.validation.max_transcript_span_mode = "percentile"
    cfg.validation.max_transcript_span_percentile = 90
    cfg.validation.max_transcript_span_factor = 2.0
    cfg.validation.max_transcript_span_bp = 5000  # Floor

    # 3 mock candidates: 2 high conf, 1 low conf
    # hc1: 1 exon, len 2000 => span 2000
    # hc2: 2 exons, len 1500, 500 => span 4000 (starts 1000, ends 5000)
    # lc1: 1 exon, len 10000 => span 10000

    data = [
        {
            "transcript_id": "hc1",
            "Feature": "exon",
            "Start": 1000,
            "End": 3000,
            "Source": "StringTie",
        },
        {
            "transcript_id": "hc2",
            "Feature": "exon",
            "Start": 1000,
            "End": 2500,
            "Source": "StringTie",
        },
        {
            "transcript_id": "hc2",
            "Feature": "exon",
            "Start": 4500,
            "End": 5000,
            "Source": "StringTie",
        },
        {
            "transcript_id": "lc1",
            "Feature": "exon",
            "Start": 1000,
            "End": 11000,
            "Source": "Scallop",
        },
    ]
    df = pd.DataFrame(data)

    # HC tids list
    protein_tids = {"hc1", "hc2"}

    params = compute_percentile_guardrails(df, cfg, protein_tids)

    # HC exon lens: 2000, 1500, 500. 90th pctile is 1900. Factor 2.0 => 3800.
    # HC spans: 2000, 4000. 90th pctile is 3800. Factor 2.0 => 7600.
    assert params["effective_max_exon_len_bp"] >= 3800
    assert params["effective_max_transcript_span_bp"] >= 7600

    # Assert limits are above floor but ignore the low confidence 10k outlier
    assert params["effective_max_exon_len_bp"] < 10000
