"""Tests for gmb.pipeline.canonical_selection.

Covers: single-transcript genes, isoforms differing on Psauron/DIAMOND
coverage/evidence breadth, incomplete-ORF penalties, missing optional
protein-validation fields, deterministic tie-breaking, domain evidence
present-but-disabled vs enabled, canonical-attribute GFF3 writing, and that
alternative isoforms are never dropped or mutated.
"""

import json
import os

import pandas as pd
import pytest

from gmb.pipeline.canonical_selection import (
    _write_annotated_gff3,
    load_transcript_records,
    run_canonical_selection,
    score_transcript,
    select_canonical_for_gene,
)
from gmb.pipeline.config import load_config


def _rec(tid, gene_id="G1", **overrides):
    base = {
        "gene_id": gene_id,
        "transcript_id": tid,
        "evidence_sources": "Scallop,StringTie",
        "exon_count": 2,
        "cds_bp": 900,
        "transcript_span_bp": 1200,
        "gmb_score": 3.0,
        "diamond_hit": "P12345.1",
        "diamond_pident": 80.0,
        "diamond_qcov": 90.0,
        "diamond_scov": 85.0,
        "diamond_bitscore": 300.0,
        "diamond_evalue": 1e-100,
        "psauron_score": 0.95,
        "protein_length": 300,
        "orf_label": "ORF:300aa ATG|STOP",
        "is_partial_5": False,
        "is_partial_3": False,
        "internal_stop_count": 0,
        "protein_coding_score": 0.8,
    }
    base.update(overrides)
    return base


@pytest.fixture
def cfg():
    return load_config().canonical_selection


class TestSingleTranscriptGenes:
    def test_single_isoform_is_canonical_high_confidence(self, cfg):
        records = [_rec("t1")]
        result = select_canonical_for_gene("G1", records, cfg, "Helixer")
        assert result["canonical_transcript_id"] == "t1"
        assert result["runner_up_transcript_id"] is None
        assert result["confidence"] == "high_confidence"
        assert result["selection_reason"] == "ONLY_ISOFORM"

    def test_single_isoform_no_protein_support_is_low_confidence(self, cfg):
        records = [_rec("t1", diamond_hit=None, psauron_score=None)]
        result = select_canonical_for_gene("G1", records, cfg, "Helixer")
        assert result["confidence"] == "low_confidence"
        assert result["selection_reason"] == "LOW_CONFIDENCE_NO_PROTEIN_SUPPORT"


class TestPsauronDifferentiation:
    def test_higher_psauron_wins_when_otherwise_tied(self, cfg):
        records = [
            _rec("t1", psauron_score=0.99, diamond_qcov=90.0, diamond_scov=90.0),
            _rec("t2", psauron_score=0.50, diamond_qcov=90.0, diamond_scov=90.0),
        ]
        result = select_canonical_for_gene("G1", records, cfg, "Helixer")
        assert result["canonical_transcript_id"] == "t1"
        assert result["selection_reason"] == "BEST_PSAURON"


class TestDiamondCoverageDifferentiation:
    def test_higher_coverage_wins_when_psauron_tied(self, cfg):
        records = [
            _rec("t1", psauron_score=0.9, diamond_qcov=95.0, diamond_scov=95.0),
            _rec("t2", psauron_score=0.9, diamond_qcov=40.0, diamond_scov=40.0),
        ]
        result = select_canonical_for_gene("G1", records, cfg, "Helixer")
        assert result["canonical_transcript_id"] == "t1"
        assert result["selection_reason"] == "BEST_DIAMOND_COVERAGE"

    def test_bitscore_alone_does_not_decide(self, cfg):
        # t2 has a much larger raw bitscore (long-protein effect) but lower
        # coverage/identity -- must not win on bitscore alone, since the
        # canonical scorer never uses raw bitscore as a component.
        records = [
            _rec(
                "t1",
                diamond_bitscore=200.0,
                diamond_pident=90.0,
                diamond_qcov=95.0,
                diamond_scov=95.0,
            ),
            _rec(
                "t2",
                diamond_bitscore=900.0,
                diamond_pident=30.0,
                diamond_qcov=20.0,
                diamond_scov=20.0,
            ),
        ]
        result = select_canonical_for_gene("G1", records, cfg, "Helixer")
        assert result["canonical_transcript_id"] == "t1"


class TestProteinVsTranscriptEvidence:
    def test_strong_protein_score_can_outweigh_more_sources(self, cfg):
        records = [
            _rec(
                "t_strong_protein",
                psauron_score=0.99,
                diamond_pident=95.0,
                diamond_qcov=98.0,
                diamond_scov=98.0,
                evidence_sources="Tiberius",
                gmb_score=1.0,
            ),
            _rec(
                "t_strong_evidence",
                psauron_score=0.3,
                diamond_hit=None,
                diamond_pident=None,
                diamond_qcov=None,
                diamond_scov=None,
                evidence_sources="Tiberius,Scallop,StringTie,GenBlast",
                gmb_score=8.0,
            ),
        ]
        result = select_canonical_for_gene("G1", records, cfg, "Tiberius")
        assert result["canonical_transcript_id"] == "t_strong_protein"


class TestIncompleteOrfPenalty:
    def test_complete_orf_beats_partial_despite_lower_raw_components(self, cfg):
        records = [
            _rec("t_complete", psauron_score=0.7, diamond_qcov=60.0, diamond_scov=60.0),
            _rec(
                "t_partial",
                psauron_score=0.9,
                diamond_qcov=90.0,
                diamond_scov=90.0,
                is_partial_5=True,
                orf_label="ORF:150aa partial5|STOP",
            ),
        ]
        result = select_canonical_for_gene("G1", records, cfg, "Helixer")
        assert result["canonical_transcript_id"] == "t_complete"
        assert result["selection_reason"] == "ONLY_COMPLETE_ORF"

    def test_internal_stop_penalised(self, cfg):
        clean = score_transcript(_rec("t1"), cfg, "Helixer", (0.0, 0.0))
        with_stop = score_transcript(_rec("t1", internal_stop_count=2), cfg, "Helixer", (0.0, 0.0))
        assert with_stop["structure_subtotal"] < clean["structure_subtotal"]
        assert with_stop["has_internal_stop"] is True

    def test_no_internal_stop_beats_internal_stop_when_orf_complete_and_tied(self, cfg):
        # Both candidates have has_complete_orf True and identical
        # protein-validation/evidence-class/gmb_score/cds_bp -- only
        # has_internal_stop differs. Rank-key tier 1 (structural/ORF
        # validity) now covers this via has_internal_stop, ahead of the
        # protein-plausibility tier.
        records = [
            _rec("t_clean"),
            _rec("t_stopped", internal_stop_count=1),
        ]
        result = select_canonical_for_gene("G1", records, cfg, "Helixer")
        assert result["canonical_transcript_id"] == "t_clean"
        assert result["selection_reason"] == "ONLY_NO_INTERNAL_STOP"


class TestBalancedDiamondCoverage:
    def test_fragment_match_not_flattered_by_averaging(self, cfg):
        # 100% query coverage / 10% target coverage is a fragment match --
        # balanced (min-based) coverage must score it near the WEAK side
        # (0.10), not the mean (0.55) an average would give.
        fragment = score_transcript(
            _rec("t1", diamond_qcov=100.0, diamond_scov=10.0), cfg, "Helixer", (0.0, 0.0)
        )
        assert fragment["diamond_balanced_coverage_component"] == 0.10

    def test_balanced_coverage_drives_protein_validation_subtotal_not_separate_sides(self, cfg):
        # Two candidates with identical (high) query coverage but very
        # different target coverage must be told apart by the balanced
        # (weaker-side) measure.
        strong = score_transcript(
            _rec("t1", diamond_qcov=95.0, diamond_scov=95.0), cfg, "Helixer", (0.0, 0.0)
        )
        fragment = score_transcript(
            _rec("t1", diamond_qcov=95.0, diamond_scov=10.0), cfg, "Helixer", (0.0, 0.0)
        )
        assert strong["protein_validation_subtotal"] > fragment["protein_validation_subtotal"]


class TestEvidenceClassFieldsInOutput:
    def test_score_transcript_reports_class_and_named_source_fields(self, cfg):
        rec = _rec("t1", evidence_sources="Scallop,StringTie")
        scored = score_transcript(rec, cfg, "Helixer", (0.0, 0.0))
        assert scored["n_evidence_classes"] >= 1
        assert "short_read_transcriptomic" in scored["evidence_classes"]
        assert scored["n_independent_sources"] == 2
        assert scored["named_evidence_sources"] == "Scallop,StringTie"

    def test_run_end_to_end_writes_evidence_class_columns(self, tmp_path):
        gff3 = tmp_path / "consensus.gff3"
        gff3.write_text(
            "##gff-version 3\n"
            "1\tGMB\tgene\t100\t500\t.\t+\t.\tID=G1\n"
            "1\tGMB\tmRNA\t100\t500\t.\t+\t.\tID=G1.t1;Parent=G1;Evidence=Scallop\n"
            "1\tGMB\texon\t100\t500\t.\t+\t.\tID=G1.t1.exon1;Parent=G1.t1\n"
            "1\tGMB\tCDS\t100\t500\t.\t+\t.\tID=G1.t1.cds1;Parent=G1.t1\n"
        )
        ev = tmp_path / "evidence_attribution.tsv"
        pd.DataFrame(
            [
                {
                    "gene_id": "G1",
                    "transcript_id": "G1.t1",
                    "evidence_sources": "Scallop,StringTie",
                    "exon_count": 1,
                    "cds_bp": 400,
                    "transcript_span_bp": 400,
                    "gmb_score": 2.0,
                }
            ]
        ).to_csv(ev, sep="\t", index=False)

        out_dir = tmp_path / "out"
        config = load_config()
        run_canonical_selection(str(gff3), str(ev), str(out_dir), config)

        canon = pd.read_csv(out_dir / "canonical_transcripts.tsv", sep="\t")
        ranking = pd.read_csv(out_dir / "transcript_ranking.tsv", sep="\t")
        for df in (canon, ranking):
            for col in ("named_sources", "n_named_sources", "evidence_classes", "n_evidence_classes"):
                assert col in df.columns
        assert ranking.loc[0, "n_named_sources"] == 2
        assert ranking.loc[0, "evidence_classes"] == "short_read_transcriptomic"


class TestMissingOptionalFields:
    def test_absent_diamond_hit_handled_safely(self, cfg):
        rec = _rec(
            "t1", diamond_hit=None, diamond_pident=None, diamond_qcov=None, diamond_scov=None
        )
        scored = score_transcript(rec, cfg, "Helixer", (0.0, 0.0))
        assert scored["diamond_hit_component"] == 0.0
        assert scored["diamond_identity_component"] == 0.0

    def test_absent_psauron_score_handled_safely(self, cfg):
        rec = _rec("t1", psauron_score=None)
        scored = score_transcript(rec, cfg, "Helixer", (0.0, 0.0))
        assert scored["psauron_component"] == 0.0
        # a real DIAMOND hit still counts as protein support even without psauron
        assert scored["has_any_protein_support"] is True

    def test_no_protein_validation_at_all(self, cfg):
        rec = _rec(
            "t1",
            diamond_hit=None,
            diamond_pident=None,
            diamond_qcov=None,
            diamond_scov=None,
            psauron_score=None,
        )
        scored = score_transcript(rec, cfg, "Helixer", (0.0, 0.0))
        assert scored["has_any_protein_support"] is False
        assert scored["protein_validation_subtotal"] == 0.0

    def test_missing_gmb_score_defaults_to_neutral_component(self, cfg):
        rec = _rec("t1", gmb_score=None)
        scored = score_transcript(rec, cfg, "Helixer", (0.0, 10.0))
        assert scored["gmb_score_component"] == 0.5

    def test_missing_evidence_sources(self, cfg):
        rec = _rec("t1", evidence_sources="")
        scored = score_transcript(rec, cfg, "Helixer", (0.0, 0.0))
        assert scored["n_independent_sources"] == 0
        assert scored["transcriptomic_support"] is False


class TestDeterministicTieBreak:
    def test_full_tie_falls_back_to_lexical_transcript_id(self, cfg):
        records = [_rec("t_zzz"), _rec("t_aaa")]
        result = select_canonical_for_gene("G1", records, cfg, "Helixer")
        assert result["canonical_transcript_id"] == "t_aaa"
        assert result["selection_reason"] == "LEXICAL_TIEBREAK"

    def test_repeated_calls_give_same_winner(self, cfg):
        records = [_rec("t2", gmb_score=5.0), _rec("t1", gmb_score=5.0), _rec("t3", gmb_score=5.0)]
        r1 = select_canonical_for_gene("G1", records, cfg, "Helixer")
        r2 = select_canonical_for_gene("G1", list(reversed(records)), cfg, "Helixer")
        assert r1["canonical_transcript_id"] == r2["canonical_transcript_id"]

    def test_broader_evidence_breaks_tie_before_gmb_score(self, cfg):
        # t_multi: Scallop+StringTie collapse to ONE short_read_transcriptomic
        # class; Tiberius is not the configured backbone here ("Helixer"),
        # so it maps to the 'other' evidence class rather than being
        # silently dropped -- giving t_multi 3 classes (short_read, other,
        # protein_validation) vs t_single's 2 (short_read,
        # protein_validation). Evidence-CLASS breadth (mandatory change --
        # see the score-provenance audit) now decides this tie before the
        # raw named-source count would have.
        records = [
            _rec("t_multi", evidence_sources="Scallop,StringTie,Tiberius", gmb_score=3.0),
            _rec("t_single", evidence_sources="Scallop", gmb_score=4.0),
        ]
        result = select_canonical_for_gene("G1", records, cfg, "Helixer")
        assert result["canonical_transcript_id"] == "t_multi"
        assert result["selection_reason"] == "BEST_EVIDENCE_CLASS_BREADTH"

    def test_named_source_count_still_breaks_ties_evidence_classes_leave_open(self, cfg):
        # Two named sources that map to the SAME evidence class (Scallop +
        # StringTie, both short_read_transcriptomic) tie on class breadth,
        # so the older, more granular named-source count still decides --
        # it was demoted, not discarded (see _rank_key's tier 4).
        records = [
            _rec("t_two_sources", evidence_sources="Scallop,StringTie", gmb_score=3.0),
            _rec("t_one_source", evidence_sources="Scallop", gmb_score=4.0),
        ]
        result = select_canonical_for_gene("G1", records, cfg, "Helixer")
        assert result["canonical_transcript_id"] == "t_two_sources"
        assert result["selection_reason"] == "BEST_MULTI_SOURCE_SUPPORT"

    def test_cds_length_tiebreak_after_evidence_and_gmb_score_tied(self, cfg):
        records = [
            _rec("t_short", cds_bp=500),
            _rec("t_long", cds_bp=900),
        ]
        result = select_canonical_for_gene("G1", records, cfg, "Helixer")
        assert result["canonical_transcript_id"] == "t_long"
        assert result["selection_reason"] == "LONGEST_VALID_CDS_TIEBREAK"


class TestDomainEvidence:
    def test_domain_evidence_present_but_disabled_has_no_effect(self, cfg):
        cfg.domains.enabled = False
        cfg.domains.domain_support_weight = 5.0  # would dominate if it were active
        domain_metrics = {"t1": {"n_significant_domains": 3, "protein_coverage_fraction": 0.9}}
        scored = score_transcript(_rec("t1"), cfg, "Helixer", (0.0, 0.0), domain_metrics)
        assert scored["domain_subtotal"] == 0.0

    def test_domain_evidence_enabled_contributes_to_score(self, cfg):
        cfg.domains.enabled = True
        cfg.domains.domain_support_weight = 1.0
        domain_metrics = {"t1": {"n_significant_domains": 3, "protein_coverage_fraction": 0.9}}
        scored = score_transcript(_rec("t1"), cfg, "Helixer", (0.0, 0.0), domain_metrics)
        assert scored["domain_subtotal"] > 0.0

    def test_no_domain_hits_is_not_penalised(self, cfg):
        cfg.domains.enabled = True
        cfg.domains.domain_support_weight = 1.0
        cfg.domains.fragmented_domain_penalty = 1.0
        # No entry at all for this transcript -- must not be treated as a
        # penalty, only as "no support" (subtotal 0, not negative).
        scored = score_transcript(_rec("t1"), cfg, "Helixer", (0.0, 0.0), domain_metrics={})
        assert scored["domain_subtotal"] == 0.0


class TestFileIOAndGff3Annotation:
    def test_load_transcript_records_from_real_files(self, tmp_path):
        gff3 = tmp_path / "consensus.gff3"
        gff3.write_text(
            "##gff-version 3\n"
            "1\tGMB\tgene\t100\t500\t.\t+\t.\tID=G1\n"
            "1\tGMB\tmRNA\t100\t500\t.\t+\t.\tID=G1.t1;Parent=G1;Evidence=Scallop\n"
            "1\tGMB\texon\t100\t500\t.\t+\t.\tID=G1.t1.exon1;Parent=G1.t1\n"
            "1\tGMB\tCDS\t100\t500\t.\t+\t.\tID=G1.t1.cds1;Parent=G1.t1\n"
        )
        ev = tmp_path / "evidence_attribution.tsv"
        pd.DataFrame(
            [
                {
                    "gene_id": "G1",
                    "transcript_id": "G1.t1",
                    "evidence_sources": "Scallop",
                    "exon_count": 1,
                    "cds_bp": 400,
                    "transcript_span_bp": 400,
                    "gmb_score": 2.0,
                }
            ]
        ).to_csv(ev, sep="\t", index=False)

        genes = load_transcript_records(str(gff3), str(ev), protein_validation_tsv=None)
        assert "G1" in genes
        assert genes["G1"][0]["transcript_id"] == "G1.t1"
        assert genes["G1"][0]["diamond_hit"] is None  # no protein_validation.tsv supplied

    def test_alternative_isoforms_unchanged_in_ranking(self, tmp_path, cfg):
        gff3 = tmp_path / "consensus.gff3"
        gff3.write_text(
            "##gff-version 3\n"
            "1\tGMB\tgene\t100\t900\t.\t+\t.\tID=G1\n"
            "1\tGMB\tmRNA\t100\t500\t.\t+\t.\tID=G1.t1;Parent=G1;Evidence=Scallop\n"
            "1\tGMB\texon\t100\t500\t.\t+\t.\tID=G1.t1.exon1;Parent=G1.t1\n"
            "1\tGMB\tCDS\t100\t500\t.\t+\t.\tID=G1.t1.cds1;Parent=G1.t1\n"
            "1\tGMB\tmRNA\t100\t900\t.\t+\t.\tID=G1.t2;Parent=G1;Evidence=StringTie\n"
            "1\tGMB\texon\t100\t900\t.\t+\t.\tID=G1.t2.exon1;Parent=G1.t2\n"
            "1\tGMB\tCDS\t100\t900\t.\t+\t.\tID=G1.t2.cds1;Parent=G1.t2\n"
        )
        ev = tmp_path / "evidence_attribution.tsv"
        pd.DataFrame(
            [
                {
                    "gene_id": "G1",
                    "transcript_id": "G1.t1",
                    "evidence_sources": "Scallop",
                    "exon_count": 1,
                    "cds_bp": 400,
                    "transcript_span_bp": 400,
                    "gmb_score": 5.0,
                },
                {
                    "gene_id": "G1",
                    "transcript_id": "G1.t2",
                    "evidence_sources": "StringTie",
                    "exon_count": 1,
                    "cds_bp": 800,
                    "transcript_span_bp": 800,
                    "gmb_score": 1.0,
                },
            ]
        ).to_csv(ev, sep="\t", index=False)

        out_dir = tmp_path / "out"
        config = load_config()
        summary = run_canonical_selection(str(gff3), str(ev), str(out_dir), config)

        assert summary["n_genes"] == 1
        ranking = pd.read_csv(out_dir / "transcript_ranking.tsv", sep="\t")
        # Both isoforms must still be present in the reporting output --
        # canonical selection is additive/annotation-only, never a filter.
        assert set(ranking["transcript_id"]) == {"G1.t1", "G1.t2"}
        # And the original GFF3 content must be completely untouched.
        original = gff3.read_text()
        assert "G1.t1" in original and "G1.t2" in original

    def test_annotate_gff3_writes_copy_without_touching_original(self, tmp_path):
        gff3 = tmp_path / "consensus.gff3"
        original_text = (
            "##gff-version 3\n"
            "1\tGMB\tgene\t100\t500\t.\t+\t.\tID=G1\n"
            "1\tGMB\tmRNA\t100\t500\t.\t+\t.\tID=G1.t1;Parent=G1;Evidence=Scallop\n"
            "1\tGMB\texon\t100\t500\t.\t+\t.\tID=G1.t1.exon1;Parent=G1.t1\n"
        )
        gff3.write_text(original_text)

        out_path = tmp_path / "annotated.gff3"
        _write_annotated_gff3(str(gff3), {"G1.t1"}, str(out_path))

        # Original untouched
        assert gff3.read_text() == original_text
        # Copy has the canonical flag
        annotated = out_path.read_text()
        assert "canonical=1" in annotated
        mrna_line = [l for l in annotated.splitlines() if "\tmRNA\t" in l][0]
        assert mrna_line.endswith("canonical=1")


class TestRunEndToEndSummaryFields:
    def test_summary_reports_expected_keys(self, tmp_path):
        gff3 = tmp_path / "consensus.gff3"
        gff3.write_text(
            "##gff-version 3\n"
            "1\tGMB\tgene\t100\t500\t.\t+\t.\tID=G1\n"
            "1\tGMB\tmRNA\t100\t500\t.\t+\t.\tID=G1.t1;Parent=G1;Evidence=Scallop\n"
            "1\tGMB\texon\t100\t500\t.\t+\t.\tID=G1.t1.exon1;Parent=G1.t1\n"
        )
        ev = tmp_path / "evidence_attribution.tsv"
        pd.DataFrame(
            [
                {
                    "gene_id": "G1",
                    "transcript_id": "G1.t1",
                    "evidence_sources": "Scallop",
                    "exon_count": 1,
                    "cds_bp": 400,
                    "transcript_span_bp": 400,
                    "gmb_score": 2.0,
                }
            ]
        ).to_csv(ev, sep="\t", index=False)

        out_dir = tmp_path / "out"
        config = load_config()
        summary = run_canonical_selection(str(gff3), str(ev), str(out_dir), config)

        for key in [
            "n_genes",
            "n_multi_isoform_genes",
            "n_single_isoform_genes",
            "n_canonical_differs_from_highest_gmb_score",
            "selection_reason_counts",
            "confidence_counts",
            "config",
        ]:
            assert key in summary

        with open(out_dir / "canonical_selection_summary.json") as fh:
            written = json.load(fh)
        assert written["n_genes"] == 1
