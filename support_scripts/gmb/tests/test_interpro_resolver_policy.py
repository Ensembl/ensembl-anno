#!/usr/bin/env python3
"""Tests for the conservative InterPro replacement policy
(gmb.pipeline.interpro_resolver.decide_canonical and its building blocks).

Covers every reason code, the fusion pre-check, and each configured
safeguard -- both when it blocks a replacement and when it deliberately
does not (the two safety-critical reason codes bypass the
protein-validation override-gap safeguard specifically; see
check_safeguards' docstring).
"""

import pytest

from gmb.pipeline.config import InterProResolverConfig
from gmb.pipeline.interpro_resolver import (
    REASON_ANTIFAM_AVOIDED,
    REASON_ARCHITECTURE_RESTORED,
    REASON_C_TERMINAL_TRUNCATION_RESOLVED,
    REASON_COMPLETE_DOMAIN_RESTORED,
    REASON_CONFLICTING,
    REASON_EQUIVALENT,
    REASON_FUSION_AVOIDED,
    REASON_N_TERMINAL_TRUNCATION_RESOLVED,
    REASON_SUPPORTS_INITIAL,
    REASON_UNINFORMATIVE,
    check_safeguards,
    classify_reason_code,
    compare_candidates,
    decide_canonical,
    detect_fusion,
    summarise_architecture,
)


def _match(acc="PF00001", ipr=None, start=1, end=100, rep=True, lib="Pfam", ftype="Domain"):
    from gmb.pipeline.interpro_resolver import classify_library as _cl

    return {
        "signature_accession": acc,
        "interpro_accession": ipr,
        "member_database": lib,
        "signature_type": ftype,
        "start": start,
        "end": end,
        "representative": rep,
        "evidence_class": _cl(lib, ftype),
    }


def _summary(matches, length=200):
    return summarise_architecture(matches, protein_length=length)


def _row(tid, **overrides):
    base = {
        "transcript_id": tid,
        "protein_length": 200,
        "has_complete_orf": True,
        "has_internal_stop": False,
        "psauron_score": None,
        "diamond_balanced_coverage": None,
    }
    base.update(overrides)
    return base


@pytest.fixture
def cfg():
    return InterProResolverConfig()


class TestReasonCodeClassification:
    def test_supports_initial(self, cfg):
        current = _summary([_match(ipr="IPR1", start=1, end=190)])
        runner_up = _summary([_match(ipr="IPR1", start=1, end=20)])
        comparison = compare_candidates(current, runner_up)
        code = classify_reason_code(current, runner_up, comparison, is_fusion=False)
        assert code == REASON_SUPPORTS_INITIAL

    def test_equivalent_when_architectures_match(self):
        current = _summary([_match(ipr="IPR1", start=1, end=100)])
        runner_up = _summary([_match(ipr="IPR1", start=1, end=100)])
        comparison = compare_candidates(current, runner_up)
        code = classify_reason_code(current, runner_up, comparison, is_fusion=False)
        assert code == REASON_EQUIVALENT

    def test_uninformative_when_neither_has_domains(self):
        current = _summary([])
        runner_up = _summary([])
        comparison = compare_candidates(current, runner_up)
        code = classify_reason_code(current, runner_up, comparison, is_fusion=False)
        assert code == REASON_UNINFORMATIVE

    def test_conflicting_when_each_side_has_exclusive_entries(self):
        current = _summary([_match(ipr="IPR1", start=1, end=50)])
        runner_up = _summary([_match(ipr="IPR2", start=1, end=50)])
        comparison = compare_candidates(current, runner_up)
        code = classify_reason_code(current, runner_up, comparison, is_fusion=False)
        assert code == REASON_CONFLICTING

    def test_antifam_avoided(self):
        current = _summary(
            [_match(acc="ANF00001", lib="AntiFam", ftype="Domain", start=1, end=50)]
        )
        runner_up = _summary([_match(ipr="IPR1", start=1, end=100)])
        comparison = compare_candidates(current, runner_up)
        code = classify_reason_code(current, runner_up, comparison, is_fusion=False)
        assert code == REASON_ANTIFAM_AVOIDED

    def test_complete_domain_restored_when_current_has_no_domains(self):
        current = _summary([])  # no strong evidence at all
        runner_up = _summary([_match(ipr="IPR1", start=1, end=100)])
        comparison = compare_candidates(current, runner_up)
        code = classify_reason_code(current, runner_up, comparison, is_fusion=False)
        assert code == REASON_COMPLETE_DOMAIN_RESTORED

    def test_n_terminal_truncation_resolved(self):
        # Runner-up's domain starts earlier (residue 1) than current's
        # (residue 80), same end -- current is missing the N-terminal part.
        current = _summary([_match(ipr="IPR1", start=80, end=190)])
        runner_up = _summary([_match(ipr="IPR1", start=1, end=190)])
        comparison = compare_candidates(current, runner_up)
        code = classify_reason_code(current, runner_up, comparison, is_fusion=False)
        assert code == REASON_N_TERMINAL_TRUNCATION_RESOLVED

    def test_c_terminal_truncation_resolved(self):
        # Runner-up's domain extends further (to residue 190) than
        # current's (to residue 20), same start -- current is missing the
        # C-terminal part.
        current = _summary([_match(ipr="IPR1", start=1, end=20)])
        runner_up = _summary([_match(ipr="IPR1", start=1, end=190)])
        comparison = compare_candidates(current, runner_up)
        code = classify_reason_code(current, runner_up, comparison, is_fusion=False)
        assert code == REASON_C_TERMINAL_TRUNCATION_RESOLVED

    def test_architecture_restored_generic_fallback(self):
        # Runner-up's span both starts earlier AND ends later (envelops
        # current's on both flanks) -- not a clean single-terminus
        # localisation, so the generic fallback applies.
        current = _summary([_match(ipr="IPR1", start=50, end=100)])
        runner_up = _summary([_match(ipr="IPR1", start=1, end=190)])
        comparison = compare_candidates(current, runner_up)
        code = classify_reason_code(current, runner_up, comparison, is_fusion=False)
        assert code == REASON_ARCHITECTURE_RESTORED

    def test_fusion_avoided_forced_regardless_of_generic_verdict(self):
        # Even if the generic entry-count comparison would otherwise favour
        # current (more integrated entries), a detected fusion always
        # reports FUSION_AVOIDED.
        current = _summary(
            [_match(ipr="IPR1", start=1, end=50), _match(ipr="IPR2", start=400, end=450)],
            length=600,
        )
        runner_up = _summary([_match(ipr="IPR1", start=1, end=50)], length=100)
        comparison = compare_candidates(current, runner_up)
        code = classify_reason_code(current, runner_up, comparison, is_fusion=True)
        assert code == REASON_FUSION_AVOIDED


class TestFusionDetection:
    def test_detects_extra_nonoverlapping_region_on_much_longer_candidate(self, cfg):
        # Current (600aa) carries runner-up's own domain (1-50) PLUS an
        # unrelated extra domain far outside it (400-450) -- the "bolted
        # on" signature of a fusion/read-through model.
        current = _summary(
            [_match(ipr="IPR1", start=1, end=50), _match(ipr="IPR2", start=400, end=450)],
            length=600,
        )
        runner_up = _summary([_match(ipr="IPR1", start=1, end=50)], length=100)
        assert (
            detect_fusion(current, runner_up, 600, 100, cfg.length_ratio_flag) is True
        )

    def test_no_fusion_when_runner_up_has_no_domains(self, cfg):
        current = _summary([_match(ipr="IPR1", start=1, end=50)], length=600)
        runner_up = _summary([], length=100)
        assert detect_fusion(current, runner_up, 600, 100, cfg.length_ratio_flag) is False

    def test_no_fusion_when_length_ratio_below_threshold(self, cfg):
        current = _summary(
            [_match(ipr="IPR1", start=1, end=50), _match(ipr="IPR2", start=60, end=90)],
            length=120,
        )
        runner_up = _summary([_match(ipr="IPR1", start=1, end=50)], length=100)
        assert detect_fusion(current, runner_up, 120, 100, cfg.length_ratio_flag) is False

    def test_no_fusion_when_longer_isoform_just_covers_more_of_same_domain(self, cfg):
        # Current is longer and has a WIDER span of the SAME region, not an
        # extra disjoint one -- this is a truncation-resolution case, not a
        # fusion, even at a large length ratio.
        current = _summary([_match(ipr="IPR1", start=1, end=190)], length=600)
        runner_up = _summary([_match(ipr="IPR1", start=1, end=50)], length=100)
        assert detect_fusion(current, runner_up, 600, 100, cfg.length_ratio_flag) is False


class TestSafeguards:
    def test_blocks_when_runner_up_has_worse_structural_class(self, cfg):
        current = _row("t_current", has_complete_orf=True)
        runner_up = _row("t_runner", has_complete_orf=False)
        result = check_safeguards(current, runner_up, REASON_ARCHITECTURE_RESTORED, cfg)
        assert result["passed"] is False
        assert result["checks"]["no_worse_structural_class"] is False

    def test_blocks_when_runner_up_introduces_new_internal_stop(self, cfg):
        current = _row("t_current", has_internal_stop=False)
        runner_up = _row("t_runner", has_internal_stop=True)
        result = check_safeguards(current, runner_up, REASON_ARCHITECTURE_RESTORED, cfg)
        assert result["passed"] is False
        assert result["checks"]["no_new_internal_stops"] is False

    def test_blocks_when_current_protein_validation_much_stronger(self, cfg):
        current = _row("t_current", psauron_score=0.95, diamond_balanced_coverage=0.95)
        runner_up = _row("t_runner", psauron_score=0.40, diamond_balanced_coverage=0.40)
        result = check_safeguards(current, runner_up, REASON_ARCHITECTURE_RESTORED, cfg)
        assert result["passed"] is False
        assert result["checks"]["no_much_stronger_protein_validation_overridden"] is False

    def test_protein_validation_gap_does_not_block_antifam_avoidance(self, cfg):
        # A spurious ORF can still score well on DIAMOND/psauron -- that is
        # exactly why AntiFam avoidance is a safety-critical override, not
        # subject to the override-gap safeguard.
        current = _row("t_current", psauron_score=0.95, diamond_balanced_coverage=0.95)
        runner_up = _row("t_runner", psauron_score=0.40, diamond_balanced_coverage=0.40)
        result = check_safeguards(current, runner_up, REASON_ANTIFAM_AVOIDED, cfg)
        assert result["passed"] is True
        assert "no_much_stronger_protein_validation_overridden" not in result["checks"]

    def test_protein_validation_gap_does_not_block_fusion_avoidance(self, cfg):
        current = _row("t_current", psauron_score=0.95, diamond_balanced_coverage=0.95)
        runner_up = _row("t_runner", psauron_score=0.40, diamond_balanced_coverage=0.40)
        result = check_safeguards(current, runner_up, REASON_FUSION_AVOIDED, cfg)
        assert result["passed"] is True

    def test_small_gap_within_threshold_passes(self, cfg):
        current = _row("t_current", psauron_score=0.80, diamond_balanced_coverage=0.80)
        runner_up = _row("t_runner", psauron_score=0.75, diamond_balanced_coverage=0.75)
        result = check_safeguards(current, runner_up, REASON_ARCHITECTURE_RESTORED, cfg)
        assert result["passed"] is True


class TestDecideCanonical:
    def test_replacement_applied_when_reason_valid_and_safeguards_pass(self, cfg):
        current_row = _row("G1.t1")
        runner_up_row = _row("G1.t2")
        cur_summary = _summary([_match(ipr="IPR1", start=1, end=20)])
        run_summary = _summary([_match(ipr="IPR1", start=1, end=190)])
        decision = decide_canonical("G1", current_row, runner_up_row, cur_summary, run_summary, cfg)
        assert decision["replaced"] is True
        assert decision["final_canonical_transcript_id"] == "G1.t2"
        assert decision["reason_code"] == REASON_C_TERMINAL_TRUNCATION_RESOLVED

    def test_no_replacement_when_apply_replacements_disabled(self, cfg):
        cfg.apply_replacements = False
        current_row = _row("G1.t1")
        runner_up_row = _row("G1.t2")
        cur_summary = _summary([_match(ipr="IPR1", start=1, end=20)])
        run_summary = _summary([_match(ipr="IPR1", start=1, end=190)])
        decision = decide_canonical("G1", current_row, runner_up_row, cur_summary, run_summary, cfg)
        assert decision["replaced"] is False
        assert decision["final_canonical_transcript_id"] == "G1.t1"
        # The recommendation is still recorded even though it was not applied.
        assert decision["interpro_recommended_transcript_id"] == "G1.t2"
        assert decision["reason_code"] == REASON_C_TERMINAL_TRUNCATION_RESOLVED

    def test_no_replacement_when_safeguard_fails(self, cfg):
        current_row = _row("G1.t1", has_complete_orf=True)
        runner_up_row = _row("G1.t2", has_complete_orf=False)  # worse structural class
        cur_summary = _summary([_match(ipr="IPR1", start=1, end=20)])
        run_summary = _summary([_match(ipr="IPR1", start=1, end=190)])
        decision = decide_canonical("G1", current_row, runner_up_row, cur_summary, run_summary, cfg)
        assert decision["replaced"] is False
        assert decision["final_canonical_transcript_id"] == "G1.t1"
        assert decision["safeguards_passed"] is False

    def test_no_replacement_candidate_never_evaluates_safeguards(self, cfg):
        # Equivalent architecture -> not even a replacement candidate, so
        # safeguards are trivially "passed" (nothing to block).
        current_row = _row("G1.t1")
        runner_up_row = _row("G1.t2")
        cur_summary = _summary([_match(ipr="IPR1", start=1, end=100)])
        run_summary = _summary([_match(ipr="IPR1", start=1, end=100)])
        decision = decide_canonical("G1", current_row, runner_up_row, cur_summary, run_summary, cfg)
        assert decision["reason_code"] == REASON_EQUIVALENT
        assert decision["replaced"] is False
        assert decision["safeguard_checks"] == {}

    def test_initial_and_recommended_always_reported_even_without_replacement(self, cfg):
        current_row = _row("G1.t1")
        runner_up_row = _row("G1.t2")
        cur_summary = _summary([])
        run_summary = _summary([])
        decision = decide_canonical("G1", current_row, runner_up_row, cur_summary, run_summary, cfg)
        assert decision["initial_canonical_transcript_id"] == "G1.t1"
        assert decision["interpro_recommended_transcript_id"] == "G1.t1"
        assert decision["final_canonical_transcript_id"] == "G1.t1"
        assert decision["reason_code"] == REASON_UNINFORMATIVE
