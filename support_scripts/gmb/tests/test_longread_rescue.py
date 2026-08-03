"""Tests for structure-aware short-read rescue (gmb.pipeline.longread.rescue).

Covers:
  - extract_junctions: ordered pairs, no snapping
  - junction_matches: explicit coordinate comparison, bin-boundary cases
  - m1_junction_complete: all-junctions policy
  - m2_full_chain: full ordered chain match
  - m3_dual_assembler: dual-assembler chain
  - s1_reciprocal_boundary: single-exon SR only, boundary tolerance, fragment detection
  - s2_dual_assembler_boundary: both assemblers S1-level
  - s3_independent_evidence: N independent SR transcripts
  - try_rescue_multi_exon / try_rescue_single_exon: dispatch
  - decide_keep_or_rescue: full decision
  - config: nested structure, validation, legacy key rejection
  - multiple short-read files: namespaced tids, provenance
"""

import os
import tempfile

import pytest

from gmb.pipeline.longread import (
    LongreadConsensusConfig,
    MultiExonRescueConfig,
    ShortreadRescueConfig,
    SingleExonRescueConfig,
    build_consensus_for_seqname,
    decide_keep_or_rescue,
    extract_junctions,
    junction_matches,
    load_longread_consensus_config,
    load_shortread_transcripts,
    m1_junction_complete,
    m2_full_chain,
    m3_dual_assembler,
    s1_reciprocal_boundary,
    s2_dual_assembler_boundary,
    s3_independent_evidence,
    try_rescue_multi_exon,
    try_rescue_single_exon,
    validate_config,
    write_rescue_attribution_tsv,
)

# ---------------------------------------------------------------------------
# Test helpers
# ---------------------------------------------------------------------------


def _tx(tid, strand, *exon_pairs, source="scallop", file_idx=0):
    """Build a transcript dict mimicking load_shortread_transcripts output."""
    ns_tid = f"{source}:{file_idx}:{tid}"
    return {
        "tid": ns_tid,
        "raw_tid": tid,
        "strand": strand,
        "exons": list(exon_pairs),
        "source": source,
        "source_file": "fake.gtf",
    }


def _sr_data(scallop=None, stringtie=None):
    s = scallop or []
    st = stringtie or []
    return {"scallop": s, "stringtie": st, "all": s + st}


def _rescue_cfg(
    enabled=True,
    me_enabled=True,
    me_policy="M1",
    me_tol=0,
    me_both=False,
    se_enabled=False,
    se_policy="S1",
    se_overlap=0.80,
    se_btol=0,
    se_req_single=True,
    se_both=False,
):
    return ShortreadRescueConfig(
        enabled=enabled,
        multi_exon=MultiExonRescueConfig(
            enabled=me_enabled,
            policy=me_policy,
            junction_tolerance_bp=me_tol,
            require_both_assemblers=me_both,
        ),
        single_exon=SingleExonRescueConfig(
            enabled=se_enabled,
            policy=se_policy,
            reciprocal_overlap_min=se_overlap,
            boundary_tolerance_bp=se_btol,
            require_single_exon_sr_model=se_req_single,
            require_both_assemblers=se_both,
        ),
    )


def _cfg(**overrides):
    rescue = overrides.pop("shortread_rescue", None)
    cfg = LongreadConsensusConfig()
    for k, v in overrides.items():
        setattr(cfg, k, v)
    if rescue is not None:
        cfg = LongreadConsensusConfig(
            **{
                f.name: getattr(cfg, f.name)
                for f in cfg.__dataclass_fields__.values()
                if f.name != "shortread_rescue"
            },
            shortread_rescue=rescue,
        )
    return cfg


# ---------------------------------------------------------------------------
# extract_junctions
# ---------------------------------------------------------------------------


class TestExtractJunctions:
    def test_single_exon_returns_empty(self):
        assert extract_junctions([(100, 200)]) == []

    def test_empty_returns_empty(self):
        assert extract_junctions([]) == []

    def test_two_exons(self):
        assert extract_junctions([(100, 200), (300, 400)]) == [(200, 300)]

    def test_three_exons(self):
        junctions = extract_junctions([(100, 200), (300, 400), (500, 600)])
        assert junctions == [(200, 300), (400, 500)]

    def test_unsorted_input_sorted_first(self):
        j_fwd = extract_junctions([(100, 200), (300, 400)])
        j_rev = extract_junctions([(300, 400), (100, 200)])
        assert j_fwd == j_rev

    def test_returns_ordered_list_not_frozenset(self):
        # Order must be preserved — not collapsed by frozenset
        junctions = extract_junctions([(100, 200), (300, 400), (500, 600)])
        assert isinstance(junctions, list)
        assert junctions[0] == (200, 300)
        assert junctions[1] == (400, 500)


# ---------------------------------------------------------------------------
# junction_matches — explicit coordinate comparison
# ---------------------------------------------------------------------------


class TestJunctionMatches:
    def test_exact_match(self):
        assert junction_matches(200, 300, 200, 300, tol_bp=10) is True

    def test_within_tolerance_donor(self):
        assert junction_matches(205, 300, 200, 300, tol_bp=10) is True

    def test_within_tolerance_acceptor(self):
        assert junction_matches(200, 308, 200, 300, tol_bp=10) is True

    def test_exactly_at_tolerance_both(self):
        assert junction_matches(210, 310, 200, 300, tol_bp=10) is True

    def test_one_outside_tolerance(self):
        assert junction_matches(211, 300, 200, 300, tol_bp=10) is False

    def test_both_outside_tolerance(self):
        assert junction_matches(215, 315, 200, 300, tol_bp=10) is False

    def test_zero_tolerance_exact_only(self):
        assert junction_matches(200, 300, 200, 300, tol_bp=0) is True
        assert junction_matches(201, 300, 200, 300, tol_bp=0) is False

    def test_bin_boundary_fix(self):
        """Old snapping-based approach fails for coordinates at bin boundaries.

        196 and 206 differ by 10 bp — exactly at tolerance.
        Old approach: 196 → round(19.6)*10 = 200; 206 → round(20.6)*10 = 210
                      These snap to different bins, so 196 vs 206 would NOT match!
        New explicit: abs(196 - 206) = 10 ≤ 10 → match.
        """
        assert junction_matches(196, 500, 206, 500, tol_bp=10) is True

    def test_near_bin_boundary_both_match(self):
        # 194 and 205 differ by 11 bp — just outside tolerance
        assert junction_matches(194, 500, 205, 500, tol_bp=10) is False
        # 195 and 205 differ by 10 bp — exactly at tolerance
        assert junction_matches(195, 500, 205, 500, tol_bp=10) is True


# ---------------------------------------------------------------------------
# m1_junction_complete
# ---------------------------------------------------------------------------


class TestM1JunctionComplete:
    def test_all_junctions_present(self):
        sr = [_tx("t1", "+", (100, 200), (300, 400), (500, 600))]
        kept, tids = m1_junction_complete([(100, 200), (300, 400), (500, 600)], "+", sr, tol_bp=0)
        assert kept is True
        assert "scallop:0:t1" in tids

    def test_one_junction_missing(self):
        sr = [_tx("t1", "+", (100, 200), (300, 400))]
        kept, tids = m1_junction_complete([(100, 200), (300, 400), (500, 600)], "+", sr, tol_bp=0)
        assert kept is False
        assert tids == []

    def test_junctions_split_across_two_sr_transcripts(self):
        # Each SR transcript covers one junction; together they cover all
        sr = [
            _tx("t1", "+", (100, 200), (300, 400)),
            _tx("t2", "+", (300, 400), (500, 600)),
        ]
        kept, tids = m1_junction_complete([(100, 200), (300, 400), (500, 600)], "+", sr, tol_bp=0)
        assert kept is True

    def test_strand_mismatch_excluded(self):
        sr = [_tx("t1", "-", (100, 200), (300, 400))]
        kept, _ = m1_junction_complete([(100, 200), (300, 400)], "+", sr, tol_bp=0)
        assert kept is False

    def test_single_exon_candidate_fails(self):
        sr = [_tx("t1", "+", (100, 400))]
        kept, _ = m1_junction_complete([(100, 400)], "+", sr, tol_bp=0)
        assert kept is False

    def test_explicit_tolerance_not_snapping(self):
        # Candidate junction (203, 298); SR junction (200, 300) — 3bp jitter
        sr = [_tx("t1", "+", (100, 200), (300, 400))]
        kept_no_tol, _ = m1_junction_complete([(100, 203), (298, 400)], "+", sr, tol_bp=0)
        kept_with_tol, _ = m1_junction_complete([(100, 203), (298, 400)], "+", sr, tol_bp=10)
        assert kept_no_tol is False
        assert kept_with_tol is True

    def test_bin_boundary_tolerance(self):
        # Candidate junction at (196, 500); SR junction at (206, 500)
        # Old snapping: 196→200, 206→210 → no match. New explicit: 10bp → match.
        sr = [_tx("t1", "+", (100, 206), (506, 600))]
        kept, _ = m1_junction_complete([(100, 196), (496, 600)], "+", sr, tol_bp=10)
        assert kept is True

    def test_plus_and_minus_strand(self):
        sr_plus = [_tx("t1", "+", (100, 200), (300, 400))]
        sr_minus = [_tx("t2", "-", (100, 200), (300, 400))]
        kept_p, _ = m1_junction_complete([(100, 200), (300, 400)], "+", sr_plus, tol_bp=0)
        kept_m, _ = m1_junction_complete([(100, 200), (300, 400)], "-", sr_minus, tol_bp=0)
        kept_cross, _ = m1_junction_complete([(100, 200), (300, 400)], "+", sr_minus, tol_bp=0)
        assert kept_p is True
        assert kept_m is True
        assert kept_cross is False


# ---------------------------------------------------------------------------
# m2_full_chain
# ---------------------------------------------------------------------------


class TestM2FullChain:
    def test_exact_match(self):
        sr = [_tx("t1", "+", (100, 200), (300, 400))]
        kept, tids = m2_full_chain([(100, 200), (300, 400)], "+", sr, tol_bp=0)
        assert kept is True
        assert "scallop:0:t1" in tids

    def test_partial_match_fails(self):
        # SR has superset of junctions — not an exact chain match
        sr = [_tx("t1", "+", (100, 200), (300, 400), (500, 600))]
        kept, _ = m2_full_chain([(100, 200), (300, 400)], "+", sr, tol_bp=0)
        assert kept is False

    def test_different_junction_fails(self):
        sr = [_tx("t1", "+", (200, 300), (400, 500))]
        kept, _ = m2_full_chain([(100, 200), (300, 400)], "+", sr, tol_bp=0)
        assert kept is False

    def test_single_exon_candidate_fails(self):
        sr = [_tx("t1", "+", (100, 400))]
        kept, _ = m2_full_chain([(100, 400)], "+", sr, tol_bp=0)
        assert kept is False

    def test_ordered_chain_not_frozenset(self):
        # Chains with junctions in different order should NOT match
        # Ordered: cand=(200,300),(400,500) vs SR=(400,500),(200,300) → different order
        # But sort() normalises the exons so order within a transcript is always ascending.
        # This tests that a different intron chain with same junctions in different order fails.
        # Exons (100,200),(400,500): junctions (200,400)
        # Exons (200,300),(500,600): junctions (300,500)
        # These are genuinely different intron chains and should not match.
        cand = [(100, 200), (400, 500)]
        sr = [_tx("t1", "+", (200, 300), (500, 600))]
        kept, _ = m2_full_chain(cand, "+", sr, tol_bp=0)
        assert kept is False

    def test_within_tolerance(self):
        # 2bp jitter on donor
        sr = [_tx("t1", "+", (100, 202), (300, 400))]
        kept_no_tol, _ = m2_full_chain([(100, 200), (300, 400)], "+", sr, tol_bp=0)
        kept_with_tol, _ = m2_full_chain([(100, 200), (300, 400)], "+", sr, tol_bp=5)
        assert kept_no_tol is False
        assert kept_with_tol is True


# ---------------------------------------------------------------------------
# m3_dual_assembler
# ---------------------------------------------------------------------------


class TestM3DualAssembler:
    def test_both_assemblers_support(self):
        scallop = [_tx("s1", "+", (100, 200), (300, 400), source="scallop")]
        stringtie = [_tx("st1", "+", (100, 200), (300, 400), source="stringtie")]
        kept, tids = m3_dual_assembler([(100, 200), (300, 400)], "+", scallop, stringtie, tol_bp=0)
        assert kept is True
        assert any("s1" in t for t in tids)
        assert any("st1" in t for t in tids)

    def test_only_scallop_fails(self):
        scallop = [_tx("s1", "+", (100, 200), (300, 400), source="scallop")]
        stringtie = [_tx("st1", "+", (200, 300), (400, 500), source="stringtie")]
        kept, _ = m3_dual_assembler([(100, 200), (300, 400)], "+", scallop, stringtie, tol_bp=0)
        assert kept is False

    def test_only_stringtie_fails(self):
        scallop = [_tx("s1", "+", (200, 300), (400, 500), source="scallop")]
        stringtie = [_tx("st1", "+", (100, 200), (300, 400), source="stringtie")]
        kept, _ = m3_dual_assembler([(100, 200), (300, 400)], "+", scallop, stringtie, tol_bp=0)
        assert kept is False

    def test_empty_assembler_fails(self):
        stringtie = [_tx("st1", "+", (100, 200), (300, 400), source="stringtie")]
        kept, _ = m3_dual_assembler([(100, 200), (300, 400)], "+", [], stringtie, tol_bp=0)
        assert kept is False


# ---------------------------------------------------------------------------
# s1_reciprocal_boundary (fixed: single-exon SR requirement, boundary tolerance)
# ---------------------------------------------------------------------------


class TestS1ReciprocalBoundary:
    def test_single_exon_sr_matches(self):
        sr = [_tx("t1", "+", (100, 900))]
        kept, tids, rej = s1_reciprocal_boundary(
            [(100, 900)], "+", sr, 0.80, boundary_tol_bp=0, require_single_exon_sr=True
        )
        assert kept is True
        assert rej is None

    def test_multi_exon_sr_rejected_when_require_single_exon(self):
        # SR transcript has two exons; its first exon (100, 900) perfectly covers
        # the candidate (100, 900) — good overlap but SR is multi-exon, so it
        # should be classified as a fragment, not rescued.
        sr = [_tx("t1", "+", (100, 900), (950, 1050))]
        kept, tids, rej = s1_reciprocal_boundary(
            [(100, 900)], "+", sr, 0.80, boundary_tol_bp=0, require_single_exon_sr=True
        )
        assert kept is False
        assert rej == "LIKELY_FRAGMENT_OF_MULTIEXON_TRANSCRIPT"

    def test_multi_exon_sr_allowed_when_not_required(self):
        # require_single_exon_sr=False — any exon of any SR transcript is OK
        sr = [_tx("t1", "+", (100, 900), (950, 1000))]
        kept, tids, rej = s1_reciprocal_boundary(
            [(100, 900)], "+", sr, 0.80, boundary_tol_bp=0, require_single_exon_sr=False
        )
        assert kept is True

    def test_low_overlap_rejected(self):
        # SR (100, 300); candidate (100, 600) → overlap 200/500 ≈ 0.40
        sr = [_tx("t1", "+", (100, 300))]
        kept, _, rej = s1_reciprocal_boundary(
            [(100, 600)], "+", sr, 0.80, boundary_tol_bp=0, require_single_exon_sr=True
        )
        assert kept is False
        assert rej is None  # not a fragment, just poor overlap

    def test_strand_mismatch_excluded(self):
        sr = [_tx("t1", "-", (100, 900))]
        kept, _, rej = s1_reciprocal_boundary(
            [(100, 900)], "+", sr, 0.80, boundary_tol_bp=0, require_single_exon_sr=True
        )
        assert kept is False

    def test_multi_exon_candidate_fails(self):
        sr = [_tx("t1", "+", (100, 900))]
        kept, _, _ = s1_reciprocal_boundary(
            [(100, 500), (600, 900)], "+", sr, 0.80, boundary_tol_bp=0, require_single_exon_sr=True
        )
        assert kept is False

    def test_boundary_tolerance_both_ends_must_agree(self):
        # Candidate (100, 900); SR (100, 920) — 3' diff = 20bp
        sr_exact = [_tx("t1", "+", (100, 900))]
        sr_off = [_tx("t2", "+", (100, 921))]
        # With boundary_tol=20: 921-900=21 > 20 → rejected
        kept_off, _, _ = s1_reciprocal_boundary(
            [(100, 900)], "+", sr_off, 0.80, boundary_tol_bp=20, require_single_exon_sr=True
        )
        kept_exact, _, _ = s1_reciprocal_boundary(
            [(100, 900)], "+", sr_exact, 0.80, boundary_tol_bp=20, require_single_exon_sr=True
        )
        assert kept_exact is True
        assert kept_off is False

    def test_boundary_tolerance_disabled_when_zero(self):
        # boundary_tol=0 → only reciprocal overlap is checked
        sr = [_tx("t1", "+", (100, 950))]
        kept, _, _ = s1_reciprocal_boundary(
            [(100, 900)], "+", sr, 0.80, boundary_tol_bp=0, require_single_exon_sr=True
        )
        # overlap 800/950 ≈ 0.84 — above threshold; boundary not checked
        assert kept is True

    def test_no_sr_transcripts(self):
        kept, tids, rej = s1_reciprocal_boundary(
            [(100, 900)], "+", [], 0.80, boundary_tol_bp=20, require_single_exon_sr=True
        )
        assert kept is False
        assert tids == []
        assert rej is None


# ---------------------------------------------------------------------------
# s2_dual_assembler_boundary
# ---------------------------------------------------------------------------


class TestS2DualAssemblerBoundary:
    def test_both_assemblers_support(self):
        scallop = [_tx("s1", "+", (100, 900), source="scallop")]
        stringtie = [_tx("st1", "+", (100, 900), source="stringtie")]
        kept, tids, rej = s2_dual_assembler_boundary(
            [(100, 900)], "+", scallop, stringtie, 0.80, 0, True
        )
        assert kept is True

    def test_only_one_assembler_fails(self):
        scallop = [_tx("s1", "+", (100, 900), source="scallop")]
        kept, _, _ = s2_dual_assembler_boundary([(100, 900)], "+", scallop, [], 0.80, 0, True)
        assert kept is False


# ---------------------------------------------------------------------------
# s3_independent_evidence
# ---------------------------------------------------------------------------


class TestS3IndependentEvidence:
    def test_two_independent_transcripts(self):
        sr = [_tx("t1", "+", (100, 900)), _tx("t2", "+", (100, 900))]
        kept, tids, rej = s3_independent_evidence(
            [(100, 900)], "+", sr, 0.80, 0, True, n_required=2
        )
        assert kept is True
        assert len(tids) >= 2

    def test_one_transcript_fails(self):
        sr = [_tx("t1", "+", (100, 900))]
        kept, _, _ = s3_independent_evidence([(100, 900)], "+", sr, 0.80, 0, True, n_required=2)
        assert kept is False

    def test_multi_exon_candidate_fails(self):
        sr = [_tx("t1", "+", (100, 900)), _tx("t2", "+", (100, 900))]
        kept, _, _ = s3_independent_evidence(
            [(100, 500), (600, 900)], "+", sr, 0.80, 0, True, n_required=2
        )
        assert kept is False

    def test_custom_n_required(self):
        sr = [_tx("t1", "+", (100, 900)), _tx("t2", "+", (100, 900))]
        kept_2, _, _ = s3_independent_evidence([(100, 900)], "+", sr, 0.80, 0, True, n_required=2)
        kept_3, _, _ = s3_independent_evidence([(100, 900)], "+", sr, 0.80, 0, True, n_required=3)
        assert kept_2 is True
        assert kept_3 is False

    def test_multiexon_sr_not_counted_when_required_single(self):
        # ONLY multi-exon SR support — no single-exon SR qualifies.
        # Rejection should be LIKELY_FRAGMENT_OF_MULTIEXON_TRANSCRIPT.
        sr = [
            _tx("t1", "+", (100, 900), (950, 1050)),  # multi-exon: first exon covers cand
        ]
        kept, _, rej = s3_independent_evidence([(100, 900)], "+", sr, 0.80, 0, True, n_required=2)
        assert kept is False
        assert rej == "LIKELY_FRAGMENT_OF_MULTIEXON_TRANSCRIPT"

    def test_one_single_exon_sr_insufficient_for_n2(self):
        # One single-exon SR matches but n_required=2 not met.
        # Rejection is None (insufficient independent evidence, not a fragment).
        sr = [
            _tx("t1", "+", (100, 900), (950, 1050)),  # multi-exon: should not count
            _tx("t2", "+", (100, 900)),  # single-exon: counts once
        ]
        kept, tids, rej = s3_independent_evidence(
            [(100, 900)], "+", sr, 0.80, 0, True, n_required=2
        )
        assert kept is False
        assert rej is None
        assert len(tids) < 2


# ---------------------------------------------------------------------------
# Dispatch: try_rescue_multi_exon / try_rescue_single_exon
# ---------------------------------------------------------------------------


class TestRescueDispatch:
    def _multi_sr(self):
        return _sr_data(
            scallop=[_tx("s1", "+", (100, 200), (300, 400), source="scallop")],
            stringtie=[_tx("st1", "+", (100, 200), (300, 400), source="stringtie")],
        )

    def test_disabled_rescue_returns_false_multi(self):
        cfg = _cfg(shortread_rescue=ShortreadRescueConfig(enabled=False))
        kept, reason, tids = try_rescue_multi_exon(
            [(100, 200), (300, 400)], "+", self._multi_sr(), cfg
        )
        assert kept is False
        assert reason is None
        assert tids == []

    def test_disabled_rescue_returns_false_single(self):
        cfg = _cfg(shortread_rescue=ShortreadRescueConfig(enabled=False))
        kept, reason, tids, rej = try_rescue_single_exon([(100, 900)], "+", _sr_data(), cfg)
        assert kept is False

    def test_m1_dispatch(self):
        cfg = _cfg(shortread_rescue=_rescue_cfg(me_policy="M1", me_tol=0))
        kept, reason, tids = try_rescue_multi_exon(
            [(100, 200), (300, 400)], "+", self._multi_sr(), cfg
        )
        assert kept is True
        assert reason == "M1_junction_complete"

    def test_m2_dispatch(self):
        cfg = _cfg(shortread_rescue=_rescue_cfg(me_policy="M2", me_tol=0))
        kept, reason, tids = try_rescue_multi_exon(
            [(100, 200), (300, 400)], "+", self._multi_sr(), cfg
        )
        assert kept is True
        assert reason == "M2_full_chain"

    def test_m3_dispatch(self):
        cfg = _cfg(shortread_rescue=_rescue_cfg(me_policy="M3", me_tol=0))
        kept, reason, tids = try_rescue_multi_exon(
            [(100, 200), (300, 400)], "+", self._multi_sr(), cfg
        )
        assert kept is True
        assert reason == "M3_dual_assembler"

    def test_s1_dispatch(self):
        sr = _sr_data(scallop=[_tx("s1", "+", (100, 900))])
        cfg = _cfg(shortread_rescue=_rescue_cfg(se_enabled=True, se_policy="S1", se_overlap=0.80))
        kept, reason, tids, rej = try_rescue_single_exon([(100, 900)], "+", sr, cfg)
        assert kept is True
        assert reason == "S1_reciprocal_boundary"

    def test_s2_dispatch(self):
        sr = _sr_data(
            scallop=[_tx("s1", "+", (100, 900))],
            stringtie=[_tx("st1", "+", (100, 900))],
        )
        cfg = _cfg(shortread_rescue=_rescue_cfg(se_enabled=True, se_policy="S2", se_overlap=0.80))
        kept, reason, tids, rej = try_rescue_single_exon([(100, 900)], "+", sr, cfg)
        assert kept is True
        assert reason == "S2_dual_assembler"

    def test_s3_dispatch(self):
        sr = _sr_data(
            scallop=[_tx("s1", "+", (100, 900)), _tx("s2", "+", (100, 900))],
        )
        cfg = _cfg(shortread_rescue=_rescue_cfg(se_enabled=True, se_policy="S3", se_overlap=0.80))
        kept, reason, tids, rej = try_rescue_single_exon([(100, 900)], "+", sr, cfg)
        assert kept is True
        assert reason == "S3_independent_evidence"

    def test_unknown_multi_exon_policy_returns_false(self):
        cfg = _cfg(shortread_rescue=_rescue_cfg(me_policy="M99"))
        kept, reason, tids = try_rescue_multi_exon(
            [(100, 200), (300, 400)], "+", self._multi_sr(), cfg
        )
        assert kept is False

    def test_unknown_single_exon_policy_returns_false(self):
        cfg = _cfg(shortread_rescue=_rescue_cfg(se_enabled=True, se_policy="S99"))
        kept, reason, tids, rej = try_rescue_single_exon([(100, 900)], "+", _sr_data(), cfg)
        assert kept is False


# ---------------------------------------------------------------------------
# decide_keep_or_rescue
# ---------------------------------------------------------------------------


class TestDecideKeepOrRescue:
    def _member(self, exons):
        return {
            "exons": sorted(exons),
            "start": min(e[0] for e in exons),
            "end": max(e[1] for e in exons),
        }

    def test_above_threshold_no_rescue(self):
        cfg = _cfg()
        kept, rescued, reason, tids, rej = decide_keep_or_rescue(3, 2, [], "+", cfg, None)
        assert kept is True
        assert rescued is False

    def test_below_threshold_no_rescue(self):
        cfg = _cfg()
        kept, rescued, reason, tids, rej = decide_keep_or_rescue(1, 2, [], "+", cfg, None)
        assert kept is False
        assert rescued is False

    def test_support_2_of_3_not_single_read(self):
        cfg = _cfg(shortread_rescue=_rescue_cfg(me_policy="M1"))
        kept, _, _, _, _ = decide_keep_or_rescue(2, 3, [], "+", cfg, None)
        assert kept is False

    def test_sr_data_none_blocks_rescue(self):
        members = [self._member([(100, 200), (300, 400)])]
        cfg = _cfg(shortread_rescue=_rescue_cfg(me_policy="M1"))
        kept, _, _, _, _ = decide_keep_or_rescue(1, 2, members, "+", cfg, None)
        assert kept is False

    def test_multi_exon_rescue_triggered(self):
        members = [self._member([(100, 200), (300, 400)])]
        sr = _sr_data(scallop=[_tx("s1", "+", (100, 200), (300, 400))])
        cfg = _cfg(shortread_rescue=_rescue_cfg(me_policy="M1", me_tol=0))
        kept, rescued, reason, tids, rej = decide_keep_or_rescue(1, 2, members, "+", cfg, sr)
        assert kept is True
        assert rescued is True
        assert reason == "M1_junction_complete"

    def test_single_exon_rescue_triggered(self):
        members = [self._member([(100, 900)])]
        sr = _sr_data(scallop=[_tx("s1", "+", (100, 900))])
        cfg = _cfg(shortread_rescue=_rescue_cfg(se_enabled=True, se_policy="S1", se_overlap=0.80))
        kept, rescued, reason, tids, rej = decide_keep_or_rescue(1, 4, members, "+", cfg, sr)
        assert kept is True
        assert rescued is True
        assert reason == "S1_reciprocal_boundary"

    def test_disabled_rescue_not_triggered(self):
        members = [self._member([(100, 900)])]
        sr = _sr_data(scallop=[_tx("s1", "+", (100, 900))])
        cfg = _cfg()  # no rescue policies enabled (shortread_rescue.enabled=False)
        kept, _, _, _, _ = decide_keep_or_rescue(1, 4, members, "+", cfg, sr)
        assert kept is False


# ---------------------------------------------------------------------------
# Configuration: nested structure, validation, legacy key rejection
# ---------------------------------------------------------------------------


class TestConfig:
    def test_default_rescue_disabled(self):
        cfg = LongreadConsensusConfig()
        assert cfg.shortread_rescue.enabled is False

    def test_nested_rescue_fields_accessible(self):
        cfg = LongreadConsensusConfig()
        assert cfg.shortread_rescue.multi_exon.policy == "M1"
        assert cfg.shortread_rescue.single_exon.reciprocal_overlap_min == 0.80

    def test_load_config_from_dict_equiv(self, tmp_path):
        import yaml as _yaml

        cfg_file = tmp_path / "cfg.yaml"
        cfg_file.write_text(
            "max_span_length: 50000\n"
            "shortread_rescue:\n"
            "  enabled: true\n"
            "  multi_exon:\n"
            "    policy: M2\n"
        )
        cfg = load_longread_consensus_config(path=str(cfg_file))
        assert cfg.max_span_length == 50000
        assert cfg.shortread_rescue.enabled is True
        assert cfg.shortread_rescue.multi_exon.policy == "M2"

    def test_unknown_top_level_key_raises(self, tmp_path):
        cfg_file = tmp_path / "cfg.yaml"
        cfg_file.write_text("totally_unknown_key: 99\n")
        with pytest.raises(ValueError, match="Unknown longread_consensus config key"):
            load_longread_consensus_config(path=str(cfg_file))

    def test_legacy_key_require_shortread_raises(self, tmp_path):
        cfg_file = tmp_path / "cfg.yaml"
        cfg_file.write_text("require_shortread_support_for_single_read_models: true\n")
        with pytest.raises(ValueError, match="removed"):
            load_longread_consensus_config(path=str(cfg_file))

    def test_legacy_key_shortread_overlap_fraction_raises(self, tmp_path):
        cfg_file = tmp_path / "cfg.yaml"
        cfg_file.write_text("shortread_overlap_fraction: 0.5\n")
        with pytest.raises(ValueError, match="removed"):
            load_longread_consensus_config(path=str(cfg_file))

    def test_legacy_key_rescue_multi_exon_policy_raises(self, tmp_path):
        cfg_file = tmp_path / "cfg.yaml"
        cfg_file.write_text("rescue_multi_exon_policy: M1\n")
        with pytest.raises(ValueError, match="Moved"):
            load_longread_consensus_config(path=str(cfg_file))

    def test_validation_passes_pure_mode(self):
        cfg = LongreadConsensusConfig()
        validate_config(cfg, [], [])  # should not raise

    def test_validation_rescue_enabled_no_files(self):
        cfg = _cfg(shortread_rescue=_rescue_cfg())
        with pytest.raises(ValueError, match="no short-read files"):
            validate_config(cfg, [], [])

    def test_validation_rescue_enabled_no_subtype(self):
        rc = ShortreadRescueConfig(
            enabled=True,
            multi_exon=MultiExonRescueConfig(enabled=False),
            single_exon=SingleExonRescueConfig(enabled=False),
        )
        cfg = _cfg(shortread_rescue=rc)
        with pytest.raises(ValueError, match="Enable at least one subtype"):
            validate_config(cfg, ["scallop.gtf"], [])

    def test_validation_invalid_multi_exon_policy(self):
        rc = _rescue_cfg(me_policy="M99")
        cfg = _cfg(shortread_rescue=rc)
        with pytest.raises(ValueError, match="Unknown multi_exon rescue policy"):
            validate_config(cfg, ["scallop.gtf"], [])

    def test_validation_invalid_single_exon_policy(self):
        rc = _rescue_cfg(se_enabled=True, se_policy="S99")
        cfg = _cfg(shortread_rescue=rc)
        with pytest.raises(ValueError, match="Unknown single_exon rescue policy"):
            validate_config(cfg, ["scallop.gtf"], [])

    def test_validation_invalid_overlap_threshold(self):
        rc = _rescue_cfg(se_enabled=True, se_overlap=1.5)
        cfg = _cfg(shortread_rescue=rc)
        with pytest.raises(ValueError, match="reciprocal_overlap_min"):
            validate_config(cfg, ["scallop.gtf"], [])

    def test_validation_negative_tolerance(self):
        rc = _rescue_cfg(me_tol=-1)
        cfg = _cfg(shortread_rescue=rc)
        with pytest.raises(ValueError, match="junction_tolerance_bp"):
            validate_config(cfg, ["scallop.gtf"], [])

    def test_validation_m3_requires_both_assemblers(self):
        rc = _rescue_cfg(me_policy="M3")
        cfg = _cfg(shortread_rescue=rc)
        with pytest.raises(ValueError, match="requires both"):
            validate_config(cfg, ["scallop.gtf"], [])  # no stringtie

    def test_preset_pfalciparum_pure_loads(self):
        cfg = load_longread_consensus_config(preset="pfalciparum_pure")
        assert cfg.shortread_rescue.enabled is False
        assert cfg.max_span_length == 45000

    def test_preset_pfalciparum_assisted_loads(self):
        cfg = load_longread_consensus_config(preset="pfalciparum_assisted")
        assert cfg.shortread_rescue.enabled is True
        assert cfg.shortread_rescue.multi_exon.policy == "M1"

    def test_unknown_preset_raises(self):
        with pytest.raises(FileNotFoundError, match="not found"):
            load_longread_consensus_config(preset="nonexistent_preset")

    def test_config_overrides_preset(self, tmp_path):
        cfg_file = tmp_path / "override.yaml"
        cfg_file.write_text("max_span_length: 99999\n")
        cfg = load_longread_consensus_config(path=str(cfg_file), preset="pfalciparum_pure")
        assert cfg.max_span_length == 99999
        assert cfg.shortread_rescue.enabled is False  # from preset


# ---------------------------------------------------------------------------
# Multiple short-read files
# ---------------------------------------------------------------------------


class TestMultipleShortReadFiles:
    def _write_gtf(self, path, transcripts):
        with open(path, "w") as fh:
            for tid, strand, exons in transcripts:
                for s, e in exons:
                    fh.write(
                        f"1\tScallop\texon\t{s}\t{e}\t.\t{strand}\t.\t" f'transcript_id "{tid}";\n'
                    )

    def test_multiple_scallop_files_merged(self, tmp_path):
        f1 = str(tmp_path / "sc1.gtf")
        f2 = str(tmp_path / "sc2.gtf")
        self._write_gtf(f1, [("t1", "+", [(100, 200), (300, 400)])])
        self._write_gtf(f2, [("t2", "+", [(500, 600), (700, 800)])])
        result = load_shortread_transcripts([f1, f2], [], "1")
        scallop_tids = {tx["raw_tid"] for tx in result["scallop"]}
        assert "t1" in scallop_tids
        assert "t2" in scallop_tids

    def test_duplicate_tid_across_files_namespaced(self, tmp_path):
        f1 = str(tmp_path / "sc1.gtf")
        f2 = str(tmp_path / "sc2.gtf")
        # Both files have transcript "t1" — different exons
        self._write_gtf(f1, [("t1", "+", [(100, 200)])])
        self._write_gtf(f2, [("t1", "+", [(500, 600)])])
        result = load_shortread_transcripts([f1, f2], [], "1")
        assert len(result["scallop"]) == 2  # both present despite same raw_tid
        ns_tids = {tx["tid"] for tx in result["scallop"]}
        assert "scallop:0:t1" in ns_tids
        assert "scallop:1:t1" in ns_tids

    def test_file_provenance_retained(self, tmp_path):
        f1 = str(tmp_path / "sc1.gtf")
        self._write_gtf(f1, [("t1", "+", [(100, 200)])])
        result = load_shortread_transcripts([f1], [], "1")
        tx = result["scallop"][0]
        assert tx["source"] == "scallop"
        assert tx["source_file"] == f1

    def test_missing_file_skipped(self, tmp_path):
        result = load_shortread_transcripts(["/does/not/exist.gtf"], [], "1")
        assert result["scallop"] == []

    def test_mixed_scallop_stringtie(self, tmp_path):
        sc = str(tmp_path / "sc.gtf")
        st = str(tmp_path / "st.gtf")
        self._write_gtf(sc, [("s1", "+", [(100, 200)])])
        self._write_gtf(st, [("st1", "+", [(300, 400)])])
        result = load_shortread_transcripts([sc], [st], "1")
        assert any(tx["raw_tid"] == "s1" for tx in result["scallop"])
        assert any(tx["raw_tid"] == "st1" for tx in result["stringtie"])
        assert len(result["all"]) == 2


# ---------------------------------------------------------------------------
# Integration: build_consensus_for_seqname with rescue
# ---------------------------------------------------------------------------


class TestBuildConsensusWithRescue:
    def _write_split_gtf(self, path, transcripts):
        with open(path, "w") as fh:
            for tid, strand, exons in transcripts:
                for s, e in exons:
                    fh.write(
                        f"1\tMinimap2\texon\t{s}\t{e}\t.\t{strand}\t.\t"
                        f'gene_id "g1"; transcript_id "{tid}";\n'
                    )

    def _write_sr_gtf(self, path, transcripts):
        with open(path, "w") as fh:
            for tid, strand, exons in transcripts:
                for s, e in exons:
                    fh.write(
                        f"1\tScallop\texon\t{s}\t{e}\t.\t{strand}\t.\t" f'transcript_id "{tid}";\n'
                    )

    def test_pure_mode_no_sr_loaded(self, tmp_path):
        split_gtf = str(tmp_path / "split.gtf")
        self._write_split_gtf(
            split_gtf,
            [
                ("lr1", "+", [(1000, 1100), (1200, 1300)]),
                ("lr2", "+", [(1000, 1100), (1200, 1300)]),
            ],
        )
        cfg = LongreadConsensusConfig()  # rescue disabled
        records, stats, dropped = build_consensus_for_seqname("1", split_gtf, cfg)
        assert stats["sr_transcripts_loaded"] == 0
        assert stats["rescued_by_shortread"] == 0

    def test_rescued_record_has_reason_and_tids(self, tmp_path):
        split_gtf = str(tmp_path / "split.gtf")
        sr_gtf = str(tmp_path / "scallop.gtf")
        self._write_split_gtf(split_gtf, [("lr1", "+", [(1000, 1100), (1200, 1300)])])
        self._write_sr_gtf(sr_gtf, [("sr1", "+", [(1000, 1100), (1200, 1300)])])
        rc = _rescue_cfg(me_policy="M1", me_tol=0)
        cfg = _cfg(min_read_support_multi_exon=2, shortread_rescue=rc)
        records, stats, dropped = build_consensus_for_seqname(
            "1", split_gtf, cfg, scallop_paths=[sr_gtf]
        )
        assert stats["rescued_by_shortread"] == 1
        assert len(records) == 1
        rec = records[0]
        assert rec["rescued_by_shortread"] is True
        assert rec["rescue_reason"] == "M1_junction_complete"
        assert any("sr1" in tid for tid in rec["rescue_sr_tids"])

    def test_no_rescue_when_policy_disabled(self, tmp_path):
        split_gtf = str(tmp_path / "split.gtf")
        sr_gtf = str(tmp_path / "scallop.gtf")
        self._write_split_gtf(split_gtf, [("lr1", "+", [(1000, 1100), (1200, 1300)])])
        self._write_sr_gtf(sr_gtf, [("sr1", "+", [(1000, 1100), (1200, 1300)])])
        cfg = LongreadConsensusConfig()  # rescue disabled
        records, stats, dropped = build_consensus_for_seqname(
            "1", split_gtf, cfg, scallop_paths=[sr_gtf]
        )
        assert stats["rescued_by_shortread"] == 0
        assert len(records) == 0

    def test_dropped_records_populated(self, tmp_path):
        split_gtf = str(tmp_path / "split.gtf")
        # One read: single-read, no SR support → will be dropped
        self._write_split_gtf(split_gtf, [("lr1", "+", [(1000, 1100), (1200, 1300)])])
        cfg = LongreadConsensusConfig()
        records, stats, dropped = build_consensus_for_seqname("1", split_gtf, cfg)
        assert len(dropped) >= 1
        reason_codes = {d["reason_code"] for d in dropped}
        assert "DROPPED_SINGLE_READ_NO_RESCUE" in reason_codes


# ---------------------------------------------------------------------------
# write_rescue_attribution_tsv
# ---------------------------------------------------------------------------


class TestWriteRescueAttributionTsv:
    def _record(self, tid, rescued, reason, sr_tids, exons=None, n_exons=1):
        exons = exons or [(100, 900)]
        return {
            "transcript_id": tid,
            "seqname": "1",
            "strand": "+",
            "exons": sorted(exons),
            "n_exons": n_exons,
            "rescued_by_shortread": rescued,
            "rescue_reason": reason,
            "rescue_sr_tids": sr_tids,
        }

    def test_only_rescued_rows_written(self, tmp_path):
        out = str(tmp_path / "rescue.tsv")
        records = [
            self._record("t1", True, "M1_junction_complete", ["scallop:0:sr1"]),
            self._record("t2", False, None, []),
            self._record("t3", True, "S1_reciprocal_boundary", ["scallop:0:sr3"]),
        ]
        n = write_rescue_attribution_tsv(records, out)
        assert n == 2
        lines = open(out).readlines()
        assert lines[0].startswith("transcript_id")
        tids_written = {line.split("\t")[0] for line in lines[1:]}
        assert "t1" in tids_written
        assert "t3" in tids_written
        assert "t2" not in tids_written

    def test_sr_tids_decoded_to_raw_tid(self, tmp_path):
        out = str(tmp_path / "rescue.tsv")
        records = [
            self._record("t1", True, "M1_junction_complete", ["scallop:0:sr1", "scallop:1:sr2"])
        ]
        write_rescue_attribution_tsv(records, out)
        lines = open(out).readlines()
        cols = lines[1].rstrip("\n").split("\t")
        rescue_tids_col = cols[7]  # rescue_sr_tids column
        assert "sr1" in rescue_tids_col
        assert "sr2" in rescue_tids_col
        # Sources column
        sources_col = cols[8]
        assert "scallop" in sources_col

    def test_no_rescued_records_writes_header_only(self, tmp_path):
        out = str(tmp_path / "rescue.tsv")
        n = write_rescue_attribution_tsv([], out)
        assert n == 0
        lines = open(out).readlines()
        assert len(lines) == 1
        assert lines[0].startswith("transcript_id")
