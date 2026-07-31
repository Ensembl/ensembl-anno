"""Unit tests for the structure-aware single-read rescue mechanism.

Tests cover:
  - _junctions_from_exons: empty, single-exon, multi-exon, snap tolerance
  - _m1_junction_complete: all-junctions policy
  - _m2_full_chain: full chain match policy
  - _m3_dual_assembler_chain: dual-assembler chain policy
  - _s1_reciprocal_boundary: reciprocal boundary policy
  - _s2_dual_assembler_boundary: dual-assembler boundary
  - _s3_independent_evidence: N-independent-transcripts policy
  - _rescue_multi_exon_structured / _rescue_single_exon_structured: dispatch
  - _keep_or_rescue: 4-tuple return, policy routing, disabled behaviour
  - _load_shortread_transcripts_for_seqname: structure, strand filtering
  - build_consensus_for_seqname: rescue attribution in records
  - write_rescue_attribution_tsv: only rescued rows, correct columns
"""

import os
import tempfile

import pytest

from gmb.pipeline.longread_consensus import (
    LongreadConsensusConfig,
    _junctions_from_exons,
    _keep_or_rescue,
    _load_shortread_transcripts_for_seqname,
    _m1_junction_complete,
    _m2_full_chain,
    _m3_dual_assembler_chain,
    _rescue_multi_exon_structured,
    _rescue_single_exon_structured,
    _s1_reciprocal_boundary,
    _s2_dual_assembler_boundary,
    _s3_independent_evidence,
    build_consensus_for_seqname,
    write_rescue_attribution_tsv,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _tx(tid, strand, *exon_pairs, source="scallop"):
    """Build a transcript dict as returned by _load_shortread_transcripts_for_seqname."""
    return {"tid": tid, "strand": strand, "exons": list(exon_pairs), "source": source}


def _sr_data(scallop=None, stringtie=None):
    s = scallop or []
    st = stringtie or []
    return {"scallop": s, "stringtie": st, "all": s + st}


def _cfg(**overrides):
    cfg = LongreadConsensusConfig()
    for k, v in overrides.items():
        setattr(cfg, k, v)
    return cfg


# ---------------------------------------------------------------------------
# _junctions_from_exons
# ---------------------------------------------------------------------------


class TestJunctionsFromExons:
    def test_single_exon_returns_empty(self):
        assert _junctions_from_exons([(100, 200)]) == frozenset()

    def test_empty_returns_empty(self):
        assert _junctions_from_exons([]) == frozenset()

    def test_two_exons(self):
        junctions = _junctions_from_exons([(100, 200), (300, 400)])
        assert junctions == frozenset({(200, 300)})

    def test_three_exons(self):
        junctions = _junctions_from_exons([(100, 200), (300, 400), (500, 600)])
        assert junctions == frozenset({(200, 300), (400, 500)})

    def test_unsorted_input_sorted_before_computing(self):
        # reversed order but same junctions
        j_fwd = _junctions_from_exons([(100, 200), (300, 400)])
        j_rev = _junctions_from_exons([(300, 400), (100, 200)])
        assert j_fwd == j_rev

    def test_snap_tolerance_rounds_boundaries(self):
        # Without snap: junction is (195, 306)
        # With snap=10: donor 195->200 (round(19.5)*10=200), acceptor 306->310 (round(30.6)*10=310)
        j_no_snap = _junctions_from_exons([(100, 195), (306, 400)], tol_bp=0)
        j_snap = _junctions_from_exons([(100, 195), (306, 400)], tol_bp=10)
        assert j_no_snap == frozenset({(195, 306)})
        assert j_snap == frozenset({(200, 310)})

    def test_snap_makes_jitter_match(self):
        # Two reads whose splice sites differ by 3bp within a 10bp snap window
        cand = _junctions_from_exons([(100, 203), (298, 400)], tol_bp=10)
        sr = _junctions_from_exons([(100, 200), (300, 400)], tol_bp=10)
        assert cand == sr


# ---------------------------------------------------------------------------
# _m1_junction_complete
# ---------------------------------------------------------------------------


class TestM1JunctionComplete:
    def test_all_junctions_present_returns_true(self):
        sr = [_tx("t1", "+", (100, 200), (300, 400), (500, 600))]
        kept, tids = _m1_junction_complete(
            [(100, 200), (300, 400), (500, 600)], "+", sr, tol_bp=0
        )
        assert kept is True
        assert "t1" in tids

    def test_one_junction_missing_returns_false(self):
        # SR only covers first junction; candidate has second junction not in SR
        sr = [_tx("t1", "+", (100, 200), (300, 400))]
        kept, tids = _m1_junction_complete(
            [(100, 200), (300, 400), (500, 600)], "+", sr, tol_bp=0
        )
        assert kept is False
        assert tids == []

    def test_junctions_split_across_two_sr_transcripts(self):
        # Each SR transcript covers one junction; together they cover all
        sr = [
            _tx("t1", "+", (100, 200), (300, 400)),
            _tx("t2", "+", (300, 400), (500, 600)),
        ]
        kept, tids = _m1_junction_complete(
            [(100, 200), (300, 400), (500, 600)], "+", sr, tol_bp=0
        )
        assert kept is True

    def test_strand_mismatch_excluded(self):
        sr = [_tx("t1", "-", (100, 200), (300, 400))]
        kept, _ = _m1_junction_complete([(100, 200), (300, 400)], "+", sr, tol_bp=0)
        assert kept is False

    def test_single_exon_candidate_returns_false(self):
        sr = [_tx("t1", "+", (100, 400))]
        kept, _ = _m1_junction_complete([(100, 400)], "+", sr, tol_bp=0)
        assert kept is False

    def test_snap_allows_jitter_match(self):
        # Candidate junction at (203, 298); SR junction at (200, 300) -- 3bp jitter
        sr = [_tx("t1", "+", (100, 200), (300, 400))]
        kept_no_snap, _ = _m1_junction_complete(
            [(100, 203), (298, 400)], "+", sr, tol_bp=0
        )
        kept_snapped, _ = _m1_junction_complete(
            [(100, 203), (298, 400)], "+", sr, tol_bp=10
        )
        assert kept_no_snap is False
        assert kept_snapped is True


# ---------------------------------------------------------------------------
# _m2_full_chain
# ---------------------------------------------------------------------------


class TestM2FullChain:
    def test_exact_match(self):
        sr = [_tx("t1", "+", (100, 200), (300, 400))]
        kept, tids = _m2_full_chain([(100, 200), (300, 400)], "+", sr, tol_bp=0)
        assert kept is True
        assert "t1" in tids

    def test_partial_match_fails(self):
        # SR has superset of junctions; not an exact chain match
        sr = [_tx("t1", "+", (100, 200), (300, 400), (500, 600))]
        kept, _ = _m2_full_chain([(100, 200), (300, 400)], "+", sr, tol_bp=0)
        assert kept is False

    def test_no_sr_match_fails(self):
        sr = [_tx("t1", "+", (200, 300), (400, 500))]
        kept, _ = _m2_full_chain([(100, 200), (300, 400)], "+", sr, tol_bp=0)
        assert kept is False

    def test_single_exon_candidate_fails(self):
        sr = [_tx("t1", "+", (100, 400))]
        kept, _ = _m2_full_chain([(100, 400)], "+", sr, tol_bp=0)
        assert kept is False


# ---------------------------------------------------------------------------
# _m3_dual_assembler_chain
# ---------------------------------------------------------------------------


class TestM3DualAssemblerChain:
    def test_both_assemblers_support(self):
        scallop = [_tx("s1", "+", (100, 200), (300, 400), source="scallop")]
        stringtie = [_tx("st1", "+", (100, 200), (300, 400), source="stringtie")]
        kept, tids = _m3_dual_assembler_chain(
            [(100, 200), (300, 400)], "+", scallop, stringtie, tol_bp=0
        )
        assert kept is True
        assert "s1" in tids
        assert "st1" in tids

    def test_only_scallop_supports_fails(self):
        scallop = [_tx("s1", "+", (100, 200), (300, 400), source="scallop")]
        stringtie = [_tx("st1", "+", (200, 300), (400, 500), source="stringtie")]
        kept, _ = _m3_dual_assembler_chain(
            [(100, 200), (300, 400)], "+", scallop, stringtie, tol_bp=0
        )
        assert kept is False

    def test_only_stringtie_supports_fails(self):
        scallop = [_tx("s1", "+", (200, 300), (400, 500), source="scallop")]
        stringtie = [_tx("st1", "+", (100, 200), (300, 400), source="stringtie")]
        kept, _ = _m3_dual_assembler_chain(
            [(100, 200), (300, 400)], "+", scallop, stringtie, tol_bp=0
        )
        assert kept is False

    def test_empty_assembler_list_fails(self):
        scallop = []
        stringtie = [_tx("st1", "+", (100, 200), (300, 400), source="stringtie")]
        kept, _ = _m3_dual_assembler_chain(
            [(100, 200), (300, 400)], "+", scallop, stringtie, tol_bp=0
        )
        assert kept is False


# ---------------------------------------------------------------------------
# _s1_reciprocal_boundary
# ---------------------------------------------------------------------------


class TestS1ReciprocalBoundary:
    def test_high_reciprocal_overlap_returns_true(self):
        # SR exon (100, 900); candidate (100, 900) => 100% overlap
        sr = [_tx("t1", "+", (100, 900))]
        kept, tids = _s1_reciprocal_boundary([(100, 900)], "+", sr, min_overlap=0.80)
        assert kept is True
        assert "t1" in tids

    def test_low_reciprocal_overlap_returns_false(self):
        # SR exon covers only 40% of candidate
        # candidate (100, 600); SR (100, 300) => overlap 200/500 = 0.40
        sr = [_tx("t1", "+", (100, 300))]
        kept, _ = _s1_reciprocal_boundary([(100, 600)], "+", sr, min_overlap=0.80)
        assert kept is False

    def test_strand_mismatch_excluded(self):
        sr = [_tx("t1", "-", (100, 900))]
        kept, _ = _s1_reciprocal_boundary([(100, 900)], "+", sr, min_overlap=0.80)
        assert kept is False

    def test_multi_exon_candidate_fails(self):
        sr = [_tx("t1", "+", (100, 900))]
        kept, _ = _s1_reciprocal_boundary([(100, 500), (600, 900)], "+", sr, min_overlap=0.80)
        assert kept is False

    def test_any_exon_of_multi_exon_sr_can_match(self):
        # SR has two exons; second exon has high overlap with single-exon candidate
        sr = [_tx("t1", "+", (50, 80), (100, 900))]
        kept, tids = _s1_reciprocal_boundary([(100, 900)], "+", sr, min_overlap=0.80)
        assert kept is True
        assert "t1" in tids


# ---------------------------------------------------------------------------
# _s2_dual_assembler_boundary
# ---------------------------------------------------------------------------


class TestS2DualAssemblerBoundary:
    def test_both_assemblers_support(self):
        scallop = [_tx("s1", "+", (100, 900), source="scallop")]
        stringtie = [_tx("st1", "+", (100, 900), source="stringtie")]
        kept, tids = _s2_dual_assembler_boundary(
            [(100, 900)], "+", scallop, stringtie, min_overlap=0.80
        )
        assert kept is True

    def test_only_one_assembler_fails(self):
        scallop = [_tx("s1", "+", (100, 900), source="scallop")]
        stringtie = []
        kept, _ = _s2_dual_assembler_boundary(
            [(100, 900)], "+", scallop, stringtie, min_overlap=0.80
        )
        assert kept is False


# ---------------------------------------------------------------------------
# _s3_independent_evidence
# ---------------------------------------------------------------------------


class TestS3IndependentEvidence:
    def test_two_independent_transcripts_return_true(self):
        sr = [
            _tx("t1", "+", (100, 900)),
            _tx("t2", "+", (100, 900)),
        ]
        kept, tids = _s3_independent_evidence([(100, 900)], "+", sr, min_overlap=0.80)
        assert kept is True
        assert len(tids) >= 2

    def test_one_transcript_only_fails(self):
        sr = [_tx("t1", "+", (100, 900))]
        kept, _ = _s3_independent_evidence([(100, 900)], "+", sr, min_overlap=0.80)
        assert kept is False

    def test_multi_exon_candidate_fails(self):
        sr = [_tx("t1", "+", (100, 900)), _tx("t2", "+", (100, 900))]
        kept, _ = _s3_independent_evidence(
            [(100, 500), (600, 900)], "+", sr, min_overlap=0.80
        )
        assert kept is False

    def test_custom_n_required(self):
        sr = [_tx("t1", "+", (100, 900)), _tx("t2", "+", (100, 900))]
        kept_2, _ = _s3_independent_evidence([(100, 900)], "+", sr, min_overlap=0.80, n_required=2)
        kept_3, _ = _s3_independent_evidence([(100, 900)], "+", sr, min_overlap=0.80, n_required=3)
        assert kept_2 is True
        assert kept_3 is False


# ---------------------------------------------------------------------------
# Dispatch: _rescue_multi_exon_structured / _rescue_single_exon_structured
# ---------------------------------------------------------------------------


class TestRescueDispatch:
    def _multi_sr(self):
        return _sr_data(
            scallop=[_tx("s1", "+", (100, 200), (300, 400), source="scallop")],
            stringtie=[_tx("st1", "+", (100, 200), (300, 400), source="stringtie")],
        )

    def test_empty_policy_always_returns_false_multi(self):
        cfg = _cfg(rescue_multi_exon_policy="")
        kept, reason, tids = _rescue_multi_exon_structured(
            [(100, 200), (300, 400)], "+", self._multi_sr(), cfg
        )
        assert kept is False
        assert reason is None
        assert tids == []

    def test_empty_policy_always_returns_false_single(self):
        cfg = _cfg(rescue_single_exon_policy="")
        kept, reason, tids = _rescue_single_exon_structured(
            [(100, 900)], "+", _sr_data(), cfg
        )
        assert kept is False

    def test_m1_dispatch(self):
        cfg = _cfg(rescue_multi_exon_policy="M1", rescue_junction_tolerance_bp=0)
        kept, reason, tids = _rescue_multi_exon_structured(
            [(100, 200), (300, 400)], "+", self._multi_sr(), cfg
        )
        assert kept is True
        assert reason == "M1_junction_complete"

    def test_m2_dispatch(self):
        cfg = _cfg(rescue_multi_exon_policy="M2", rescue_junction_tolerance_bp=0)
        kept, reason, tids = _rescue_multi_exon_structured(
            [(100, 200), (300, 400)], "+", self._multi_sr(), cfg
        )
        assert kept is True
        assert reason == "M2_full_chain"

    def test_m3_dispatch(self):
        cfg = _cfg(rescue_multi_exon_policy="M3", rescue_junction_tolerance_bp=0)
        kept, reason, tids = _rescue_multi_exon_structured(
            [(100, 200), (300, 400)], "+", self._multi_sr(), cfg
        )
        assert kept is True
        assert reason == "M3_dual_assembler"

    def test_s1_dispatch(self):
        sr = _sr_data(scallop=[_tx("s1", "+", (100, 900))])
        cfg = _cfg(rescue_single_exon_policy="S1", rescue_reciprocal_overlap_min=0.80)
        kept, reason, tids = _rescue_single_exon_structured([(100, 900)], "+", sr, cfg)
        assert kept is True
        assert reason == "S1_reciprocal_boundary"

    def test_s2_dispatch(self):
        sr = _sr_data(
            scallop=[_tx("s1", "+", (100, 900))],
            stringtie=[_tx("st1", "+", (100, 900))],
        )
        cfg = _cfg(rescue_single_exon_policy="S2", rescue_reciprocal_overlap_min=0.80)
        kept, reason, tids = _rescue_single_exon_structured([(100, 900)], "+", sr, cfg)
        assert kept is True
        assert reason == "S2_dual_assembler"

    def test_s3_dispatch(self):
        sr = _sr_data(
            scallop=[_tx("s1", "+", (100, 900)), _tx("s2", "+", (100, 900))],
        )
        cfg = _cfg(rescue_single_exon_policy="S3", rescue_reciprocal_overlap_min=0.80)
        kept, reason, tids = _rescue_single_exon_structured([(100, 900)], "+", sr, cfg)
        assert kept is True
        assert reason == "S3_independent_evidence"

    def test_unknown_policy_returns_false(self):
        cfg = _cfg(rescue_multi_exon_policy="M99")
        kept, reason, tids = _rescue_multi_exon_structured(
            [(100, 200), (300, 400)], "+", self._multi_sr(), cfg
        )
        assert kept is False
        cfg2 = _cfg(rescue_single_exon_policy="S99")
        kept2, reason2, tids2 = _rescue_single_exon_structured(
            [(100, 900)], "+", _sr_data(), cfg2
        )
        assert kept2 is False


# ---------------------------------------------------------------------------
# _keep_or_rescue (4-tuple return)
# ---------------------------------------------------------------------------


class TestKeepOrRescue:
    def _member(self, exons):
        return {"exons": sorted(exons), "start": min(e[0] for e in exons), "end": max(e[1] for e in exons)}

    def test_above_threshold_no_rescue(self):
        cfg = _cfg()
        kept, rescued, reason, tids = _keep_or_rescue(3, 2, [], "+", None, cfg, None)
        assert kept is True
        assert rescued is False
        assert reason is None
        assert tids == []

    def test_below_threshold_no_rescue_returns_false(self):
        cfg = _cfg()
        kept, rescued, reason, tids = _keep_or_rescue(1, 2, [], "+", None, cfg, None)
        assert kept is False
        assert rescued is False

    def test_support_2_of_3_not_single_read_no_rescue_attempted(self):
        cfg = _cfg(rescue_multi_exon_policy="M1")
        kept, _, _, _ = _keep_or_rescue(2, 3, [], "+", None, cfg, None)
        assert kept is False

    def test_new_rescue_multi_exon_triggered(self):
        members = [self._member([(100, 200), (300, 400)])]
        sr = _sr_data(scallop=[_tx("s1", "+", (100, 200), (300, 400))])
        cfg = _cfg(rescue_multi_exon_policy="M1", rescue_junction_tolerance_bp=0)
        kept, rescued, reason, tids = _keep_or_rescue(1, 2, members, "+", None, cfg, sr)
        assert kept is True
        assert rescued is True
        assert reason == "M1_junction_complete"
        assert "s1" in tids

    def test_new_rescue_single_exon_triggered(self):
        members = [self._member([(100, 900)])]
        sr = _sr_data(scallop=[_tx("s1", "+", (100, 900))])
        cfg = _cfg(rescue_single_exon_policy="S1", rescue_reciprocal_overlap_min=0.80)
        kept, rescued, reason, tids = _keep_or_rescue(1, 4, members, "+", None, cfg, sr)
        assert kept is True
        assert rescued is True
        assert reason == "S1_reciprocal_boundary"

    def test_old_legacy_rescue_still_works(self):
        members = [self._member([(100, 900)])]
        shortread_exons = [(100, 900, "+")]
        cfg = _cfg(
            require_shortread_support_for_single_read_models=True,
            shortread_overlap_fraction=0.5,
        )
        kept, rescued, reason, tids = _keep_or_rescue(1, 4, members, "+", shortread_exons, cfg, None)
        assert kept is True
        assert rescued is True
        assert reason == "legacy_overlap"
        assert tids == []

    def test_no_rescue_when_sr_data_none(self):
        members = [self._member([(100, 200), (300, 400)])]
        cfg = _cfg(rescue_multi_exon_policy="M1")
        kept, _, _, _ = _keep_or_rescue(1, 2, members, "+", None, cfg, None)
        assert kept is False

    def test_new_rescue_not_triggered_when_policy_empty(self):
        members = [self._member([(100, 900)])]
        sr = _sr_data(scallop=[_tx("s1", "+", (100, 900))])
        cfg = _cfg()  # no rescue policies set
        kept, _, _, _ = _keep_or_rescue(1, 4, members, "+", None, cfg, sr)
        assert kept is False


# ---------------------------------------------------------------------------
# _load_shortread_transcripts_for_seqname
# ---------------------------------------------------------------------------


class TestLoadShorteadTranscripts:
    def _write_gtf(self, path, rows):
        with open(path, "w") as fh:
            for row in rows:
                fh.write("\t".join(str(x) for x in row) + "\n")

    def test_loads_exons_grouped_by_transcript(self, tmp_path):
        gtf = str(tmp_path / "scallop.gtf")
        self._write_gtf(
            gtf,
            [
                # seqname source feat start end score strand frame attrs
                ["1", "Scallop", "exon", 100, 200, ".", "+", ".", 'gene_id "g1"; transcript_id "t1";'],
                ["1", "Scallop", "exon", 300, 400, ".", "+", ".", 'gene_id "g1"; transcript_id "t1";'],
                ["1", "Scallop", "exon", 100, 500, ".", "+", ".", 'gene_id "g2"; transcript_id "t2";'],
                ["2", "Scallop", "exon", 100, 200, ".", "+", ".", 'gene_id "g3"; transcript_id "t3";'],
            ],
        )
        result = _load_shortread_transcripts_for_seqname([gtf], [], "1")
        scallop = result["scallop"]
        tids = {tx["tid"] for tx in scallop}
        assert "t1" in tids
        assert "t2" in tids
        assert "t3" not in tids  # different seqname

        t1 = next(tx for tx in scallop if tx["tid"] == "t1")
        assert sorted(t1["exons"]) == [(100, 200), (300, 400)]
        assert t1["strand"] == "+"

    def test_separate_scallop_stringtie_sources(self, tmp_path):
        scallop_gtf = str(tmp_path / "scallop.gtf")
        stringtie_gtf = str(tmp_path / "stringtie.gtf")
        self._write_gtf(
            scallop_gtf,
            [["1", "Scallop", "exon", 100, 200, ".", "+", ".", 'transcript_id "s1";']],
        )
        self._write_gtf(
            stringtie_gtf,
            [["1", "StringTie", "exon", 300, 400, ".", "-", ".", 'transcript_id "st1";']],
        )
        result = _load_shortread_transcripts_for_seqname([scallop_gtf], [stringtie_gtf], "1")
        assert any(tx["tid"] == "s1" for tx in result["scallop"])
        assert any(tx["tid"] == "st1" for tx in result["stringtie"])
        all_tids = {tx["tid"] for tx in result["all"]}
        assert "s1" in all_tids and "st1" in all_tids

    def test_missing_path_handled_gracefully(self, tmp_path):
        result = _load_shortread_transcripts_for_seqname(
            ["/does/not/exist.gtf"], [], "1"
        )
        assert result["scallop"] == []
        assert result["all"] == []


# ---------------------------------------------------------------------------
# Integration: build_consensus_for_seqname with rescue
# ---------------------------------------------------------------------------


class TestBuildConsensusWithRescue:
    def _write_split_gtf(self, path, transcripts):
        """Write a minimal per-seqname split GTF (exon rows only)."""
        with open(path, "w") as fh:
            for tid, strand, exons in transcripts:
                for s, e in exons:
                    fh.write(
                        f'1\tMinimap2\texon\t{s}\t{e}\t.\t{strand}\t.\t'
                        f'gene_id "g1"; transcript_id "{tid}";\n'
                    )

    def _write_sr_gtf(self, path, transcripts):
        with open(path, "w") as fh:
            for tid, strand, exons in transcripts:
                for s, e in exons:
                    fh.write(
                        f'1\tScallop\texon\t{s}\t{e}\t.\t{strand}\t.\t'
                        f'transcript_id "{tid}";\n'
                    )

    def test_rescued_record_has_reason_and_tids(self, tmp_path):
        split_gtf = str(tmp_path / "split.gtf")
        sr_gtf = str(tmp_path / "scallop.gtf")

        # One single-read multi-exon candidate whose junctions are in SR
        self._write_split_gtf(
            split_gtf,
            [("lr1", "+", [(1000, 1100), (1200, 1300)])],
        )
        self._write_sr_gtf(
            sr_gtf,
            [("sr1", "+", [(1000, 1100), (1200, 1300)])],
        )

        cfg = _cfg(
            min_read_support_multi_exon=2,
            rescue_multi_exon_policy="M1",
            rescue_junction_tolerance_bp=0,
        )
        records, stats = build_consensus_for_seqname(
            "1",
            split_gtf,
            cfg,
            scallop_paths=[sr_gtf],
            stringtie_paths=[],
        )
        assert stats["rescued_by_shortread"] == 1
        assert len(records) == 1
        rec = records[0]
        assert rec["rescued_by_shortread"] is True
        assert rec["rescue_reason"] == "M1_junction_complete"
        assert "sr1" in rec["rescue_sr_tids"]

    def test_unrescued_record_has_no_reason(self, tmp_path):
        split_gtf = str(tmp_path / "split.gtf")
        # Two reads with same chain → support=2 → no rescue needed
        self._write_split_gtf(
            split_gtf,
            [
                ("lr1", "+", [(1000, 1100), (1200, 1300)]),
                ("lr2", "+", [(1000, 1100), (1200, 1300)]),
            ],
        )
        cfg = _cfg(min_read_support_multi_exon=2)
        records, _ = build_consensus_for_seqname("1", split_gtf, cfg)
        assert len(records) == 1
        rec = records[0]
        assert rec["rescued_by_shortread"] is False
        assert rec["rescue_reason"] is None
        assert rec["rescue_sr_tids"] == []

    def test_no_rescue_when_policy_empty(self, tmp_path):
        split_gtf = str(tmp_path / "split.gtf")
        sr_gtf = str(tmp_path / "scallop.gtf")
        self._write_split_gtf(
            split_gtf,
            [("lr1", "+", [(1000, 1100), (1200, 1300)])],
        )
        self._write_sr_gtf(
            sr_gtf,
            [("sr1", "+", [(1000, 1100), (1200, 1300)])],
        )
        cfg = _cfg(min_read_support_multi_exon=2)  # no rescue policy
        records, stats = build_consensus_for_seqname(
            "1", split_gtf, cfg, scallop_paths=[sr_gtf]
        )
        assert stats["rescued_by_shortread"] == 0
        assert len(records) == 0


# ---------------------------------------------------------------------------
# write_rescue_attribution_tsv
# ---------------------------------------------------------------------------


class TestWriteRescueAttributionTsv:
    def _record(self, tid, rescued, reason, tids, exons=None, n_exons=1):
        exons = exons or [(100, 900)]
        return {
            "transcript_id": tid,
            "seqname": "1",
            "strand": "+",
            "exons": sorted(exons),
            "n_exons": n_exons,
            "rescued_by_shortread": rescued,
            "rescue_reason": reason,
            "rescue_sr_tids": tids,
        }

    def test_only_rescued_rows_written(self, tmp_path):
        out = str(tmp_path / "rescue.tsv")
        records = [
            self._record("t1", True, "M1_junction_complete", ["sr1", "sr2"]),
            self._record("t2", False, None, []),
            self._record("t3", True, "S1_reciprocal_boundary", ["sr3"]),
        ]
        n = write_rescue_attribution_tsv(records, out)
        assert n == 2
        lines = open(out).readlines()
        assert lines[0].startswith("transcript_id")  # header
        tids_written = {l.split("\t")[0] for l in lines[1:]}
        assert "t1" in tids_written
        assert "t3" in tids_written
        assert "t2" not in tids_written

    def test_sr_tids_pipe_separated(self, tmp_path):
        out = str(tmp_path / "rescue.tsv")
        records = [self._record("t1", True, "M1_junction_complete", ["sr1", "sr2"])]
        write_rescue_attribution_tsv(records, out)
        lines = open(out).readlines()
        data_line = lines[1]
        cols = data_line.rstrip("\n").split("\t")
        assert cols[-1] == "sr1|sr2"

    def test_empty_sr_tids_written_as_empty_string(self, tmp_path):
        out = str(tmp_path / "rescue.tsv")
        records = [self._record("t1", True, "legacy_overlap", [])]
        write_rescue_attribution_tsv(records, out)
        lines = open(out).readlines()
        cols = lines[1].rstrip("\n").split("\t")
        assert cols[-1] == ""

    def test_no_rescued_records_writes_header_only(self, tmp_path):
        out = str(tmp_path / "rescue.tsv")
        n = write_rescue_attribution_tsv([], out)
        assert n == 0
        lines = open(out).readlines()
        assert len(lines) == 1
        assert lines[0].startswith("transcript_id")


# ---------------------------------------------------------------------------
# Regression: disabled config (variant A) matches original behaviour exactly
# ---------------------------------------------------------------------------


class TestRegressionDisabledConfig:
    def test_no_new_fields_changes_output_count(self, tmp_path):
        """With all rescue disabled (default config), rescued_by_shortread must be 0."""
        split_gtf = str(tmp_path / "split.gtf")
        # Write 1 single-read multi-exon transcript (would need rescue)
        with open(split_gtf, "w") as fh:
            fh.write(
                '1\tMM\texon\t100\t200\t.\t+\t.\tgene_id "g"; transcript_id "t1";\n'
                '1\tMM\texon\t300\t400\t.\t+\t.\tgene_id "g"; transcript_id "t1";\n'
            )
        cfg = LongreadConsensusConfig()  # all defaults, no rescue
        records, stats = build_consensus_for_seqname("1", split_gtf, cfg)
        assert stats["rescued_by_shortread"] == 0
        assert all(not r["rescued_by_shortread"] for r in records)
        assert all(r["rescue_reason"] is None for r in records)
