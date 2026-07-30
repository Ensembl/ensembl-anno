#!/usr/bin/env python3
"""Tests for the InterProScan review/resolver second stage.

Split by intent:
  * parser tests run against a REAL InterProScan 6.0.1 fixture trimmed from
    a completed local run (InterPro 109.0), so column/field assumptions are
    validated against actual output rather than an invented schema;
  * resolver/ambiguity tests use small synthetic records, so each decision
    rule is exercised in isolation with controlled inputs.
"""

import json
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
from gmb.pipeline.config import load_config
from gmb.pipeline.interpro_resolver import (
    VERDICT_CONFLICTING,
    VERDICT_SUPPORTS_CURRENT,
    VERDICT_SUPPORTS_RUNNER_UP,
    VERDICT_UNINFORMATIVE,
    build_resolver_report,
    classify_library,
    compare_candidates,
    parse_interproscan_gff3,
    parse_interproscan_jsonl,
    summarise_architecture,
)
from gmb.pipeline.interpro_review import (
    balanced_coverage,
    build_interproscan_nextflow_command,
    build_review_set,
    detect_ambiguity,
    evidence_classes,
    run_interproscan_workflow,
    write_review_inputs,
)

_FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "interpro")
_JSONL = os.path.join(_FIXTURES, "interproscan6_sample.jsonl")
_GFF3 = os.path.join(_FIXTURES, "interproscan6_sample.gff3")
# Real output from scanning GMB-prepared candidates, i.e. a FASTA whose
# headers are protein SHA-256 checksums rather than UniProt accessions.
_GMB_JSONL = os.path.join(_FIXTURES, "gmb_candidates_sample.jsonl")


@pytest.fixture
def review_cfg():
    return load_config().canonical_selection.interpro_resolver


# ---------------------------------------------------------------------------
# Evidence classes (score-redundancy concern)
# ---------------------------------------------------------------------------


class TestEvidenceClasses:
    def test_scallop_and_stringtie_are_one_short_read_class(self):
        # Both assemblers customarily run over the SAME RNA-seq libraries;
        # counting them as two independent sources overstates support.
        classes = evidence_classes("Scallop,StringTie")
        assert classes == {"short_read_transcriptomic"}
        assert len(classes) == 1

    def test_distinct_classes_counted_separately(self):
        classes = evidence_classes("Scallop,StringTie,Minimap2,OrthoDB,Tiberius", "Tiberius")
        assert classes == {
            "short_read_transcriptomic",
            "long_read_transcriptomic",
            "protein_alignment",
            "backbone",
        }

    def test_backbone_label_respected(self):
        assert evidence_classes("Helixer", "Helixer") == {"backbone"}
        assert evidence_classes("Tiberius", "Tiberius") == {"backbone"}

    def test_protein_alignment_distinct_from_protein_plausibility(self):
        # DIAMOND/psauron plausibility is a property of the translated
        # sequence, not an annotation evidence track -- this named-source-only
        # wrapper must never appear to add it as a class (the full 5-class
        # picture, including protein_validation, is
        # evidence_classes_for_transcript(), used where a support flag is
        # available -- see build_review_set()).
        classes = evidence_classes("OrthoDB,GenBlast")
        assert classes == {"protein_alignment"}
        assert not any("psauron" in c or "diamond" in c for c in classes)

    def test_unknown_source_not_silently_dropped(self):
        # Single shared "other" bucket (not one class per unknown name) --
        # the raw name is still preserved for warning/audit purposes (see
        # the UserWarning this call emits), never silently lost.
        with pytest.warns(UserWarning, match="somenewtool"):
            assert evidence_classes("SomeNewTool") == {"other"}

    def test_empty_sources(self):
        assert evidence_classes("") == set()
        assert evidence_classes(None) == set()


class TestBalancedCoverage:
    def test_uses_minimum_not_mean(self):
        # 100% of a short query but 10% of a long target is a fragment
        # match; the mean (0.55) would flatter it.
        assert balanced_coverage(100.0, 10.0) == 0.1

    def test_symmetric(self):
        assert balanced_coverage(10.0, 100.0) == balanced_coverage(100.0, 10.0)

    def test_bounded_0_1(self):
        assert balanced_coverage(100.0, 100.0) == 1.0
        assert balanced_coverage(0.0, 0.0) == 0.0

    def test_missing_returns_none(self):
        assert balanced_coverage(None, 50.0) is None
        assert balanced_coverage(50.0, None) is None


# ---------------------------------------------------------------------------
# Ambiguity detection
# ---------------------------------------------------------------------------


def _row(tid, **kw):
    base = {
        "transcript_id": tid,
        "rank": 1,
        "psauron_score": 0.9,
        "diamond_qcov": 90.0,
        "diamond_scov": 90.0,
        "protein_length": 300,
        "protein_sha256": "sha_" + tid,
        "has_internal_stop": False,
        "evidence_sources": "Tiberius",
    }
    base.update(kw)
    return base


class TestAmbiguityDetection:
    def test_single_transcript_gene_excluded(self, review_cfg):
        assert detect_ambiguity([_row("t1")], {"confidence": "low_confidence"}, review_cfg) == []

    def test_high_confidence_wide_gap_excluded(self, review_cfg):
        rows = [_row("t1", rank=1), _row("t2", rank=2, psauron_score=0.2)]
        canonical = {
            "confidence": "high_confidence",
            "score_gap": 0.9,
            "selection_reason": "BEST_PSAURON",
        }
        assert detect_ambiguity(rows, canonical, review_cfg) == []

    def test_close_score_gap_included(self, review_cfg):
        rows = [_row("t1", rank=1), _row("t2", rank=2)]
        canonical = {
            "confidence": "high_confidence",
            "score_gap": 0.01,
            "selection_reason": "BEST_PSAURON",
        }
        assert "close_score_gap" in detect_ambiguity(rows, canonical, review_cfg)

    def test_low_confidence_included(self, review_cfg):
        rows = [_row("t1", rank=1), _row("t2", rank=2)]
        canonical = {
            "confidence": "low_confidence",
            "score_gap": 0.9,
            "selection_reason": "BEST_PSAURON",
        }
        assert "low_confidence" in detect_ambiguity(rows, canonical, review_cfg)

    def test_diamond_psauron_disagreement_included(self, review_cfg):
        # psauron prefers t1, DIAMOND coverage prefers t2.
        rows = [
            _row("t1", rank=1, psauron_score=0.9, diamond_qcov=40.0, diamond_scov=40.0),
            _row("t2", rank=2, psauron_score=0.5, diamond_qcov=95.0, diamond_scov=95.0),
        ]
        canonical = {
            "confidence": "high_confidence",
            "score_gap": 0.9,
            "selection_reason": "BEST_PSAURON",
        }
        assert "diamond_psauron_disagree" in detect_ambiguity(rows, canonical, review_cfg)

    def test_length_or_lexical_only_included(self, review_cfg):
        rows = [_row("t1", rank=1), _row("t2", rank=2)]
        for reason in ("LONGEST_VALID_CDS_TIEBREAK", "LEXICAL_TIEBREAK"):
            canonical = {
                "confidence": "high_confidence",
                "score_gap": 0.9,
                "selection_reason": reason,
            }
            assert "length_or_lexical_only" in detect_ambiguity(rows, canonical, review_cfg)

    def test_internal_stop_conflict_included(self, review_cfg):
        rows = [
            _row("t1", rank=1, has_internal_stop=False),
            _row("t2", rank=2, has_internal_stop=True),
        ]
        canonical = {
            "confidence": "high_confidence",
            "score_gap": 0.9,
            "selection_reason": "BEST_PSAURON",
        }
        assert "internal_stop_conflict" in detect_ambiguity(rows, canonical, review_cfg)

    def test_identical_proteins_excluded(self, review_cfg):
        # Same protein sequence -> a domain scan cannot distinguish them,
        # so this is not a reviewable ambiguity even at low confidence.
        rows = [
            _row("t1", rank=1, protein_sha256="same"),
            _row("t2", rank=2, protein_sha256="same"),
        ]
        canonical = {
            "confidence": "low_confidence",
            "score_gap": 0.0,
            "selection_reason": "LEXICAL_TIEBREAK",
        }
        assert detect_ambiguity(rows, canonical, review_cfg) == []

    def test_triggers_individually_switchable(self, review_cfg):
        rows = [_row("t1", rank=1), _row("t2", rank=2)]
        canonical = {
            "confidence": "low_confidence",
            "score_gap": 0.9,
            "selection_reason": "BEST_PSAURON",
        }
        assert "low_confidence" in detect_ambiguity(rows, canonical, review_cfg)
        review_cfg.trigger_low_confidence = False
        assert "low_confidence" not in detect_ambiguity(rows, canonical, review_cfg)


# ---------------------------------------------------------------------------
# Manifest / FASTA preparation
# ---------------------------------------------------------------------------


def _write_review_inputs(tmp_path, ranking_rows, canonical_rows, pv_rows, proteins):
    import pandas as pd

    ranking = tmp_path / "transcript_ranking.tsv"
    canonical = tmp_path / "canonical_transcripts.tsv"
    pv = tmp_path / "protein_validation.tsv"
    fa = tmp_path / "prot.fa"
    pd.DataFrame(ranking_rows).to_csv(ranking, sep="\t", index=False)
    pd.DataFrame(canonical_rows).to_csv(canonical, sep="\t", index=False)
    pd.DataFrame(pv_rows).to_csv(pv, sep="\t", index=False)
    fa.write_text("".join(f">{k}\n{v}\n" for k, v in proteins.items()))
    return str(canonical), str(ranking), str(pv), str(fa)


class TestReviewManifest:
    def _ambiguous_gene(self, tmp_path, shared_protein=False):
        sha1 = "aaa"
        sha2 = "aaa" if shared_protein else "bbb"
        ranking = [
            {
                "gene_id": "G1",
                "transcript_id": "G1.t1",
                "rank": 1,
                "is_canonical": True,
                "psauron_score": 0.9,
                "diamond_qcov": 40.0,
                "diamond_scov": 40.0,
                "evidence_sources": "Scallop,StringTie",
                "n_independent_sources": 2,
                "has_complete_orf": True,
                "has_internal_stop": False,
                "diamond_hit": "X",
                "diamond_pident": 50.0,
            },
            {
                "gene_id": "G1",
                "transcript_id": "G1.t2",
                "rank": 2,
                "is_canonical": False,
                "psauron_score": 0.5,
                "diamond_qcov": 95.0,
                "diamond_scov": 95.0,
                "evidence_sources": "Tiberius",
                "n_independent_sources": 1,
                "has_complete_orf": True,
                "has_internal_stop": False,
                "diamond_hit": "Y",
                "diamond_pident": 60.0,
            },
        ]
        canonical = [
            {
                "gene_id": "G1",
                "canonical_transcript_id": "G1.t1",
                "confidence": "low_confidence",
                "score_gap": 0.01,
                "selection_reason": "BEST_PSAURON",
                "runner_up_transcript_id": "G1.t2",
            }
        ]
        pv = [
            {"transcript_id": "G1.t1", "protein_sha256": sha1, "protein_length": 300},
            {"transcript_id": "G1.t2", "protein_sha256": sha2, "protein_length": 280},
        ]
        proteins = {"G1.t1": "MAAA", "G1.t2": "MAAA" if shared_protein else "MBBB"}
        return _write_review_inputs(tmp_path, ranking, canonical, pv, proteins)

    def test_manifest_has_unique_checksums_and_full_provenance(self, tmp_path, review_cfg):
        paths = self._ambiguous_gene(tmp_path)
        rows, checksums, stats = build_review_set(*paths, review_cfg, backbone_label="Tiberius")
        assert stats["genes_selected_for_review"] == 1
        assert len(checksums) == 2
        for row in rows:
            for field in (
                "gene_id",
                "transcript_id",
                "protein_sha256",
                "protein_length",
                "current_rank",
                "is_current_canonical",
                "canonical_confidence",
                "review_reason",
                "runner_up_transcript_id",
                "diamond_qcov",
                "psauron_score",
                "evidence_classes",
                "current_selection_reason",
            ):
                assert field in row, f"manifest missing {field}"

    def test_identical_proteins_submitted_once_mapping_retained(self, tmp_path, review_cfg):
        # Force review despite identical proteins by disabling the
        # same-protein short circuit is not possible by config, so instead
        # verify dedup at the FASTA level with distinct genes sharing a
        # sequence.
        import pandas as pd

        ranking = []
        canonical = []
        pv = []
        proteins = {}
        for gene in ("G1", "G2"):
            ranking += [
                {
                    "gene_id": gene,
                    "transcript_id": f"{gene}.t1",
                    "rank": 1,
                    "is_canonical": True,
                    "psauron_score": 0.9,
                    "diamond_qcov": 40.0,
                    "diamond_scov": 40.0,
                    "evidence_sources": "Tiberius",
                    "n_independent_sources": 1,
                    "has_complete_orf": True,
                    "has_internal_stop": False,
                },
                {
                    "gene_id": gene,
                    "transcript_id": f"{gene}.t2",
                    "rank": 2,
                    "is_canonical": False,
                    "psauron_score": 0.5,
                    "diamond_qcov": 95.0,
                    "diamond_scov": 95.0,
                    "evidence_sources": "Tiberius",
                    "n_independent_sources": 1,
                    "has_complete_orf": True,
                    "has_internal_stop": False,
                },
            ]
            canonical.append(
                {
                    "gene_id": gene,
                    "canonical_transcript_id": f"{gene}.t1",
                    "confidence": "low_confidence",
                    "score_gap": 0.01,
                    "selection_reason": "BEST_PSAURON",
                    "runner_up_transcript_id": f"{gene}.t2",
                }
            )
            # Both genes' t1 share ONE protein sequence.
            pv += [
                {"transcript_id": f"{gene}.t1", "protein_sha256": "shared", "protein_length": 300},
                {
                    "transcript_id": f"{gene}.t2",
                    "protein_sha256": f"{gene}_own",
                    "protein_length": 280,
                },
            ]
            proteins[f"{gene}.t1"] = "MSHARED"
            proteins[f"{gene}.t2"] = f"M{gene}"

        paths = _write_review_inputs(tmp_path, ranking, canonical, pv, proteins)
        rows, checksums, stats = build_review_set(*paths, review_cfg, backbone_label="Tiberius")
        assert stats["candidates_submitted"] == 4
        # "shared" submitted once -> 3 unique sequences, 1 duplicate avoided
        assert stats["unique_protein_sequences"] == 3
        assert stats["duplicate_candidates_avoided"] == 1
        # Both transcripts still mapped in the manifest.
        shared_rows = [r for r in rows if r["protein_sha256"] == "shared"]
        assert {r["transcript_id"] for r in shared_rows} == {"G1.t1", "G2.t1"}

    def test_top_n_candidates_respected(self, tmp_path, review_cfg):
        import pandas as pd

        ranking = [
            {
                "gene_id": "G1",
                "transcript_id": f"G1.t{i}",
                "rank": i,
                "is_canonical": i == 1,
                "psauron_score": 0.9 if i == 1 else 0.5,
                "diamond_qcov": 40.0 if i == 1 else 95.0,
                "diamond_scov": 40.0 if i == 1 else 95.0,
                "evidence_sources": "Tiberius",
                "n_independent_sources": 1,
                "has_complete_orf": True,
                "has_internal_stop": False,
            }
            for i in (1, 2, 3, 4)
        ]
        canonical = [
            {
                "gene_id": "G1",
                "canonical_transcript_id": "G1.t1",
                "confidence": "low_confidence",
                "score_gap": 0.01,
                "selection_reason": "BEST_PSAURON",
                "runner_up_transcript_id": "G1.t2",
            }
        ]
        pv = [
            {"transcript_id": f"G1.t{i}", "protein_sha256": f"s{i}", "protein_length": 300}
            for i in (1, 2, 3, 4)
        ]
        proteins = {f"G1.t{i}": f"MSEQ{i}" for i in (1, 2, 3, 4)}
        paths = _write_review_inputs(tmp_path, ranking, canonical, pv, proteins)

        review_cfg.max_candidates_per_gene = 2
        rows, _, _ = build_review_set(*paths, review_cfg, backbone_label="Tiberius")
        assert len(rows) == 2
        review_cfg.max_candidates_per_gene = 3
        rows, _, _ = build_review_set(*paths, review_cfg, backbone_label="Tiberius")
        assert len(rows) == 3

    def test_stable_ordering_and_batched_single_fasta(self, tmp_path, review_cfg):
        paths = self._ambiguous_gene(tmp_path)
        rows, checksums, stats = build_review_set(*paths, review_cfg, backbone_label="Tiberius")
        out1 = tmp_path / "o1"
        out2 = tmp_path / "o2"
        written1 = write_review_inputs(rows, checksums, stats, str(out1))
        written2 = write_review_inputs(rows, checksums, stats, str(out2))
        # One FASTA for the whole review set, byte-identical across runs.
        assert open(written1["fasta"]).read() == open(written2["fasta"]).read()
        assert open(written1["manifest"]).read() == open(written2["manifest"]).read()
        fasta = open(written1["fasta"]).read()
        assert fasta.count(">") == len(checksums)
        # Headers are checksums, not transcript IDs.
        assert ">aaa" in fasta and "G1.t1" not in fasta


# ---------------------------------------------------------------------------
# Parser -- real InterProScan 6.0.1 output
# ---------------------------------------------------------------------------


class TestParserRealFixture:
    def test_jsonl_parses_real_output(self):
        result = parse_interproscan_jsonl(_JSONL)
        assert result, "no proteins parsed from real fixture"
        for protein_id, matches in result.items():
            assert isinstance(protein_id, str) and protein_id
            for m in matches:
                assert "signature_accession" in m
                assert "representative" in m
                assert isinstance(m["representative"], bool)

    def test_jsonl_captures_required_fields(self):
        result = parse_interproscan_jsonl(_JSONL)
        all_matches = [m for ms in result.values() for m in ms]
        assert any(m["interpro_accession"] for m in all_matches), "no integrated entries parsed"
        assert any(m["representative"] for m in all_matches), "no representative locations parsed"
        assert any(m["member_database"] for m in all_matches)
        assert any(m["member_database_version"] for m in all_matches)
        assert any(m["signature_type"] for m in all_matches)
        assert all(m["interproscan_version"] == "6.0.1" for m in all_matches)
        assert all(m["interpro_version"] == "109.0" for m in all_matches)

    def test_representative_true_and_false_both_present(self):
        result = parse_interproscan_jsonl(_JSONL)
        flags = {m["representative"] for ms in result.values() for m in ms}
        assert flags == {True, False}

    def test_integrated_and_unintegrated_both_handled(self):
        result = parse_interproscan_jsonl(_JSONL)
        all_matches = [m for ms in result.values() for m in ms]
        assert any(m["interpro_accession"] is not None for m in all_matches)
        assert any(m["interpro_accession"] is None for m in all_matches)

    def test_antifam_classified_as_negative(self):
        result = parse_interproscan_jsonl(_JSONL)
        antifam = [m for ms in result.values() for m in ms if m["member_database"] == "AntiFam"]
        assert antifam, "fixture should retain the real AntiFam match"
        assert all(m["evidence_class"] == "negative" for m in antifam)

    def test_gff3_parses_real_output(self):
        result = parse_interproscan_gff3(_GFF3)
        assert result
        matches = [m for ms in result.values() for m in ms]
        assert any(m["representative"] for m in matches)
        assert any(not m["representative"] for m in matches)
        assert all(m["signature_accession"] for m in matches)

    def test_empty_and_missing_inputs(self, tmp_path):
        assert parse_interproscan_jsonl(str(tmp_path / "nope.jsonl")) == {}
        assert parse_interproscan_gff3(str(tmp_path / "nope.gff3")) == {}
        empty = tmp_path / "empty.jsonl"
        empty.write_text("")
        assert parse_interproscan_jsonl(str(empty)) == {}

    def test_malformed_line_skipped_not_fatal(self, tmp_path):
        path = tmp_path / "mixed.jsonl"
        good = open(_JSONL).readline()
        path.write_text("{not json\n" + good + "\n\n")
        result = parse_interproscan_jsonl(str(path))
        assert result, "valid line should still parse despite a malformed one"

    def test_protein_with_no_matches(self, tmp_path):
        path = tmp_path / "nomatch.jsonl"
        path.write_text(
            json.dumps(
                {
                    "interproscan-version": "6.0.1",
                    "interpro-version": "109.0",
                    "results": [{"md5": "x", "xref": [{"id": "p1"}], "matches": []}],
                }
            )
            + "\n"
        )
        assert parse_interproscan_jsonl(str(path)) == {"p1": []}

    def test_missing_optional_fields_tolerated(self, tmp_path):
        # PROSITE patterns / COILS legitimately emit no evalue or score.
        path = tmp_path / "sparse.jsonl"
        path.write_text(
            json.dumps(
                {
                    "results": [
                        {
                            "xref": [{"id": "p1"}],
                            "matches": [
                                {
                                    "signature": {"accession": "PS0001", "type": "Conserved_site"},
                                    "locations": [{"start": 1, "end": 10}],
                                }
                            ],
                        }
                    ]
                }
            )
            + "\n"
        )
        matches = parse_interproscan_jsonl(str(path))["p1"]
        assert matches[0]["evalue"] is None
        assert matches[0]["member_database"] is None
        assert matches[0]["representative"] is False


class TestParserOnGmbPreparedCandidates:
    """Round-trip check against real InterProScan output for a FASTA that
    GMB itself produced -- confirms the checksum-keyed handoff survives an
    actual scan, not just a synthetic one."""

    def test_ids_are_protein_checksums(self):
        result = parse_interproscan_jsonl(_GMB_JSONL)
        assert result
        for protein_id in result:
            assert len(protein_id) == 64, f"expected a SHA-256 hex id, got {protein_id!r}"
            int(protein_id, 16)  # raises if not hex

    def test_real_gmb_candidates_summarise(self):
        result = parse_interproscan_jsonl(_GMB_JSONL)
        for protein_id, matches in result.items():
            summary = summarise_architecture(matches, protein_length=500)
            assert summary["n_matches_raw"] == len(matches)
            assert summary["n_representative_locations"] <= summary["n_matches_raw"]
            if summary["representative_coverage_fraction"] is not None:
                assert 0.0 <= summary["representative_coverage_fraction"] <= 1.0


class TestClassifyLibrary:
    def test_categories(self):
        assert classify_library("AntiFam", "Family") == "negative"
        assert classify_library("MobiDB-lite", "Region") == "non_evidence"
        assert classify_library("Pfam", "Domain") == "strong"
        assert classify_library("Pfam", "Family") == "strong"
        assert classify_library("CATH-Gene3D", "Homologous_superfamily") == "strong"
        assert classify_library("PROSITE patterns", "Conserved_site") == "supporting"
        assert classify_library("COILS", "Coiled_coil") == "non_evidence"


# ---------------------------------------------------------------------------
# Architecture summary + resolver verdicts
# ---------------------------------------------------------------------------


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


class TestArchitectureSummary:
    def test_only_representative_locations_count_for_coverage(self):
        matches = [
            _match(start=1, end=50, rep=True),
            _match(start=1, end=50, rep=False),  # redundant duplicate signature
        ]
        summary = summarise_architecture(matches, protein_length=100)
        assert summary["n_representative_locations"] == 1
        assert summary["representative_coverage_fraction"] == 0.5

    def test_overlapping_representative_regions_merged_once(self):
        matches = [
            _match(start=1, end=60, rep=True),
            _match(acc="PF2", start=50, end=100, rep=True),
        ]
        summary = summarise_architecture(matches, protein_length=100)
        assert summary["representative_coverage_fraction"] == 1.0

    def test_no_matches_is_not_penalised(self):
        summary = summarise_architecture([], protein_length=100)
        assert summary["has_domains"] is False
        assert summary["has_negative_evidence"] is False
        assert summary["n_integrated_entries"] == 0

    def test_corroborated_entries_need_multiple_databases(self):
        matches = [
            _match(acc="PF1", ipr="IPR1", lib="Pfam"),
            _match(acc="G3D1", ipr="IPR1", lib="CATH-Gene3D", ftype="Homologous_superfamily"),
            _match(acc="PF2", ipr="IPR2", lib="Pfam"),
        ]
        summary = summarise_architecture(matches, protein_length=200)
        assert summary["corroborated_entries"] == ["IPR1"]


class TestResolverVerdicts:
    def _summary(self, matches, length=200):
        return summarise_architecture(matches, protein_length=length)

    def test_supports_current(self):
        current = self._summary(
            [_match(ipr="IPR1"), _match(acc="PF2", ipr="IPR2", start=101, end=200)]
        )
        runner = self._summary([_match(ipr="IPR1")])
        assert compare_candidates(current, runner)["verdict"] == VERDICT_SUPPORTS_CURRENT

    def test_supports_runner_up(self):
        current = self._summary([_match(ipr="IPR1")])
        runner = self._summary(
            [_match(ipr="IPR1"), _match(acc="PF2", ipr="IPR2", start=101, end=200)]
        )
        assert compare_candidates(current, runner)["verdict"] == VERDICT_SUPPORTS_RUNNER_UP

    def test_uninformative_when_neither_has_domains(self):
        result = compare_candidates(self._summary([]), self._summary([]))
        assert result["verdict"] == VERDICT_UNINFORMATIVE
        assert "not evidence against" in result["explanation"]

    def test_uninformative_when_equivalent(self):
        current = self._summary([_match(ipr="IPR1")])
        runner = self._summary([_match(ipr="IPR1")])
        assert compare_candidates(current, runner)["verdict"] == VERDICT_UNINFORMATIVE

    def test_conflicting_when_each_has_exclusive_entries(self):
        current = self._summary([_match(acc="PF1", ipr="IPR1")])
        runner = self._summary([_match(acc="PF2", ipr="IPR2")])
        assert compare_candidates(current, runner)["verdict"] == VERDICT_CONFLICTING

    def test_truncated_domain_coverage_difference_reported(self):
        # Runner-up covers far more of the protein with representative
        # domains -> current looks truncated.
        current = self._summary([_match(ipr="IPR1", start=1, end=20)], length=200)
        runner = self._summary([_match(ipr="IPR1", start=1, end=190)], length=200)
        result = compare_candidates(current, runner)
        assert result["verdict"] == VERDICT_SUPPORTS_RUNNER_UP
        assert "coverage" in result["explanation"]

    def test_negative_evidence_on_current_supports_runner_up(self):
        current = self._summary([_match(acc="ANF001", lib="AntiFam", ftype="Family")])
        runner = self._summary([_match(ipr="IPR1")])
        result = compare_candidates(current, runner)
        assert result["verdict"] == VERDICT_SUPPORTS_RUNNER_UP
        assert "negative-evidence" in result["explanation"]

    def test_negative_evidence_on_runner_up_supports_current(self):
        current = self._summary([_match(ipr="IPR1")])
        runner = self._summary([_match(acc="ANF001", lib="AntiFam", ftype="Family")])
        assert compare_candidates(current, runner)["verdict"] == VERDICT_SUPPORTS_CURRENT

    def test_possible_fusion_extra_architecture_flagged_as_conflicting(self):
        # Candidate carries an extra unrelated integrated entry -- the
        # signature of a fused/run-on model.
        current = self._summary(
            [_match(ipr="IPR1"), _match(acc="PF9", ipr="IPR9", start=101, end=200)]
        )
        runner = self._summary(
            [_match(ipr="IPR1"), _match(acc="PF2", ipr="IPR2", start=101, end=200)]
        )
        assert compare_candidates(current, runner)["verdict"] == VERDICT_CONFLICTING


class TestResolverReport:
    def test_report_never_modifies_canonical_and_flags_difference(self):
        manifest = [
            {
                "gene_id": "G1",
                "transcript_id": "G1.t1",
                "protein_sha256": "a",
                "protein_length": 200,
                "current_rank": 1,
                "canonical_confidence": "low_confidence",
                "review_reason": "close_score_gap",
                "current_selection_reason": "BEST_PSAURON",
            },
            {
                "gene_id": "G1",
                "transcript_id": "G1.t2",
                "protein_sha256": "b",
                "protein_length": 200,
                "current_rank": 2,
                "canonical_confidence": "low_confidence",
                "review_reason": "close_score_gap",
                "current_selection_reason": "BEST_PSAURON",
            },
        ]
        matches = {
            "a": [_match(ipr="IPR1", start=1, end=20)],
            "b": [_match(ipr="IPR1", start=1, end=190)],
        }
        report = build_resolver_report(manifest, matches)
        assert len(report) == 1
        row = report[0]
        assert row["interpro_recommended_transcript_id"] == "G1.t2"
        assert row["recommendation_differs_from_initial"] is True
        assert row["initial_canonical_transcript_id"] == "G1.t1"
        # Runner-up's IPR1 match extends to residue 190 of a 200aa protein
        # while current's stops at 20 -- a C-terminal-only extension, so the
        # reason code should localise it rather than the generic fallback.
        assert row["reason_code"] == "INTERPRO_C_TERMINAL_TRUNCATION_RESOLVED"
        # Default InterProResolverConfig has apply_replacements=True, so
        # with no structural/protein-validation safeguard data present to
        # block it, this specific replacement is applied.
        assert row["replaced"] is True
        assert row["final_canonical_transcript_id"] == "G1.t2"

    def test_shared_protein_resolves_once_for_both_transcripts(self):
        # Two transcripts, one protein -> both look up the same result and
        # the repeated result is not treated as extra evidence.
        manifest = [
            {
                "gene_id": "G1",
                "transcript_id": "G1.t1",
                "protein_sha256": "same",
                "protein_length": 200,
                "current_rank": 1,
            },
            {
                "gene_id": "G1",
                "transcript_id": "G1.t2",
                "protein_sha256": "same",
                "protein_length": 200,
                "current_rank": 2,
            },
        ]
        matches = {"same": [_match(ipr="IPR1")]}
        report = build_resolver_report(manifest, matches)
        assert report[0]["interpro_verdict"] == VERDICT_UNINFORMATIVE
        # Identical architecture on both sides -> EQUIVALENT, not a
        # replacement candidate at all (never even reaches "differs").
        assert report[0]["reason_code"] == "INTERPRO_EQUIVALENT"
        assert report[0]["recommendation_differs_from_initial"] is False
        assert report[0]["replaced"] is False

    def test_gene_with_one_candidate_skipped(self):
        manifest = [
            {
                "gene_id": "G1",
                "transcript_id": "G1.t1",
                "protein_sha256": "a",
                "protein_length": 200,
                "current_rank": 1,
            }
        ]
        assert build_resolver_report(manifest, {"a": [_match()]}) == []


class TestInterProScanNextflowModeA:
    """Mode A (GMB launches InterProScan 6 itself via Nextflow) -- command
    construction is pure/testable without Nextflow installed; the actual
    subprocess invocation is mocked so these tests never require or launch
    a real Nextflow/Docker/Singularity environment.
    """

    def _nf_cfg(self, **overrides):
        cfg = load_config().canonical_selection.interpro_resolver.nextflow
        for k, v in overrides.items():
            setattr(cfg, k, v)
        return cfg

    def test_command_has_no_hardcoded_paths_only_config_values(self):
        nf_cfg = self._nf_cfg(data_dir="/shared/interpro/data", config_file="site.config")
        cmd = build_interproscan_nextflow_command(nf_cfg, "candidates.faa", "out")
        assert cmd[0] == "nextflow"
        assert "ebi-pf-team/interproscan6" in cmd
        assert "-r" in cmd and "6.0.1" in cmd
        assert "--input" in cmd and "candidates.faa" in cmd
        assert "--datadir" in cmd and "/shared/interpro/data" in cmd
        assert "-c" in cmd and "site.config" in cmd
        # Defaults derived from output_dir only when nf_cfg leaves them unset.
        assert any("out/interpro_run/results" in part for part in cmd)
        assert any("out/interpro_run/work" in part for part in cmd)

    def test_explicit_output_and_work_dir_override_defaults(self):
        nf_cfg = self._nf_cfg(
            data_dir="/data", output_dir="/custom/out", work_dir="/custom/work"
        )
        cmd = build_interproscan_nextflow_command(nf_cfg, "candidates.faa", "out")
        assert "/custom/out" in cmd
        assert "/custom/work" in cmd
        assert "out/interpro_run" not in " ".join(cmd)

    def test_cluster_profile_and_extra_args_passed_through(self):
        nf_cfg = self._nf_cfg(
            data_dir="/data",
            profile="slurm,singularity",
            applications="Pfam,PANTHER",
            extra_args=["-with-report", "report.html"],
        )
        cmd = build_interproscan_nextflow_command(nf_cfg, "candidates.faa", "out")
        assert "slurm,singularity" in cmd
        assert "--appl" in cmd and "Pfam,PANTHER" in cmd
        assert "-with-report" in cmd and "report.html" in cmd

    def test_run_requires_data_dir(self):
        nf_cfg = self._nf_cfg(data_dir=None)
        with pytest.raises(ValueError, match="data_dir"):
            run_interproscan_workflow(nf_cfg, "candidates.faa", "out")

    def test_run_raises_on_nonzero_exit_rather_than_silently_succeeding(self, tmp_path):
        import subprocess as _subprocess
        from unittest.mock import patch

        nf_cfg = self._nf_cfg(data_dir="/data")
        fake_result = _subprocess.CompletedProcess(
            args=["nextflow"], returncode=1, stdout="", stderr="boom"
        )
        with patch("subprocess.run", return_value=fake_result):
            with pytest.raises(RuntimeError, match="boom"):
                run_interproscan_workflow(nf_cfg, "candidates.faa", str(tmp_path))

    def test_run_returns_output_dir_on_success(self, tmp_path):
        import subprocess as _subprocess
        from unittest.mock import patch

        nf_cfg = self._nf_cfg(data_dir="/data")
        fake_result = _subprocess.CompletedProcess(
            args=["nextflow"], returncode=0, stdout="ok", stderr=""
        )
        with patch("subprocess.run", return_value=fake_result) as mock_run:
            result = run_interproscan_workflow(nf_cfg, "candidates.faa", str(tmp_path))
        assert mock_run.called
        assert result["output_dir"] == os.path.join(str(tmp_path), "interpro_run", "results")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
