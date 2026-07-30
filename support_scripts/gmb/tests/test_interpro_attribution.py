#!/usr/bin/env python3
"""Tests for InterPro attribution outputs (Part 12/13):
canonical_decisions.tsv, the final-canonical annotated GFF3, and
cross-output consistency (exactly one final canonical per gene, no
transcripts ever deleted, non-reviewed genes trivially initial==final).
"""

import os

import pandas as pd
import pytest

from gmb.pipeline.interpro_resolver import (
    build_canonical_decisions,
    write_canonical_decisions,
    write_final_annotated_gff3,
)


def _canonical_row(gene_id, tid, reason="ONLY_ISOFORM"):
    return {
        "gene_id": gene_id,
        "canonical_transcript_id": tid,
        "selection_reason": reason,
    }


def _resolver_decision(
    gene_id,
    initial_id,
    recommended_id,
    final_id,
    reason_code,
    replaced,
    verdict="supports_runner_up",
):
    return {
        "gene_id": gene_id,
        "initial_canonical_transcript_id": initial_id,
        "interpro_recommended_transcript_id": recommended_id,
        "final_canonical_transcript_id": final_id,
        "replaced": replaced,
        "interpro_verdict": verdict,
        "reason_code": reason_code,
        "explanation": "test",
        "safeguards_passed": True,
        "safeguard_checks": {},
        "safeguard_notes": "",
    }


class TestBuildCanonicalDecisions:
    def test_non_reviewed_gene_is_trivial_initial_equals_final(self):
        canonical_rows = [_canonical_row("G1", "G1.t1")]
        decisions = build_canonical_decisions(canonical_rows, resolver_report=[])
        assert len(decisions) == 1
        row = decisions[0]
        assert row["canonical_selection_stage"] == "initial"
        assert row["initial_canonical_transcript_id"] == "G1.t1"
        assert row["final_canonical_transcript_id"] == "G1.t1"
        assert row["replaced"] is False
        assert row["interpro_verdict"] is None

    def test_reviewed_gene_with_replacement_carries_interpro_stage(self):
        canonical_rows = [_canonical_row("G1", "G1.t1", reason="BEST_DIAMOND_COVERAGE")]
        resolver_report = [
            _resolver_decision(
                "G1", "G1.t1", "G1.t2", "G1.t2", "INTERPRO_C_TERMINAL_TRUNCATION_RESOLVED", True
            )
        ]
        decisions = build_canonical_decisions(canonical_rows, resolver_report)
        row = decisions[0]
        assert row["canonical_selection_stage"] == "interpro"
        assert row["initial_canonical_transcript_id"] == "G1.t1"
        assert row["final_canonical_transcript_id"] == "G1.t2"
        assert row["replaced"] is True
        assert row["reason_code"] == "INTERPRO_C_TERMINAL_TRUNCATION_RESOLVED"

    def test_reviewed_gene_without_replacement_still_reports_stage_interpro(self):
        # InterPro was consulted (gene was reviewed) even though it did not
        # change the outcome -- canonical_selection_stage records
        # "consulted", not "changed anything".
        canonical_rows = [_canonical_row("G1", "G1.t1")]
        resolver_report = [
            _resolver_decision(
                "G1",
                "G1.t1",
                "G1.t1",
                "G1.t1",
                "INTERPRO_SUPPORTS_INITIAL",
                False,
                verdict="supports_current",
            )
        ]
        decisions = build_canonical_decisions(canonical_rows, resolver_report)
        row = decisions[0]
        assert row["canonical_selection_stage"] == "interpro"
        assert row["final_canonical_transcript_id"] == "G1.t1"
        assert row["replaced"] is False

    def test_exactly_one_row_per_gene_across_mixed_reviewed_and_not(self):
        canonical_rows = [
            _canonical_row("G1", "G1.t1"),
            _canonical_row("G2", "G2.t1"),
            _canonical_row("G3", "G3.t1"),
        ]
        resolver_report = [
            _resolver_decision(
                "G2", "G2.t1", "G2.t2", "G2.t2", "INTERPRO_COMPLETE_DOMAIN_RESTORED", True
            )
        ]
        decisions = build_canonical_decisions(canonical_rows, resolver_report)
        assert len(decisions) == 3
        assert {d["gene_id"] for d in decisions} == {"G1", "G2", "G3"}
        # Exactly one final canonical transcript recorded per gene.
        final_by_gene = {d["gene_id"]: d["final_canonical_transcript_id"] for d in decisions}
        assert final_by_gene == {"G1": "G1.t1", "G2": "G2.t2", "G3": "G3.t1"}

    def test_mismatched_initial_canonical_raises_rather_than_silently_reconciling(self):
        # Resolver report built from a stale/different canonical_transcripts.tsv
        # than the one passed in here -- must fail loudly, not guess.
        canonical_rows = [_canonical_row("G1", "G1.t1")]
        resolver_report = [
            _resolver_decision(
                "G1", "G1.t_STALE", "G1.t2", "G1.t2", "INTERPRO_ARCHITECTURE_RESTORED", True
            )
        ]
        with pytest.raises(ValueError, match="mismatch"):
            build_canonical_decisions(canonical_rows, resolver_report)


class TestWriteCanonicalDecisions:
    def test_writes_tsv_with_expected_columns(self, tmp_path):
        rows = build_canonical_decisions([_canonical_row("G1", "G1.t1")], [])
        path = write_canonical_decisions(rows, str(tmp_path))
        assert os.path.basename(path) == "canonical_decisions.tsv"
        df = pd.read_csv(path, sep="\t")
        for col in (
            "gene_id",
            "initial_canonical_transcript_id",
            "canonical_selection_stage",
            "final_canonical_transcript_id",
            "replaced",
            "reason_code",
        ):
            assert col in df.columns


class TestWriteFinalAnnotatedGff3:
    def _gff3_text(self):
        return (
            "##gff-version 3\n"
            "1\tGMB\tgene\t100\t900\t.\t+\t.\tID=G1\n"
            "1\tGMB\tmRNA\t100\t500\t.\t+\t.\tID=G1.t1;Parent=G1\n"
            "1\tGMB\texon\t100\t500\t.\t+\t.\tID=G1.t1.exon1;Parent=G1.t1\n"
            "1\tGMB\tmRNA\t100\t900\t.\t+\t.\tID=G1.t2;Parent=G1\n"
            "1\tGMB\texon\t100\t900\t.\t+\t.\tID=G1.t2.exon1;Parent=G1.t2\n"
        )

    def test_final_canonical_flag_reflects_replacement_not_initial(self, tmp_path):
        gff3 = tmp_path / "consensus.gff3"
        gff3.write_text(self._gff3_text())
        out_path = tmp_path / "annotated.gff3"

        decisions = [
            {
                "gene_id": "G1",
                "initial_canonical_transcript_id": "G1.t1",
                "final_canonical_transcript_id": "G1.t2",  # replaced
                "canonical_selection_stage": "interpro",
                "reason_code": "INTERPRO_C_TERMINAL_TRUNCATION_RESOLVED",
            }
        ]
        write_final_annotated_gff3(str(gff3), decisions, str(out_path))

        # Source file must remain byte-identical -- never modified in place.
        assert gff3.read_text() == self._gff3_text()

        lines = out_path.read_text().splitlines()
        mrna_lines = {
            line.split(";")[0].split("ID=")[1]: line for line in lines if "\tmRNA\t" in line
        }
        assert "canonical=0" in mrna_lines["G1.t1"]
        assert "canonical=1" in mrna_lines["G1.t2"]
        assert "canonical_selection_stage=interpro" in mrna_lines["G1.t2"]
        assert "initial_canonical=G1.t1" in mrna_lines["G1.t2"]
        assert "canonical_reason=INTERPRO_C_TERMINAL_TRUNCATION_RESOLVED" in mrna_lines["G1.t2"]
        # Both isoforms still present -- annotation never deletes a transcript.
        assert "G1.t1" in mrna_lines and "G1.t2" in mrna_lines

    def test_exactly_one_canonical_per_gene_in_output(self, tmp_path):
        gff3 = tmp_path / "consensus.gff3"
        gff3.write_text(self._gff3_text())
        out_path = tmp_path / "annotated.gff3"
        decisions = [
            {
                "gene_id": "G1",
                "initial_canonical_transcript_id": "G1.t1",
                "final_canonical_transcript_id": "G1.t1",
                "canonical_selection_stage": "initial",
                "reason_code": "ONLY_ISOFORM",
            }
        ]
        write_final_annotated_gff3(str(gff3), decisions, str(out_path))
        lines = [l for l in out_path.read_text().splitlines() if "\tmRNA\t" in l]
        n_canonical = sum(1 for l in lines if "canonical=1" in l)
        assert n_canonical == 1

    def test_gene_missing_from_decisions_left_untouched_not_guessed(self, tmp_path):
        gff3 = tmp_path / "consensus.gff3"
        gff3.write_text(self._gff3_text())
        out_path = tmp_path / "annotated.gff3"
        write_final_annotated_gff3(str(gff3), [], str(out_path))
        out_text = out_path.read_text()
        assert "canonical=" not in out_text
