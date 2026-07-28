"""Tests for gmb.pipeline.duplicate_isoform_collapse.

Covers the 16 scenarios from the duplicate-isoform-collapse task: exact
collapse of identical structures from different sources, evidence
unioning without double-counting, scalar score recalculation, deterministic
retained-ID selection order-independent of input order, and the various
"must never collapse" guards (different UTRs, different ends, different
CDS start/stop, different proteins, opposite strand, different seqname,
genuine alternative isoforms) -- plus provenance-mapping output and
compatibility with canonical selection afterwards.
"""

import pandas as pd
import pytest

from gmb.pipeline.canonical_selection import select_canonical_for_gene
from gmb.pipeline.config import load_config
from gmb.pipeline.duplicate_isoform_collapse import (
    CATEGORY_GENUINE_ALTERNATIVE,
    CATEGORY_IDENTICAL_STRUCTURE_DIFFERENT_EVIDENCE,
    CATEGORY_IDENTICAL_STRUCTURE_SAME_EVIDENCE,
    CATEGORY_NOT_COMPARABLE,
    CATEGORY_SAME_CDS_DIFFERENT_UTR,
    CATEGORY_SAME_INTRON_CHAIN_DIFFERENT_ENDS,
    classify_pair,
    collapse_exact_duplicate_transcripts,
    group_exact_duplicates,
    is_exact_duplicate_pair,
    merge_evidence_sources,
    select_retained_transcript,
)


def _t(tid, exons, cds=None, strand="+", chrom="1", protein="MKV", evidence="Tiberius", **overrides):
    rec = {
        "transcript_id": tid,
        "gene_id": "G1",
        "chrom": chrom,
        "strand": strand,
        "exons": exons,
        "cds": cds if cds is not None else exons,
        "cds_phases": None,
        "protein": protein,
        "evidence_sources": evidence,
        "gmb_score": 3.0,
        "protein_coding_score": 0.5,
        "orf_label": "ORF:3aa ATG|STOP",
        "is_partial_5": False,
        "is_partial_3": False,
    }
    rec.update(overrides)
    return rec


@pytest.fixture
def cfg():
    return load_config().duplicate_isoform_handling


# ---------------------------------------------------------------------------
# 1, 5. identical exon/CDS from two sources collapse deterministically
# ---------------------------------------------------------------------------


class TestExactDuplicateCollapse:
    def test_identical_structure_different_source_is_collapse_eligible(self, cfg):
        t1 = _t("t1", [(100, 200), (300, 400)], evidence="Tiberius")
        t2 = _t("t2", [(100, 200), (300, 400)], evidence="Scallop")
        eligible, category, _reason = is_exact_duplicate_pair(t1, t2, cfg)
        assert eligible is True
        assert category == CATEGORY_IDENTICAL_STRUCTURE_DIFFERENT_EVIDENCE

    def test_identical_structure_identical_evidence(self, cfg):
        t1 = _t("t1", [(100, 200)], evidence="Tiberius")
        t2 = _t("t2", [(100, 200)], evidence="Tiberius")
        _eligible, category, _ = is_exact_duplicate_pair(t1, t2, cfg)
        assert category == CATEGORY_IDENTICAL_STRUCTURE_SAME_EVIDENCE

    def test_grouping_collects_all_exact_duplicates_together(self, cfg):
        transcripts = [
            _t("t1", [(100, 200)]),
            _t("t2", [(100, 200)]),
            _t("t3", [(100, 200)]),
        ]
        groups = group_exact_duplicates(transcripts, cfg)
        assert len(groups) == 1
        assert len(groups[0]) == 3


# ---------------------------------------------------------------------------
# 2, 3. evidence provenance unioned, not double-counted
# ---------------------------------------------------------------------------


class TestEvidenceUnion:
    def test_evidence_sources_unioned(self):
        group = [
            _t("t1", [(100, 200)], evidence="Tiberius"),
            _t("t2", [(100, 200)], evidence="Scallop,StringTie"),
        ]
        merged = merge_evidence_sources(group)
        assert set(merged.split(",")) == {"Tiberius", "Scallop", "StringTie"}

    def test_duplicate_source_not_double_counted(self):
        group = [
            _t("t1", [(100, 200)], evidence="Tiberius,Scallop"),
            _t("t2", [(100, 200)], evidence="Scallop"),
        ]
        merged = merge_evidence_sources(group)
        # "Scallop" appears in both inputs but must appear once in the union.
        assert merged.split(",").count("Scallop") == 1
        assert set(merged.split(",")) == {"Tiberius", "Scallop"}


# ---------------------------------------------------------------------------
# 4. scalar scores recalculated/aggregated as documented (never summed)
# ---------------------------------------------------------------------------


class TestScoreHandling:
    def test_collapse_does_not_sum_gmb_scores(self):
        gff_rows = _make_gff_rows(
            [
                ("t1", [(100, 200), (300, 400)], "Tiberius", 3.0),
                ("t2", [(100, 200), (300, 400)], "Tiberius", 3.0),
            ]
        )
        config = load_config()
        new_rows, log_rows, stats = collapse_exact_duplicate_transcripts(gff_rows, config)
        assert stats["transcripts_removed"] == 1
        retained_mrna = [r for r in new_rows if r["Feature"] == "mRNA"][0]
        # Recalculated/kept score must not be anywhere near a naive sum
        # (3.0 + 3.0 = 6.0) -- duplicated upstream evidence must not grant
        # an artificial score advantage.
        assert retained_mrna["gmb_score"] < 5.0

    def test_score_gap_documented_in_log(self):
        gff_rows = _make_gff_rows(
            [
                ("t1", [(100, 200)], "Tiberius", 2.0),
                ("t2", [(100, 200)], "Tiberius", 2.0),
            ]
        )
        config = load_config()
        _new_rows, log_rows, _stats = collapse_exact_duplicate_transcripts(gff_rows, config)
        removed_row = [r for r in log_rows if r["removed_transcript_id"] is not None][0]
        assert "gmb_score_before" in removed_row
        assert "gmb_score_recalculated" in removed_row


# ---------------------------------------------------------------------------
# 6-9, 11. must NOT collapse: different UTR/ends/CDS bounds/proteins/seqname
# ---------------------------------------------------------------------------


class TestNeverCollapse:
    def test_same_cds_different_utr_retained_by_default(self, cfg):
        # Same CDS, but exon boundaries differ (implies different UTR extent).
        t1 = _t("t1", [(90, 200), (300, 410)], cds=[(100, 200), (300, 400)])
        t2 = _t("t2", [(80, 200), (300, 420)], cds=[(100, 200), (300, 400)])
        eligible, category, reason = is_exact_duplicate_pair(t1, t2, cfg)
        assert category == CATEGORY_SAME_CDS_DIFFERENT_UTR
        assert eligible is False
        assert "preserve_distinct_utrs" in reason

    def test_same_intron_chain_different_ends_retained(self, cfg):
        t1 = _t("t1", [(100, 200), (300, 400)])
        t2 = _t("t2", [(90, 200), (300, 410)])  # same intron junction (200-300), different outer ends
        eligible, category, _reason = is_exact_duplicate_pair(t1, t2, cfg)
        assert category == CATEGORY_SAME_INTRON_CHAIN_DIFFERENT_ENDS
        assert eligible is False

    def test_different_cds_start_stop_retained(self, cfg):
        t1 = _t("t1", [(100, 400)], cds=[(100, 400)])
        t2 = _t("t2", [(100, 400)], cds=[(130, 400)])  # different CDS start
        eligible, category, _reason = is_exact_duplicate_pair(t1, t2, cfg)
        assert eligible is False
        assert category == CATEGORY_GENUINE_ALTERNATIVE

    def test_different_proteins_retained(self, cfg):
        t1 = _t("t1", [(100, 200)], cds=[(100, 200)], protein="MKV")
        t2 = _t("t2", [(100, 200)], cds=[(100, 200)], protein="MKX")
        eligible, category, reason = is_exact_duplicate_pair(t1, t2, cfg)
        assert eligible is False
        assert "protein" in reason

    def test_opposite_strand_never_collapsed(self, cfg):
        t1 = _t("t1", [(100, 200)], strand="+")
        t2 = _t("t2", [(100, 200)], strand="-")
        eligible, category, _reason = is_exact_duplicate_pair(t1, t2, cfg)
        assert eligible is False
        assert category == CATEGORY_NOT_COMPARABLE

    def test_different_seqname_never_collapsed(self, cfg):
        t1 = _t("t1", [(100, 200)], chrom="1")
        t2 = _t("t2", [(100, 200)], chrom="2")
        eligible, category, _reason = is_exact_duplicate_pair(t1, t2, cfg)
        assert eligible is False
        assert category == CATEGORY_NOT_COMPARABLE

    def test_genuine_alternative_splice_isoforms_remain_separate(self, cfg):
        t1 = _t("t1", [(100, 200), (300, 400)])
        t2 = _t("t2", [(100, 250), (300, 400)])  # genuinely different splice site
        eligible, category, _reason = is_exact_duplicate_pair(t1, t2, cfg)
        assert eligible is False
        assert category == CATEGORY_GENUINE_ALTERNATIVE


# ---------------------------------------------------------------------------
# 13. deterministic retained-ID selection, independent of input order
# ---------------------------------------------------------------------------


class TestDeterministicRetainedId:
    def test_retained_id_independent_of_order(self):
        group_a = [
            _t("t_zzz", [(100, 200)], gmb_score=5.0),
            _t("t_aaa", [(100, 200)], gmb_score=5.0),
        ]
        group_b = list(reversed(group_a))
        assert (
            select_retained_transcript(group_a)["transcript_id"]
            == select_retained_transcript(group_b)["transcript_id"]
            == "t_aaa"  # full tie -> lexical
        )

    def test_complete_orf_preferred_over_partial(self):
        group = [
            _t("t_partial", [(100, 200)], is_partial_5=True, gmb_score=10.0),
            _t("t_complete", [(100, 200)], is_partial_5=False, gmb_score=1.0),
        ]
        assert select_retained_transcript(group)["transcript_id"] == "t_complete"

    def test_no_internal_stop_preferred(self):
        group = [
            _t("t_stop", [(100, 200)], protein="MK*V", gmb_score=10.0),
            _t("t_clean", [(100, 200)], protein="MKV", gmb_score=1.0),
        ]
        assert select_retained_transcript(group)["transcript_id"] == "t_clean"

    def test_richer_evidence_preferred(self):
        group = [
            _t("t_single", [(100, 200)], evidence="Tiberius", gmb_score=10.0),
            _t("t_multi", [(100, 200)], evidence="Tiberius,Scallop", gmb_score=1.0),
        ]
        assert select_retained_transcript(group)["transcript_id"] == "t_multi"


# ---------------------------------------------------------------------------
# 14. removed-to-retained provenance mapping is written
# ---------------------------------------------------------------------------


def _make_gff_rows(transcripts, gene_id="G1", chrom="1", strand="+"):
    """transcripts: list of (tid, exons, evidence, gmb_score)."""
    all_start = min(e[0] for _t, exs, _ev, _sc in transcripts for e in exs)
    all_end = max(e[1] for _t, exs, _ev, _sc in transcripts for e in exs)
    rows = [
        {
            "Chromosome": chrom, "Source": "GMB", "Feature": "gene",
            "Start": all_start, "End": all_end, "Score": ".", "Strand": strand,
            "Frame": ".", "ID": gene_id, "Parent": "",
        }
    ]
    for tid, exons, evidence, gmb_score in transcripts:
        start = min(e[0] for e in exons)
        end = max(e[1] for e in exons)
        rows.append(
            {
                "Chromosome": chrom, "Source": "GMB", "Feature": "mRNA",
                "Start": start, "End": end, "Score": ".", "Strand": strand,
                "Frame": ".", "ID": tid, "Parent": gene_id, "Evidence": evidence,
                "gmb_score": gmb_score, "protein_coding_score": 0.5,
                "orf_label": "ORF:3aa ATG|STOP", "is_partial_5": False, "is_partial_3": False,
            }
        )
        for i, (s, e) in enumerate(exons, start=1):
            rows.append(
                {
                    "Chromosome": chrom, "Source": "GMB", "Feature": "exon",
                    "Start": s, "End": e, "Score": ".", "Strand": strand,
                    "Frame": ".", "ID": f"{tid}.exon{i}", "Parent": tid,
                }
            )
            rows.append(
                {
                    "Chromosome": chrom, "Source": "GMB", "Feature": "CDS",
                    "Start": s, "End": e, "Score": ".", "Strand": strand,
                    "Frame": ".", "ID": f"{tid}.cds{i}", "Parent": tid,
                }
            )
    return rows


class TestCollapsePipelineIntegration:
    def test_provenance_mapping_written(self):
        gff_rows = _make_gff_rows(
            [
                ("t1", [(100, 200)], "Tiberius", 3.0),
                ("t2", [(100, 200)], "Scallop", 3.0),
            ]
        )
        config = load_config()
        protein_by_tid = {"t1": "MKV", "t2": "MKV"}
        new_rows, log_rows, stats = collapse_exact_duplicate_transcripts(
            gff_rows, config, protein_by_tid=protein_by_tid
        )
        assert stats["transcripts_removed"] == 1
        assert len(log_rows) == 2  # one row for retained, one for removed
        removed = [r for r in log_rows if r["removed_transcript_id"] is not None]
        assert len(removed) == 1
        assert removed[0]["retained_transcript_id"] in ("t1", "t2")
        # The GFF3 output itself must carry a traceability attribute.
        mrna_rows = [r for r in new_rows if r["Feature"] == "mRNA"]
        assert len(mrna_rows) == 1
        assert "CollapsedFrom" in mrna_rows[0]

    def test_removed_transcript_and_children_dropped_from_gff(self):
        gff_rows = _make_gff_rows(
            [
                ("t1", [(100, 200), (300, 400)], "Tiberius", 3.0),
                ("t2", [(100, 200), (300, 400)], "Tiberius", 3.0),
            ]
        )
        config = load_config()
        protein_by_tid = {"t1": "MKV", "t2": "MKV"}
        new_rows, _log, stats = collapse_exact_duplicate_transcripts(
            gff_rows, config, protein_by_tid=protein_by_tid
        )
        ids_present = {r.get("ID") for r in new_rows}
        parents_present = {r.get("Parent") for r in new_rows}
        removed_tid = "t2" if stats["transcripts_removed"] and "t1" in ids_present else "t1"
        # Neither the removed mRNA nor any of its exon/CDS children survive.
        assert removed_tid not in ids_present
        assert removed_tid not in parents_present

    def test_disabled_config_is_a_no_op(self):
        gff_rows = _make_gff_rows(
            [
                ("t1", [(100, 200)], "Tiberius", 3.0),
                ("t2", [(100, 200)], "Tiberius", 3.0),
            ]
        )
        config = load_config()
        config.duplicate_isoform_handling.collapse_exact_duplicates = False
        new_rows, log_rows, stats = collapse_exact_duplicate_transcripts(gff_rows, config)
        assert new_rows == gff_rows
        assert log_rows == []
        assert stats["transcripts_removed"] == 0


# ---------------------------------------------------------------------------
# 15. canonical selection still works after duplicate collapse
# ---------------------------------------------------------------------------


class TestCanonicalSelectionAfterCollapse:
    def test_canonical_selection_runs_cleanly_post_collapse(self):
        gff_rows = _make_gff_rows(
            [
                ("t1", [(100, 200), (300, 400)], "Tiberius", 3.0),
                ("t2", [(100, 200), (300, 400)], "Tiberius", 3.0),
                ("t3", [(100, 250), (300, 400)], "Scallop", 1.0),  # genuine alt isoform
            ]
        )
        config = load_config()
        protein_by_tid = {"t1": "MKV", "t2": "MKV", "t3": "MKQ"}
        new_rows, _log, stats = collapse_exact_duplicate_transcripts(
            gff_rows, config, protein_by_tid=protein_by_tid
        )
        assert stats["transcripts_removed"] == 1

        mrna_rows = [r for r in new_rows if r["Feature"] == "mRNA"]
        assert len(mrna_rows) == 2  # one retained duplicate winner + the genuine alternative

        records = []
        for m in mrna_rows:
            records.append(
                {
                    "transcript_id": m["ID"],
                    "gene_id": m["Parent"],
                    "evidence_sources": m.get("Evidence", ""),
                    "gmb_score": m.get("gmb_score"),
                    "cds_bp": 100,
                    "diamond_hit": None,
                    "psauron_score": None,
                    "protein_length": 3,
                    "orf_label": "ORF:3aa ATG|STOP",
                    "is_partial_5": False,
                    "is_partial_3": False,
                    "internal_stop_count": 0,
                    "protein_coding_score": None,
                }
            )
        result = select_canonical_for_gene("G1", records, config.canonical_selection, "Tiberius")
        assert result["canonical_transcript_id"] in {r["transcript_id"] for r in records}
        assert result["n_isoforms"] == 2


# ---------------------------------------------------------------------------
# 16. unrelated/non-duplicate transcripts remain unchanged
# ---------------------------------------------------------------------------


class TestNonDuplicateTranscriptsUnchanged:
    def test_gene_with_no_duplicates_untouched(self):
        gff_rows = _make_gff_rows(
            [
                ("t1", [(100, 200), (300, 400)], "Tiberius", 3.0),
                ("t2", [(100, 260), (300, 400)], "Scallop", 2.0),  # genuinely different
            ]
        )
        config = load_config()
        new_rows, log_rows, stats = collapse_exact_duplicate_transcripts(gff_rows, config)
        assert stats["transcripts_removed"] == 0
        assert log_rows == []
        assert new_rows == gff_rows

    def test_single_isoform_genes_untouched(self):
        gff_rows = _make_gff_rows([("t1", [(100, 200)], "Tiberius", 3.0)])
        config = load_config()
        new_rows, log_rows, stats = collapse_exact_duplicate_transcripts(gff_rows, config)
        assert stats["transcripts_removed"] == 0
        assert new_rows == gff_rows
