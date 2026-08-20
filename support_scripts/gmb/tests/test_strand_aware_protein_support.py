#!/usr/bin/env python3
"""Tests for same-strand protein-alignment support.

Regression cover for the defect where `PyRanges.overlap` silently matched
antisense protein alignments: every track loaded by `pr.read_gtf` carries
`CategoricalDtype(['.', '-', '+'])`, and pyranges' `check_strandedness`
rejects any categorical Strand whose *levels* include something outside
`+`/`-` even when no row uses it. Both operands therefore reported
`stranded is False` and the overlap degraded to unstranded.
"""

import os
import sys

import pandas as pd
import pyranges as pr
import pytest
from pandas.api.types import CategoricalDtype

sys.path.insert(0, os.path.dirname(__file__))
from gmb.utils.intervals import VALID_STRANDS, same_strand_overlap_ids

# The exact dtype pr.read_gtf produces.
GTF_STRAND_DTYPE = CategoricalDtype([".", "-", "+"], ordered=True)


def _df(rows, strand_dtype=GTF_STRAND_DTYPE):
    """rows: (chrom, start, end, strand, tid)"""
    df = pd.DataFrame(
        {
            "Chromosome": [r[0] for r in rows],
            "Start": [r[1] for r in rows],
            "End": [r[2] for r in rows],
            "Strand": [r[3] for r in rows],
            "transcript_id": [r[4] for r in rows],
        }
    )
    if strand_dtype is not None:
        df["Strand"] = df["Strand"].astype(strand_dtype)
    return df


class TestPyRangesStrandednessQuirk:
    """Pin the upstream behaviour this helper exists to work around."""

    def test_read_gtf_dtype_is_reported_unstranded(self):
        df = _df([("1", 10, 20, "+", "a")])
        assert pr.PyRanges(df).stranded is False

    def test_plain_string_strand_is_stranded(self):
        df = _df([("1", 10, 20, "+", "a")], strand_dtype=None)
        assert pr.PyRanges(df).stranded is True

    def test_plain_overlap_matches_antisense(self):
        """The defect itself: opposite strands still overlap."""
        q = _df([("1", 10, 100, "+", "q1")])
        s = _df([("1", 20, 80, "-", "s1")])
        ovl = pr.PyRanges(q).overlap(pr.PyRanges(s))
        assert not ovl.df.empty, "expected the unstranded fallback to match"


class TestSameStrandOverlapIds:
    def test_same_strand_overlap_is_found(self):
        q = _df([("1", 10, 100, "+", "q1")])
        s = _df([("1", 20, 80, "+", "s1")])
        assert same_strand_overlap_ids(q, s) == {"q1"}

    def test_antisense_overlap_is_rejected(self):
        q = _df([("1", 10, 100, "+", "q1")])
        s = _df([("1", 20, 80, "-", "s1")])
        assert same_strand_overlap_ids(q, s) == set()

    def test_mixed_strands_partitioned_correctly(self):
        q = _df(
            [
                ("1", 10, 100, "+", "plus_hit"),
                ("1", 10, 100, "-", "minus_hit"),
                ("1", 500, 600, "+", "plus_miss"),
            ]
        )
        s = _df([("1", 20, 80, "+", "sp"), ("1", 30, 90, "-", "sm")])
        assert same_strand_overlap_ids(q, s) == {"plus_hit", "minus_hit"}

    def test_unstranded_query_receives_no_support(self):
        """A '.' candidate cannot be strand-matched, so it gets nothing."""
        q = _df([("1", 10, 100, ".", "dot")])
        s = _df([("1", 20, 80, "+", "s1"), ("1", 20, 80, "-", "s2")])
        assert same_strand_overlap_ids(q, s) == set()

    def test_different_chromosome_does_not_match(self):
        q = _df([("1", 10, 100, "+", "q1")])
        s = _df([("2", 10, 100, "+", "s1")])
        assert same_strand_overlap_ids(q, s) == set()

    def test_non_overlapping_same_strand_does_not_match(self):
        q = _df([("1", 10, 20, "+", "q1")])
        s = _df([("1", 100, 200, "+", "s1")])
        assert same_strand_overlap_ids(q, s) == set()

    def test_result_is_independent_of_strand_dtype(self):
        """Behaviour must not depend on the categorical-vs-object quirk."""
        rows_q = [("1", 10, 100, "+", "q1"), ("1", 10, 100, "-", "q2")]
        rows_s = [("1", 20, 80, "+", "s1")]
        cat = same_strand_overlap_ids(_df(rows_q), _df(rows_s))
        obj = same_strand_overlap_ids(
            _df(rows_q, strand_dtype=None), _df(rows_s, strand_dtype=None)
        )
        assert cat == obj == {"q1"}

    @pytest.mark.parametrize("empty_side", ["query", "subject", "both"])
    def test_empty_frames_return_empty_set(self, empty_side):
        q = _df([("1", 10, 100, "+", "q1")])
        s = _df([("1", 20, 80, "+", "s1")])
        if empty_side in ("query", "both"):
            q = q.iloc[0:0]
        if empty_side in ("subject", "both"):
            s = s.iloc[0:0]
        assert same_strand_overlap_ids(q, s) == set()

    def test_none_inputs_return_empty_set(self):
        q = _df([("1", 10, 100, "+", "q1")])
        assert same_strand_overlap_ids(None, q) == set()
        assert same_strand_overlap_ids(q, None) == set()

    def test_custom_id_column(self):
        q = _df([("1", 10, 100, "+", "q1")]).rename(columns={"transcript_id": "gene_id"})
        s = _df([("1", 20, 80, "+", "s1")])
        assert same_strand_overlap_ids(q, s, id_col="gene_id") == {"q1"}

    def test_valid_strands_constant(self):
        assert VALID_STRANDS == ("+", "-")


class TestProteinSupportSemantics:
    """The property that matters for the build: antisense grants no support."""

    def test_antisense_only_candidate_is_unsupported_but_sense_one_is_not(self):
        candidates = _df(
            [
                ("1", 1000, 2000, "+", "sense_candidate"),
                ("1", 1000, 2000, "-", "antisense_candidate"),
            ]
        )
        protein = _df([("1", 1200, 1800, "+", "orthodb_aln")])
        supported = same_strand_overlap_ids(candidates, protein)
        assert supported == {"sense_candidate"}
        # Under the old unstranded behaviour both would have been supported.
        legacy = pr.PyRanges(candidates).overlap(pr.PyRanges(protein))
        assert set(legacy.df["transcript_id"]) == {
            "sense_candidate",
            "antisense_candidate",
        }
