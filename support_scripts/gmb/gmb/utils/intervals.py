"""Interval arithmetic shared across the pipeline."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd


# Strand values that denote a real, usable orientation.  Anything else
# ("." in GTF/GFF3, NaN) is treated as "orientation unknown".
VALID_STRANDS = ("+", "-")


def same_strand_overlap_ids(
    query_df: "pd.DataFrame",
    subject_df: "pd.DataFrame",
    id_col: str = "transcript_id",
) -> set:
    """IDs in *query_df* overlapping *subject_df* on the **same strand**.

    Why this exists rather than ``PyRanges.overlap(..., strandedness="same")``:
    ``pyranges`` refuses stranded operations unless *both* operands report
    ``PyRanges.stranded is True``, and its ``check_strandedness`` rejects any
    ``Strand`` column whose *categorical dtype* lists a level outside ``+``/``-``
    -- **regardless of whether that level is actually used**.  Every track loaded
    via ``pr.read_gtf`` carries ``CategoricalDtype(['.', '-', '+'])``, so all of
    them report ``stranded is False`` and a plain ``.overlap()`` silently
    degrades to unstranded matching.  Normalising the dtype globally is not a
    safe fix either: real candidate sets legitimately contain ``"."`` rows
    (models whose strand could not be resolved), which would keep the operand
    unstranded anyway, and re-typing ``Strand`` pipeline-wide would additionally
    make ``PyRanges.cluster`` strand-aware and silently redefine every locus.

    This helper therefore partitions both frames by strand and intersects within
    each orientation, which is dtype-agnostic and changes nothing else.  Rows
    whose strand is not ``+`` or ``-`` are excluded on both sides: an
    unoriented interval cannot be strand-matched.

    Parameters
    ----------
    query_df, subject_df : pandas.DataFrame
        Interval frames with ``Chromosome``, ``Start``, ``End``, ``Strand``.
        *query_df* must also carry *id_col*.
    id_col : str
        Column of *query_df* holding the identifier to return.

    Returns
    -------
    set
        Distinct *id_col* values with at least one same-strand overlap.
    """
    import pyranges as pr

    if query_df is None or subject_df is None or query_df.empty or subject_df.empty:
        return set()

    hits: set = set()
    q_strand = query_df["Strand"].astype(str)
    s_strand = subject_df["Strand"].astype(str)

    for strand in VALID_STRANDS:
        q = query_df[q_strand == strand]
        s = subject_df[s_strand == strand]
        if q.empty or s.empty:
            continue
        # Within a single-strand partition, orientation is already guaranteed,
        # so an unstranded overlap is exactly a same-strand overlap.
        ovl = pr.PyRanges(q).overlap(pr.PyRanges(s))
        if not ovl.df.empty:
            hits.update(ovl.df[id_col].unique())

    return hits


def cds_span_compatible_ids(
    candidate_spans: dict,
    protein_spans: dict,
    candidates_with_cds: set,
) -> set:
    """Candidates whose protein support is **CDS-span compatible**.

    A candidate qualifies when at least one same-strand protein alignment lies
    entirely inside the candidate's transcript span, and the candidate itself
    has a coding sequence.  This is a stronger statement than bare positional
    overlap: the aligned protein is contained by, rather than merely touching,
    the model it is being counted as evidence for.

    Parameters
    ----------
    candidate_spans : dict
        ``{id: (chrom, strand, start, end)}`` for candidate transcripts.
    protein_spans : dict
        ``{id: (chrom, strand, start, end)}`` for filtered protein alignments.
    candidates_with_cds : set
        Ids of candidates that have a CDS; candidates outside this set can never
        qualify, since there is no coding structure for the alignment to support.

    Returns
    -------
    set
        Ids from *candidate_spans* meeting the criterion.
    """
    import bisect
    from collections import defaultdict

    by_cs = defaultdict(list)
    for chrom, strand, start, end in protein_spans.values():
        if strand in VALID_STRANDS:
            by_cs[(chrom, strand)].append((start, end))

    starts_index = {}
    for key, spans in by_cs.items():
        spans.sort()
        starts_index[key] = ([s for s, _e in spans], spans)

    hits = set()
    for cid, (chrom, strand, c_start, c_end) in candidate_spans.items():
        if cid not in candidates_with_cds or strand not in VALID_STRANDS:
            continue
        idx = starts_index.get((chrom, strand))
        if idx is None:
            continue
        starts, spans = idx
        # Any protein alignment contained in [c_start, c_end] must start at or
        # after c_start; scan forward until its start passes c_end.
        i = bisect.bisect_left(starts, c_start)
        while i < len(spans) and spans[i][0] <= c_end:
            if spans[i][1] <= c_end:
                hits.add(cid)
                break
            i += 1
    return hits


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping intervals into a disjoint set."""
    if not intervals:
        return []
    sorted_ivs = sorted(intervals)
    merged = [sorted_ivs[0]]
    for current in sorted_ivs[1:]:
        last = merged[-1]
        if current[0] <= last[1]:
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)
    return merged


def reciprocal_overlap(s1: int, e1: int, s2: int, e2: int) -> float:
    """Symmetric reciprocal overlap: min of forward and reverse fractions.

    Returns 0.0 if either interval has zero length.
    Used by evidence_filter for redundancy collapsing.
    """
    overlap = max(0, min(e1, e2) - max(s1, s2))
    len1 = e1 - s1
    len2 = e2 - s2
    if len1 == 0 or len2 == 0:
        return 0.0
    return min(overlap / len1, overlap / len2)


def containment_overlap(s1: int, e1: int, s2: int, e2: int) -> float:
    """Containment-style overlap: overlap / min(len1, len2).

    Returns 0.0 if either interval has zero length.
    Used by dedup_genes for structural deduplication.
    """
    overlap = max(0, min(e1, e2) - max(s1, s2))
    len1 = e1 - s1
    len2 = e2 - s2
    if min(len1, len2) == 0:
        return 0.0
    return overlap / min(len1, len2)
