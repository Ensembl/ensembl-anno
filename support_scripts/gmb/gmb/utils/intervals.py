"""Interval arithmetic shared across the pipeline."""

from __future__ import annotations


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
