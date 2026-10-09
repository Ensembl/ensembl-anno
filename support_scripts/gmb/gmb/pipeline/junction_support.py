"""Junction-level evidence for candidate structures.

Positional evidence ("an alignment overlaps this model") is near-saturated in
compact genomes: almost every candidate at a locus overlaps some transcript and
some protein, so it cannot tell competing structures apart. The measures here
compare *splice structure*:

* ``unsupported_intron_counts`` -- how many of a candidate's introns no
  transcriptomic model observes exactly. Junctions are taken from every loaded
  transcriptomic transcript, including those later removed as read-through
  chimeras: a chimera's individual splice junctions are still observed splices.
* ``protein_compatibility_counts`` -- how many spliced protein alignments
  agree with a candidate's CDS intron chain where they overlap, and how many
  contradict it. Unspliced overlaps are uninformative and counted as neither.

Both are reported for every selected transcript; they change selection only
through ``scoring.primary_selection: junction_supported``.
"""

from __future__ import annotations

import bisect
from collections import defaultdict

import pandas as pd

from gmb.utils.intervals import VALID_STRANDS

# Shortest protein/CDS overlap (bp) that is allowed to count as agreement or
# contradiction; below this an alignment only grazes the CDS.
MIN_PROTEIN_CDS_OVERLAP_BP = 30


def _introns(intervals) -> list[tuple[int, int]]:
    ivs = sorted(intervals)
    return [(ivs[i][1], ivs[i + 1][0]) for i in range(len(ivs) - 1)]


def _transcripts(exon_df: pd.DataFrame):
    """Yield (transcript_id, chrom, strand, sorted exons) for each transcript."""
    if exon_df is None or exon_df.empty:
        return
    df = exon_df[["transcript_id", "Chromosome", "Strand", "Start", "End"]]
    for tid, g in df.groupby("transcript_id", observed=True, sort=False):
        yield (tid, str(g["Chromosome"].iloc[0]), str(g["Strand"].iloc[0]),
               sorted(zip(g["Start"].astype(int), g["End"].astype(int))))


def observed_introns(exon_df: pd.DataFrame) -> set[tuple[str, str, int, int]]:
    """Every (chrom, strand, start, end) intron of the transcripts in *exon_df*."""
    out = set()
    for _tid, chrom, strand, exons in _transcripts(exon_df):
        if strand in VALID_STRANDS:
            out.update((chrom, strand, s, e) for s, e in _introns(exons))
    return out


def unsupported_intron_counts(
    candidate_exons: pd.DataFrame, observed: set
) -> dict[str, tuple[int, int]]:
    """``{transcript_id: (n_introns, n_introns_not_in_observed)}``."""
    out = {}
    for tid, chrom, strand, exons in _transcripts(candidate_exons):
        ints = _introns(exons)
        out[tid] = (len(ints), sum((chrom, strand, s, e) not in observed for s, e in ints))
    return out


class ProteinIntronIndex:
    """Spliced protein alignments indexed by (chrom, strand) and start."""

    def __init__(self, protein_exons: pd.DataFrame):
        by_key = defaultdict(list)
        self.max_span = 0
        for _tid, chrom, strand, exons in _transcripts(protein_exons):
            if strand not in VALID_STRANDS:
                continue
            start, end = exons[0][0], exons[-1][1]
            by_key[(chrom, strand)].append((start, end, tuple(_introns(exons))))
            self.max_span = max(self.max_span, end - start)
        self._index = {}
        for key, rows in by_key.items():
            rows.sort()
            self._index[key] = ([r[0] for r in rows], rows)

    def overlapping(self, chrom, strand, lo, hi):
        entry = self._index.get((chrom, strand))
        if entry is None:
            return
        starts, rows = entry
        i = bisect.bisect_left(starts, lo - self.max_span)
        while i < len(rows) and rows[i][0] < hi:
            if rows[i][1] > lo:
                yield rows[i]
            i += 1


def protein_compatibility(cds, chrom, strand, index: ProteinIntronIndex) -> tuple[int, int]:
    """(compatible, incompatible) spliced protein alignments over one CDS."""
    if not cds or strand not in VALID_STRANDS:
        return 0, 0
    cds = sorted(cds)
    lo, hi = cds[0][0], cds[-1][1]
    cds_introns = _introns(cds)
    compatible = incompatible = 0
    for p_start, p_end, p_introns in index.overlapping(chrom, strand, lo, hi):
        a, b = max(lo, p_start), min(hi, p_end)
        if b - a < MIN_PROTEIN_CDS_OVERLAP_BP:
            continue
        p_in = {i for i in p_introns if i[0] >= a and i[1] <= b}
        c_in = {i for i in cds_introns if i[0] >= a and i[1] <= b}
        if not p_in and not c_in:
            continue
        if p_in == c_in:
            compatible += 1
        else:
            incompatible += 1
    return compatible, incompatible


def protein_compatibility_counts(
    candidate_exons: pd.DataFrame,
    candidate_cds: dict[str, list],
    protein_exons: pd.DataFrame,
) -> dict[str, tuple[int, int]]:
    """``{transcript_id: (compatible, incompatible)}`` for candidates with a CDS."""
    if protein_exons is None or protein_exons.empty or candidate_exons.empty:
        return {}
    index = ProteinIntronIndex(protein_exons)
    out = {}
    for tid, chrom, strand, _exons in _transcripts(candidate_exons):
        cds = candidate_cds.get(tid)
        if cds:
            out[tid] = protein_compatibility(cds, chrom, strand, index)
    return out
