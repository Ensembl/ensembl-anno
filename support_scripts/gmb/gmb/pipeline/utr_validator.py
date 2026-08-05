"""Maintained UTR structural validator for GMB GFF3 output.

Enforces the per-end UTR-CDS adjacency invariant:

    For every coding mRNA that carries an explicit five_prime_UTR or
    three_prime_UTR feature, the exonic bases between the UTR proximal
    end and the CDS must be zero.

This module is independent of the external acceptance script (utr_qc.py)
and may be imported and called from the pipeline or from tests.

Public API
----------
validate_utr_partition(gff3_rows) → list[UtrViolation]
    Check a flat list of GFF3 row dicts and return one UtrViolation per
    transcript that violates the invariant.

summarise_violations(violations, n_applicable) → dict
    Aggregate violation list into a summary dict suitable for JSON output.
"""

from __future__ import annotations

from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Interval helpers (no external dependencies)
# ---------------------------------------------------------------------------


def _merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    ivs = sorted(intervals)
    merged = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [tuple(x) for x in merged]  # type: ignore[return-value]


def _total(intervals: list[tuple[int, int]]) -> int:
    return sum(e - s for s, e in intervals)


def _intersect_exonic(gap_start: int, gap_end: int, exon_union: list[tuple[int, int]]) -> int:
    """Return the bp of [gap_start, gap_end) that are covered by exon_union."""
    total = 0
    for es, ee in exon_union:
        s, e = max(gap_start, es), min(gap_end, ee)
        if s < e:
            total += e - s
    return total


# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------


@dataclass
class UtrViolation:
    transcript_id: str
    seqname: str
    strand: str
    utr5_cds_gap_bp: int
    utr3_cds_gap_bp: int
    utr5_len: int
    utr3_len: int
    cds_len: int

    def total_gap_bp(self) -> int:
        return self.utr5_cds_gap_bp + self.utr3_cds_gap_bp


# ---------------------------------------------------------------------------
# Core per-transcript check
# ---------------------------------------------------------------------------


def _check_transcript(
    mrna_id: str,
    seqname: str,
    strand: str,
    exon_ivs: list[tuple[int, int]],
    cds_ivs: list[tuple[int, int]],
    utr5_ivs: list[tuple[int, int]],
    utr3_ivs: list[tuple[int, int]],
) -> UtrViolation | None:
    """Return a UtrViolation if the per-end invariant is violated, else None."""
    if not cds_ivs or not (utr5_ivs or utr3_ivs):
        return None

    exon_union = _merge(exon_ivs)
    cds_union = _merge(cds_ivs)
    utr5_union = _merge(utr5_ivs)
    utr3_union = _merge(utr3_ivs)

    cds_min = cds_union[0][0]
    cds_max = cds_union[-1][1]

    def _gap(utr: list[tuple[int, int]], proximal_is_high: bool) -> int:
        if not utr:
            return 0
        if proximal_is_high:
            utr_prox = utr[-1][1]
            if utr_prox >= cds_min:
                return 0
            return _intersect_exonic(utr_prox, cds_min, exon_union)
        else:
            utr_prox = utr[0][0]
            if cds_max >= utr_prox:
                return 0
            return _intersect_exonic(cds_max, utr_prox, exon_union)

    # 5' UTR: at LOW coords → proximal is HIGH end (max UTR coord <= cds_min)
    utr5_prox_high = bool(utr5_union and utr5_union[-1][1] <= cds_min)
    gap5 = _gap(utr5_union, proximal_is_high=utr5_prox_high)

    # 3' UTR: at LOW coords → proximal is HIGH end (max UTR coord <= cds_min)
    utr3_prox_high = bool(utr3_union and utr3_union[-1][1] <= cds_min)
    gap3 = _gap(utr3_union, proximal_is_high=utr3_prox_high)

    if gap5 == 0 and gap3 == 0:
        return None

    return UtrViolation(
        transcript_id=mrna_id,
        seqname=seqname,
        strand=strand,
        utr5_cds_gap_bp=gap5,
        utr3_cds_gap_bp=gap3,
        utr5_len=_total(utr5_union),
        utr3_len=_total(utr3_union),
        cds_len=_total(cds_union),
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def validate_utr_partition(gff3_rows: list[dict]) -> list[UtrViolation]:
    """Check per-end UTR-CDS adjacency for all transcripts in `gff3_rows`.

    Parameters
    ----------
    gff3_rows:
        Flat list of GFF3 row dicts.  Each dict must have at minimum:
        ``Feature``, ``Start``, ``End``, ``ID``, ``Parent``,
        ``Chromosome`` (or ``Seqname``), ``Strand``.
        Coordinates must be 0-based half-open (GMB internal convention).

    Returns
    -------
    list[UtrViolation]
        One entry per transcript that violates the invariant.
        Empty list means the invariant holds for all applicable transcripts.
    """
    from collections import defaultdict

    mrnas: dict[str, dict] = {}
    children: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(
        lambda: {"exon": [], "CDS": [], "five_prime_UTR": [], "three_prime_UTR": []}
    )

    for row in gff3_rows:
        feat = row.get("Feature", "")
        fid = row.get("ID", "")
        parent = row.get("Parent", "")
        start = int(row["Start"])
        end = int(row["End"])
        seqname = row.get("Chromosome") or row.get("Seqname", "")
        strand = row.get("Strand", "")

        if feat == "mRNA":
            mrnas[fid] = {"seqname": seqname, "strand": strand}
        elif feat in ("exon", "CDS", "five_prime_UTR", "three_prime_UTR"):
            if parent:
                children[parent][feat].append((start, end))

    violations: list[UtrViolation] = []
    for mrna_id, info in mrnas.items():
        ch = children.get(mrna_id, {})
        v = _check_transcript(
            mrna_id=mrna_id,
            seqname=info["seqname"],
            strand=info["strand"],
            exon_ivs=ch.get("exon", []),
            cds_ivs=ch.get("CDS", []),
            utr5_ivs=ch.get("five_prime_UTR", []),
            utr3_ivs=ch.get("three_prime_UTR", []),
        )
        if v is not None:
            violations.append(v)

    return violations


def summarise_violations(violations: list[UtrViolation], n_applicable: int) -> dict:
    """Aggregate violations into a summary suitable for JSON or log output."""
    n_violations = len(violations)
    pct = round(100 * n_violations / n_applicable, 2) if n_applicable else None
    return {
        "applicable_transcripts": n_applicable,
        "violations": n_violations,
        "violation_pct": pct,
        "invariant_ok": n_violations == 0,
        "total_gap_bp": sum(v.total_gap_bp() for v in violations),
    }
