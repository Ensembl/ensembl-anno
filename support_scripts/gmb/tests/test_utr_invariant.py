"""Tests for the GFF3 UTR structural invariant:

    union(exons) == union(CDS u five_prime_UTR u three_prime_UTR)

for every coding transcript that carries explicit UTR features.

Covers two root-cause bugs fixed in gff3_validate.py:
  (A) _hard_cap direction: must keep the CDS-proximal portion of the UTR,
      not the distal portion.  On + strand the 5' UTR sits at LOW genomic
      coordinates; its proximal end is the HIGH end.  On - strand the 3' UTR
      sits at LOW coords; same applies.  Before the fix, trimming kept the
      distal rows and left a gap between the trimmed UTR and the CDS.
  (B) Exon synchronisation: after any UTR trimming the exon features were
      never updated.  The fix rebuilds exon_rows as the intersection of
      (CDS + trimmed UTRs) with the original exon union.
"""

import pytest

from gmb.pipeline.config import PipelineConfig, UtrConfig, ValidationConfig
from gmb.pipeline.gff3_validate import trim_utrs, validate_and_fix_gff3

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _row(feature, start, end, **kw):
    return {"Feature": feature, "Start": start, "End": end, **kw}


def _utr_cfg(**kw):
    defaults = dict(
        trim_policy="hard_cap",
        max_5p_bp=3000,
        max_3p_bp=3000,
        max_total_bp=6000,
        max_utr_to_cds_ratio=10.0,
    )
    defaults.update(kw)
    return UtrConfig(**defaults)


def _total_bp(rows):
    return sum(r["End"] - r["Start"] for r in rows)


# ---------------------------------------------------------------------------
# trim_utrs -- direction correction (Bug A)
# ---------------------------------------------------------------------------


class TestHardCapDirection:
    """_hard_cap must keep the CDS-proximal portion and trim the distal end."""

    def test_5p_plus_strand_keeps_proximal(self):
        """5' UTR at LOW coords (+ strand): proximal = HIGH end.  After cap,
        the kept rows must be the HIGH portion (adjacent to the CDS)."""
        # Exons:    [1000,2000]  [2500,4000]  [5000,8000]
        # CDS:                               [6000,8000]
        # 5' UTR from derive_utrs: [1000,2000] + [2500,4000] + [5000,6000] = 4500 bp
        cds = [_row("CDS", 6000, 8000)]
        exon_rows = [
            _row("exon", 1000, 2000),
            _row("exon", 2500, 4000),
            _row("exon", 5000, 8000),
        ]
        utr_5p = [
            _row("five_prime_UTR", 1000, 2000),
            _row("five_prime_UTR", 2500, 4000),
            _row("five_prime_UTR", 5000, 6000),
        ]
        cfg = _utr_cfg(max_5p_bp=3000)

        out_5p, out_3p, was_trimmed = trim_utrs(utr_5p, [], cds, exon_rows, cfg)

        assert was_trimmed
        assert _total_bp(out_5p) == 3000
        # The proximal portion (HIGH end) must touch the CDS: max UTR coord == CDS start
        assert max(r["End"] for r in out_5p) == 6000

    def test_3p_minus_strand_keeps_proximal(self):
        """3' UTR at LOW coords (- strand): proximal = HIGH end.  After cap,
        the kept rows must be the HIGH portion (adjacent to the CDS).

        On - strand the CDS sits at HIGH genomic coords; the 3' UTR is at LOW
        genomic coords.  Proximal end of the 3' UTR is the HIGH end (nearest
        to cds_min).
        """
        # CDS at HIGH coords [9000,12000]; 3' UTR at LOW coords
        cds = [_row("CDS", 9000, 12000)]
        exon_rows = [
            _row("exon", 1000, 4000),
            _row("exon", 4500, 6000),
            _row("exon", 7000, 12000),
        ]
        utr_3p = [
            _row("three_prime_UTR", 1000, 4000),
            _row("three_prime_UTR", 4500, 6000),
            _row("three_prime_UTR", 7000, 9000),
        ]
        cfg = _utr_cfg(max_3p_bp=3000)

        out_5p, out_3p, was_trimmed = trim_utrs([], utr_3p, cds, exon_rows, cfg)

        assert was_trimmed
        assert _total_bp(out_3p) == 3000
        # Proximal end of 3' UTR (for - strand at LOW coords) is the HIGH end
        assert max(r["End"] for r in out_3p) == 9000

    def test_3p_plus_strand_keeps_proximal(self):
        """3' UTR at HIGH coords (+ strand): proximal = LOW end.  Original
        behaviour was already correct; must remain so after the fix."""
        # Exons:  [1000,3000]  [5000,12000]
        # CDS:    [1000,3000]  [5000,6000]
        # 3' UTR:              [6000,12000] = 6000 bp -> cap to 3000
        cds = [_row("CDS", 1000, 3000), _row("CDS", 5000, 6000)]
        exon_rows = [_row("exon", 1000, 3000), _row("exon", 5000, 12000)]
        utr_3p = [_row("three_prime_UTR", 6000, 12000)]
        cfg = _utr_cfg(max_3p_bp=3000)

        out_5p, out_3p, was_trimmed = trim_utrs([], utr_3p, cds, exon_rows, cfg)

        assert was_trimmed
        assert _total_bp(out_3p) == 3000
        # Proximal end at LOW coords of 3' UTR == CDS end
        cds_max = max(r["End"] for r in cds)
        assert min(r["Start"] for r in out_3p) == cds_max

    def test_5p_minus_strand_keeps_proximal(self):
        """5' UTR at HIGH coords (- strand): proximal = LOW end.  Original
        behaviour was already correct; must remain so.

        The CDS is at LOW coords [1000,5000] and the 5' UTR is at HIGH coords
        [8000,14000].  The proximal end is the LOW end of the UTR (8000) --
        the end closest to the CDS (5000).  After the cap the kept 3000 bp
        must include position 8000, not 11000.
        """
        cds = [_row("CDS", 1000, 5000)]
        exon_rows = [_row("exon", 1000, 5000), _row("exon", 8000, 14000)]
        utr_5p = [_row("five_prime_UTR", 8000, 14000)]
        cfg = _utr_cfg(max_5p_bp=3000)

        out_5p, out_3p, was_trimmed = trim_utrs(utr_5p, [], cds, exon_rows, cfg)

        assert was_trimmed
        assert _total_bp(out_5p) == 3000
        # Kept from the LOW (proximal) end: must start at 8000 and end at 11000
        assert min(r["Start"] for r in out_5p) == 8000
        assert max(r["End"] for r in out_5p) == 11000

    def test_no_gap_after_trim_both_ends(self):
        """Both UTRs trimmed: the proximal (CDS-adjacent) portion is retained.

        The 5' UTR [0,5000] is at LOW coords (+ strand); proximal = HIGH end.
        The 3' UTR [11000,16000] is at HIGH coords; proximal = LOW end.
        Both UTRs are separated from the CDS exon by introns, so the cap
        preserves the innermost portions.
        """
        cds = [_row("CDS", 6000, 10000)]
        exon_rows = [
            _row("exon", 0, 5000),
            _row("exon", 6000, 10000),
            _row("exon", 11000, 16000),
        ]
        utr_5p = [_row("five_prime_UTR", 0, 5000)]
        utr_3p = [_row("three_prime_UTR", 11000, 16000)]
        cfg = _utr_cfg(max_5p_bp=3000, max_3p_bp=3000)

        out_5p, out_3p, was_trimmed = trim_utrs(utr_5p, utr_3p, cds, exon_rows, cfg)

        assert was_trimmed
        assert _total_bp(out_5p) == 3000
        assert _total_bp(out_3p) == 3000
        # 5' proximal = HIGH end of [0,5000]; kept portion ends at 5000 (exon boundary)
        assert max(r["End"] for r in out_5p) == 5000
        # 3' proximal = LOW end of [11000,16000]; kept portion starts at 11000
        assert min(r["Start"] for r in out_3p) == 11000


# ---------------------------------------------------------------------------
# validate_and_fix_gff3 -- exon synchronisation (Bug B)
# ---------------------------------------------------------------------------


def _make_config(max_5p=3000, max_3p=3000):
    cfg = PipelineConfig()
    cfg.utr = UtrConfig(
        trim_policy="hard_cap",
        max_5p_bp=max_5p,
        max_3p_bp=max_3p,
        max_total_bp=max_5p + max_3p,
        max_utr_to_cds_ratio=50.0,  # high so ratio cap doesn't interfere
    )
    cfg.validation = ValidationConfig(
        mode="drop_transcript",
        feature_outside_exons_policy="trim",
        log_violations=False,
    )
    return cfg


def _build_gff3(gene_id, tid, strand, exons, cds, utr_5p=None, utr_3p=None):
    """Construct a minimal GFF3 row list for one gene/transcript."""
    gene_s = min(s for s, e in exons)
    gene_e = max(e for s, e in exons)
    rows = [
        {
            "Feature": "gene",
            "Chromosome": "1",
            "Start": gene_s,
            "End": gene_e,
            "Strand": strand,
            "ID": gene_id,
            "Parent": "",
            "Source": "GMB",
            "Score": ".",
            "Frame": ".",
        },
        {
            "Feature": "mRNA",
            "Chromosome": "1",
            "Start": gene_s,
            "End": gene_e,
            "Strand": strand,
            "ID": tid,
            "Parent": gene_id,
            "Source": "GMB",
            "Score": ".",
            "Frame": ".",
        },
    ]
    for i, (s, e) in enumerate(sorted(exons)):
        rows.append(
            {
                "Feature": "exon",
                "Chromosome": "1",
                "Start": s,
                "End": e,
                "Strand": strand,
                "ID": f"{tid}.exon{i+1}",
                "Parent": tid,
                "Source": "GMB",
                "Score": ".",
                "Frame": ".",
            }
        )
    for i, (s, e) in enumerate(sorted(cds)):
        rows.append(
            {
                "Feature": "CDS",
                "Chromosome": "1",
                "Start": s,
                "End": e,
                "Strand": strand,
                "ID": f"{tid}.cds{i+1}",
                "Parent": tid,
                "Source": "GMB",
                "Score": ".",
                "Frame": "0",
            }
        )
    for i, (s, e) in enumerate(sorted(utr_5p or [])):
        rows.append(
            {
                "Feature": "five_prime_UTR",
                "Chromosome": "1",
                "Start": s,
                "End": e,
                "Strand": strand,
                "ID": f"{tid}.5utr{i+1}",
                "Parent": tid,
                "Source": "GMB",
                "Score": ".",
                "Frame": ".",
            }
        )
    for i, (s, e) in enumerate(sorted(utr_3p or [])):
        rows.append(
            {
                "Feature": "three_prime_UTR",
                "Chromosome": "1",
                "Start": s,
                "End": e,
                "Strand": strand,
                "ID": f"{tid}.3utr{i+1}",
                "Parent": tid,
                "Source": "GMB",
                "Score": ".",
                "Frame": ".",
            }
        )
    return rows


def _exon_union(out_rows, tid):
    exons = [r for r in out_rows if r["Feature"] == "exon" and r.get("Parent") == tid]
    if not exons:
        return set()
    intervals = sorted((r["Start"], r["End"]) for r in exons)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return {tuple(iv) for iv in merged}


def _feature_union(out_rows, tid, features):
    rows = [r for r in out_rows if r["Feature"] in features and r.get("Parent") == tid]
    if not rows:
        return set()
    intervals = sorted((r["Start"], r["End"]) for r in rows)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return {tuple(iv) for iv in merged}


class TestExonInvariant:
    """After validate_and_fix_gff3, union(exons) must equal union(CDS+UTRs)."""

    def test_at_cap_plus_strand_5p_gap_resolved(self):
        """+ strand: 5' UTR too long; after cap the distal exon must be trimmed."""
        # Exons:  [0,4000]  [5000,8000]  [9000,12000]
        # CDS:                           [9000,12000]
        # 5' UTR: [0,4000]  [5000,8000]  -- total 7000 bp; cap at 3000
        # After cap (keep HIGH = proximal): [5000,8000] = 3000 bp
        # Exon [0,4000] is no longer covered -> removed by exon rebuild
        tid = "t1"
        gff = _build_gff3(
            "g1",
            tid,
            "+",
            exons=[(0, 4000), (5000, 8000), (9000, 12000)],
            cds=[(9000, 12000)],
            utr_5p=[(0, 4000), (5000, 8000)],
        )
        cfg = _make_config(max_5p=3000)
        out, stats = validate_and_fix_gff3(gff, cfg)

        eu = _exon_union(out, tid)
        fu = _feature_union(out, tid, {"CDS", "five_prime_UTR", "three_prime_UTR"})
        assert eu == fu, f"exon union {eu} != feature union {fu}"

        # After cap, the retained UTR occupies [5000,8000]
        utr_rows = [r for r in out if r["Feature"] == "five_prime_UTR" and r.get("Parent") == tid]
        if utr_rows:
            assert max(r["End"] for r in utr_rows) == 8000

    def test_at_cap_minus_strand_3p_gap_resolved(self):
        """- strand: 3' UTR at LOW coords too long; after cap the distal exon is trimmed."""
        # CDS at HIGH coords [9000,12000] (typical for - strand)
        # 3' UTR at LOW coords [0,4000] + [5000,8000] = 7000 bp; cap at 3000
        # After cap (keep HIGH = proximal): [5000,8000] = 3000 bp
        # Exon [0,4000] is no longer covered -> removed
        tid = "t2"
        gff = _build_gff3(
            "g2",
            tid,
            "-",
            exons=[(0, 4000), (5000, 8000), (9000, 12000)],
            cds=[(9000, 12000)],
            utr_3p=[(0, 4000), (5000, 8000)],
        )
        cfg = _make_config(max_3p=3000)
        out, stats = validate_and_fix_gff3(gff, cfg)

        eu = _exon_union(out, tid)
        fu = _feature_union(out, tid, {"CDS", "five_prime_UTR", "three_prime_UTR"})
        assert eu == fu, f"exon union {eu} != feature union {fu}"

        utr_rows = [r for r in out if r["Feature"] == "three_prime_UTR" and r.get("Parent") == tid]
        if utr_rows:
            assert max(r["End"] for r in utr_rows) == 8000

    def test_no_gap_no_exon_change(self):
        """When no UTR trimming occurs exon_rows must be unchanged."""
        # Exons: [0,1000]  [2000,5000]
        # CDS:             [2000,5000]
        # 5' UTR: [0,1000] = 1000 bp (well below cap of 3000)
        tid = "t3"
        gff = _build_gff3(
            "g3",
            tid,
            "+",
            exons=[(0, 1000), (2000, 5000)],
            cds=[(2000, 5000)],
            utr_5p=[(0, 1000)],
        )
        cfg = _make_config(max_5p=3000)
        out, stats = validate_and_fix_gff3(gff, cfg)

        eu = _exon_union(out, tid)
        fu = _feature_union(out, tid, {"CDS", "five_prime_UTR", "three_prime_UTR"})
        assert eu == fu

    def test_ratio_cap_triggers_same_fix(self):
        """When max_utr_to_cds_ratio triggers a secondary cap, invariant still holds."""
        # CDS: 500 bp.  max_utr_to_cds_ratio=2.0 -> max total UTR = 1000 bp.
        # 3' UTR: [1000,5000] = 4000 bp (above hard cap of 3000 too).
        tid = "t4"
        gff = _build_gff3(
            "g4",
            tid,
            "+",
            exons=[(0, 500), (1000, 5000)],
            cds=[(0, 500)],
            utr_3p=[(1000, 5000)],
        )
        cfg = _make_config(max_3p=3000)
        cfg.utr.max_utr_to_cds_ratio = 2.0
        out, stats = validate_and_fix_gff3(gff, cfg)

        eu = _exon_union(out, tid)
        fu = _feature_union(out, tid, {"CDS", "five_prime_UTR", "three_prime_UTR"})
        assert eu == fu

    def test_multi_exon_trim_preserves_introns(self):
        """After trimming, exon boundaries must respect the original intron structure."""
        # + strand.  3' UTR spans three exons; cap keeps only the first 3000 bp from CDS.
        # Exons: [0,2000]  [3000,8000]  [9000,12000]  [13000,16000]
        # CDS:   [0,2000]
        # 3' UTR:[3000,8000] [9000,12000] [13000,16000] = 9000 bp -> cap to 3000
        # After cap (proximal = LOW end of 3' UTR): [3000,6000] = 3000 bp
        # Exons should become: [0,2000] [3000,6000]
        tid = "t5"
        gff = _build_gff3(
            "g5",
            tid,
            "+",
            exons=[(0, 2000), (3000, 8000), (9000, 12000), (13000, 16000)],
            cds=[(0, 2000)],
            utr_3p=[(3000, 8000), (9000, 12000), (13000, 16000)],
        )
        cfg = _make_config(max_3p=3000)
        out, stats = validate_and_fix_gff3(gff, cfg)

        eu = _exon_union(out, tid)
        fu = _feature_union(out, tid, {"CDS", "five_prime_UTR", "three_prime_UTR"})
        assert eu == fu

        # Exon intervals should be exactly {(0,2000), (3000,6000)}
        exon_intervals = sorted(
            (r["Start"], r["End"])
            for r in out
            if r["Feature"] == "exon" and r.get("Parent") == tid
        )
        assert exon_intervals == [(0, 2000), (3000, 6000)]
        # Intron between (2000,3000) is preserved (no exon spans it)
        for s, e in exon_intervals:
            assert not (s < 2000 and e > 3000)
