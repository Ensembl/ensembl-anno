"""Unit tests for STAR junction loading, indexing and diagnostics (synthetic data)."""

import gzip
import os

import numpy as np
import pytest

from gmb.pipeline.junction_diagnostics import (
    JunctionAnnotator,
    ModelSpec,
    genome_motif,
    introns_of,
    read_models,
)
from gmb.pipeline.star_junctions import (
    FAILED_EMPTY,
    LOADED,
    MALFORMED,
    MISSING_FILE,
    NO_KNOWN_SEQIDS,
    NOT_SELECTED,
    JunctionIndex,
    JunctionManifestError,
    parse_sj_file,
    preflight_checks,
    read_manifest,
)

HEADER = [
    "library_id",
    "path",
    "role",
    "study",
    "strain",
    "condition",
    "count_unit",
    "max_intron",
    "min_intron",
    "noncanonical_filtered",
]


def sj(*rows):
    """STAR rows: (chrom, s1, e1, strand, motif, annot, unique, multi, overhang)."""
    return "".join("\t".join(str(x) for x in r) + "\n" for r in rows)


def write(path, text, gz=False):
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        with open(path, "w") as fh:
            fh.write(text)
    return str(path)


def manifest(tmp_path, libs, header=HEADER, name="manifest.tsv"):
    """libs: list of dicts; missing columns default."""
    defaults = dict(
        role="construction",
        study="S1",
        strain="X",
        condition="c1",
        count_unit="mate",
        max_intron="1000",
        min_intron="21",
        noncanonical_filtered="true",
    )
    lines = ["\t".join(header)]
    for lib in libs:
        row = {**defaults, **lib}
        lines.append("\t".join(str(row.get(c, "")) for c in header))
    return write(tmp_path / name, "\n".join(lines) + "\n")


# ----------------------------------------------------------------- parsing
class TestParse:
    def test_coordinates_counts_and_strand(self, tmp_path):
        p = write(
            tmp_path / "a.sj.tab",
            sj(("1", 201, 300, 1, 1, 1, 7, 2, 40), ("1", 501, 650, 2, 2, 0, 0, 3, 12)),
        )
        status, msg, df, unknown = parse_sj_file(p, {"1": 1000})
        assert status == LOADED and unknown == 0
        cols = ["start0", "end", "strand", "motif", "unique", "multi", "overhang"]
        # STAR 1-based inclusive 201-300 -> half-open 0-based [200, 300)
        assert [int(df[c].iloc[0]) for c in cols] == [200, 300, 1, 1, 7, 2, 40]
        assert [int(df[c].iloc[1]) for c in ("strand", "unique", "multi")] == [2, 0, 3]

    def test_star_coordinates_match_gmb_exon_introns(self, tmp_path):
        """GTF exons 101-200 and 301-400 (1-based) bound the STAR intron 201-300."""
        gtf = write(
            tmp_path / "t.gtf",
            '1\tX\tCDS\t101\t200\t.\t+\t0\ttranscript_id "t1"; gene_id "g1";\n'
            '1\tX\tCDS\t301\t400\t.\t+\t0\ttranscript_id "t1"; gene_id "g1";\n',
        )
        (m,) = read_models(gtf)
        assert introns_of(m.segments) == [(200, 300)]
        p = write(tmp_path / "a.sj.tab", sj(("1", 201, 300, 1, 1, 1, 5, 0, 30)))
        _, _, df, _ = parse_sj_file(p)
        assert (int(df.start0[0]), int(df.end[0])) == (200, 300)

    def test_pyranges_intron_convention(self, tmp_path):
        """read_models and the builder's pyranges loading give the same intron."""
        pr = pytest.importorskip("pyranges")
        gtf = write(
            tmp_path / "t.gtf",
            '1\tX\texon\t101\t200\t.\t+\t.\ttranscript_id "t1"; gene_id "g1";\n'
            '1\tX\texon\t301\t400\t.\t+\t.\ttranscript_id "t1"; gene_id "g1";\n',
        )
        df = pr.read_gtf(gtf).df.sort_values("Start")
        builder_intron = (int(df.End.iloc[0]), int(df.Start.iloc[1]))
        (m,) = read_models(gtf, feature="exon")
        assert introns_of(m.segments) == [builder_intron] == [(200, 300)]

    def test_gzip_detected_by_magic_not_suffix(self, tmp_path):
        p = write(tmp_path / "a.sj.tab", sj(("1", 201, 300, 1, 1, 1, 5, 0, 30)), gz=True)
        assert parse_sj_file(p)[0] == LOADED

    @pytest.mark.parametrize(
        "text,fragment",
        [
            ("1\t201\t300\t1\t1\t1\t5\t0\n", "columns"),
            ("1\t201\t300\t1\t1\t1\tfive\t0\t30\n", "non-integer"),
            ("1\t201\t300\t3\t1\t1\t5\t0\t30\n", "out of range"),
            ("1\t301\t300\t1\t1\t1\t5\t0\t30\n", "out of range"),
            ("1\t201\t300\t1\t9\t1\t5\t0\t30\n", "out of range"),
            ("1\t201\t300\t1\t1\t1\t5\t0\t30", "newline"),
            ("1\t201\t300\t1\t1\t1\t5\t0\t30\n1\t201\t300\t1\t1\t1\t6\t0\t30\n", "duplicate"),
            ("1\t201\t300\t1\t1\t1\t5\t0\t30\t9\n", ""),
        ],
    )
    def test_malformed_rejects_whole_file(self, tmp_path, text, fragment):
        status, msg, df, _ = parse_sj_file(write(tmp_path / "bad.sj.tab", text))
        assert status == MALFORMED and df.empty
        assert fragment in msg

    def test_zero_byte_is_failed_not_zero_junctions(self, tmp_path):
        status, msg, df, _ = parse_sj_file(write(tmp_path / "e.sj.tab", ""))
        assert status == FAILED_EMPTY and "not zero junctions" in msg

    def test_missing_file(self, tmp_path):
        assert parse_sj_file(str(tmp_path / "nope.sj.tab"))[0] == MISSING_FILE

    def test_unknown_seqids_skipped_and_counted(self, tmp_path):
        p = write(
            tmp_path / "a.sj.tab",
            sj(("1", 201, 300, 1, 1, 1, 5, 0, 30), ("chrUn", 201, 300, 1, 1, 1, 5, 0, 30)),
        )
        status, _, df, unknown = parse_sj_file(p, {"1": 1000})
        assert status == LOADED and len(df) == 1 and unknown == 1

    def test_only_unknown_seqids(self, tmp_path):
        p = write(tmp_path / "a.sj.tab", sj(("chrUn", 201, 300, 1, 1, 1, 5, 0, 30)))
        assert parse_sj_file(p, {"1": 1000})[0] == NO_KNOWN_SEQIDS

    def test_beyond_sequence_end(self, tmp_path):
        p = write(tmp_path / "a.sj.tab", sj(("1", 201, 3000, 1, 1, 1, 5, 0, 30)))
        assert parse_sj_file(p, {"1": 1000})[0] == MALFORMED


# ---------------------------------------------------------------- manifest
class TestManifest:
    def test_required_columns(self, tmp_path):
        m = write(tmp_path / "m.tsv", "library_id\trole\nL1\tconstruction\n")
        with pytest.raises(JunctionManifestError, match="path"):
            read_manifest(m)

    def test_duplicate_library(self, tmp_path):
        m = manifest(tmp_path, [dict(library_id="L1", path="a"), dict(library_id="L1", path="b")])
        with pytest.raises(JunctionManifestError, match="duplicate"):
            read_manifest(m)

    def test_no_role_column_needs_explicit_all(self, tmp_path):
        write(tmp_path / "a.sj.tab", sj(("1", 201, 300, 1, 1, 1, 5, 0, 30)))
        m = write(tmp_path / "m.tsv", "library_id\tpath\nL1\ta.sj.tab\n")
        with pytest.raises(JunctionManifestError, match="role"):
            read_manifest(m)
        (lib,) = read_manifest(m, roles=("all",))
        assert lib.status == "" and lib.path == str(tmp_path / "a.sj.tab")  # relative to manifest

    def test_group_defaults(self, tmp_path):
        m = manifest(
            tmp_path,
            [
                dict(library_id="L1", path="a"),
                dict(library_id="L2", path="b", study="", strain="", condition=""),
            ],
        )
        a, b = read_manifest(m)
        assert a.group == "S1|X|c1" and b.group == "L2" and b.study == "L2"

    def test_bad_integer(self, tmp_path):
        m = manifest(tmp_path, [dict(library_id="L1", path="a", max_intron="lots")])
        with pytest.raises(JunctionManifestError, match="max_intron"):
            read_manifest(m)

    def test_holdout_never_read(self, tmp_path):
        """A held-out library is not opened: even a malformed file stays NOT_SELECTED."""
        write(tmp_path / "c.sj.tab", sj(("1", 201, 300, 1, 1, 1, 5, 0, 30)))
        write(tmp_path / "h.sj.tab", "garbage\n")
        m = manifest(
            tmp_path,
            [
                dict(library_id="C", path="c.sj.tab"),
                dict(library_id="H", path="h.sj.tab", role="holdout"),
            ],
        )
        ix = JunctionIndex.from_manifest(m)
        status = {lib.library_id: lib.status for lib in ix.libraries}
        assert status == {"C": LOADED, "H": NOT_SELECTED}
        assert [lib.library_id for lib in ix.loaded] == ["C"]


# ------------------------------------------------------------- aggregation
@pytest.fixture
def three_libs(tmp_path):
    # L1, L2 same condition (group A); L3 another study (group B).
    write(
        tmp_path / "l1.sj.tab",
        sj(
            ("1", 201, 300, 1, 1, 1, 10, 1, 40),
            ("1", 201, 320, 1, 1, 1, 2, 0, 20),
            ("1", 901, 950, 1, 1, 1, 0, 4, 10),
        ),
    )
    write(
        tmp_path / "l2.sj.tab",
        sj(("1", 201, 300, 1, 1, 1, 5, 0, 50), ("2", 101, 200, 2, 2, 1, 3, 0, 30)),
    )
    write(
        tmp_path / "l3.sj.tab",
        sj(("1", 201, 300, 1, 1, 1, 1, 0, 25), ("1", 1601, 1700, 0, 0, 0, 2, 0, 30)),
    )
    write(tmp_path / "e.sj.tab", "")
    libs = [
        dict(library_id="L1", path="l1.sj.tab"),
        dict(library_id="L2", path="l2.sj.tab"),
        dict(library_id="L3", path="l3.sj.tab", study="S2", condition="c2"),
        dict(library_id="E", path="e.sj.tab"),
    ]
    return manifest(tmp_path, libs)


class TestIndex:
    def test_aggregates(self, three_libs):
        ix = JunctionIndex.from_manifest(three_libs, seq_lengths={"1": 5000, "2": 5000})
        s = ix.lookup("1", 200, 300)
        assert (s.unique_reads, s.multi_reads, s.n_libraries, s.n_libraries_any) == (16, 1, 3, 3)
        assert (s.n_groups, s.n_studies, s.max_overhang, s.strand, s.motif) == (
            2,
            2,
            50,
            "+",
            "GT/AG",
        )
        m = ix.lookup("1", 900, 950)
        assert (m.unique_reads, m.multi_reads, m.n_libraries, m.n_libraries_any) == (0, 4, 0, 1)
        assert ix.lookup("2", 100, 200).strand == "-"
        assert ix.lookup("1", 100, 300) is None
        assert sorted(ix.library_counts(ix.find("1", 200, 300))) == [
            ("L1", 10, 1),
            ("L2", 5, 0),
            ("L3", 1, 0),
        ]
        status = {lib.library_id: lib.status for lib in ix.libraries}
        assert status["E"] == FAILED_EMPTY and ix.stats["libraries_loaded"] == 3

    def test_deterministic_and_order_independent(self, tmp_path, three_libs):
        a = JunctionIndex.from_manifest(three_libs)
        lines = open(three_libs).read().splitlines()
        rev = write(tmp_path / "rev.tsv", "\n".join([lines[0]] + lines[1:][::-1]) + "\n")
        b = JunctionIndex.from_manifest(rev)
        for attr in ("j_start", "j_end", "j_unique", "j_multi", "j_nlib", "j_ngroup", "j_nstudy"):
            assert np.array_equal(getattr(a, attr), getattr(b, attr))

    def test_site_partners_plus_and_minus(self, tmp_path):
        write(
            tmp_path / "l.sj.tab",
            sj(
                ("1", 201, 300, 1, 1, 1, 10, 0, 40),  # + intron
                ("1", 201, 320, 1, 1, 1, 3, 0, 40),  # + same donor (start), other acceptor
                ("1", 181, 300, 1, 1, 1, 2, 0, 40),  # + other donor, same acceptor
                ("1", 501, 600, 2, 2, 1, 9, 0, 40),  # - intron: donor at 600, acceptor at 501
                ("1", 504, 600, 2, 2, 1, 4, 0, 40),  # - shares donor (high end)
                ("1", 801, 900, 0, 0, 0, 1, 0, 40),
            ),
        )  # unstranded (non-canonical)
        ix = JunctionIndex.from_manifest(
            manifest(tmp_path, [dict(library_id="L", path="l.sj.tab")])
        )
        p = ix.site_partners("1", "+", 200, 300)
        assert [(ix.support(j).start1, ix.support(j).end) for j in p["donor"]] == [(201, 320)]
        assert [(ix.support(j).start1, ix.support(j).end) for j in p["acceptor"]] == [(181, 300)]
        q = ix.site_partners("1", "-", 500, 600)
        assert [(ix.support(j).start1, ix.support(j).end) for j in q["donor"]] == [(504, 600)]
        assert q["acceptor"] == []
        # an unstranded junction is offered to queries on either strand
        assert (
            ix.site_partners("1", "+", 800, 950)["donor"]
            == ix.site_partners("1", "-", 700, 900)["donor"]
        )

    def test_overlapping(self, three_libs):
        ix = JunctionIndex.from_manifest(three_libs)
        got = sorted(
            (ix.support(j).start1, ix.support(j).end) for j in ix.overlapping("1", 250, 260)
        )
        assert got == [(201, 300), (201, 320)]
        assert ix.overlapping("1", 2000, 3000) == [] and ix.overlapping("9", 0, 10) == []


# ------------------------------------------------------------- diagnostics
GENOME = {
    "1": ("A" * 200 + "GT" + "A" * 96 + "AG" + "A" * 300 + "CC" + "A" * 96 + "TT" + "A" * 2400)
}


@pytest.fixture
def annot(tmp_path):
    write(
        tmp_path / "l1.sj.tab",
        sj(
            ("1", 201, 300, 1, 1, 1, 10, 0, 40),
            ("1", 1201, 1300, 1, 1, 1, 0, 3, 20),  # multimap only
            ("1", 1401, 1500, 2, 2, 1, 6, 0, 20),  # minus-strand motif
            ("1", 1601, 1700, 1, 1, 1, 8, 0, 20),  # alternative to 1601-1720
            ("1", 1801, 1900, 1, 1, 1, 2, 0, 20),  # minor vs 1801-1950
            ("1", 1801, 1950, 1, 1, 1, 9, 0, 20),
            ("1", 2101, 2200, 1, 1, 1, 5, 0, 20),
        ),
    )  # somewhere inside a later model
    write(
        tmp_path / "l2.sj.tab",
        sj(("1", 201, 300, 1, 1, 1, 4, 0, 40), ("1", 1801, 1900, 1, 1, 1, 1, 0, 20)),
    )
    m = manifest(
        tmp_path,
        [
            dict(library_id="L1", path="l1.sj.tab"),
            dict(library_id="L2", path="l2.sj.tab", condition="c2"),
        ],
    )
    return JunctionAnnotator(JunctionIndex.from_manifest(m))  # no genome: motif unknown


def one(annot, start1, end, strand="+"):
    return annot.intron("1", strand, start1 - 1, end, span=(start1 - 50, end + 50))


class TestDiagnostics:
    def test_supported(self, annot):
        annot.genome = GENOME
        r = one(annot, 201, 300)
        assert r["status"] == "SUPPORTED" and r["evidence_class"] == "supported"
        assert (r["unique_reads"], r["n_libraries"], r["n_groups"]) == (14, 2, 2)
        assert r["genome_motif"] == "GT/AG" and r["donor_usage"] == 1.0

    def test_multimap_only(self, annot):
        assert one(annot, 1201, 1300)["status"] == "MULTIMAP_ONLY"

    def test_opposite_strand_is_contradiction(self, annot):
        r = one(annot, 1401, 1500, "+")
        assert r["status"] == "OPPOSITE_STRAND" and r["evidence_class"] == "contradicted"
        assert one(annot, 1401, 1500, "-")["status"] == "SUPPORTED"

    def test_alternative_site_is_contradiction(self, annot):
        r = one(annot, 1601, 1720)
        assert r["status"] == "ALTERNATIVE_SITE_OBSERVED" and r["evidence_class"] == "contradicted"
        assert r["best_alternative"] == "1601-1700" and r["best_alternative_side"] == "donor"
        assert r["donor_usage"] == 0.0

    def test_minor_site_and_usage(self, annot):
        r = one(annot, 1801, 1900)
        assert r["status"] == "SUPPORTED" and r["site_conflict"] == "minor_site"
        assert r["donor_usage"] == pytest.approx(3 / 12)
        r2 = one(annot, 1801, 1950)
        assert r2["site_conflict"] == "" and r2["donor_usage"] == pytest.approx(9 / 12)

    def test_long_intron_unobservable_not_penalised(self, annot):
        r = annot.intron("1", "+", 3000, 4500)
        assert r["status"] == "UNOBSERVABLE_LONG" and r["evidence_class"] == "unobservable"
        rows, s = annot.model(ModelSpec("m", "t", "1", "+", [(2900, 3000), (4500, 4600)]))
        assert s["n_contradicted"] == 0 and s["model_class"] == "no_junction_evidence"

    def test_short_intron_unobservable(self, annot):
        assert annot.intron("1", "+", 3000, 3010)["status"] == "UNOBSERVABLE_SHORT"

    def test_noncanonical_unobservable(self, annot):
        annot.genome = GENOME
        r = annot.intron("1", "+", 600, 700)  # CC...TT in GENOME
        assert r["genome_motif"] == "CC/TT" and r["status"] == "UNOBSERVABLE_NONCANONICAL"

    def test_observed_evidence_outranks_noncanonical(self, annot):
        """With a genome, these introns read non-canonical, yet the observed opposite
        strand / competing site is still reported as contradiction."""
        annot.genome = GENOME
        assert one(annot, 1401, 1500, "+")["status"] == "OPPOSITE_STRAND"
        r = one(annot, 1601, 1720)
        assert r["status"] == "ALTERNATIVE_SITE_OBSERVED" and "non-canonical" in r["reasons"]

    def test_unobserved_expressed_vs_no_evidence(self, annot):
        assert (
            annot.intron("1", "+", 2300, 2400, span=(2000, 2500))["status"]
            == "UNOBSERVED_EXPRESSED_LOCUS"
        )
        r = annot.intron("1", "+", 2700, 2800, span=(2650, 2850))
        assert r["status"] == "UNOBSERVED_NO_LOCAL_EVIDENCE" and r["evidence_class"] == "unknown"

    def test_model_summary(self, annot):
        _, s = annot.model(ModelSpec("m1", "Helixer", "1", "+", [(100, 200), (300, 400)]))
        assert s["model_class"] == "all_supported" and s["libraries_supporting_all_introns"] == 2
        assert s["groups_supporting_all_introns"] == 2
        rows, s = annot.model(
            ModelSpec("m2", "x", "1", "+", [(1500, 1600), (1720, 1800), (1900, 2000)])
        )
        assert [r["status"] for r in rows] == ["ALTERNATIVE_SITE_OBSERVED", "SUPPORTED"]
        assert s["model_class"] == "contradicted" and s["n_minor_site"] == 1
        _, s = annot.model(ModelSpec("m3", "x", "1", "+", [(100, 500)]))
        assert s["model_class"] == "no_cds_introns" and s["n_cds_introns"] == 0

    def test_minus_strand_numbering_5_to_3(self, annot):
        rows, _ = annot.model(
            ModelSpec("m", "x", "1", "-", [(1000, 1400), (1500, 1600), (1700, 1800)])
        )
        assert [r["intron_number"] for r in rows] == [1, 2]
        assert rows[0]["intron_start"] == 1601  # 5'-most intron of a minus-strand model

    def test_genome_motif_minus(self):
        seq = "A" * 10 + "CT" + "A" * 10 + "AC" + "A" * 10
        assert genome_motif(seq, 10, 24, "-") == "GT/AG"

    def test_read_models_gff3(self, tmp_path):
        g = write(
            tmp_path / "a.gff3",
            "##gff-version 3\n1\tH\tgene\t101\t400\t.\t+\t.\tID=g1\n"
            "1\tH\tmRNA\t101\t400\t.\t+\t.\tID=t1;Parent=g1\n"
            "1\tH\tCDS\t101\t200\t.\t+\t0\tID=c1;Parent=t1\n"
            "1\tH\tCDS\t301\t400\t.\t+\t2\tID=c2;Parent=t1\n",
        )
        (m,) = read_models(g, source="H")
        assert (m.model_id, m.gene_id, m.segments, m.source) == (
            "t1",
            "g1",
            [(100, 200), (300, 400)],
            "H",
        )


class TestPreflight:
    def test_verdicts(self, tmp_path, three_libs):
        checks = {
            c[0]: c[1]
            for c in preflight_checks(three_libs, ("construction",), {"1": 5000, "2": 5000})
        }
        assert checks["star_junction_libraries"] == "WARN"  # the empty library
        assert checks["star_junction_count_unit"] == "WARN"
        assert preflight_checks(three_libs, ("nonexistent",), None)[0][1] == "FAIL"
        assert (
            preflight_checks(str(tmp_path / "missing.tsv"), ("construction",), None)[0][1]
            == "FAIL"
        )
