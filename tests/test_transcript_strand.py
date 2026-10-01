"""Tests for Minimap2 transcript-strand resolution and SAM/BAM -> BED12 conversion."""

import shutil
import tempfile
import unittest
from pathlib import Path

from src.python.ensembl.tools.anno.transcriptomic_annotation.transcript_strand import (
    MissingTranscriptStrandError,
    MissingTsPolicy,
    transcript_strand,
)
from src.python.ensembl.tools.anno.transcriptomic_annotation.splice_to_bed import (
    iter_alignments,
    splice_to_bed,
)
from src.python.ensembl.tools.anno.transcriptomic_annotation.minimap import _bed_to_gtf

DATA = Path(__file__).parent / "data"
FIXTURE_SAM = DATA / "minimap_ts_fixture.sam"
FIXTURE_BAM = DATA / "minimap_ts_fixture.bam"

#: Transcript strand each fixture read must end up with.
EXPECTED_STRAND = {
    "read_fwd_ts_plus": "+",
    "read_rev_ts_plus": "-",
    "read_fwd_ts_minus": "-",
    "read_rev_ts_minus": "+",
    "read_no_ts": ".",
    "read_bad_ts": ".",
    "read_two_primary": "+",
}


def _bed_rows(path):
    """Read a BED12 file into a list of field lists."""
    return [line.rstrip("\n").split("\t") for line in path.read_text(encoding="utf8").splitlines()]


def _by_name(path):
    """Map BED name -> list of rows for that name."""
    out = {}
    for row in _bed_rows(path):
        out.setdefault(row[3], []).append(row)
    return out


class TestTranscriptStrandHelper(unittest.TestCase):
    """The strand rule itself: ts combined with alignment orientation."""

    def test_forward_alignment_ts_plus(self):
        """Case 1: forward alignment, ts:+ -> transcript on the forward strand."""
        self.assertEqual(transcript_strand("+", is_reverse=False), "+")

    def test_reverse_alignment_ts_plus(self):
        """Case 2: reverse alignment, ts:+ -> transcript on the reverse strand."""
        self.assertEqual(transcript_strand("+", is_reverse=True), "-")

    def test_forward_alignment_ts_minus(self):
        """Case 3: forward alignment, ts:- -> transcript on the reverse strand."""
        self.assertEqual(transcript_strand("-", is_reverse=False), "-")

    def test_reverse_alignment_ts_minus(self):
        """Case 4: reverse alignment, ts:- -> transcript on the forward strand."""
        self.assertEqual(transcript_strand("-", is_reverse=True), "+")

    def test_missing_ts_defaults_to_unstranded(self):
        """Case 5: absent ts yields '.' rather than a fabricated strand."""
        self.assertEqual(transcript_strand(None, is_reverse=False), ".")
        self.assertEqual(transcript_strand(None, is_reverse=True), ".")

    def test_invalid_ts_is_treated_as_missing(self):
        """Case 6: ts:A:* or junk is not a strand and must not be guessed."""
        for bad in ("*", "?", "", "x", "+-"):
            self.assertEqual(transcript_strand(bad, is_reverse=False), ".")

    def test_missing_ts_error_policy(self):
        """The error policy refuses to invent a strand."""
        with self.assertRaises(MissingTranscriptStrandError):
            transcript_strand(None, is_reverse=False, on_missing=MissingTsPolicy.ERROR)

    def test_missing_ts_flag_policy_is_opt_in(self):
        """The historical FLAG behaviour remains reachable, but only explicitly."""
        self.assertEqual(transcript_strand(None, True, on_missing=MissingTsPolicy.FLAG), "-")
        self.assertEqual(transcript_strand(None, False, on_missing=MissingTsPolicy.FLAG), "+")

    def test_unknown_policy_rejected(self):
        """An unrecognised policy is a programming error, not a silent default."""
        with self.assertRaises(ValueError):
            transcript_strand("+", False, on_missing="whatever")

    def test_ts_is_authoritative_regardless_of_policy(self):
        """When ts is present the policy never applies."""
        for policy in MissingTsPolicy.ALL:
            self.assertEqual(transcript_strand("-", is_reverse=True, on_missing=policy), "+")


class TestBedConversion(unittest.TestCase):
    """SAM/BAM -> BED12 behaviour."""

    def setUp(self):
        """Create a temporary working directory."""
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        """Remove temporary files."""
        shutil.rmtree(self.temp_dir)

    def test_bed_strand_emitted_correctly(self):
        """Case 7: every fixture read gets its ts-derived strand in BED column 6."""
        bed = self.temp_dir / "out.bed"
        splice_to_bed(FIXTURE_SAM, bed)
        rows = _by_name(bed)
        for name, expected in EXPECTED_STRAND.items():
            for row in rows[name]:
                self.assertEqual(row[5], expected, f"{name} strand")

    def test_secondary_alignments_skipped_by_default(self):
        """Filtering policy is unchanged from paftools splice2bed: secondaries dropped."""
        bed = self.temp_dir / "out.bed"
        stats = splice_to_bed(FIXTURE_SAM, bed)
        self.assertNotIn("read_secondary", _by_name(bed))
        self.assertEqual(stats.skipped_secondary, 1)

    def test_secondary_alignments_kept_on_request(self):
        """keep_secondary is a separate knob from the strand fix."""
        bed = self.temp_dir / "out.bed"
        splice_to_bed(FIXTURE_SAM, bed, keep_secondary=True)
        rows = _by_name(bed)
        self.assertIn("read_secondary", rows)
        self.assertEqual(rows["read_secondary"][0][8], "0,192,0")

    def test_supplementary_retained_as_separate_record(self):
        """Supplementary alignments still become their own BED record, as before."""
        bed = self.temp_dir / "out.bed"
        splice_to_bed(FIXTURE_SAM, bed)
        self.assertEqual(len(_by_name(bed)["read_two_primary"]), 2)

    def test_non_strand_fields_unchanged(self):
        """Case 10: name, score, thick range, blocks and colour keep their old meaning."""
        bed = self.temp_dir / "out.bed"
        splice_to_bed(FIXTURE_SAM, bed)
        row = _by_name(bed)["read_fwd_ts_plus"][0]
        self.assertEqual(row[0], "chr1")
        self.assertEqual(row[1], "100")  # POS 101, 1-based -> 100, 0-based
        self.assertEqual(row[2], "400")  # 50 + 200 + 50 reference bases
        self.assertEqual(row[4], "1000")
        self.assertEqual((row[6], row[7]), ("100", "400"))  # thick == chrom range
        self.assertEqual(row[8], "0,128,255")  # single primary
        self.assertEqual(row[9], "2")
        self.assertEqual(row[10], "50,50,")
        self.assertEqual(row[11], "0,250,")

    def test_multi_primary_colour(self):
        """A read with primary + supplementary is recoloured, as splice2bed does."""
        bed = self.temp_dir / "out.bed"
        splice_to_bed(FIXTURE_SAM, bed)
        for row in _by_name(bed)["read_two_primary"]:
            self.assertEqual(row[8], "255,0,0")

    def test_bam_and_sam_give_identical_bed(self):
        """Generic SAM/BAM support: the BAM path must agree byte for byte."""
        from_sam = self.temp_dir / "sam.bed"
        from_bam = self.temp_dir / "bam.bed"
        splice_to_bed(FIXTURE_SAM, from_sam)
        splice_to_bed(FIXTURE_BAM, from_bam)
        self.assertEqual(
            from_sam.read_text(encoding="utf8"),
            from_bam.read_text(encoding="utf8"),
        )

    def test_bam_reader_recovers_ts_tags(self):
        """The native BAM tag walker finds ts among other tags."""
        tags = {rec.query_name: rec.ts_tag for rec in iter_alignments(FIXTURE_BAM)}
        self.assertEqual(tags["read_fwd_ts_plus"], "+")
        self.assertEqual(tags["read_rev_ts_minus"], "-")
        self.assertIsNone(tags["read_no_ts"])

    def test_error_policy_propagates_from_converter(self):
        """on_missing_ts=error surfaces the unstranded read instead of hiding it."""
        with self.assertRaises(MissingTranscriptStrandError):
            splice_to_bed(FIXTURE_SAM, self.temp_dir / "out.bed", on_missing_ts=MissingTsPolicy.ERROR)

    def test_stats_report_strand_provenance(self):
        """The caller can see how many strands came from ts and how many did not."""
        stats = splice_to_bed(FIXTURE_SAM, self.temp_dir / "out.bed")
        self.assertEqual(stats.strand_from_ts, 6)
        self.assertEqual(stats.strand_unstranded, 2)
        self.assertEqual(stats.spliced_without_ts, 1)  # read_bad_ts is spliced; read_no_ts is not


class TestCoordinateSemantics(unittest.TestCase):
    """BED is 0-based half-open; GTF is 1-based closed. That must not drift."""

    def setUp(self):
        """Create a temporary working directory."""
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        """Remove temporary files."""
        shutil.rmtree(self.temp_dir)

    def test_bed_is_zero_based_half_open(self):
        """Case 9a: a 50M200N50M alignment at POS 101 spans [100, 400)."""
        bed = self.temp_dir / "out.bed"
        splice_to_bed(FIXTURE_SAM, bed)
        row = _by_name(bed)["read_fwd_ts_plus"][0]
        start, end = int(row[1]), int(row[2])
        self.assertEqual(start, 100)
        self.assertEqual(end - start, 300)

    def test_gtf_is_one_based_closed_and_strand_propagates(self):
        """Cases 8 and 9b: BED -> GTF keeps coordinates and carries the strand through."""
        splice_to_bed(FIXTURE_SAM, self.temp_dir / "reads.bed")
        _bed_to_gtf(self.temp_dir)
        gtf = (self.temp_dir / "annotation.gtf").read_text(encoding="utf8")

        exons, transcripts = [], []
        for line in gtf.splitlines():
            parts = line.split("\t")
            (exons if parts[2] == "exon" else transcripts).append(parts)

        # The first fixture read: exons 101-150 and 351-400 on the forward strand.
        first = [e for e in exons if e[3] == "101"]
        self.assertEqual(len(first), 1)
        self.assertEqual((first[0][3], first[0][4]), ("101", "150"))
        self.assertEqual(first[0][6], "+")
        second = [e for e in exons if e[3] == "351"]
        self.assertEqual((second[0][3], second[0][4]), ("351", "400"))

        # Case 8: GTF strand matches the BED strand for every read.
        bed_strands = {row[3]: row[5] for row in _bed_rows(self.temp_dir / "reads.bed")}
        self.assertEqual(set(bed_strands.values()), {"+", "-", "."})
        self.assertIn("\t.\t", gtf)  # unstranded reads stay unstranded in GTF

        # A GTF exon start is its BED block start + 1; the end is inclusive.
        for row in _bed_rows(self.temp_dir / "reads.bed"):
            offset = int(row[1])
            sizes = [int(x) for x in row[10].split(",") if x]
            starts = [int(x) for x in row[11].split(",") if x]
            for size, rel in zip(sizes, starts):
                want = (str(offset + rel + 1), str(offset + rel + size))
                self.assertTrue(
                    any((e[3], e[4]) == want and e[6] == row[5] for e in exons),
                    f"no GTF exon for BED block {want}",
                )

        # Transcript spans the outermost exon bounds.
        tx = [t for t in transcripts if t[3] == "101"][0]
        self.assertEqual((tx[3], tx[4]), ("101", "400"))


class TestHistoricalFlagRegression(unittest.TestCase):
    """The bug: splice2bed wrote the SAM FLAG orientation as the transcript strand."""

    def setUp(self):
        """Create a temporary working directory."""
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        """Remove temporary files."""
        shutil.rmtree(self.temp_dir)

    @staticmethod
    def _historical_flag_strand(record):
        """Reproduce paftools.js splice2bed: a1[5] = (flag & 16) ? '-' : '+'."""
        return "-" if record.is_reverse else "+"

    def test_flag_only_handling_is_wrong_for_discriminating_reads(self):
        """Two fixture reads have a transcript strand opposite to their FLAG."""
        wrong = {}
        for record in iter_alignments(FIXTURE_SAM):
            if record.is_secondary:
                continue
            old = self._historical_flag_strand(record)
            new = transcript_strand(record.ts_tag, record.is_reverse)
            if new in ("+", "-") and old != new:
                wrong[record.query_name] = (old, new)

        self.assertEqual(
            wrong,
            {
                "read_fwd_ts_minus": ("+", "-"),
                "read_rev_ts_minus": ("-", "+"),
            },
        )

    def test_new_implementation_emits_the_expected_strands(self):
        """End to end: the BED produced now carries the ts-derived strand."""
        bed = self.temp_dir / "out.bed"
        splice_to_bed(FIXTURE_SAM, bed)
        rows = _by_name(bed)
        self.assertEqual(rows["read_fwd_ts_minus"][0][5], "-")
        self.assertEqual(rows["read_rev_ts_minus"][0][5], "+")

    def test_flag_policy_reproduces_the_old_output_for_unstranded_reads(self):
        """Opting in to the FLAG policy restores the historical strand for ts-less reads."""
        bed = self.temp_dir / "out.bed"
        splice_to_bed(FIXTURE_SAM, bed, on_missing_ts=MissingTsPolicy.FLAG)
        rows = _by_name(bed)
        self.assertEqual(rows["read_no_ts"][0][5], "+")
        # ...but a read with ts is still resolved correctly, not from the FLAG.
        self.assertEqual(rows["read_rev_ts_minus"][0][5], "+")


if __name__ == "__main__":
    unittest.main()
