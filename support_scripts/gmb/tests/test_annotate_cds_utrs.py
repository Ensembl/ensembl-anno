#!/usr/bin/env python3
"""Unit tests for annotate_cds_utrs module.

Uses a tiny region of Chromosome 1 from the Candida auris genome
and Helixer predictions as ground truth.
"""

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from annotate_cds_utrs import (
    _make_orf_label,
    annotate_transcript,
    build_spliced_seq,
    check_frame_continuity,
    check_splice_sites,
    derive_utrs,
    find_best_orf,
    get_start_stop_positions,
    load_genome,
    map_cds_to_genomic,
    reverse_complement,
    translate,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

GENOME_PATH = os.path.join(os.path.dirname(__file__), "candida_auris_softmasked_toplevel.fa")

# Helixer transcript 1: Chr 1, + strand, single exon
# Exon:  6847-11220   CDS: 6990-11118   5'UTR: 6847-6990   3'UTR: 11118-11220
HELIXER_TX1 = {
    "chrom": "1",
    "strand": "+",
    "exons": [(6847, 11220)],
    "cds": [(6990, 11118)],
    "five_utr": [(6847, 6990)],
    "three_utr": [(11118, 11220)],
}

# Helixer transcript 2: Chr 1, + strand, single exon
# Exon: 11338-14171   CDS: 11535-13908   5'UTR: 11338-11535   3'UTR: 13908-14171
HELIXER_TX2 = {
    "chrom": "1",
    "strand": "+",
    "exons": [(11338, 14171)],
    "cds": [(11535, 13908)],
    "five_utr": [(11338, 11535)],
    "three_utr": [(13908, 14171)],
}


@pytest.fixture(scope="session")
def genome():
    """Load the genome once for all tests."""
    if not os.path.exists(GENOME_PATH):
        pytest.skip(f"Genome file not found: {GENOME_PATH}")
    return load_genome(GENOME_PATH)


# ---------------------------------------------------------------------------
# Test: FASTA loading
# ---------------------------------------------------------------------------


class TestLoadGenome:
    def test_loads_all_chromosomes(self, genome):
        assert len(genome) == 7
        for i in range(1, 8):
            assert str(i) in genome

    def test_sequence_is_uppercase(self, genome):
        # Softmasked FASTA has lowercase for repeats; we uppercase everything
        assert genome["1"][:100] == genome["1"][:100].upper()

    def test_small_fasta(self, tmp_path):
        fa = tmp_path / "test.fa"
        fa.write_text(">chrA\nACGT\nAAAA\n>chrB\nTTTT\n")
        g = load_genome(str(fa))
        assert g == {"chrA": "ACGTAAAA", "chrB": "TTTT"}


# ---------------------------------------------------------------------------
# Test: reverse_complement
# ---------------------------------------------------------------------------


class TestReverseComplement:
    def test_basic(self):
        assert reverse_complement("ACGT") == "ACGT"
        assert reverse_complement("AAAA") == "TTTT"
        assert reverse_complement("ATCG") == "CGAT"

    def test_with_n(self):
        assert reverse_complement("ACNGT") == "ACNGT"


# ---------------------------------------------------------------------------
# Test: build_spliced_seq
# ---------------------------------------------------------------------------


class TestBuildSplicedSeq:
    def test_forward_single_exon(self, genome):
        """Spliced cDNA for + strand is just the genomic subsequence."""
        exons = HELIXER_TX1["exons"]
        seq = build_spliced_seq(exons, "+", genome["1"])
        expected_len = exons[0][1] - exons[0][0]
        assert len(seq) == expected_len
        assert seq == genome["1"][exons[0][0] : exons[0][1]]

    def test_reverse_single_exon(self, genome):
        """Spliced cDNA for - strand is rev-comp of genomic subsequence."""
        exons = HELIXER_TX1["exons"]
        seq = build_spliced_seq(exons, "-", genome["1"])
        expected = reverse_complement(genome["1"][exons[0][0] : exons[0][1]])
        assert seq == expected

    def test_multi_exon_forward(self, genome):
        """Multi-exon + strand: concatenate exons in genomic order."""
        exons = [(100, 200), (300, 400)]
        seq = build_spliced_seq(exons, "+", genome["1"])
        assert len(seq) == 200
        assert seq == genome["1"][100:200] + genome["1"][300:400]

    def test_multi_exon_reverse(self, genome):
        """Multi-exon - strand: concatenate then rev-comp."""
        exons = [(100, 200), (300, 400)]
        seq = build_spliced_seq(exons, "-", genome["1"])
        genomic = genome["1"][100:200] + genome["1"][300:400]
        assert seq == reverse_complement(genomic)


# ---------------------------------------------------------------------------
# Test: find_best_orf
# ---------------------------------------------------------------------------


class TestFindBestOrf:
    def test_simple_orf(self):
        """ATG ... TAA in frame 0."""
        # 10 codons: ATG + 8 codons + TAA
        orf = "ATG" + "AAA" * 8 + "TAA"
        result = find_best_orf(orf, min_codons=1)
        assert result is not None
        start, end, p5, p3 = result
        assert start == 0
        assert end == len(orf)
        assert p5 is False
        assert p3 is False

    def test_prefers_atg(self):
        """Should prefer ATG-initiated ORF over 5′-partial."""
        # Frame 0: partial (no ATG), frame 1: has ATG
        seq = "A" + "ATG" + "GGG" * 10 + "TAA" + "CC"
        result = find_best_orf(seq, min_codons=1)
        assert result is not None
        start, end, p5, p3 = result
        assert start == 1  # ATG at position 1
        assert p5 is False

    def test_partial_5prime(self):
        """ORF starting at base 0 without ATG = 5′ partial."""
        # No ATG anywhere, just coding in frame 0
        seq = "GGG" * 20 + "TAA"
        result = find_best_orf(seq, min_codons=1)
        assert result is not None
        start, end, p5, p3 = result
        assert p5 is True

    def test_partial_3prime(self):
        """ATG but no stop codon = 3′ partial."""
        seq = "ATG" + "GGG" * 20
        result = find_best_orf(seq, min_codons=1)
        assert result is not None
        start, end, p5, p3 = result
        assert start == 0
        assert p3 is True

    def test_min_codons_filter(self):
        """Short ORF below threshold falls back to best available."""
        seq = "ATG" + "GGG" * 5 + "TAA"  # 7 codons
        # With min_codons=100, should still return it as fallback
        result = find_best_orf(seq, min_codons=100)
        assert result is not None
        assert result[0] == 0

    def test_no_orf_in_tiny_seq(self):
        result = find_best_orf("AA", min_codons=1)
        assert result is None

    def test_longest_atg_orf_wins(self):
        """If multiple ATG ORFs, pick the longest."""
        # Short ATG ORF in frame 0, long ATG ORF in frame 1
        short = "ATG" + "AAA" * 3 + "TAA"  # 5 codons
        long_orf = "ATG" + "GGG" * 20 + "TAG"  # 22 codons
        seq = short + "A" + long_orf
        result = find_best_orf(seq, min_codons=1)
        assert result is not None
        _, _, p5, _ = result
        # The long ORF should win
        assert (result[1] - result[0]) // 3 >= 20


# ---------------------------------------------------------------------------
# Test: translate
# ---------------------------------------------------------------------------


class TestTranslate:
    def test_simple(self):
        assert translate("ATGGGG") == "MG"

    def test_stop_codon(self):
        assert translate("ATGTAA") == "M*"
        assert translate("ATGTAG") == "M*"
        assert translate("ATGTGA") == "M*"

    def test_all_amino_acids(self):
        # Just check a few
        assert translate("TTT") == "F"
        assert translate("CTG") == "L"
        assert translate("GAT") == "D"
        assert translate("TGG") == "W"

    def test_unknown_codon(self):
        assert translate("NNN") == "X"


# ---------------------------------------------------------------------------
# Test: map_cds_to_genomic
# ---------------------------------------------------------------------------


class TestMapCdsToGenomic:
    def test_single_exon_forward(self):
        """CDS within a single + exon."""
        exons = [(1000, 2000)]
        # CDS at transcript positions 100-800
        result = map_cds_to_genomic(100, 800, exons, "+")
        assert result == [(1100, 1800)]

    def test_single_exon_reverse(self):
        """CDS within a single - exon."""
        exons = [(1000, 2000)]
        # For - strand, transcript pos 0 = genomic 2000,
        # so CDS at tx:100-800 → genomic: 1200-1900
        result = map_cds_to_genomic(100, 800, exons, "-")
        assert result == [(1200, 1900)]

    def test_multi_exon_forward(self):
        """CDS spanning two + exons."""
        exons = [(1000, 1100), (2000, 2200)]  # 100bp + 200bp = 300bp
        # CDS at tx:50-250 → exon1: 50-100 (genomic 1050-1100)
        #                  → exon2: 0-150  (genomic 2000-2150)
        result = map_cds_to_genomic(50, 250, exons, "+")
        assert result == [(1050, 1100), (2000, 2150)]

    def test_multi_exon_reverse(self):
        """CDS spanning two - exons."""
        exons = [(1000, 1100), (2000, 2200)]  # 100bp + 200bp = 300bp
        # For - strand, tx order is: exon2(2000-2200) then exon1(1000-1100)
        # CDS at tx:50-250
        # exon2 covers tx:0-200 → CDS overlap tx:50-200 →
        #   genomic: 2200-(200)=2000, offset 50 from right →
        #   g_cds_end = 2200-50=2150, g_cds_start = 2200-200=2000
        #   so genomic (2000, 2150)
        # exon1 covers tx:200-300 → CDS overlap tx:200-250 →
        #   offset 0-50 from tx start of this exon →
        #   g_cds_end = 1100-0=1100, g_cds_start = 1100-50=1050
        #   so genomic (1050, 1100)
        result = map_cds_to_genomic(50, 250, exons, "-")
        assert result == [(1050, 1100), (2000, 2150)]

    def test_full_exon_coverage(self):
        """CDS covers entire transcript."""
        exons = [(500, 600), (700, 800)]
        result = map_cds_to_genomic(0, 200, exons, "+")
        assert result == [(500, 600), (700, 800)]


# ---------------------------------------------------------------------------
# Test: derive_utrs
# ---------------------------------------------------------------------------


class TestDeriveUtrs:
    def test_forward_strand(self):
        exons = [(100, 500)]
        cds = [(200, 400)]
        five, three = derive_utrs(exons, cds, "+")
        assert five == [(100, 200)]
        assert three == [(400, 500)]

    def test_reverse_strand(self):
        """On - strand, 5′ UTR is rightmost, 3′ UTR is leftmost."""
        exons = [(100, 500)]
        cds = [(200, 400)]
        five, three = derive_utrs(exons, cds, "-")
        assert five == [(400, 500)]  # 5′ = downstream genomically
        assert three == [(100, 200)]  # 3′ = upstream genomically

    def test_no_utr(self):
        """CDS covers entire exon."""
        exons = [(100, 500)]
        cds = [(100, 500)]
        five, three = derive_utrs(exons, cds, "+")
        assert five == []
        assert three == []

    def test_multi_exon_utr(self):
        """UTR in first exon, CDS in second and third."""
        exons = [(100, 200), (300, 400), (500, 600)]
        cds = [(150, 200), (300, 400), (500, 550)]
        five, three = derive_utrs(exons, cds, "+")
        assert five == [(100, 150)]
        assert three == [(550, 600)]


# ---------------------------------------------------------------------------
# Test: annotate_transcript (integration)
# ---------------------------------------------------------------------------


class TestAnnotateTranscript:
    def test_with_existing_cds(self, genome):
        """Path A: CDS provided, derive UTR."""
        tx = HELIXER_TX1
        exon_df = pd.DataFrame(
            {
                "Start": [tx["exons"][0][0]],
                "End": [tx["exons"][0][1]],
            }
        )
        cds_df = pd.DataFrame(
            {
                "Start": [tx["cds"][0][0]],
                "End": [tx["cds"][0][1]],
            }
        )
        result = annotate_transcript(exon_df, tx["chrom"], tx["strand"], genome, cds_df=cds_df)

        assert result["cds"] == tx["cds"]
        assert result["five_prime_utr"] == tx["five_utr"]
        assert result["three_prime_utr"] == tx["three_utr"]
        assert len(result["cdna"]) == tx["exons"][0][1] - tx["exons"][0][0]
        assert result["protein"] is not None
        assert len(result["protein"]) > 0

    def test_cds_matches_helixer(self, genome):
        """Path B: ORF prediction should recover Helixer CDS boundaries."""
        tx = HELIXER_TX1
        exon_df = pd.DataFrame(
            {
                "Start": [tx["exons"][0][0]],
                "End": [tx["exons"][0][1]],
            }
        )
        # Don't provide CDS — let it predict
        result = annotate_transcript(exon_df, tx["chrom"], tx["strand"], genome, cds_df=None)

        assert result["cds"], "Should find an ORF"
        predicted_cds = result["cds"][0]
        expected_cds = tx["cds"][0]

        # The ORF prediction may not exactly match Helixer boundaries
        # (Helixer uses neural network, we use longest-ORF), but it should
        # be close. Check that the predicted CDS is within a small tolerance.
        assert (
            abs(predicted_cds[0] - expected_cds[0]) <= 30
        ), f"CDS start off by {abs(predicted_cds[0] - expected_cds[0])}bp"
        assert (
            abs(predicted_cds[1] - expected_cds[1]) <= 30
        ), f"CDS end off by {abs(predicted_cds[1] - expected_cds[1])}bp"

    def test_protein_starts_with_m(self, genome):
        """Predicted protein should start with M (ATG start)."""
        tx = HELIXER_TX2
        exon_df = pd.DataFrame(
            {
                "Start": [tx["exons"][0][0]],
                "End": [tx["exons"][0][1]],
            }
        )
        result = annotate_transcript(exon_df, tx["chrom"], tx["strand"], genome, cds_df=None)
        if result["protein"] and not result["is_partial_5"]:
            assert (
                result["protein"][0] == "M"
            ), f"Protein should start with M, got {result['protein'][:5]}"

    def test_utr_segmentation(self, genome):
        """UTR + CDS should cover entire exon span."""
        tx = HELIXER_TX1
        exon_df = pd.DataFrame(
            {
                "Start": [tx["exons"][0][0]],
                "End": [tx["exons"][0][1]],
            }
        )
        cds_df = pd.DataFrame(
            {
                "Start": [tx["cds"][0][0]],
                "End": [tx["cds"][0][1]],
            }
        )
        result = annotate_transcript(exon_df, tx["chrom"], tx["strand"], genome, cds_df=cds_df)

        # Total coverage = UTR + CDS should equal exon length
        total = 0
        for s, e in result["five_prime_utr"]:
            total += e - s
        for s, e in result["cds"]:
            total += e - s
        for s, e in result["three_prime_utr"]:
            total += e - s
        exon_len = tx["exons"][0][1] - tx["exons"][0][0]
        assert total == exon_len, f"UTR+CDS={total}, exon={exon_len}"


class TestPartialOrf:
    def test_partial_orf(self):
        """Handles transcripts without ATG start or stop codon."""
        # Sequence with no ATG or stop
        seq = "GGG" * 50
        result = find_best_orf(seq, min_codons=1)
        assert result is not None
        start, end, p5, p3 = result
        assert p5 is True
        assert p3 is True


# ---------------------------------------------------------------------------
# Test: get_start_stop_positions
# ---------------------------------------------------------------------------


class TestGetStartStopPositions:
    def test_forward_strand(self):
        cds = [(1000, 1300), (2000, 2200)]
        start, stop = get_start_stop_positions(cds, "+")
        assert start == 1000  # first base of CDS
        assert stop == 2200 - 3  # first base of stop codon

    def test_reverse_strand(self):
        cds = [(1000, 1300), (2000, 2200)]
        start, stop = get_start_stop_positions(cds, "-")
        assert start == 2200 - 3  # start codon at rightmost end (- strand)
        assert stop == 1000  # stop codon at leftmost end (- strand)

    def test_partial_none(self):
        start, stop = get_start_stop_positions([], "+")
        assert start is None
        assert stop is None


# ---------------------------------------------------------------------------
# Test: check_splice_sites
# ---------------------------------------------------------------------------


class TestCheckSpliceSites:
    def test_canonical_forward(self):
        """GT-AG splice on + strand."""
        # Exon1: 0-10, Exon2: 20-30
        # Intron: 10..20, donor at pos 10-11 = 'GT', acceptor at 18-19 = 'AG'
        seq = "A" * 10 + "GT" + "N" * 6 + "AG" + "A" * 10
        result = check_splice_sites([(0, 10), (20, 30)], "+", seq)
        assert len(result) == 1
        assert result[0]["class"] == "canonical"
        assert result[0]["donor"] == "GT"
        assert result[0]["acceptor"] == "AG"

    def test_noncanonical(self):
        """Non-standard splice dinucleotides."""
        seq = "A" * 10 + "AA" + "N" * 6 + "TT" + "A" * 10
        result = check_splice_sites([(0, 10), (20, 30)], "+", seq)
        assert len(result) == 1
        assert result[0]["class"] == "noncanonical"

    def test_known_noncanonical_gc_ag(self):
        """GC-AG is a known non-canonical splice."""
        seq = "A" * 10 + "GC" + "N" * 6 + "AG" + "A" * 10
        result = check_splice_sites([(0, 10), (20, 30)], "+", seq)
        assert result[0]["class"] == "known_noncanonical"


# ---------------------------------------------------------------------------
# Test: check_frame_continuity
# ---------------------------------------------------------------------------


class TestCheckFrameContinuity:
    def test_frame_ok(self):
        """Total CDS = 300bp (divisible by 3)."""
        cds = [(100, 250), (300, 450)]  # 150 + 150 = 300
        assert check_frame_continuity(cds, "+") is True

    def test_frame_break(self):
        """Total CDS = 301bp (not divisible by 3)."""
        cds = [(100, 251), (300, 450)]  # 151 + 150 = 301
        assert check_frame_continuity(cds, "+") is False

    def test_empty(self):
        assert check_frame_continuity([], "+") is True


# ---------------------------------------------------------------------------
# Test: _make_orf_label
# ---------------------------------------------------------------------------


class TestMakeOrfLabel:
    def test_complete(self):
        label = _make_orf_label("M" * 100, False, False)
        assert "100aa" in label
        assert "ATG" in label
        assert "STOP" in label

    def test_partial_5(self):
        label = _make_orf_label("G" * 50, True, False)
        assert "partial5" in label
        assert "STOP" in label

    def test_no_orf(self):
        assert _make_orf_label(None, False, False) == "no ORF"


# ---------------------------------------------------------------------------
# Test: synthetic full transcripts (+ and - strand)
# ---------------------------------------------------------------------------


class TestSyntheticTranscripts:
    """Integration tests on tiny synthetic transcripts.

    Constructs a 2-exon transcript with known sequence containing an ATG
    ORF and verifies CDS, UTR, start/stop, splice, and frame results.
    """

    def _make_genome_and_df(self, strand):
        """Build a small genome with a 2-exon transcript.

        Genome layout (0-based):
            0-100:   exon1
            100-120: intron (GT...AG)
            120-220: exon2

        For + strand: embed ATG at exon1 pos 30, stop TAA ending at exon2 pos 80
        cDNA = exon1[0:100] + exon2[0:100] = 200bp
        ORF at cDNA 30..183 (ATG at 30, TAA at 180-182, end=183) = 51 codons
        """
        # Build the genome seq ensuring GT-AG splice
        if strand == "+":
            # exon1: 30bp UTR + ATG + coding
            utr5_exon1 = "AAA" * 10  # 30bp
            coding_exon1 = "ATG" + "GCT" * 22 + "G"  # 3 + 66 + 1 = 70bp
            exon1 = utr5_exon1 + coding_exon1  # 100bp
            intron = "GT" + "A" * 16 + "AG"  # 20bp
            # coding continues 83bp into exon2,
            # cDNA so far: 30bp utr + 70bp coding_exon1 = offset 100.
            # Need stop at cDNA pos 183, which is exon2 offset 83.
            # 83 - 1 (carry G from exon1) = 82bp more coding in exon2
            # then TAA at 82..85, but we need it frame-aligned:
            # From exon1 we have 70bp coding = 23 codons + 1bp leftover (G)
            # In exon2: GCT*27 = 81bp, then first codon is (G from exon1 +
            # GC from exon2 = GGC). Let me be precise:
            # cDNA[30:33] = ATG (codon 1)
            # cDNA[33:36] = GCT (codon 2)
            # ... continues in frame
            # cDNA[30 + 51*3] = cDNA[183] should be first base after stop
            # Coding region in exon2: need 100 - 70 = 30bp of coding + stop
            # From cDNA 100..100+80 = exon2[0:80]
            # cDNA[30] = ATG, cDNA[30 + 51*3 - 3] = cDNA[180] = start of TAA
            # So stop at cDNA 180..183. That's exon2[80:83].
            # exon2: coding[0:80] + TAA + UTR3
            coding_exon2 = "GCT" * 26 + "GC"  # 80bp
            stop_codon = "TAA"
            utr3_exon2 = "CCC" * 5 + "CC"  # 17bp
            exon2 = coding_exon2 + stop_codon + utr3_exon2  # 100bp
            genome_seq = exon1 + intron + exon2
        else:
            # For - strand, we reverse-complement the design.
            # We'll construct a + strand sequence and the transcript reads
            # from right to left.
            # Keep same physical layout: exon1=0-100, intron=100-120, exon2=120-220
            # Transcript order (5'→3'): exon2 rev-comp, then exon1 rev-comp
            # Build desired cDNA first, then lay it out as rev-comp
            utr5 = "AAA" * 10  # 30bp
            coding = "ATG" + "GCT" * 49 + "TAA"  # 3 + 147 + 3 = 153bp = 51 codons
            utr3 = "CCC" * 5 + "CC"  # 17bp
            cdna = utr5 + coding + utr3  # 200bp

            # cDNA [0:100] maps to exon2 rev-comp, [100:200] to exon1 rev-comp
            exon2_revcomp = cdna[:100]
            exon1_revcomp = cdna[100:200]
            exon1 = reverse_complement(exon1_revcomp)
            exon2 = reverse_complement(exon2_revcomp)
            # Need GT-AG intron on + strand
            intron = "CT" + "A" * 16 + "AC"  # complementary: AG...GT on - strand
            genome_seq = exon1 + intron + exon2

        genome = {"chrT": genome_seq}
        exon_df = pd.DataFrame(
            {
                "Start": [0, 120],
                "End": [100, 220],
            }
        )
        return genome, exon_df, strand

    def test_forward_full(self):
        """+ strand synthetic transcript: UTR + CDS + start/stop + splice."""
        genome, exon_df, strand = self._make_genome_and_df("+")
        result = annotate_transcript(exon_df, "chrT", strand, genome, min_codons=1)

        # Should find CDS
        assert result["cds"], "Should find an ORF"
        assert result["protein"] is not None
        assert result["protein"][0] == "M", "Should start with M"
        assert "*" not in result["protein"], "Protein should not contain stop"

        # UTRs should exist
        assert result["five_prime_utr"], "Should have 5-UTR"
        assert result["three_prime_utr"], "Should have 3-UTR"

        # Start/stop positions
        assert result["start_pos"] is not None
        assert result["stop_pos"] is not None
        assert result["start_pos"] < result["stop_pos"]  # + strand

        # Splice sites — should be canonical GT-AG
        assert len(result["splice_sites"]) == 1
        assert result["splice_sites"][0]["class"] == "canonical"

        # Frame continuity
        assert result["frame_ok"]

        # ORF label
        assert "ATG" in result["orf_label"]
        assert "STOP" in result["orf_label"]

    def test_reverse_full(self):
        """- strand synthetic transcript: UTR + CDS + start/stop + splice."""
        genome, exon_df, strand = self._make_genome_and_df("-")
        result = annotate_transcript(exon_df, "chrT", strand, genome, min_codons=1)

        # Should find CDS
        assert result["cds"], "Should find an ORF"
        assert result["protein"] is not None
        assert result["protein"][0] == "M", "Should start with M"
        assert "*" not in result["protein"], "No internal stops"

        # UTRs
        assert result["five_prime_utr"], "Should have 5-UTR"
        assert result["three_prime_utr"], "Should have 3-UTR"

        # Start/stop: on - strand, start_pos > stop_pos (genomically)
        assert result["start_pos"] is not None
        assert result["stop_pos"] is not None
        assert result["start_pos"] > result["stop_pos"]  # - strand

        # Splice
        assert len(result["splice_sites"]) == 1

        # Frame
        assert result["frame_ok"]


def test_compute_utr_end_support():
    """Test utr end support logic."""
    import pandas as pd

    from config import PipelineConfig
    from gene_model_builder import compute_utr_end_support

    cfg = PipelineConfig()
    cfg.utr.require_end_support = True
    cfg.utr.end_support_mode = "multisource_end_agreement"
    cfg.utr.end_support_sources = ["StringTie", "Scallop"]
    cfg.utr.end_tolerance_bp = 10
    cfg.utr.require_multisource_for_utr_5p = True
    cfg.utr.require_multisource_for_utr_3p = True
    cfg.utr.fallback_policy_when_unsupported = "drop_utr"

    # Locus df with other models
    data = [
        # Scallop model (matching 5p, mismatched 3p)
        {
            "transcript_id": "sc_1",
            "Feature": "exon",
            "Start": 1005,
            "End": 1500,
            "Source": "Scallop",
            "Strand": "+",
        },
        {
            "transcript_id": "sc_1",
            "Feature": "exon",
            "Start": 1600,
            "End": 2500,
            "Source": "Scallop",
            "Strand": "+",
        },
    ]
    df = pd.DataFrame(data)

    # The model we are testing
    model = {
        "id": "st_1",
        "strand": "+",
        "start": 1000,
        "end": 2900,
        # 'protein_coding_score' not set
    }

    res = compute_utr_end_support(model, df, cfg)

    # 5p should be supported (1000 vs 1005 with tol 10)
    assert res["supported_5p"] is True
    # 3p should not be supported (2900 vs 2500, diff 400 > tol)
    assert res["supported_3p"] is False
    assert res["action_3p"] == "dropped"

    # Try with drop_transcript fallback
    cfg.utr.fallback_policy_when_unsupported = "drop_transcript"
    res = compute_utr_end_support(model, df, cfg)
    assert res["action_3p"] == "drop_transcript"


# ---------------------------------------------------------------------------
# Test: multi-exon CDS translation correctness (plus- and minus-strand)
#
# These tests pin the exact peptide produced from a synthetic genome where
# the correct protein is derivable by inspection.  They serve as regression
# tests for the RC(X+Y)==RC(Y)+RC(X) ordering bug on minus-strand CDS.
# ---------------------------------------------------------------------------


class TestCdsTranslationMultiExon:
    """Verify annotate_transcript produces the right protein for multi-exon genes.

    Genome layout (0-based, half-open):
        pos  0-9 : padding
        pos 10-15: exon1_low  (used as-is for + ; RC'd for - strand 3' segment)
        pos 16-49: padding
        pos 50-55: exon2_high (used as-is for + ; RC'd for - strand 5' segment)
        pos 56+  : padding

    Plus strand CDS = genome[10:16] + genome[50:56]
    Minus strand CDS (5'→3' mRNA) = RC(genome[50:56]) + RC(genome[10:16])

    We choose sequences so that both orientations yield protein "MPG"
    when assembled correctly:
        + : "ATGCCC" + "GGGTAA" → ATG-CCC-GGG-TAA → MPG(stop)
        - : RC("GGGTAA") + RC("ATGCCC") = "TTACCC"+"GGGCAT" in genome;
            ascending-concat-then-RC → RC("TTACCC"+"GGGCAT")
                                      = RC("GGGCAT")+RC("TTACCC")    [RC algebra]
                                      = "ATGCCC" + "GGGTAA"         ✓
            → same mRNA CDS → MPG(stop)
    """

    _PADDING = "AAAAAAAAAA"  # 10 bp
    # Exon low [10:16]
    _EXON1_PLUS = "ATGCCC"   # → codons ATG-CCC (M-P)
    _EXON1_MINUS = "TTACCC"  # RC = GGGTAA (provides stop + 3' tail of - mRNA)
    # Exon high [50:56]
    _EXON2_PLUS = "GGGTAA"   # → codons GGG-TAA (G-stop)
    _EXON2_MINUS = "GGGCAT"  # RC = ATGCCC (provides ATG + M-P of - mRNA)

    _GAP = "N" * 34  # fills [16:50]
    _TAIL = "TTTTTTTTTT"

    @staticmethod
    def _make_genome(exon1_seq: str, exon2_seq: str) -> dict:
        """Build {chrom: seq} with exon1 at [10:16) and exon2 at [50:56)."""
        padding = "AAAAAAAAAA"
        gap = "N" * 34
        tail = "TTTTTTTTTT"
        seq = padding + exon1_seq + gap + exon2_seq + tail
        assert seq[10:16] == exon1_seq
        assert seq[50:56] == exon2_seq
        return {"chr1": seq}

    @staticmethod
    def _cds_df(intervals):
        return pd.DataFrame({"Start": [s for s, e in intervals],
                             "End":   [e for s, e in intervals]})

    @staticmethod
    def _exon_df(intervals):
        return pd.DataFrame({"Start": [s for s, e in intervals],
                             "End":   [e for s, e in intervals]})

    # ── plus strand ──────────────────────────────────────────────────────────

    def test_plus_strand_multiexon_expected_peptide(self):
        """+ strand multi-exon CDS: two exons → known peptide MPG, no internal stops."""
        genome = self._make_genome(self._EXON1_PLUS, self._EXON2_PLUS)
        exon_df = self._exon_df([(10, 16), (50, 56)])
        cds_df  = self._cds_df([(10, 16), (50, 56)])

        result = annotate_transcript(exon_df, "chr1", "+", genome, cds_df=cds_df)

        protein = result["protein"]
        assert protein is not None, "Should produce a protein"
        assert "*" not in protein, f"Internal stop in + strand protein: {protein!r}"
        assert protein == "MPG", f"Expected MPG, got: {protein!r}"

    def test_plus_strand_multiexon_phase_continuity(self):
        """+ strand: CDS total length must be multiple of 3 (frame preserved across exons)."""
        genome = self._make_genome(self._EXON1_PLUS, self._EXON2_PLUS)
        exon_df = self._exon_df([(10, 16), (50, 56)])
        cds_df  = self._cds_df([(10, 16), (50, 56)])

        result = annotate_transcript(exon_df, "chr1", "+", genome, cds_df=cds_df)

        total_cds_bp = sum(e - s for s, e in result["cds"])
        assert total_cds_bp % 3 == 0, (
            f"CDS length {total_cds_bp} is not a multiple of 3 — reading frame broken"
        )
        assert result.get("frame_ok"), "frame_ok should be True"

    # ── minus strand ─────────────────────────────────────────────────────────

    def test_minus_strand_multiexon_expected_peptide(self):
        """- strand multi-exon CDS: ascending-then-RC gives correct protein MPG.

        If exons were concatenated in *descending* order before RC (the old bug),
        the mRNA sequence would be reversed at the junction and produce internal
        stop codons instead.
        """
        genome = self._make_genome(self._EXON1_MINUS, self._EXON2_MINUS)
        # Exon intervals: low [10,16), high [50,56)  — sorted ascending in GFF3
        exon_df = self._exon_df([(10, 16), (50, 56)])
        cds_df  = self._cds_df([(10, 16), (50, 56)])

        result = annotate_transcript(exon_df, "chr1", "-", genome, cds_df=cds_df)

        protein = result["protein"]
        assert protein is not None, "Should produce a protein"
        assert "*" not in protein, (
            f"Internal stop codon in - strand protein: {protein!r}\n"
            "(Hint: exons are being concatenated in wrong order before RC)"
        )
        assert protein[0] == "M", f"Protein should start with Met, got: {protein!r}"
        assert protein == "MPG", f"Expected MPG, got: {protein!r}"

    def test_minus_strand_multiexon_phase_continuity(self):
        """- strand: CDS total length must be multiple of 3 across exon boundary."""
        genome = self._make_genome(self._EXON1_MINUS, self._EXON2_MINUS)
        exon_df = self._exon_df([(10, 16), (50, 56)])
        cds_df  = self._cds_df([(10, 16), (50, 56)])

        result = annotate_transcript(exon_df, "chr1", "-", genome, cds_df=cds_df)

        total_cds_bp = sum(e - s for s, e in result["cds"])
        assert total_cds_bp % 3 == 0, (
            f"CDS length {total_cds_bp} is not a multiple of 3 — reading frame broken"
        )
        assert result.get("frame_ok"), "frame_ok should be True"

    def test_plus_minus_symmetry(self):
        """+ and - strand CDS built from mirror sequences yield the same protein."""
        # Plus genome
        gp = self._make_genome(self._EXON1_PLUS, self._EXON2_PLUS)
        # Minus genome (mirror): swap sequences so RC ordering produces same mRNA
        gm = self._make_genome(self._EXON1_MINUS, self._EXON2_MINUS)

        exon_df = self._exon_df([(10, 16), (50, 56)])
        cds_df  = self._cds_df([(10, 16), (50, 56)])

        rp = annotate_transcript(exon_df, "chr1", "+", gp, cds_df=cds_df)
        rm = annotate_transcript(exon_df, "chr1", "-", gm, cds_df=cds_df)

        assert rp["protein"] == rm["protein"], (
            f"Plus protein {rp['protein']!r} != minus protein {rm['protein']!r}"
        )

    def test_minus_strand_three_exons_no_internal_stop(self):
        """- strand with three CDS exons must produce no internal stop codons.

        Three exons stress-test the ordering fix more thoroughly than two.
        """
        # Build genome with 3 CDS exons (each 6 bp = 2 codons) at [10:16), [50:56), [90:96)
        # mRNA CDS (5'→3') = RC([90:96]) + RC([50:56]) + RC([10:16])
        # Choose exon sequences so mRNA = ATG-CCC-GGG-TTT-AAA-TAA → MPGFKstop
        # RC([90:96])=ATGCCC → [90:96]=GGGCAT
        # RC([50:56])=GGGTTТ → [50:56]=AAACCC (wait: RC("AAACCC")=GGGTTT ✓)
        # RC([10:16])=TTTAAA → [10:16]=TTTAAA (self-complementary for A/T pairs... RC("TTTAAA")=TTTAAA ✓)
        # Actually: RC("TTTAAA") = complement(rev("TTTAAA")) = complement("AAATTT") = "TTTAAA" ✓ (palindrome)

        # Verify: RC("GGGCAT")=ATGCCC, RC("AAACCC")=GGGTTT, RC("TTTAAA")=TTTAAA
        # mRNA = ATGCCC + GGGTTT + TTTAAA → ATG-CCC-GGG-TTT-TTT-AAA
        # protein = MPGFFK  (no stop in the middle, stop not included = partial? let's add one)
        # Adjust: [10:16] = RC of "TAAYYY" for stop — let's use "TTAGGG" → RC="CCCTAA"
        # mRNA = ATGCCC + GGGTTT + CCCTAA = ATG-CCC-GGG-TTT-CCC-TAA → MPGFP(stop)
        # [10:16] = "TTAGGG"

        seq = ("A" * 10       # [0:10]  padding
               + "TTAGGG"     # [10:16] exon1_low  (RC=CCCTAA → 3' mRNA end)
               + "N" * 34     # [16:50] gap
               + "AAACCC"     # [50:56] exon2_mid  (RC=GGGTTT → middle of mRNA)
               + "N" * 34     # [56:90] gap
               + "GGGCAT"     # [90:96] exon3_high (RC=ATGCCC → 5' mRNA start)
               + "A" * 10)    # padding

        genome = {"chr1": seq}
        exon_df = self._exon_df([(10, 16), (50, 56), (90, 96)])
        cds_df  = self._cds_df([(10, 16), (50, 56), (90, 96)])

        result = annotate_transcript(exon_df, "chr1", "-", genome, cds_df=cds_df)

        protein = result["protein"]
        assert protein is not None, "Should produce a protein"
        assert "*" not in protein, (
            f"Internal stop in 3-exon - strand protein: {protein!r}"
        )
        assert protein[0] == "M", f"Protein should start with Met, got: {protein!r}"
        assert protein == "MPGFP", f"Expected MPGFP, got: {protein!r}"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
