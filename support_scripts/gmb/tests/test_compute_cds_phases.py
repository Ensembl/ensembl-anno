"""Tests for compute_cds_phases in gmb.pipeline.builder.

GFF3 CDS phase: the number of bases at the biological 5' start of each CDS
segment to skip to reach the first complete codon boundary (0, 1, or 2).

For + strand, biological order == ascending genomic order.
For - strand, biological order == descending genomic order, but GFF3 always
reports features in ascending genomic order with phase mapped back accordingly.
"""

from gmb.pipeline.builder import compute_cds_phases


class TestComputeCdsPhasesEmpty:
    def test_empty_intervals_returns_empty(self):
        assert compute_cds_phases([], "+") == []
        assert compute_cds_phases([], "-") == []

    def test_unresolved_strand_returns_empty(self):
        assert compute_cds_phases([(100, 109)], ".") == []
        assert compute_cds_phases([(100, 109)], "") == []


class TestComputeCdsPhasSingleExon:
    def test_single_exon_forward_phase_is_zero(self):
        # Biological 5' start is the first base -- always phase 0.
        assert compute_cds_phases([(100, 109)], "+") == [0]

    def test_single_exon_reverse_phase_is_zero(self):
        assert compute_cds_phases([(100, 109)], "-") == [0]

    def test_single_exon_length_divisible_by_3(self):
        assert compute_cds_phases([(0, 12)], "+") == [0]
        assert compute_cds_phases([(0, 12)], "-") == [0]


class TestComputeCdsPhasesMultiExonForward:
    def test_two_exons_lengths_divisible_by_3(self):
        # 6 + 6 bp: each exon starts a clean codon boundary.
        assert compute_cds_phases([(100, 106), (200, 206)], "+") == [0, 0]

    def test_two_exons_with_phase_carry(self):
        # First exon: 7 bp (not divisible by 3).
        # Second exon phase = (3 - 7%3) % 3 = (3-1)%3 = 2.
        assert compute_cds_phases([(100, 107), (200, 209)], "+") == [0, 2]

    def test_two_exons_phase_carry_1(self):
        # First exon: 8 bp. Phase carry = (3 - 8%3) % 3 = (3-2)%3 = 1.
        assert compute_cds_phases([(100, 108), (200, 209)], "+") == [0, 1]

    def test_three_exons_accumulating_phase(self):
        # Exon lengths: 7, 5, 9
        # Start: phase=0, after 7 bp carry=(3-1)%3=2
        # Second: phase=2, after 7+5=12 bp carry=(3-0)%3=0
        # Third: phase=0
        result = compute_cds_phases([(0, 7), (10, 15), (20, 29)], "+")
        assert result == [0, 2, 0]

    def test_order_is_ascending_genomic(self):
        # Unsorted input should produce the same result as sorted input.
        unsorted = [(200, 207), (100, 107)]
        sorted_in = [(100, 107), (200, 207)]
        assert compute_cds_phases(unsorted, "+") == compute_cds_phases(sorted_in, "+")


class TestComputeCdsPhasesMultiExonReverse:
    def test_two_exons_lengths_divisible_by_3(self):
        # Biological order is right-to-left: (200,206) then (100,106).
        # Both 6 bp: phases are [0, 0] in bio order.
        # GFF3 ascending = reversed: [0, 0].
        assert compute_cds_phases([(100, 106), (200, 206)], "-") == [0, 0]

    def test_two_exons_with_phase_carry(self):
        # Bio order: (200,207) first (7 bp), then (100,109).
        # Bio phases: [0, 2].
        # GFF3 ascending order of [(100,109),(200,207)] → reversed: [2, 0].
        result = compute_cds_phases([(100, 109), (200, 207)], "-")
        assert result == [2, 0]

    def test_three_exons_reverse(self):
        # Exons (ascending): (10,13), (20,27), (30,33) → lengths 3, 7, 3
        # Bio order (descending): (30,33)→3 bp, (20,27)→7 bp, (10,13)→3 bp
        # Bio phases:
        #   (30,33): cum=0 → phase 0; after: cum=3
        #   (20,27): cum=3 → (3-0)%3=0; after: cum=10
        #   (10,13): cum=10 → (3-1)%3=2
        # Bio phases = [0, 0, 2]; reversed for GFF3 = [2, 0, 0]
        result = compute_cds_phases([(10, 13), (20, 27), (30, 33)], "-")
        assert result == [2, 0, 0]

    def test_reverse_phase_count_matches_intervals(self):
        cds = [(0, 5), (10, 17), (20, 26)]
        result = compute_cds_phases(cds, "-")
        assert len(result) == 3
        assert all(p in (0, 1, 2) for p in result)

    def test_unsorted_input_reverse(self):
        # Provide intervals out of ascending order.
        unsorted = [(200, 207), (100, 109)]
        sorted_in = [(100, 109), (200, 207)]
        assert compute_cds_phases(unsorted, "-") == compute_cds_phases(sorted_in, "-")


class TestComputeCdsPhasesBoundary:
    def test_all_phase_values_are_0_1_or_2(self):
        # Stress test with lengths that cycle through all residues mod 3.
        # Lengths: 1, 2, 3, 4, 5, 6 → carries: 0→2→0→2→0→0
        cds = [(i * 10, i * 10 + (i + 1)) for i in range(6)]
        result = compute_cds_phases(cds, "+")
        assert len(result) == 6
        assert all(p in (0, 1, 2) for p in result)

    def test_first_phase_always_zero(self):
        # The biological 5' start of any transcript is always phase 0.
        for strand in ("+", "-"):
            result = compute_cds_phases([(0, 7), (10, 17)], strand)
            if strand == "+":
                assert result[0] == 0, f"+ strand first phase must be 0, got {result}"
            else:
                # For -, GFF3 output has the rightmost interval first (bio 5').
                # The last element in ascending order is the first in bio order.
                assert result[-1] == 0, f"- strand rightmost phase must be 0, got {result}"
