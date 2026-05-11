#!/usr/bin/env python3
"""Tests for FASTA QC validation and the FASTA export reconciliation fix.

Covers:
  C1 - Minimal FASTA export integrity test (multi-exon +/- strand, phase, no-CDS)
  C2 - Regression test for ghost proteins (validate_and_fix / dedup dropping transcripts)
  C3 - fasta_qc module checks
"""

import os
import json
import pytest

from gmb.pipeline.annotate_cds_utrs import reverse_complement, translate
from gmb.pipeline.fasta_qc import (
    parse_fasta_ids,
    parse_fasta_records,
    parse_gff3,
    run_coverage_checks,
    run_sequence_checks,
    validate_fasta,
)


# ---------------------------------------------------------------------------
# Fixtures: tiny genome + GFF for C1
# ---------------------------------------------------------------------------

# Contig sequence: 500bp of known content with embedded ORFs.
# We place two genes:
#   gene_plus  on + strand:  exons [50..120), [200..350) → CDS [60..120), [200..330)
#   gene_minus on - strand:  exons [370..430), [440..490) → CDS [370..430), [440..480)
#   gene_nocds: exon only, no CDS → should NOT appear in prot.fa
_CONTIG_SEQ = (
    "A" * 50  # 0..50: padding
    + "AAAAAAAAAA"  # 50..60: 5' UTR of plus gene (exon1 start)
    + "ATGGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT"  # 60..120: CDS exon1 (60bp, 20 codons incl ATG)
    + "A" * 80  # 120..200: intron
    + "GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTTAA"  # 200..330: CDS exon2 (130bp)
    + "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"  # 330..370: 3' UTR + padding (40bp to reach 370)
    + "ATGGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTTAA"  # 370..430: minus CDS exon1 (60bp)
    + "A" * 10  # 430..440: intron
    + "GCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCT"  # 440..480: minus CDS exon2 (40bp)
    + "CCCCCCCCCC"  # 480..490: minus UTR (10bp)
    + "T" * 10  # 490..500: padding
)


def _write_genome(tmp_path):
    """Write tiny genome FASTA, return path and dict."""
    fa_path = tmp_path / "genome.fa"
    fa_path.write_text(f">chr1\n{_CONTIG_SEQ}\n")
    return str(fa_path), {"chr1": _CONTIG_SEQ.upper()}


def _write_gff3(tmp_path, include_nocds=True, extra_ghost_tid=None):
    """Write a tiny consensus GFF3.

    If extra_ghost_tid is set, we do NOT include it in GFF (simulates a dropped transcript),
    but the caller should add it to FASTA to simulate the bug.
    """
    lines = ["##gff-version 3"]

    # Gene plus (+ strand, multi-exon with CDS)
    lines.append("chr1\tGMB\tgene\t51\t350\t.\t+\t.\tID=GENE_001")
    lines.append(
        "chr1\tGMB\tmRNA\t51\t350\t.\t+\t.\tID=GENE_001.t1;Parent=GENE_001"
    )
    lines.append(
        "chr1\tGMB\texon\t51\t120\t.\t+\t.\tID=GENE_001.t1.exon1;Parent=GENE_001.t1"
    )
    lines.append(
        "chr1\tGMB\texon\t201\t350\t.\t+\t.\tID=GENE_001.t1.exon2;Parent=GENE_001.t1"
    )
    lines.append(
        "chr1\tGMB\tCDS\t61\t120\t.\t+\t0\tID=GENE_001.t1.cds1;Parent=GENE_001.t1"
    )
    lines.append(
        "chr1\tGMB\tCDS\t201\t330\t.\t+\t0\tID=GENE_001.t1.cds2;Parent=GENE_001.t1"
    )

    # Gene minus (- strand, multi-exon with CDS)
    lines.append("chr1\tGMB\tgene\t371\t490\t.\t-\t.\tID=GENE_002")
    lines.append(
        "chr1\tGMB\tmRNA\t371\t490\t.\t-\t.\tID=GENE_002.t1;Parent=GENE_002"
    )
    lines.append(
        "chr1\tGMB\texon\t371\t430\t.\t-\t.\tID=GENE_002.t1.exon1;Parent=GENE_002.t1"
    )
    lines.append(
        "chr1\tGMB\texon\t441\t490\t.\t-\t.\tID=GENE_002.t1.exon2;Parent=GENE_002.t1"
    )
    lines.append(
        "chr1\tGMB\tCDS\t371\t430\t.\t-\t0\tID=GENE_002.t1.cds1;Parent=GENE_002.t1"
    )
    lines.append(
        "chr1\tGMB\tCDS\t441\t480\t.\t-\t0\tID=GENE_002.t1.cds2;Parent=GENE_002.t1"
    )

    # Gene without CDS (exon only)
    if include_nocds:
        lines.append("chr1\tGMB\tgene\t1\t30\t.\t+\t.\tID=GENE_003")
        lines.append(
            "chr1\tGMB\tmRNA\t1\t30\t.\t+\t.\tID=GENE_003.t1;Parent=GENE_003"
        )
        lines.append(
            "chr1\tGMB\texon\t1\t30\t.\t+\t.\tID=GENE_003.t1.exon1;Parent=GENE_003.t1"
        )

    gff_path = tmp_path / "consensus.gff3"
    gff_path.write_text("\n".join(lines) + "\n")
    return str(gff_path)


def _write_prot_fa(tmp_path, records: dict[str, str]):
    """Write prot.fa from {id: sequence} dict."""
    path = tmp_path / "prot.fa"
    lines = []
    for tid, seq in records.items():
        lines.append(f">{tid}")
        lines.append(seq)
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def _write_cdna_fa(tmp_path, records: dict[str, str]):
    """Write cdna.fa from {id: sequence} dict."""
    path = tmp_path / "cdna.fa"
    lines = []
    for tid, seq in records.items():
        lines.append(f">{tid}")
        lines.append(seq)
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def _reconstruct_protein_plus():
    """Reconstruct expected protein for GENE_001.t1 (+ strand)."""
    # CDS exon1: [60..120) = 60bp, CDS exon2: [200..330) = 130bp → 190bp total
    cds_nuc = _CONTIG_SEQ[60:120] + _CONTIG_SEQ[200:330]
    prot = translate(cds_nuc)
    if prot.endswith("*"):
        prot = prot[:-1]
    return prot


def _reconstruct_protein_minus():
    """Reconstruct expected protein for GENE_002.t1 (- strand).

    CDS intervals [370..430) and [440..480) must be concatenated in *ascending*
    genomic order before reverse-complementing.  RC(X+Y)==RC(Y)+RC(X), so
    ascending-then-RC places the highest-coordinate (5') exon first in the
    resulting mRNA-orientation CDS — matching the corrected fasta_qc logic.
    """
    # CDS: [370..430) = 60bp, [440..480) = 40bp → 100bp total
    # Ascending genomic order then revcomp (NOT reversed order then revcomp)
    cds_nuc = _CONTIG_SEQ[370:430] + _CONTIG_SEQ[440:480]
    cds_nuc = reverse_complement(cds_nuc)
    prot = translate(cds_nuc)
    if prot.endswith("*"):
        prot = prot[:-1]
    return prot


def _reconstruct_cdna_plus():
    """Reconstruct expected cDNA for GENE_001.t1."""
    return _CONTIG_SEQ[50:120] + _CONTIG_SEQ[200:350]


def _reconstruct_cdna_minus():
    """Reconstruct expected cDNA for GENE_002.t1."""
    cdna = _CONTIG_SEQ[370:430] + _CONTIG_SEQ[440:490]
    return reverse_complement(cdna)


def _reconstruct_cdna_nocds():
    """Reconstruct expected cDNA for GENE_003.t1 (no CDS)."""
    return _CONTIG_SEQ[0:30]


# ---------------------------------------------------------------------------
# C1: Minimal FASTA export integrity test
# ---------------------------------------------------------------------------

class TestFastaExportIntegrity:
    """C1: Verify FASTA QC passes for correct outputs."""

    def test_correct_output_passes_qc(self, tmp_path):
        """When FASTA matches GFF exactly, QC should pass."""
        _write_genome(tmp_path)
        out_dir = tmp_path / "output"
        out_dir.mkdir()

        _write_gff3(out_dir)

        prot_plus = _reconstruct_protein_plus()
        prot_minus = _reconstruct_protein_minus()
        _write_prot_fa(
            out_dir,
            {"GENE_001.t1": prot_plus, "GENE_002.t1": prot_minus},
        )
        _write_cdna_fa(
            out_dir,
            {
                "GENE_001.t1": _reconstruct_cdna_plus(),
                "GENE_002.t1": _reconstruct_cdna_minus(),
                "GENE_003.t1": _reconstruct_cdna_nocds(),
            },
        )

        report = validate_fasta(str(out_dir))
        assert report["pass"], f"QC should pass but got: {json.dumps(report, indent=2)}"
        assert report["n_prot_records"] == 2
        assert report["n_cdna_records"] == 3
        assert report["n_cds_transcripts"] == 2

    def test_nocds_transcript_excluded_from_prot(self, tmp_path):
        """Transcript without CDS should be in cdna.fa but NOT prot.fa."""
        out_dir = tmp_path / "output"
        out_dir.mkdir()
        _write_gff3(out_dir)

        prot_plus = _reconstruct_protein_plus()
        prot_minus = _reconstruct_protein_minus()

        # Correct: only CDS-bearing transcripts in prot.fa
        _write_prot_fa(
            out_dir,
            {"GENE_001.t1": prot_plus, "GENE_002.t1": prot_minus},
        )
        # Include GENE_003.t1 only in cdna
        _write_cdna_fa(
            out_dir,
            {
                "GENE_001.t1": _reconstruct_cdna_plus(),
                "GENE_002.t1": _reconstruct_cdna_minus(),
                "GENE_003.t1": _reconstruct_cdna_nocds(),
            },
        )

        report = validate_fasta(str(out_dir))
        assert report["pass"]
        assert report["missing_proteins_total"] == 0

    def test_no_duplicate_headers(self, tmp_path):
        """Duplicate FASTA headers should be detected."""
        out_dir = tmp_path / "output"
        out_dir.mkdir()
        _write_gff3(out_dir, include_nocds=False)

        prot_plus = _reconstruct_protein_plus()
        # Manually write duplicate
        prot_path = out_dir / "prot.fa"
        prot_path.write_text(
            f">GENE_001.t1\n{prot_plus}\n>GENE_001.t1\n{prot_plus}\n>GENE_002.t1\nMAA\n"
        )
        _write_cdna_fa(
            out_dir,
            {"GENE_001.t1": "ATGAAA", "GENE_002.t1": "ATGCCC"},
        )

        report = validate_fasta(str(out_dir))
        assert not report["pass"]
        assert len(report["duplicate_prot_headers"]) > 0

    def test_sequence_check_plus_strand(self, tmp_path):
        """Reconstructed protein should match for + strand transcript."""
        genome_path, genome = _write_genome(tmp_path)
        out_dir = tmp_path / "output"
        out_dir.mkdir()
        _write_gff3(out_dir, include_nocds=False)

        prot_plus = _reconstruct_protein_plus()
        prot_minus = _reconstruct_protein_minus()
        _write_prot_fa(
            out_dir,
            {"GENE_001.t1": prot_plus, "GENE_002.t1": prot_minus},
        )
        _write_cdna_fa(
            out_dir,
            {
                "GENE_001.t1": _reconstruct_cdna_plus(),
                "GENE_002.t1": _reconstruct_cdna_minus(),
            },
        )

        report = validate_fasta(str(out_dir), genome_path)
        assert report["pass"]
        sc = report["sequence_checks"]
        assert sc["protein_mismatches_total"] == 0
        assert sc["cdna_mismatches_total"] == 0


# ---------------------------------------------------------------------------
# C2: Regression test for ghost proteins (the bug we found)
# ---------------------------------------------------------------------------

class TestGhostProteinRegression:
    """C2: FASTA records for transcripts dropped by post-processing must be detected."""

    def test_extra_proteins_detected(self, tmp_path):
        """Ghost protein (in FASTA but not GFF) should cause QC failure."""
        out_dir = tmp_path / "output"
        out_dir.mkdir()

        # GFF has GENE_001.t1 and GENE_002.t1 only
        _write_gff3(out_dir, include_nocds=False)

        prot_plus = _reconstruct_protein_plus()
        prot_minus = _reconstruct_protein_minus()
        # Add a ghost protein: GENE_999.t1 not in GFF
        _write_prot_fa(
            out_dir,
            {
                "GENE_001.t1": prot_plus,
                "GENE_002.t1": prot_minus,
                "GENE_999.t1": "MAAAA",
            },
        )
        _write_cdna_fa(
            out_dir,
            {
                "GENE_001.t1": _reconstruct_cdna_plus(),
                "GENE_002.t1": _reconstruct_cdna_minus(),
                "GENE_999.t1": "ATGGCTGCTGCTGCTTAA",
            },
        )

        report = validate_fasta(str(out_dir))
        assert not report["pass"]
        assert report["extra_proteins_total"] == 1
        assert "GENE_999.t1" in report["extra_proteins"]
        assert report["extra_cdna_total"] == 1

    def test_missing_proteins_detected(self, tmp_path):
        """Missing protein (in GFF but not FASTA) should cause QC failure."""
        out_dir = tmp_path / "output"
        out_dir.mkdir()
        _write_gff3(out_dir, include_nocds=False)

        # Only include one of two CDS-bearing transcripts
        prot_plus = _reconstruct_protein_plus()
        _write_prot_fa(out_dir, {"GENE_001.t1": prot_plus})
        _write_cdna_fa(
            out_dir,
            {"GENE_001.t1": _reconstruct_cdna_plus()},
        )

        report = validate_fasta(str(out_dir))
        assert not report["pass"]
        assert report["missing_proteins_total"] == 1
        assert "GENE_002.t1" in report["missing_proteins"]

    def test_reconciliation_filter(self, tmp_path):
        """Simulate the fix: filtering FASTA to surviving IDs should pass QC."""
        out_dir = tmp_path / "output"
        out_dir.mkdir()
        _write_gff3(out_dir, include_nocds=False)

        # Start with ghost entries
        all_prot = {
            "GENE_001.t1": _reconstruct_protein_plus(),
            "GENE_002.t1": _reconstruct_protein_minus(),
            "GHOST.t1": "MAAAA",
        }
        all_cdna = {
            "GENE_001.t1": _reconstruct_cdna_plus(),
            "GENE_002.t1": _reconstruct_cdna_minus(),
            "GHOST.t1": "ATGGCTGCTTAA",
        }

        # Simulate the fix: filter to surviving GFF mRNA IDs
        gff_data = parse_gff3(str(out_dir / "consensus.gff3"))
        surviving = set(gff_data["transcripts"].keys())

        filtered_prot = {k: v for k, v in all_prot.items() if k in surviving}
        filtered_cdna = {k: v for k, v in all_cdna.items() if k in surviving}

        _write_prot_fa(out_dir, filtered_prot)
        _write_cdna_fa(out_dir, filtered_cdna)

        report = validate_fasta(str(out_dir))
        assert report["pass"]


# ---------------------------------------------------------------------------
# C3: fasta_qc module unit tests
# ---------------------------------------------------------------------------

class TestFastaQcHelpers:
    def test_parse_fasta_ids(self, tmp_path):
        fa = tmp_path / "test.fa"
        fa.write_text(">seq1 desc\nATG\n>seq2\nGCT\n")
        ids = parse_fasta_ids(str(fa))
        assert ids == ["seq1", "seq2"]

    def test_parse_fasta_ids_missing_file(self, tmp_path):
        ids = parse_fasta_ids(str(tmp_path / "missing.fa"))
        assert ids == []

    def test_parse_fasta_records(self, tmp_path):
        fa = tmp_path / "test.fa"
        fa.write_text(">s1\nATG\nGCT\n>s2\nAAA\n")
        records = parse_fasta_records(str(fa))
        assert records == {"s1": "ATGGCT", "s2": "AAA"}

    def test_parse_gff3_counts(self, tmp_path):
        gff_path = _write_gff3(tmp_path)
        data = parse_gff3(gff_path)
        assert data["n_genes"] == 3
        assert data["n_transcripts"] == 3
        assert data["n_cds_transcripts"] == 2

    def test_validate_fasta_no_gff(self, tmp_path):
        report = validate_fasta(str(tmp_path))
        assert "error" in report


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
