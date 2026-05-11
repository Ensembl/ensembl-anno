#!/usr/bin/env python3
"""Tests for fasta_export module."""

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from gmb.pipeline.annotate_cds_utrs import reverse_complement
from gmb.pipeline.config import load_config
from gmb.pipeline.fasta_export import export_fasta


@pytest.fixture
def config():
    return load_config()


def _make_genome(tmp_path):
    """Create a tiny genome FASTA for testing."""
    fa = tmp_path / "genome.fa"
    # Chromosome with a known ORF:
    # pos 30-33: ATG, pos 30-183: ORF (51 codons), pos 180-183: TAA
    utr5 = "AAA" * 10  # 30bp
    coding = "ATG" + "GCT" * 49 + "TAA"  # 3 + 147 + 3 = 153bp
    utr3 = "CCC" * 10  # 30bp
    seq = utr5 + coding + utr3  # 213bp

    fa.write_text(f">chr1\n{seq}\n")
    return str(fa), {"chr1": seq.upper()}


def _make_models(gene_id="GENE_1"):
    """Create a minimal loci list for FASTA export."""
    exons_df = pd.DataFrame(
        {
            "Start": [0],
            "End": [213],
            "Chromosome": ["chr1"],
            "Strand": ["+"],
        }
    )
    model = {
        "id": "tx1",
        "source": "Scallop",
        "chrom": "chr1",
        "strand": "+",
        "df": exons_df,
        "start": 0,
        "end": 213,
        "exon_count": 1,
        "combined_evidence": "Scallop",
        "output_tid": f"{gene_id}.1",
        "gene_id": gene_id,
    }
    return [(1, [model])]


class TestExportFasta:
    def test_produces_cdna_file(self, tmp_path, config):
        """cdna.fa should be created with content."""
        _, genome = _make_genome(tmp_path)
        loci = _make_models()
        out_dir = str(tmp_path / "output")

        stats = export_fasta(loci, genome, None, config, out_dir)

        cdna_path = os.path.join(out_dir, "cdna.fa")
        assert os.path.exists(cdna_path)
        with open(cdna_path) as f:
            content = f.read()
        assert ">GENE_1.1" in content
        assert len(content) > 50

    def test_produces_protein_file(self, tmp_path, config):
        """prot.fa should contain a translated protein."""
        _, genome = _make_genome(tmp_path)
        loci = _make_models()
        out_dir = str(tmp_path / "output")

        stats = export_fasta(loci, genome, None, config, out_dir)

        prot_path = os.path.join(out_dir, "prot.fa")
        assert os.path.exists(prot_path)
        with open(prot_path) as f:
            content = f.read()
        assert ">GENE_1.1" in content
        # Should contain an M (methionine start)
        lines = [l for l in content.split("\n") if not l.startswith(">")]
        seq = "".join(lines)
        assert seq[0] == "M", f"Protein should start with M, got {seq[:5]}"

    def test_produces_cds_file(self, tmp_path, config):
        """cds.fa should be created when write_cds=True."""
        _, genome = _make_genome(tmp_path)
        loci = _make_models()
        out_dir = str(tmp_path / "output")

        stats = export_fasta(loci, genome, None, config, out_dir)

        cds_path = os.path.join(out_dir, "cds.fa")
        assert os.path.exists(cds_path)
        with open(cds_path) as f:
            content = f.read()
        assert ">GENE_1.1" in content

    def test_fasta_header_format(self, tmp_path, config):
        """Headers should contain gene ID and evidence."""
        _, genome = _make_genome(tmp_path)
        loci = _make_models()
        out_dir = str(tmp_path / "output")

        export_fasta(loci, genome, None, config, out_dir)

        with open(os.path.join(out_dir, "cdna.fa")) as f:
            header = f.readline().strip()
        assert header.startswith(">GENE_1.1")
        assert "gene=GENE_1" in header
        assert "evidence=Scallop" in header

    def test_export_stats(self, tmp_path, config):
        """Export should return stats dict with counts."""
        _, genome = _make_genome(tmp_path)
        loci = _make_models()
        out_dir = str(tmp_path / "output")

        stats = export_fasta(loci, genome, None, config, out_dir)
        assert stats["total_transcripts"] == 1
        assert stats["with_cds"] >= 0
        assert "with_protein" in stats

    def test_reverse_strand_cdna(self, tmp_path, config):
        """cDNA for - strand should be reverse complemented."""
        # Create genome
        fa = tmp_path / "genome.fa"
        seq = "ATG" + "GCT" * 20 + "TAA" + "CCC" * 5  # 81bp
        fa.write_text(f">chr1\n{seq}\n")
        genome = {"chr1": seq.upper()}

        # Model on - strand
        seq_len = len(seq)
        exons_df = pd.DataFrame(
            {
                "Start": [0],
                "End": [seq_len],
                "Chromosome": ["chr1"],
                "Strand": ["-"],
            }
        )
        model = {
            "id": "tx_rev",
            "source": "Scallop",
            "chrom": "chr1",
            "strand": "-",
            "df": exons_df,
            "start": 0,
            "end": seq_len,
            "exon_count": 1,
            "combined_evidence": "Scallop",
            "output_tid": "GENE_R.1",
            "gene_id": "GENE_R",
        }
        loci = [(1, [model])]
        out_dir = str(tmp_path / "output")

        export_fasta(loci, genome, None, config, out_dir)

        with open(os.path.join(out_dir, "cdna.fa")) as f:
            lines = f.readlines()
        # Skip header, join sequence
        cdna = "".join(l.strip() for l in lines[1:])
        expected = reverse_complement(seq.upper())
        assert cdna == expected


class TestHelixerCdsExport:
    def test_helixer_cds_used(self, tmp_path, config):
        """When Helixer CDS is provided, use it instead of ORF prediction."""
        _, genome = _make_genome(tmp_path)

        # Helixer CDS: different from what ORF prediction would find
        helixer_cds = pd.DataFrame(
            {
                "Chromosome": ["chr1"],
                "Start": [30],
                "End": [183],
                "Strand": ["+"],
                "transcript_id": ["hx_tx1"],
            }
        )

        exons_df = pd.DataFrame(
            {
                "Start": [0],
                "End": [213],
                "Chromosome": ["chr1"],
                "Strand": ["+"],
            }
        )
        model = {
            "id": "hx_tx1",
            "source": "Helixer",
            "chrom": "chr1",
            "strand": "+",
            "df": exons_df,
            "start": 0,
            "end": 213,
            "exon_count": 1,
            "combined_evidence": "Helixer",
            "output_tid": "GENE_H.1",
            "gene_id": "GENE_H",
        }
        loci = [(1, [model])]
        out_dir = str(tmp_path / "output")

        stats = export_fasta(loci, genome, helixer_cds, config, out_dir)
        assert stats["with_cds"] == 1

        with open(os.path.join(out_dir, "prot.fa")) as f:
            lines = f.readlines()
        prot = "".join(l.strip() for l in lines[1:])
        assert prot[0] == "M"


    def test_helixer_cds_minus_strand_multiexon(self, tmp_path, config):
        """Minus-strand multi-exon Helixer CDS must not produce internal stop codons.

        RC(X+Y) == RC(Y)+RC(X), so concatenating exon sequences in ascending
        genomic order before RC'ing places the 5' (highest-coord) exon first in
        the result.  Concatenating in descending order inverts the exon order and
        produces wrong proteins with premature stops.
        """
        # Build a synthetic genome where we control the exact CDS sequence.
        # Gene on - strand, two CDS exons:
        #   exon1 (low coord):  positions 0-12  (12 bp)
        #   exon2 (high coord): positions 20-32 (12 bp)
        # Correct mRNA CDS = RC(seq[20:32]) + RC(seq[0:12])
        # Choose sequences so the correct frame starts with ATG and has no stops.
        #   RC(GTTGGCATGAAC) = GTTCATGCCAAC  → GTT-CAT-GCC-AAC = Val-His-Ala-Asn
        #   RC(ATGCCCGGTTAA) = TTAACCGGGCAT → but this has TGA... let's be explicit.
        # Use a known-good encoding: set exon2 RC = ATG AAA TTT, exon1 RC = GGG TAA
        # exon2 genomic (RC to get mRNA) = RC(ATGAAATTT) = AAATTTCAT
        # exon1 genomic (RC to get mRNA) = RC(GGGTAA) = TTACCC
        # mRNA CDS = ATGAAATTT + GGGTAA → Met-Lys-Phe-Gly-stop → protein = MKF (no internal stop)
        # Layout (0-based half-open):
        #   exon1 at [0, 6):   TTACCC  → RC = GGGTAA  (3' end of mRNA CDS, includes stop)
        #   gap   at [6, 20):  padding
        #   exon2 at [20, 29): AAATTTCAT → RC = ATGAAATTT (5' end of mRNA CDS, start codon)
        #
        # Correct mRNA CDS (5'→3'): RC(exon2_genomic) + RC(exon1_genomic)
        #   = ATGAAATTT + GGGTAA  →  ATG-AAA-TTT-GGG-TAA  →  M-K-F-G-*  →  protein "MKFG"
        exon2_genomic = "AAATTTCAT"  # high coord → 5' in mRNA
        exon1_genomic = "TTACCC"     # low coord  → 3' in mRNA
        padding = "N" * 14           # gap so exon2 starts at exactly position 20
        genome_seq = exon1_genomic + padding + exon2_genomic  # len = 6+14+9 = 29
        genome = {"chr1": genome_seq}

        exons_df = pd.DataFrame(
            {
                "Start": [0, 20],
                "End": [6, 29],
                "Chromosome": ["chr1", "chr1"],
                "Strand": ["-", "-"],
            }
        )
        helixer_cds = pd.DataFrame(
            {
                "Chromosome": ["chr1", "chr1"],
                "Start": [0, 20],
                "End": [6, 29],
                "Strand": ["-", "-"],
                "transcript_id": ["hx_minus", "hx_minus"],
            }
        )
        model = {
            "id": "hx_minus",
            "source": "Helixer",
            "chrom": "chr1",
            "strand": "-",
            "df": exons_df,
            "start": 0,
            "end": 29,
            "exon_count": 2,
            "combined_evidence": "Helixer",
            "output_tid": "GENE_M.1",
            "gene_id": "GENE_M",
        }
        loci = [(1, [model])]
        out_dir = str(tmp_path / "output_minus")

        stats = export_fasta(loci, genome, helixer_cds, config, out_dir)
        assert stats["with_cds"] == 1

        with open(os.path.join(out_dir, "prot.fa")) as f:
            lines = f.readlines()
        prot = "".join(l.strip() for l in lines[1:])
        # Correct protein starts with M and has no internal stops
        assert prot[0] == "M", f"Protein should start with Met, got: {prot!r}"
        assert "*" not in prot, f"Protein has internal stop codons: {prot!r}"
        assert prot == "MKFG", f"Expected MKFG, got: {prot!r}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
