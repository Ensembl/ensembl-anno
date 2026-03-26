#!/usr/bin/env python3
"""Tests for fasta_export module."""

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from annotate_cds_utrs import reverse_complement
from config import load_config
from fasta_export import export_fasta


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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
