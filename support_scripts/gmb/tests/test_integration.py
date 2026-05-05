#!/usr/bin/env python3
"""Integration test: run the pipeline on a tiny synthetic dataset.

Creates:
  - A toy genome (1 chromosome, ~500bp)
  - A Scallop GTF (1 transcript)
  - A StringTie GTF (1 transcript)
  - A Helixer GFF3 (1 gene with CDS)
  - A protein evidence GTF (1 mapping)

Runs the pipeline and verifies outputs exist and are well-formed.
"""

import os
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))

# Check if pyranges is available
try:
    import pyranges
except ImportError:
    pytest.skip("pyranges required", allow_module_level=True)


def _create_test_data(tmpdir):
    """Create minimal test input files."""
    # Genome: 1 chromosome, 500 bp
    utr5 = "AAA" * 10  # 30bp
    coding = "ATG" + "GCT" * 49 + "TAA"  # 153bp
    utr3 = "CCC" * 20  # 60bp
    padding = "NNN" * 86  # 258bp (padding to fill to 501bp)
    seq = utr5 + coding + utr3 + padding  # 501bp

    genome_fa = os.path.join(tmpdir, "genome.fa")
    with open(genome_fa, "w") as f:
        f.write(f">1\n{seq}\n")

    # Scallop GTF: 1 transcript
    scallop_gtf = os.path.join(tmpdir, "scallop.gtf")
    with open(scallop_gtf, "w") as f:
        f.write(
            "1\tScallop\ttranscript\t1\t243\t.\t+\t.\t"
            'gene_id "sc_gene1"; transcript_id "sc_tx1";\n'
        )
        f.write(
            "1\tScallop\texon\t1\t243\t.\t+\t.\t" 'gene_id "sc_gene1"; transcript_id "sc_tx1";\n'
        )

    # StringTie GTF: 1 transcript (identical locus)
    stringtie_gtf = os.path.join(tmpdir, "stringtie.gtf")
    with open(stringtie_gtf, "w") as f:
        f.write(
            "1\tStringTie\ttranscript\t1\t243\t.\t+\t.\t"
            'gene_id "st_gene1"; transcript_id "st_tx1";\n'
        )
        f.write(
            "1\tStringTie\texon\t1\t243\t.\t+\t.\t" 'gene_id "st_gene1"; transcript_id "st_tx1";\n'
        )

    # Helixer GFF3: 1 gene with exon + CDS
    helixer_gff3 = os.path.join(tmpdir, "helixer.gff3")
    with open(helixer_gff3, "w") as f:
        f.write("##gff-version 3\n")
        f.write("1\tHelixer\tgene\t1\t243\t.\t+\t.\tID=hx_gene1\n")
        f.write("1\tHelixer\tmRNA\t1\t243\t.\t+\t.\t" "ID=hx_tx1;Parent=hx_gene1\n")
        f.write("1\tHelixer\texon\t1\t243\t.\t+\t.\t" "ID=hx_tx1.exon1;Parent=hx_tx1\n")
        # CDS at positions 31-183 (1-based)
        f.write("1\tHelixer\tCDS\t31\t183\t.\t+\t0\t" "ID=hx_tx1.cds1;Parent=hx_tx1\n")

    # Protein evidence GTF (OrthoDB): 1 mapping
    orthodb_gtf = os.path.join(tmpdir, "orthodb.gtf")
    with open(orthodb_gtf, "w") as f:
        f.write(
            "1\tOrthoDB\ttranscript\t31\t183\t.\t+\t.\t"
            'gene_id "odb_gene1"; transcript_id "odb_tx1";\n'
        )
        f.write(
            "1\tOrthoDB\texon\t31\t183\t.\t+\t.\t"
            'gene_id "odb_gene1"; transcript_id "odb_tx1";\n'
        )

    # UniProt GTF: 1 mapping
    uniprot_gtf = os.path.join(tmpdir, "uniprot.gtf")
    with open(uniprot_gtf, "w") as f:
        f.write(
            "1\tUniProt\ttranscript\t31\t183\t.\t+\t.\t"
            'gene_id "up_gene1"; transcript_id "up_tx1";\n'
        )
        f.write(
            "1\tUniProt\texon\t31\t183\t.\t+\t.\t" 'gene_id "up_gene1"; transcript_id "up_tx1";\n'
        )

    return {
        "genome": genome_fa,
        "scallop": scallop_gtf,
        "stringtie": stringtie_gtf,
        "helixer": helixer_gff3,
        "orthodb": orthodb_gtf,
        "uniprot": uniprot_gtf,
    }


# Config override that disables protein validation so the integration tests
# don't require diamond/psauron to be installed on the test machine.
_TEST_CONFIG = os.path.join(os.path.dirname(__file__), "fixtures", "test_no_pv.yaml")


def _add_test_config(cmd: list) -> list:
    """Append --config test_no_pv.yaml to a pipeline command list."""
    return cmd + ["--config", _TEST_CONFIG]


class TestIntegration:
    """Run the full pipeline on a tiny synthetic dataset."""

    def test_pipeline_runs(self, tmp_path):
        """Pipeline should complete without errors."""
        tmpdir = str(tmp_path)
        files = _create_test_data(tmpdir)
        output_dir = os.path.join(tmpdir, "output")

        script = os.path.join(os.path.dirname(os.path.dirname(__file__)), "gene_model_builder.py")
        cmd = [
            sys.executable,
            script,
            "--scallop",
            files["scallop"],
            "--stringtie",
            files["stringtie"],
            "--helixer",
            files["helixer"],
            "--orthodb",
            files["orthodb"],
            "--uniprot",
            files["uniprot"],
            "--genome",
            files["genome"],
            "--output-dir",
            output_dir,
            "--gene-prefix",
            "TEST",
        ]
        result = subprocess.run(_add_test_config(cmd), capture_output=True, text=True, cwd=os.path.dirname(__file__))

        assert result.returncode == 0, (
            f"Pipeline failed:\nstdout:\n{result.stdout}\n" f"stderr:\n{result.stderr}"
        )

    def test_output_files_exist(self, tmp_path):
        """Expected output files should be created."""
        tmpdir = str(tmp_path)
        files = _create_test_data(tmpdir)
        output_dir = os.path.join(tmpdir, "output")

        script = os.path.join(os.path.dirname(os.path.dirname(__file__)), "gene_model_builder.py")
        cmd = [
            sys.executable,
            script,
            "--scallop",
            files["scallop"],
            "--stringtie",
            files["stringtie"],
            "--helixer",
            files["helixer"],
            "--orthodb",
            files["orthodb"],
            "--uniprot",
            files["uniprot"],
            "--genome",
            files["genome"],
            "--output-dir",
            output_dir,
            "--gene-prefix",
            "TEST",
        ]
        subprocess.run(_add_test_config(cmd), capture_output=True, text=True, cwd=os.path.dirname(__file__))

        assert os.path.exists(os.path.join(output_dir, "consensus.gff3"))
        assert os.path.exists(os.path.join(output_dir, "cdna.fa"))
        assert os.path.exists(os.path.join(output_dir, "prot.fa"))
        assert os.path.exists(os.path.join(output_dir, "summary.json"))
        assert os.path.exists(os.path.join(output_dir, "summary.tsv"))

    def test_gff3_has_genes(self, tmp_path):
        """GFF3 should contain gene features."""
        tmpdir = str(tmp_path)
        files = _create_test_data(tmpdir)
        output_dir = os.path.join(tmpdir, "output")

        script = os.path.join(os.path.dirname(os.path.dirname(__file__)), "gene_model_builder.py")
        cmd = [
            sys.executable,
            script,
            "--scallop",
            files["scallop"],
            "--stringtie",
            files["stringtie"],
            "--helixer",
            files["helixer"],
            "--orthodb",
            files["orthodb"],
            "--uniprot",
            files["uniprot"],
            "--genome",
            files["genome"],
            "--output-dir",
            output_dir,
            "--gene-prefix",
            "TEST",
        ]
        script_dir = os.path.dirname(__file__)
        env = os.environ.copy()
        env["PYTHONPATH"] = os.path.dirname(script_dir)
        subprocess.run(_add_test_config(cmd), capture_output=True, text=True, cwd=script_dir, env=env)

        gff3_path = os.path.join(output_dir, "consensus.gff3")
        with open(gff3_path) as f:
            content = f.read()

        assert "gene" in content
        assert "mRNA" in content
        assert "exon" in content

    def test_protein_fasta_has_protein(self, tmp_path):
        """prot.fa should contain at least one protein sequence."""
        tmpdir = str(tmp_path)
        files = _create_test_data(tmpdir)
        output_dir = os.path.join(tmpdir, "output")

        script = os.path.join(os.path.dirname(os.path.dirname(__file__)), "gene_model_builder.py")
        cmd = [
            sys.executable,
            script,
            "--scallop",
            files["scallop"],
            "--stringtie",
            files["stringtie"],
            "--helixer",
            files["helixer"],
            "--orthodb",
            files["orthodb"],
            "--uniprot",
            files["uniprot"],
            "--genome",
            files["genome"],
            "--output-dir",
            output_dir,
            "--gene-prefix",
            "TEST",
        ]
        script_dir = os.path.dirname(__file__)
        env = os.environ.copy()
        env["PYTHONPATH"] = os.path.dirname(script_dir)
        subprocess.run(_add_test_config(cmd), capture_output=True, text=True, cwd=script_dir, env=env)

        prot_path = os.path.join(output_dir, "prot.fa")
        with open(prot_path) as f:
            content = f.read()

        # Should have at least one header and one sequence
        assert content.count(">") >= 1
        lines = [l for l in content.split("\n") if l and not l.startswith(">")]
        protein_seq = "".join(lines)
        assert len(protein_seq) > 0

    def test_summary_report(self, tmp_path):
        """summary.json should have expected structure."""
        import json

        tmpdir = str(tmp_path)
        files = _create_test_data(tmpdir)
        output_dir = os.path.join(tmpdir, "output")

        script = os.path.join(os.path.dirname(os.path.dirname(__file__)), "gene_model_builder.py")
        cmd = [
            sys.executable,
            script,
            "--scallop",
            files["scallop"],
            "--stringtie",
            files["stringtie"],
            "--helixer",
            files["helixer"],
            "--orthodb",
            files["orthodb"],
            "--uniprot",
            files["uniprot"],
            "--genome",
            files["genome"],
            "--output-dir",
            output_dir,
            "--gene-prefix",
            "TEST",
        ]
        script_dir = os.path.dirname(__file__)
        env = os.environ.copy()
        env["PYTHONPATH"] = os.path.dirname(script_dir)
        subprocess.run(_add_test_config(cmd), capture_output=True, text=True, cwd=script_dir, env=env)

        with open(os.path.join(output_dir, "summary.json")) as f:
            report = json.load(f)

        assert "summary" in report
        assert "total_loci" in report["summary"]
        assert report["summary"]["total_loci"] >= 1
        assert "filtering" in report


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
