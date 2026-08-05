"""Tests for finalisation workflow, FASTA regeneration, FASTA QC, and canonical GFF3."""

from __future__ import annotations

import json
import os
import tempfile

import pytest


# ---------------------------------------------------------------------------
# FASTA regeneration (Part 4)
# ---------------------------------------------------------------------------
class TestRegenerateFinalFasta:
    """Test that regenerate_final_fasta derives correct sequences from GFF3 rows."""

    @staticmethod
    def _make_rows_and_genome():
        genome = {"chr1": "AAAAACCCCCGGGGGTTTTTNNNNN" * 10}
        rows = [
            {
                "Feature": "gene",
                "Start": 0,
                "End": 50,
                "Chromosome": "chr1",
                "Strand": "+",
                "ID": "G001",
                "Parent": "",
                "Source": "GMB",
                "Score": ".",
                "Frame": ".",
            },
            {
                "Feature": "mRNA",
                "Start": 0,
                "End": 50,
                "Chromosome": "chr1",
                "Strand": "+",
                "ID": "G001.t1",
                "Parent": "G001",
                "Source": "GMB",
                "Score": ".",
                "Frame": ".",
            },
            {
                "Feature": "exon",
                "Start": 0,
                "End": 30,
                "Chromosome": "chr1",
                "Strand": "+",
                "ID": "G001.t1.exon1",
                "Parent": "G001.t1",
                "Source": "GMB",
                "Score": ".",
                "Frame": ".",
            },
            {
                "Feature": "exon",
                "Start": 40,
                "End": 50,
                "Chromosome": "chr1",
                "Strand": "+",
                "ID": "G001.t1.exon2",
                "Parent": "G001.t1",
                "Source": "GMB",
                "Score": ".",
                "Frame": ".",
            },
            {
                "Feature": "CDS",
                "Start": 5,
                "End": 30,
                "Chromosome": "chr1",
                "Strand": "+",
                "ID": "G001.t1.cds1",
                "Parent": "G001.t1",
                "Source": "GMB",
                "Score": ".",
                "Frame": "0",
            },
            {
                "Feature": "CDS",
                "Start": 40,
                "End": 49,
                "Chromosome": "chr1",
                "Strand": "+",
                "ID": "G001.t1.cds2",
                "Parent": "G001.t1",
                "Source": "GMB",
                "Score": ".",
                "Frame": "1",
            },
        ]
        return rows, genome

    def test_produces_all_three_files(self):
        from gmb.pipeline.builder import regenerate_final_fasta

        rows, genome = self._make_rows_and_genome()
        with tempfile.TemporaryDirectory() as d:
            stats = regenerate_final_fasta(rows, genome, d)
            assert os.path.exists(os.path.join(d, "cdna.fa"))
            assert os.path.exists(os.path.join(d, "cds.fa"))
            assert os.path.exists(os.path.join(d, "prot.fa"))
            assert stats["cdna"] == 1
            assert stats["cds"] == 1
            assert stats["prot"] == 1

    def test_cdna_matches_exon_coords(self):
        from gmb.pipeline.annotate_cds_utrs import build_spliced_seq
        from gmb.pipeline.builder import regenerate_final_fasta

        rows, genome = self._make_rows_and_genome()
        with tempfile.TemporaryDirectory() as d:
            regenerate_final_fasta(rows, genome, d)
            from gmb.utils.fasta import parse_fasta_records

            cdna = parse_fasta_records(os.path.join(d, "cdna.fa"))
            expected = build_spliced_seq([(0, 30), (40, 50)], "+", genome["chr1"])
            assert cdna["G001.t1"] == expected

    def test_cds_matches_cds_coords(self):
        from gmb.pipeline.builder import regenerate_final_fasta

        rows, genome = self._make_rows_and_genome()
        with tempfile.TemporaryDirectory() as d:
            regenerate_final_fasta(rows, genome, d)
            from gmb.utils.fasta import parse_fasta_records

            cds = parse_fasta_records(os.path.join(d, "cds.fa"))
            expected = genome["chr1"][5:30] + genome["chr1"][40:49]
            assert cds["G001.t1"] == expected

    def test_minus_strand(self):
        from gmb.pipeline.annotate_cds_utrs import reverse_complement
        from gmb.pipeline.builder import regenerate_final_fasta

        genome = {"chr1": "ATCGATCGATCGATCGATCG" * 5}
        rows = [
            {
                "Feature": "gene",
                "Start": 0,
                "End": 40,
                "Chromosome": "chr1",
                "Strand": "-",
                "ID": "G002",
                "Parent": "",
                "Source": "GMB",
                "Score": ".",
                "Frame": ".",
            },
            {
                "Feature": "mRNA",
                "Start": 0,
                "End": 40,
                "Chromosome": "chr1",
                "Strand": "-",
                "ID": "G002.t1",
                "Parent": "G002",
                "Source": "GMB",
                "Score": ".",
                "Frame": ".",
            },
            {
                "Feature": "exon",
                "Start": 0,
                "End": 20,
                "Chromosome": "chr1",
                "Strand": "-",
                "ID": "G002.t1.exon1",
                "Parent": "G002.t1",
                "Source": "GMB",
                "Score": ".",
                "Frame": ".",
            },
            {
                "Feature": "exon",
                "Start": 30,
                "End": 40,
                "Chromosome": "chr1",
                "Strand": "-",
                "ID": "G002.t1.exon2",
                "Parent": "G002.t1",
                "Source": "GMB",
                "Score": ".",
                "Frame": ".",
            },
            {
                "Feature": "CDS",
                "Start": 5,
                "End": 20,
                "Chromosome": "chr1",
                "Strand": "-",
                "ID": "G002.t1.cds1",
                "Parent": "G002.t1",
                "Source": "GMB",
                "Score": ".",
                "Frame": "0",
            },
            {
                "Feature": "CDS",
                "Start": 30,
                "End": 39,
                "Chromosome": "chr1",
                "Strand": "-",
                "ID": "G002.t1.cds2",
                "Parent": "G002.t1",
                "Source": "GMB",
                "Score": ".",
                "Frame": "0",
            },
        ]
        with tempfile.TemporaryDirectory() as d:
            regenerate_final_fasta(rows, genome, d)
            from gmb.utils.fasta import parse_fasta_records

            cds = parse_fasta_records(os.path.join(d, "cds.fa"))
            expected = reverse_complement(genome["chr1"][5:20] + genome["chr1"][30:39])
            assert cds["G002.t1"] == expected

    def test_no_removed_transcripts(self):
        from gmb.pipeline.builder import regenerate_final_fasta

        genome = {"chr1": "A" * 100}
        rows = [
            {
                "Feature": "gene",
                "Start": 0,
                "End": 50,
                "Chromosome": "chr1",
                "Strand": "+",
                "ID": "G001",
                "Parent": "",
                "Source": "GMB",
                "Score": ".",
                "Frame": ".",
            },
            {
                "Feature": "mRNA",
                "Start": 0,
                "End": 50,
                "Chromosome": "chr1",
                "Strand": "+",
                "ID": "G001.t1",
                "Parent": "G001",
                "Source": "GMB",
                "Score": ".",
                "Frame": ".",
            },
            {
                "Feature": "exon",
                "Start": 0,
                "End": 50,
                "Chromosome": "chr1",
                "Strand": "+",
                "ID": "G001.t1.exon1",
                "Parent": "G001.t1",
                "Source": "GMB",
                "Score": ".",
                "Frame": ".",
            },
        ]
        with tempfile.TemporaryDirectory() as d:
            stats = regenerate_final_fasta(rows, genome, d)
            assert stats["cdna"] == 1
            assert stats["cds"] == 0
            assert stats["prot"] == 0


# ---------------------------------------------------------------------------
# FASTA QC pass/fail (Part 5)
# ---------------------------------------------------------------------------
class TestFastaQcPassFail:

    def test_zero_mismatches_pass(self):
        from gmb.pipeline.fasta_qc import run_coverage_checks, run_sequence_checks

        gff_data = {
            "transcripts": {
                "t1": {
                    "chrom": "c",
                    "strand": "+",
                    "exons": [(0, 10)],
                    "cds": [(0, 9)],
                    "cds_phases": [0],
                }
            },
            "n_genes": 1,
            "n_transcripts": 1,
            "n_cds_transcripts": 1,
            "cds_transcript_ids": {"t1"},
            "genes": {"g1": {}},
        }
        genome = {"c": "ATGATGATGA" + "A" * 90}
        from gmb.pipeline.annotate_cds_utrs import build_spliced_seq, reverse_complement, translate

        cdna = build_spliced_seq([(0, 10)], "+", genome["c"])
        cds_nuc = genome["c"][0:9]
        prot = translate(cds_nuc)
        if prot.endswith("*"):
            prot = prot[:-1]

        report = run_coverage_checks(gff_data, ["t1"], ["t1"])
        sc = run_sequence_checks(
            gff_data, {"t1": prot}, {"t1": cdna}, genome, cds_records={"t1": cds_nuc}
        )
        report["sequence_checks"] = sc

        assert sc["cdna_mismatches_total"] == 0
        assert sc["protein_mismatches_total"] == 0
        assert sc["cds_mismatches_total"] == 0

    def test_cdna_mismatch_fails(self):
        report = {
            "missing_proteins_total": 0,
            "extra_proteins_total": 0,
            "missing_cdna_total": 0,
            "extra_cdna_total": 0,
            "duplicate_prot_headers": {},
            "duplicate_cdna_headers": {},
            "sequence_checks": {
                "cdna_mismatches_total": 1,
                "protein_mismatches_total": 0,
                "cds_mismatches_total": 0,
                "internal_stops_total": 0,
            },
        }
        failed = []
        if report["sequence_checks"]["cdna_mismatches_total"] > 0:
            failed.append("cdna_mismatches")
        report["pass"] = len(failed) == 0
        assert not report["pass"]

    def test_protein_mismatch_fails(self):
        report = {
            "missing_proteins_total": 0,
            "extra_proteins_total": 0,
            "missing_cdna_total": 0,
            "extra_cdna_total": 0,
            "duplicate_prot_headers": {},
            "duplicate_cdna_headers": {},
            "sequence_checks": {
                "cdna_mismatches_total": 0,
                "protein_mismatches_total": 1,
                "cds_mismatches_total": 0,
                "internal_stops_total": 0,
            },
        }
        failed = []
        if report["sequence_checks"]["protein_mismatches_total"] > 0:
            failed.append("protein_mismatches")
        report["pass"] = len(failed) == 0
        assert not report["pass"]

    def test_internal_stop_fails(self):
        report = {
            "missing_proteins_total": 0,
            "extra_proteins_total": 0,
            "missing_cdna_total": 0,
            "extra_cdna_total": 0,
            "duplicate_prot_headers": {},
            "duplicate_cdna_headers": {},
            "sequence_checks": {
                "cdna_mismatches_total": 0,
                "protein_mismatches_total": 0,
                "cds_mismatches_total": 0,
                "internal_stops_total": 1,
            },
        }
        failed = []
        if report["sequence_checks"]["internal_stops_total"] > 0:
            failed.append("internal_stops")
        report["pass"] = len(failed) == 0
        assert not report["pass"]

    def test_missing_fasta_id_fails(self):
        report = {
            "missing_proteins_total": 1,
            "extra_proteins_total": 0,
            "missing_cdna_total": 0,
            "extra_cdna_total": 0,
            "duplicate_prot_headers": {},
            "duplicate_cdna_headers": {},
        }
        failed = []
        if report["missing_proteins_total"] > 0:
            failed.append("missing_proteins")
        report["pass"] = len(failed) == 0
        assert not report["pass"]

    def test_duplicate_header_fails(self):
        report = {
            "missing_proteins_total": 0,
            "extra_proteins_total": 0,
            "missing_cdna_total": 0,
            "extra_cdna_total": 0,
            "duplicate_prot_headers": {"t1": 2},
            "duplicate_cdna_headers": {},
        }
        failed = []
        if report["duplicate_prot_headers"]:
            failed.append("duplicate_prot_headers")
        report["pass"] = len(failed) == 0
        assert not report["pass"]


# ---------------------------------------------------------------------------
# Canonical GFF3 annotation (Part 8)
# ---------------------------------------------------------------------------
class TestCanonicalGff3Annotation:

    def _write_gff3(self, path, lines):
        with open(path, "w") as fh:
            for line in lines:
                fh.write(line + "\n")

    def test_one_tagged_mrna_per_gene(self):
        from gmb.pipeline.canonical_selection import _write_annotated_gff3

        gff_lines = [
            "##gff-version 3",
            "chr1\tGMB\tgene\t1\t1000\t.\t+\t.\tID=G1",
            "chr1\tGMB\tmRNA\t1\t1000\t.\t+\t.\tID=G1.t1;Parent=G1",
            "chr1\tGMB\tmRNA\t1\t1000\t.\t+\t.\tID=G1.t2;Parent=G1",
        ]
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "in.gff3")
            out = os.path.join(d, "out.gff3")
            self._write_gff3(src, gff_lines)
            _write_annotated_gff3(src, {"G1.t1"}, out)

            with open(out) as fh:
                lines = fh.readlines()
            tagged = [l for l in lines if "Ensembl_canonical" in l]
            assert len(tagged) == 1
            assert "G1.t1" in tagged[0]
            alt = [l for l in lines if "G1.t2" in l]
            assert len(alt) == 1
            assert "Ensembl_canonical" not in alt[0]

    def test_existing_tag_preserved(self):
        from gmb.pipeline.canonical_selection import _write_annotated_gff3

        gff_lines = [
            "##gff-version 3",
            "chr1\tGMB\tmRNA\t1\t100\t.\t+\t.\tID=T1;Parent=G1;tag=basic",
        ]
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "in.gff3")
            out = os.path.join(d, "out.gff3")
            self._write_gff3(src, gff_lines)
            _write_annotated_gff3(src, {"T1"}, out)

            with open(out) as fh:
                content = fh.read()
            assert "tag=basic,Ensembl_canonical" in content

    def test_no_duplicate_ensembl_canonical(self):
        from gmb.pipeline.canonical_selection import _write_annotated_gff3

        gff_lines = [
            "##gff-version 3",
            "chr1\tGMB\tmRNA\t1\t100\t.\t+\t.\tID=T1;Parent=G1;tag=Ensembl_canonical",
        ]
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "in.gff3")
            out = os.path.join(d, "out.gff3")
            self._write_gff3(src, gff_lines)
            _write_annotated_gff3(src, {"T1"}, out)

            with open(out) as fh:
                content = fh.read()
            assert content.count("Ensembl_canonical") == 1

    def test_non_mrna_lines_unchanged(self):
        from gmb.pipeline.canonical_selection import _write_annotated_gff3

        gff_lines = [
            "##gff-version 3",
            "chr1\tGMB\tgene\t1\t100\t.\t+\t.\tID=G1",
            "chr1\tGMB\texon\t1\t100\t.\t+\t.\tID=E1;Parent=T1",
        ]
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "in.gff3")
            out = os.path.join(d, "out.gff3")
            self._write_gff3(src, gff_lines)
            _write_annotated_gff3(src, set(), out)

            with open(src) as f1, open(out) as f2:
                assert f1.read() == f2.read()
