"""Regression tests for cross-chromosome gene_id/transcript_id collisions.

Tiberius restarts gene numbering ("g1", "g2", ...) independently on every
chromosome, so the same gene_id/transcript_id string can refer to entirely
different genes. The loader must namespace these internally so distinct
per-chromosome genes/transcripts are not collapsed together, while leaving
genome-wide-unique IDs (the normal case for Ensembl GTF/GFF3) untouched.
"""

import os

import pytest

from gmb.compare.annotation_loader import load_annotation


TIBERIUS_STYLE_GTF = """\
chr1\tTiberius\tgene\t100\t500\t.\t+\t.\tgene_id "g1";
chr1\tTiberius\ttranscript\t100\t500\t.\t+\t.\tgene_id "g1"; transcript_id "g1.t1";
chr1\tTiberius\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "g1.t1";
chr1\tTiberius\texon\t300\t500\t.\t+\t.\tgene_id "g1"; transcript_id "g1.t1";
chr1\tTiberius\tCDS\t100\t200\t.\t+\t0\tgene_id "g1"; transcript_id "g1.t1";
chr1\tTiberius\tCDS\t300\t500\t.\t+\t0\tgene_id "g1"; transcript_id "g1.t1";
chr2\tTiberius\tgene\t9000\t9400\t.\t-\t.\tgene_id "g1";
chr2\tTiberius\ttranscript\t9000\t9400\t.\t-\t.\tgene_id "g1"; transcript_id "g1.t1";
chr2\tTiberius\texon\t9000\t9400\t.\t-\t.\tgene_id "g1"; transcript_id "g1.t1";
chr2\tTiberius\tCDS\t9000\t9400\t.\t-\t0\tgene_id "g1"; transcript_id "g1.t1";
"""

NORMAL_GTF = """\
chr1\tEnsembl\tgene\t100\t500\t.\t+\t.\tgene_id "ENSG001";
chr1\tEnsembl\ttranscript\t100\t500\t.\t+\t.\tgene_id "ENSG001"; transcript_id "ENST001";
chr1\tEnsembl\texon\t100\t200\t.\t+\t.\tgene_id "ENSG001"; transcript_id "ENST001";
chr1\tEnsembl\texon\t300\t500\t.\t+\t.\tgene_id "ENSG001"; transcript_id "ENST001";
chr2\tEnsembl\tgene\t9000\t9400\t.\t-\t.\tgene_id "ENSG002";
chr2\tEnsembl\ttranscript\t9000\t9400\t.\t-\t.\tgene_id "ENSG002"; transcript_id "ENST002";
chr2\tEnsembl\texon\t9000\t9400\t.\t-\t.\tgene_id "ENSG002"; transcript_id "ENST002";
"""


@pytest.fixture
def tiberius_style_gtf(tmp_path):
    path = tmp_path / "tiberius_style.gtf"
    path.write_text(TIBERIUS_STYLE_GTF)
    return str(path)


@pytest.fixture
def normal_gtf(tmp_path):
    path = tmp_path / "normal.gtf"
    path.write_text(NORMAL_GTF)
    return str(path)


def test_duplicate_gene_ids_across_chromosomes_become_distinct(tiberius_style_gtf):
    genes, mrna, exons, cds = load_annotation(tiberius_style_gtf, "Tiberius")

    # Feature counts are preserved exactly (no rows dropped or merged).
    assert len(genes) == 2
    assert len(mrna) == 2
    assert len(exons) == 3
    assert len(cds) == 3

    # The two "g1" genes on chr1/chr2 must be recognised as distinct genes.
    assert genes["gene_id"].nunique() == 2
    assert mrna["transcript_id"].nunique() == 2

    # Original IDs are preserved for debugging.
    assert set(genes["original_gene_id"]) == {"g1"}
    assert set(mrna["original_transcript_id"]) == {"g1.t1"}

    # Each gene's internal ID is namespaced by its chromosome.
    chr1_gene = genes[genes["Chromosome"] == "chr1"].iloc[0]
    chr2_gene = genes[genes["Chromosome"] == "chr2"].iloc[0]
    assert chr1_gene["gene_id"] != chr2_gene["gene_id"]
    assert chr1_gene["gene_id"].endswith("g1")
    assert chr2_gene["gene_id"].endswith("g1")

    # Exons/CDS remain correctly linked to their own chromosome's gene/transcript.
    chr1_exons = exons[exons["Chromosome"] == "chr1"]
    chr2_exons = exons[exons["Chromosome"] == "chr2"]
    assert set(chr1_exons["gene_id"]) == {chr1_gene["gene_id"]}
    assert set(chr2_exons["gene_id"]) == {chr2_gene["gene_id"]}
    assert chr1_exons["transcript_id"].nunique() == 1
    assert chr2_exons["transcript_id"].nunique() == 1
    assert set(chr1_exons["transcript_id"]) != set(chr2_exons["transcript_id"])

    chr1_cds = cds[cds["Chromosome"] == "chr1"]
    chr2_cds = cds[cds["Chromosome"] == "chr2"]
    assert len(chr1_cds) == 2
    assert len(chr2_cds) == 1
    assert set(chr1_cds["gene_id"]) == {chr1_gene["gene_id"]}
    assert set(chr2_cds["gene_id"]) == {chr2_gene["gene_id"]}


def test_globally_unique_ids_are_left_untouched(normal_gtf):
    """Normal Ensembl-style GTF with genome-wide unique IDs must not be rewritten."""
    genes, mrna, exons, cds = load_annotation(normal_gtf, "Reference")

    assert len(genes) == 2
    assert set(genes["gene_id"]) == {"ENSG001", "ENSG002"}
    assert (genes["gene_id"] == genes["original_gene_id"]).all()

    assert set(mrna["transcript_id"]) == {"ENST001", "ENST002"}
    assert (mrna["transcript_id"] == mrna["original_transcript_id"]).all()


def test_real_tiberius_fixture_has_collisions_if_present():
    """Sanity-check against the real Tiberius output, when available locally.

    This isn't a hard dependency of the suite (the large GTF lives outside
    the repo), but if it's present we confirm the loader collapses the
    ~20x per-chromosome ID reuse into fully-unique internal IDs while
    preserving every gene/transcript/exon/CDS row.
    """
    candidate = os.path.join(
        os.path.dirname(__file__), "..", "..", "..", "..",
        "tiberius_annotations", "GCA_009914755.4_tiberius.gtf",
    )
    candidate = os.path.abspath(candidate)
    if not os.path.exists(candidate):
        pytest.skip("Real Tiberius GTF fixture not available in this environment")

    genes, mrna, exons, cds = load_annotation(candidate, "Tiberius")

    assert len(genes) == 22941
    assert len(mrna) == 22941
    assert len(exons) == 205646
    assert len(cds) == 205646

    # Internal IDs are fully unique even though original IDs repeat heavily.
    assert genes["gene_id"].nunique() == len(genes)
    assert mrna["transcript_id"].nunique() == len(mrna)

    # Original IDs reused across chromosomes, confirming the fixture still
    # exhibits the collision this fix targets.
    assert genes["original_gene_id"].nunique() < len(genes)
