import os
import sys
import tempfile

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from gmb.compare.compare_annotations import classify_locus_pairs
from gmb.pipeline.subset_utils import load_seqname_map, remap_df_seqnames


def test_load_seqname_map():
    with tempfile.NamedTemporaryFile("w", delete=False) as f:
        f.write("# comment\n")
        f.write("from_seqname\tto_seqname\n")
        f.write("scaffold_1\t1\n")
        f.write("scaffold_2,2\n")  # Handle CSV line too
        f.write("  scaffold_3 \t 3  \n")
        f.write("\n")
        path = f.name

    try:
        mapping = load_seqname_map(path)
        assert mapping["scaffold_1"] == "1"
        assert mapping["scaffold_2"] == "2"
        assert mapping["scaffold_3"] == "3"
        assert "from_seqname" not in mapping
    finally:
        os.unlink(path)


def test_mapping_precedence():
    assembly_map = {"scaffold_1": "chr1", "scaffold_2": "chr2"}
    seq_map = {"scaffold_1": "1", "scaffold_3": "3"}

    combined = {}
    combined.update(assembly_map)
    combined.update(seq_map)

    assert combined["scaffold_1"] == "1"
    assert combined["scaffold_2"] == "chr2"
    assert combined["scaffold_3"] == "3"


def test_synthetic_classification_with_remapping():
    ref_genes = pd.DataFrame(
        [
            {
                "Chromosome": "Ca22chr1A_C_Albicans",
                "Start": 1000,
                "End": 2000,
                "Strand": "+",
                "ID": "gene1",
            }
        ]
    )
    ref_exons = pd.DataFrame(
        [
            {
                "Chromosome": "Ca22chr1A_C_Albicans",
                "Start": 1000,
                "End": 2000,
                "Strand": "+",
                "transcript_id": "tx1",
                "Parent": "gene1",
            }
        ]
    )
    ref_cds = pd.DataFrame()

    cons_genes = pd.DataFrame(
        [{"Chromosome": "1", "Start": 1000, "End": 2000, "Strand": "+", "ID": "cons_gene1"}]
    )
    cons_exons = pd.DataFrame(
        [
            {
                "Chromosome": "1",
                "Start": 1000,
                "End": 2000,
                "Strand": "+",
                "transcript_id": "cons_tx1",
                "Parent": "cons_gene1",
            }
        ]
    )
    cons_cds = pd.DataFrame()
    cons_mrna = pd.DataFrame(
        [
            {
                "Chromosome": "1",
                "Start": 1000,
                "End": 2000,
                "Strand": "+",
                "transcript_id": "cons_tx1",
                "Parent": "cons_gene1",
                "Evidence": "StringTie",
            }
        ]
    )

    mapping = {"Ca22chr1A_C_Albicans": "1"}
    mapped_ref_genes = remap_df_seqnames(ref_genes, mapping, label="Reference")
    mapped_ref_exons = remap_df_seqnames(ref_exons, mapping)

    ref_res, cons_res = classify_locus_pairs(
        mapped_ref_genes, mapped_ref_exons, ref_cds, cons_genes, cons_exons, cons_cds, cons_mrna
    )

    assert len(ref_res) == 1
    assert ref_res[0]["classification"] == "Exact_Match"
    assert ref_res[0]["matched_id"] == "cons_gene1"

    assert len(cons_res) == 1
    assert cons_res[0]["classification"] == "Exact_Match"


def test_isoform_aware_and_cds_classification():
    # Ref gene with 2 transcripts
    # tx1 is 1000-2000
    # tx2 is 1000-1500, 1800-2000 (len 700)
    ref_genes = pd.DataFrame(
        [{"Chromosome": "1", "Start": 1000, "End": 2000, "Strand": "+", "ID": "gene1"}]
    )
    ref_exons = pd.DataFrame(
        [
            {
                "Chromosome": "1",
                "Start": 1000,
                "End": 2000,
                "Strand": "+",
                "transcript_id": "tx1",
                "Parent": "gene1",
            },
            {
                "Chromosome": "1",
                "Start": 1000,
                "End": 1500,
                "Strand": "+",
                "transcript_id": "tx2",
                "Parent": "gene1",
            },
            {
                "Chromosome": "1",
                "Start": 1800,
                "End": 2000,
                "Strand": "+",
                "transcript_id": "tx2",
                "Parent": "gene1",
            },
        ]
    )
    # Ref CDS: tx1 1100-1900, tx2 1100-1500 and 1800-1900
    ref_cds = pd.DataFrame(
        [
            {
                "Chromosome": "1",
                "Start": 1100,
                "End": 1900,
                "Strand": "+",
                "transcript_id": "tx1",
                "Parent": "tx1",
            },
            {
                "Chromosome": "1",
                "Start": 1100,
                "End": 1500,
                "Strand": "+",
                "transcript_id": "tx2",
                "Parent": "tx2",
            },
            {
                "Chromosome": "1",
                "Start": 1800,
                "End": 1900,
                "Strand": "+",
                "transcript_id": "tx2",
                "Parent": "tx2",
            },
        ]
    )

    # Cons gene matching tx2 exactly in CDS, but different UTR boundaries (extra UTR intron)
    # cons_exons: 900-920, 950-1500, 1800-2050 -> intron chain: 920-950, 1500-1800 (differs from ref 1500-1800)
    # Overlap: 950-1500 (len 500) + 1800-2000 (len 200). Total 700.
    # Cons length: 20 + 550 + 250 = 820. recip = 700/820 = 0.853 > 0.8
    # CDS: 1100-1500, 1800-1900 (matches perfectly)
    cons_genes = pd.DataFrame(
        [{"Chromosome": "1", "Start": 900, "End": 2050, "Strand": "+", "ID": "cons_gene1"}]
    )
    cons_exons = pd.DataFrame(
        [
            {
                "Chromosome": "1",
                "Start": 900,
                "End": 920,
                "Strand": "+",
                "transcript_id": "cons_tx1",
                "Parent": "cons_gene1",
            },
            {
                "Chromosome": "1",
                "Start": 950,
                "End": 1500,
                "Strand": "+",
                "transcript_id": "cons_tx1",
                "Parent": "cons_gene1",
            },
            {
                "Chromosome": "1",
                "Start": 1800,
                "End": 2050,
                "Strand": "+",
                "transcript_id": "cons_tx1",
                "Parent": "cons_gene1",
            },
        ]
    )
    cons_cds = pd.DataFrame(
        [
            {
                "Chromosome": "1",
                "Start": 1100,
                "End": 1500,
                "Strand": "+",
                "transcript_id": "cons_tx1",
                "Parent": "cons_tx1",
            },
            {
                "Chromosome": "1",
                "Start": 1800,
                "End": 1900,
                "Strand": "+",
                "transcript_id": "cons_tx1",
                "Parent": "cons_tx1",
            },
        ]
    )
    cons_mrna = pd.DataFrame(
        [
            {
                "Chromosome": "1",
                "Start": 900,
                "End": 2050,
                "Strand": "+",
                "transcript_id": "cons_tx1",
                "Parent": "cons_gene1",
                "Evidence": "Helixer",
            }
        ]
    )

    ref_res, cons_res = classify_locus_pairs(
        ref_genes, ref_exons, ref_cds, cons_genes, cons_exons, cons_cds, cons_mrna
    )

    assert len(ref_res) == 1
    # Exons should be Structural_Mismatch because intron chains differ (due to extra UTR intron), but overlap > 0.8
    assert ref_res[0]["classification"] == "Structural_Mismatch"
    # CDS should be Exact_Match
    assert ref_res[0]["classification_cds"] == "Exact_Match"

    # Needs to have selected tx2, since it matches the intron chain (1500-1800)
    assert ref_res[0]["best_ref_transcript_id"] == "cons_tx1"  # The consensus ID mapped to it
    assert ref_res[0]["cds_intron_chain_match"] is True
