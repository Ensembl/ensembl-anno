import os
import tempfile
import pandas as pd
import pytest

from compare_annotations import (
    load_seqname_map,
    remap_df_seqnames,
    classify_locus_pairs
)

def test_load_seqname_map():
    with tempfile.NamedTemporaryFile('w', delete=False) as f:
        f.write("# comment\n")
        f.write("from_seqname\tto_seqname\n")
        f.write("scaffold_1\t1\n")
        f.write("scaffold_2,2\n") # Handle CSV line too
        f.write("  scaffold_3 \t 3  \n")
        f.write("\n")
        path = f.name
    
    try:
        mapping = load_seqname_map(path)
        assert mapping['scaffold_1'] == '1'
        assert mapping['scaffold_2'] == '2'
        assert mapping['scaffold_3'] == '3'
        assert 'from_seqname' not in mapping
    finally:
        os.unlink(path)

def test_mapping_precedence():
    assembly_map = {'scaffold_1': 'chr1', 'scaffold_2': 'chr2'}
    seq_map = {'scaffold_1': '1', 'scaffold_3': '3'}
    
    combined = {}
    combined.update(assembly_map)
    combined.update(seq_map)
    
    assert combined['scaffold_1'] == '1'
    assert combined['scaffold_2'] == 'chr2'
    assert combined['scaffold_3'] == '3'

def test_synthetic_classification_with_remapping():
    ref_genes = pd.DataFrame([
        {'Chromosome': 'Ca22chr1A_C_Albicans', 'Start': 1000, 'End': 2000, 'Strand': '+', 'ID': 'gene1'}
    ])
    ref_exons = pd.DataFrame([
        {'Chromosome': 'Ca22chr1A_C_Albicans', 'Start': 1000, 'End': 2000, 'Strand': '+', 'transcript_id': 'tx1', 'Parent': 'gene1'}
    ])
    
    cons_genes = pd.DataFrame([
        {'Chromosome': '1', 'Start': 1000, 'End': 2000, 'Strand': '+', 'ID': 'cons_gene1'}
    ])
    cons_exons = pd.DataFrame([
        {'Chromosome': '1', 'Start': 1000, 'End': 2000, 'Strand': '+', 'transcript_id': 'cons_tx1', 'Parent': 'cons_gene1'}
    ])
    cons_mrna = pd.DataFrame([
        {'Chromosome': '1', 'Start': 1000, 'End': 2000, 'Strand': '+', 'transcript_id': 'cons_tx1', 'Parent': 'cons_gene1', 'Evidence': 'StringTie'}
    ])
    
    mapping = {'Ca22chr1A_C_Albicans': '1'}
    mapped_ref_genes = remap_df_seqnames(ref_genes, mapping, label="Reference")
    mapped_ref_exons = remap_df_seqnames(ref_exons, mapping)
    
    ref_res, cons_res = classify_locus_pairs(
        mapped_ref_genes, mapped_ref_exons,
        cons_genes, cons_exons, cons_mrna
    )
    
    assert len(ref_res) == 1
    assert ref_res[0]['classification'] == 'Exact_Match'
    assert ref_res[0]['matched_id'] == 'cons_gene1'
    
    assert len(cons_res) == 1
    assert cons_res[0]['classification'] == 'Matched'
