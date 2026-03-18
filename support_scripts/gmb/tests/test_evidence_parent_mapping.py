import pytest
import pandas as pd
import os
import tempfile
from collections import defaultdict
import subprocess

def test_evidence_parent_mapping():
    """Test that the evidence_attribution.tsv correctly maps transcript_id to the actual
    final Parent gene_id found in the consensus.gff3 file."""
    
    # Normally we'd rely on a full run, but let's mock the output directories
    # In integration test context, we'll run the actual gene_model_builder
    # on a small test.
    with tempfile.TemporaryDirectory() as tmpdir:
        # We can just write a dummy testing script here or assert based on 
        # actual integration outputs. Since we don't have test framework mock setup entirely,
        # let's assume we read from an actual output directory if it exists, or mock it if not.
        pass

    # A better test: Parse the real outputs from the integration run if available
    # Or simulate parsing logic directly.
    
    gff3_path = 'output/consensus.gff3'
    ev_path = 'output/evidence_attribution.tsv'
    
    if not os.path.exists(gff3_path) or not os.path.exists(ev_path):
        pytest.skip("Output files not found, skipping evidence parent mapping test.")
        
    # 1. Parse GFF3
    mrna_parent = {}
    with open(gff3_path, 'r') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] == 'mRNA':
                attrs = dict(pair.split('=') for pair in parts[8].split(';') if '=' in pair)
                tid = attrs.get('ID')
                parent = attrs.get('Parent')
                if tid and parent:
                    mrna_parent[tid] = parent
                    
    # 2. Parse TSV
    ev_df = pd.read_csv(ev_path, sep='\t')
    
    # 3. Assert matches 100%
    for _, row in ev_df.iterrows():
        tid = row['transcript_id']
        expected_gene_id = mrna_parent.get(tid)
        assert expected_gene_id is not None, f"Transcript {tid} not found in GFF3."
        assert row['gene_id'] == expected_gene_id, f"Mismatch for {tid}: {row['gene_id']} != {expected_gene_id}"
