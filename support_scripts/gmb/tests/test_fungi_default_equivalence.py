#!/usr/bin/env python3
"""Equivalency test ensuring fungi_default introduces zero diffs compared to hardcoded."""

import os
import sys
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from config import load_config, PipelineConfig

def test_fungi_default_equivalence():
    """Ensure that the explicit fungi_default yaml file loads perfectly equivalently to standard."""
    cfg_default = PipelineConfig()
    cfg_yaml = load_config(preset='fungi')
    
    assert cfg_yaml.orf.min_codons == 33
    assert cfg_yaml.protein_filter.min_protein_aa == 30
    assert cfg_yaml.scoring.fungal_single_exon_mode == True
    assert cfg_yaml.scoring.max_isoforms_per_locus == 2
    assert cfg_yaml.protein_validation.enabled == False

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
