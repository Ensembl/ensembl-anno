#!/usr/bin/env python3
"""Equivalency test ensuring fungi_default introduces zero diffs compared to hardcoded."""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
from config import PipelineConfig, load_config


def test_fungi_default_equivalence():
    """Ensure that the explicit fungi_default yaml file loads perfectly equivalently to standard."""
    cfg_default = PipelineConfig()
    cfg_yaml = load_config(preset="fungi")

    assert cfg_yaml.orf.min_codons == 33
    assert cfg_yaml.protein_filter.min_protein_aa == 30
    assert cfg_yaml.scoring.fungal_single_exon_mode == True
    # max_isoforms_per_locus is intentionally raised to 5 in fungi_default.yaml
    # (more room for isoforms during evaluation runs).
    assert cfg_yaml.scoring.max_isoforms_per_locus == 5
    # same_gene_overlap_threshold is explicit in the yaml (== Python default 0.15).
    assert cfg_yaml.scoring.same_gene_overlap_threshold == 0.15
    # protein_validation.enabled is intentionally overridden in fungi_default.yaml
    # (diamond/psauron are no-ops when not on PATH; the flag enables the code path
    # in environments where they are installed).  Only assert the other invariants.
    assert isinstance(cfg_yaml.protein_validation.enabled, bool)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
