#!/usr/bin/env python3
"""Tests for config module."""

import os
import sys
import tempfile

import pytest

sys.path.insert(0, os.path.dirname(__file__))
from config import load_config, PipelineConfig


class TestDefaultConfig:
    def test_loads_without_yaml(self):
        """Default config should load even with no YAML file."""
        cfg = PipelineConfig()
        assert cfg.preset == 'fungi'
        assert cfg.orf.min_codons == 33
        assert cfg.protein_filter.min_protein_aa == 30

    def test_load_config_function(self):
        cfg = load_config()
        assert isinstance(cfg, PipelineConfig)
        assert cfg.orf.min_codons == 33

    def test_fungal_preset_values(self):
        cfg = load_config()
        assert cfg.orf.min_codons == 33
        assert cfg.orf.allow_partial_5 is True
        assert cfg.protein_filter.top_n_per_locus == 3
        assert cfg.scoring.max_isoforms_per_locus == 2
        assert cfg.scoring.fungal_single_exon_mode is True
        assert cfg.transcriptomic_filter.allow_single_exon is True


class TestYamlOverride:
    def test_override_min_codons(self, tmp_path):
        yaml_file = tmp_path / 'test.yaml'
        yaml_file.write_text('orf:\n  min_codons: 50\n')
        cfg = load_config(str(yaml_file))
        assert cfg.orf.min_codons == 50
        # Other values should remain default
        assert cfg.protein_filter.min_protein_aa == 30

    def test_override_nested(self, tmp_path):
        yaml_file = tmp_path / 'test.yaml'
        yaml_file.write_text(
            'scoring:\n'
            '  max_isoforms_per_locus: 5\n'
            '  weights:\n'
            '    helixer: 3.0\n'
        )
        cfg = load_config(str(yaml_file))
        assert cfg.scoring.max_isoforms_per_locus == 5
        assert cfg.scoring.weights.helixer == 3.0
        assert cfg.scoring.weights.scallop == 1.0  # unchanged

    def test_missing_file_returns_defaults(self):
        cfg = load_config('/nonexistent/path.yaml')
        assert cfg.orf.min_codons == 33


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
