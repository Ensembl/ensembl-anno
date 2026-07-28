#!/usr/bin/env python3
"""Tests for config module."""

import glob
import os
import sys

import pytest
import yaml

sys.path.insert(0, os.path.dirname(__file__))
from gmb.pipeline.config import PipelineConfig, load_config

_CONFIGS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "configs")
_FUNGI_DEFAULT_YAML = os.path.join(_CONFIGS_DIR, "fungi_default.yaml")


def _fungi_default_raw() -> dict:
    with open(_FUNGI_DEFAULT_YAML) as fh:
        return yaml.safe_load(fh)


class TestDefaultConfig:
    def test_loads_without_yaml(self):
        """Default config should load even with no YAML file."""
        cfg = PipelineConfig()
        assert cfg.preset == "fungi"
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
        # max_isoforms_per_locus's tuned value has drifted once already (5 -> 3,
        # 2026-05-12) without this assertion being updated -- see
        # fungi_default.yaml's comment on the key for the history. Read the
        # tracked value directly rather than re-hard-coding a number here, so
        # this test validates that config loading is faithful to the YAML
        # (and fails if the key is ever silently removed), not a specific
        # untracked business decision.
        raw = _fungi_default_raw()
        assert "max_isoforms_per_locus" in raw["scoring"]
        assert cfg.scoring.max_isoforms_per_locus == raw["scoring"]["max_isoforms_per_locus"]
        assert cfg.scoring.fungal_single_exon_mode is True
        assert cfg.transcriptomic_filter.allow_single_exon is True


class TestYamlOverride:
    def test_override_min_codons(self, tmp_path):
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text("orf:\n  min_codons: 50\n")
        cfg = load_config(str(yaml_file))
        assert cfg.orf.min_codons == 50
        # Other values should remain default
        assert cfg.protein_filter.min_protein_aa == 30

    def test_override_nested(self, tmp_path):
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text(
            "scoring:\n" "  max_isoforms_per_locus: 5\n" "  weights:\n" "    helixer: 3.0\n"
        )
        cfg = load_config(str(yaml_file))
        assert cfg.scoring.max_isoforms_per_locus == 5
        assert cfg.scoring.weights.helixer == 3.0
        assert cfg.scoring.weights.scallop == 1.0  # unchanged

    def test_override_list(self, tmp_path):
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text("qc:\n" "  skip_orf_inference_tracks:\n" "    - UniProt\n")
        cfg = load_config(str(yaml_file))
        assert list(cfg.qc.skip_orf_inference_tracks) == ["UniProt"]

    def test_unknown_key_raises(self, tmp_path):
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text("nonexistent_key: true\n")
        with pytest.raises(ValueError, match="Unknown configuration key"):
            cfg = load_config(str(yaml_file))

    def test_missing_file_returns_defaults(self):
        cfg = load_config("/nonexistent/path.yaml")
        assert cfg.orf.min_codons == 33


class TestUtrSupportConfig:
    def test_utr_support_defaults(self):
        cfg = PipelineConfig()
        assert cfg.utr.require_end_support is True
        assert cfg.utr.end_support_mode == "multisource_end_agreement"
        assert cfg.utr.end_support_sources == ["Scallop", "StringTie"]
        assert cfg.utr.end_tolerance_bp == 50
        assert cfg.utr.require_multisource_for_utr_5p is True
        assert cfg.utr.require_multisource_for_utr_3p is True
        assert cfg.utr.fallback_policy_when_unsupported == "drop_utr"
        assert cfg.utr.min_protein_coding_score_for_utr is None
        assert cfg.utr.max_end_extension_bp is None

    def test_utr_support_yaml_override(self, tmp_path):
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text(
            "utr:\n"
            "  require_end_support: false\n"
            '  end_support_mode: "protein_validated"\n'
            "  end_support_sources:\n"
            "    - Helixer\n"
            '  fallback_policy_when_unsupported: "drop_transcript"\n'
        )
        cfg = load_config(str(yaml_file))
        assert cfg.utr.require_end_support is False
        assert cfg.utr.end_support_mode == "protein_validated"
        assert cfg.utr.end_support_sources == ["Helixer"]
        assert cfg.utr.fallback_policy_when_unsupported == "drop_transcript"

    def test_invalid_end_support_mode_raises(self, tmp_path):
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text('utr:\n  end_support_mode: "invalid_mode"\n')
        with pytest.raises(ValueError, match="Invalid end_support_mode: invalid_mode"):
            load_config(str(yaml_file))

    def test_invalid_fallback_policy_raises(self, tmp_path):
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text('utr:\n  fallback_policy_when_unsupported: "magic"\n')
        with pytest.raises(ValueError, match="Invalid fallback_policy_when_unsupported: magic"):
            load_config(str(yaml_file))


class TestPolicyAlias:
    def test_penalize_accepted(self, tmp_path):
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text('protein_validation:\n  policy: "penalize"\n')
        cfg = load_config(str(yaml_file))
        assert cfg.protein_validation.policy == "penalize"

    def test_penalise_accepted(self, tmp_path):
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text('protein_validation:\n  policy: "penalise"\n')
        cfg = load_config(str(yaml_file))
        assert cfg.protein_validation.policy == "penalise"

    def test_default_preset_uses_penalize(self):
        cfg = load_config()
        assert cfg.protein_validation.policy == "penalize"


class TestLayeredConfig:
    """gmb-build (and other CLIs) accept a repeated --config; load_config()
    accepts path as a str (backward compatible) or a list of str, applied
    in order on top of the preset default.
    """

    def _write(self, tmp_path, name, text):
        p = tmp_path / name
        p.write_text(text)
        return str(p)

    def test_single_config_string_backward_compatible(self, tmp_path):
        one = self._write(tmp_path, "one.yaml", "orf:\n  min_codons: 50\n")
        cfg_str = load_config(one)
        cfg_list = load_config([one])
        assert cfg_str.orf.min_codons == cfg_list.orf.min_codons == 50

    def test_two_configs_later_wins_on_shared_key(self, tmp_path):
        one = self._write(tmp_path, "one.yaml", "orf:\n  min_codons: 50\n")
        two = self._write(tmp_path, "two.yaml", "orf:\n  min_codons: 60\n")
        cfg = load_config([one, two])
        assert cfg.orf.min_codons == 60
        # Order matters: reversed application gives the other result.
        cfg_rev = load_config([two, one])
        assert cfg_rev.orf.min_codons == 50

    def test_three_configs_each_contributes(self, tmp_path):
        one = self._write(tmp_path, "one.yaml", "orf:\n  min_codons: 50\n")
        two = self._write(tmp_path, "two.yaml", "protein_filter:\n  min_protein_aa: 40\n")
        three = self._write(tmp_path, "three.yaml", "scoring:\n  max_isoforms_per_locus: 7\n")
        cfg = load_config([one, two, three])
        assert cfg.orf.min_codons == 50
        assert cfg.protein_filter.min_protein_aa == 40
        assert cfg.scoring.max_isoforms_per_locus == 7

    def test_nested_override_across_files_preserves_sibling_keys(self, tmp_path):
        # First file enables protein_validation and sets one field; second
        # file (a small overlay, matching the real chr1 use case) sets a
        # different field in the *same* nested section -- both must survive.
        base = self._write(
            tmp_path, "base.yaml", "protein_validation:\n  enabled: true\n  min_score: 0.9\n"
        )
        overlay = self._write(
            tmp_path, "overlay.yaml", "protein_validation:\n  diamond_db: /tmp/x.dmnd\n"
        )
        cfg = load_config([base, overlay])
        assert cfg.protein_validation.enabled is True
        assert cfg.protein_validation.min_score == 0.9
        assert cfg.protein_validation.diamond_db == "/tmp/x.dmnd"

    def test_list_valued_key_is_replaced_not_concatenated_across_files(self, tmp_path):
        one = self._write(
            tmp_path, "one.yaml", "qc:\n  skip_orf_inference_tracks:\n    - UniProt\n"
        )
        two = self._write(
            tmp_path, "two.yaml", "qc:\n  skip_orf_inference_tracks:\n    - GenBlast\n"
        )
        cfg = load_config([one, two])
        # Complete replacement: GenBlast only, not UniProt+GenBlast.
        assert list(cfg.qc.skip_orf_inference_tracks) == ["GenBlast"]

    def test_unknown_key_in_second_file_raises(self, tmp_path):
        one = self._write(tmp_path, "one.yaml", "orf:\n  min_codons: 50\n")
        two = self._write(tmp_path, "two.yaml", "not_a_real_key: 1\n")
        with pytest.raises(ValueError, match="Unknown configuration key"):
            load_config([one, two])

    def test_missing_file_in_list_is_skipped_others_still_applied(self, tmp_path):
        one = self._write(tmp_path, "one.yaml", "orf:\n  min_codons: 50\n")
        missing = str(tmp_path / "does_not_exist.yaml")
        cfg = load_config([one, missing])
        assert cfg.orf.min_codons == 50  # first file's override still applied

    def test_empty_list_and_none_both_use_defaults(self, tmp_path):
        assert load_config([]).orf.min_codons == load_config(None).orf.min_codons == 33

    def test_apicomplexa_layered_example_matches_documented_use_case(self):
        # The real motivating case: a small protein-validation overlay
        # layered on top of the established apicomplexa tuning, instead of
        # duplicating that tuning's deltas inside the overlay file.
        first_pass = os.path.join(_CONFIGS_DIR, "apicomplexa_first_pass.yaml")
        overlay = os.path.join(_CONFIGS_DIR, "apicomplexa_chr1_protein_validation.example.yaml")
        cfg = load_config([first_pass, overlay])
        # From apicomplexa_first_pass.yaml:
        assert cfg.transcriptomic_filter.max_transcript_length == 35000
        # From the overlay, layered on top:
        assert cfg.protein_validation.enabled is True


class TestAllTrackedConfigsLoad:
    """Every tracked configs/*.yaml must load cleanly as a fungi-preset delta.

    Catches the class of bug this is meant to prevent: a config key that was
    renamed/removed in the dataclasses but left stale in a tracked YAML file
    (load_config() raises ValueError on any unknown key), or a value that
    fails a dataclass __post_init__ validator.
    """

    @pytest.mark.parametrize(
        "yaml_path",
        sorted(glob.glob(os.path.join(_CONFIGS_DIR, "*.yaml"))),
        ids=lambda p: os.path.basename(p),
    )
    def test_config_loads(self, yaml_path):
        if os.path.basename(yaml_path) == "fungi_default.yaml":
            # This is the base preset itself, loaded automatically by every
            # other case here -- exercised directly via load_config().
            cfg = load_config()
        else:
            cfg = load_config(yaml_path)
        assert isinstance(cfg, PipelineConfig)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
