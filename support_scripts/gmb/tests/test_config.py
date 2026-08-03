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
            "scoring:\n" "  max_isoforms_per_locus: 5\n" "  weights:\n" "    backbone: 3.0\n"
        )
        cfg = load_config(str(yaml_file))
        assert cfg.scoring.max_isoforms_per_locus == 5
        assert cfg.scoring.weights.backbone == 3.0
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

    def test_missing_file_raises_error(self):
        with pytest.raises(FileNotFoundError, match="Config file not found"):
            load_config("/nonexistent/path.yaml")


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


class TestProteinValidationConfig:
    def test_enabled_without_diamond_db_raises(self, tmp_path):
        yaml_file = tmp_path / "pv.yaml"
        yaml_file.write_text("protein_validation:\n  enabled: true\n  diamond_db: null\n")
        with pytest.raises(ValueError, match="diamond_db must be set"):
            load_config(str(yaml_file))

    def test_enabled_with_diamond_db_is_valid(self, tmp_path):
        yaml_file = tmp_path / "pv.yaml"
        yaml_file.write_text(
            "protein_validation:\n  enabled: true\n  diamond_db: /data/swissprot.dmnd\n"
        )
        cfg = load_config(str(yaml_file))
        assert cfg.protein_validation.enabled is True
        assert cfg.protein_validation.diamond_db == "/data/swissprot.dmnd"

    def test_disabled_with_no_diamond_db_is_valid(self):
        cfg = load_config()
        assert cfg.protein_validation.enabled is False
        assert cfg.protein_validation.diamond_db is None

    def test_enabled_with_zero_diamond_weight_no_db_is_valid(self, tmp_path):
        yaml_file = tmp_path / "pv.yaml"
        yaml_file.write_text(
            "protein_validation:\n  enabled: true\n  diamond_weight: 0\n" "  psauron_weight: 1.0\n"
        )
        cfg = load_config(str(yaml_file))
        assert cfg.protein_validation.enabled is True
        assert cfg.protein_validation.diamond_db is None


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

    def test_missing_file_in_list_raises_error(self, tmp_path):
        one = self._write(tmp_path, "one.yaml", "orf:\n  min_codons: 50\n")
        missing = str(tmp_path / "does_not_exist.yaml")
        with pytest.raises(FileNotFoundError, match="Config file not found"):
            load_config([one, missing])

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


class TestBackboneNamingAliases:
    """scoring.keep_helixer_without_support/helixer_filter/weights.helixer
    were renamed to their generic backbone_* equivalents (config.py's
    _DEPRECATED_KEY_ALIASES). Legacy keys must keep working with one
    DeprecationWarning each; the new key always wins if both are set.
    """

    def _write(self, tmp_path, name, text):
        p = tmp_path / name
        p.write_text(text)
        return str(p)

    def test_legacy_keep_helixer_without_support_still_works(self, tmp_path):
        legacy = self._write(
            tmp_path, "legacy.yaml", "scoring:\n  keep_helixer_without_support: false\n"
        )
        with pytest.warns(DeprecationWarning, match="keep_helixer_without_support"):
            cfg = load_config(legacy)
        assert cfg.scoring.keep_backbone_without_support is False

    def test_new_keep_backbone_without_support_no_warning(self, tmp_path, recwarn):
        new = self._write(
            tmp_path, "new.yaml", "scoring:\n  keep_backbone_without_support: false\n"
        )
        cfg = load_config(new)
        assert cfg.scoring.keep_backbone_without_support is False
        assert not any(issubclass(w.category, DeprecationWarning) for w in recwarn.list)

    def test_both_keep_helixer_and_backbone_new_wins_with_warning(self, tmp_path):
        both = self._write(
            tmp_path,
            "both.yaml",
            "scoring:\n"
            "  keep_helixer_without_support: false\n"
            "  keep_backbone_without_support: true\n",
        )
        with pytest.warns(DeprecationWarning, match="ignored because"):
            cfg = load_config(both)
        # New key wins; legacy value is not silently combined or averaged.
        assert cfg.scoring.keep_backbone_without_support is True

    def test_legacy_helixer_filter_section_still_works(self, tmp_path):
        legacy = self._write(
            tmp_path, "legacy.yaml", "helixer_filter:\n  min_cds_bp: 123\n  max_exons: 7\n"
        )
        with pytest.warns(DeprecationWarning, match="helixer_filter"):
            cfg = load_config(legacy)
        assert cfg.backbone_filter.min_cds_bp == 123
        assert cfg.backbone_filter.max_exons == 7

    def test_new_backbone_filter_section_no_warning(self, tmp_path, recwarn):
        new = self._write(tmp_path, "new.yaml", "backbone_filter:\n  min_cds_bp: 123\n")
        cfg = load_config(new)
        assert cfg.backbone_filter.min_cds_bp == 123
        assert not any(issubclass(w.category, DeprecationWarning) for w in recwarn.list)

    def test_legacy_weights_helixer_still_works(self, tmp_path):
        legacy = self._write(tmp_path, "legacy.yaml", "scoring:\n  weights:\n    helixer: 9.0\n")
        with pytest.warns(DeprecationWarning, match="weights.helixer"):
            cfg = load_config(legacy)
        assert cfg.scoring.weights.backbone == 9.0

    def test_legacy_and_new_keys_produce_identical_config(self, tmp_path):
        legacy = self._write(
            tmp_path,
            "legacy.yaml",
            "scoring:\n"
            "  keep_helixer_without_support: false\n"
            "  weights:\n"
            "    helixer: 9.0\n"
            "helixer_filter:\n"
            "  min_cds_bp: 123\n",
        )
        new = self._write(
            tmp_path,
            "new.yaml",
            "scoring:\n"
            "  keep_backbone_without_support: false\n"
            "  weights:\n"
            "    backbone: 9.0\n"
            "backbone_filter:\n"
            "  min_cds_bp: 123\n",
        )
        cfg_legacy = load_config(legacy)
        cfg_new = load_config(new)
        assert (
            cfg_legacy.scoring.keep_backbone_without_support
            == cfg_new.scoring.keep_backbone_without_support
        )
        assert cfg_legacy.scoring.weights.backbone == cfg_new.scoring.weights.backbone
        assert cfg_legacy.backbone_filter.min_cds_bp == cfg_new.backbone_filter.min_cds_bp

    def test_backbone_label_independent_of_weight_field_rename(self, tmp_path):
        # backbone_label (which evidence-source string the generic "backbone"
        # weight/gate actually applies to for this run) is untouched by this
        # rename -- it's set by the CLI from --helixer vs --tiberius, not by
        # scoring.weights.backbone/keep_backbone_without_support themselves.
        cfg_helixer = load_config()
        cfg_helixer.scoring.backbone_label = "Helixer"
        cfg_tiberius = load_config()
        cfg_tiberius.scoring.backbone_label = "Tiberius"
        assert cfg_helixer.scoring.weights.backbone == cfg_tiberius.scoring.weights.backbone
        assert (
            cfg_helixer.scoring.keep_backbone_without_support
            == cfg_tiberius.scoring.keep_backbone_without_support
        )

    def test_tracked_configs_use_current_names_no_deprecation_warning(self, recwarn):
        # The repo's own tracked configs should already be migrated -- if
        # this starts warning, one of them regressed to a legacy key.
        load_config()
        load_config(os.path.join(_CONFIGS_DIR, "apicomplexa_first_pass.yaml"))
        assert not any(issubclass(w.category, DeprecationWarning) for w in recwarn.list)


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


class TestInterProResolverConfig:
    """canonical_selection.interpro_resolver -- nesting, defaults, override shape."""

    def test_lives_under_canonical_selection_not_top_level(self):
        cfg = load_config()
        assert hasattr(cfg.canonical_selection, "interpro_resolver")
        assert not hasattr(cfg, "interpro_review")

    def test_master_switch_defaults_to_false(self):
        # InterPro involvement must be opt-in: a first-pass GMB build must
        # complete without it, and accidentally enabling it must never
        # happen just by loading defaults.
        cfg = load_config()
        resolver = cfg.canonical_selection.interpro_resolver
        assert resolver.enabled is False
        assert resolver.run_interproscan is False
        # apply_replacements defaults True, but is inert while enabled=False.
        assert resolver.apply_replacements is True

    def test_nextflow_exec_settings_have_no_hardcoded_site_paths(self):
        cfg = load_config()
        nf = cfg.canonical_selection.interpro_resolver.nextflow
        # Only the workflow identity/version and a generic default profile
        # are meaningful defaults; anything site-specific (data dir, work
        # dir, output dir, extra config file) must default to None/empty --
        # never a guessed laptop or cluster path.
        assert nf.data_dir is None
        assert nf.work_dir is None
        assert nf.output_dir is None
        assert nf.config_file is None
        assert nf.extra_args == []
        assert nf.workflow == "ebi-pf-team/interproscan6"

    def test_yaml_override_reaches_nested_resolver_config(self, tmp_path):
        override = tmp_path / "override.yaml"
        override.write_text(
            "canonical_selection:\n"
            "  interpro_resolver:\n"
            "    enabled: true\n"
            "    run_interproscan: true\n"
            "    apply_replacements: false\n"
            "    nextflow:\n"
            "      profile: slurm,singularity\n"
            "      data_dir: /shared/interpro/data\n"
        )
        cfg = load_config(str(override))
        resolver = cfg.canonical_selection.interpro_resolver
        assert resolver.enabled is True
        assert resolver.run_interproscan is True
        assert resolver.apply_replacements is False
        assert resolver.nextflow.profile == "slurm,singularity"
        assert resolver.nextflow.data_dir == "/shared/interpro/data"
        # Untouched nested defaults must survive the partial override.
        assert resolver.min_isoforms == 2

    def test_old_top_level_interpro_review_key_now_rejected(self, tmp_path):
        # The old top-level `interpro_review:` section moved under
        # canonical_selection -- a config still using the old location
        # should fail loudly (unknown key), not be silently ignored.
        override = tmp_path / "old_shape.yaml"
        override.write_text("interpro_review:\n  enabled: true\n")
        with pytest.raises(ValueError, match="interpro_review"):
            load_config(str(override))


class TestLayeredPresetArchitecture:
    """Tests for the standard → preset → user-override hierarchy.

    Parts 2, 3, 12, 13, 14, and 16 of the architecture refactor spec.
    """

    # ------------------------------------------------------------------
    # Part 13a: preset resolution
    # ------------------------------------------------------------------

    def test_list_build_presets_returns_nonempty_list(self):
        from gmb.pipeline.config import list_build_presets

        presets = list_build_presets()
        assert isinstance(presets, list)
        assert len(presets) >= 2, f"Expected at least fungi and apicomplexa, got {presets}"
        assert "fungi" in presets
        assert "apicomplexa" in presets

    def test_list_build_presets_excludes_standard(self):
        from gmb.pipeline.config import list_build_presets

        presets = list_build_presets()
        assert "standard" not in presets, "standard is the internal base, not a user preset"

    def test_unknown_preset_raises_file_not_found(self):
        with pytest.raises(FileNotFoundError, match="Unknown preset"):
            load_config(preset="does_not_exist")

    # ------------------------------------------------------------------
    # Part 13b: standard-only preset (preset=None)
    # ------------------------------------------------------------------

    def test_preset_none_loads_standard_only(self):
        # standard.yaml sets fungal_single_exon_mode: false (neutral).
        # fungi.yaml sets it true.  preset=None must not load fungi.
        cfg = load_config(preset=None)
        assert cfg.scoring.fungal_single_exon_mode is False

    def test_preset_string_none_also_loads_standard_only(self):
        cfg = load_config(preset="none")
        assert cfg.scoring.fungal_single_exon_mode is False

    def test_standard_has_permissive_transcript_length(self):
        cfg = load_config(preset=None)
        # standard sets 2,000,000 (effectively unlimited); fungi tightens to 20,000
        assert cfg.transcriptomic_filter.max_transcript_length >= 100_000

    def test_standard_has_permissive_intron_length(self):
        cfg = load_config(preset=None)
        assert cfg.transcriptomic_filter.max_intron_length >= 100_000

    def test_standard_has_no_transcript_splitting(self):
        cfg = load_config(preset=None)
        assert cfg.transcript_splitting.split_enabled is False

    def test_standard_has_more_isoforms_per_locus(self):
        cfg_std = load_config(preset=None)
        cfg_fungi = load_config(preset="fungi")
        # standard is permissive (5); fungi is compact (3)
        assert cfg_std.scoring.max_isoforms_per_locus > cfg_fungi.scoring.max_isoforms_per_locus

    # ------------------------------------------------------------------
    # Part 13c: fungi preset regression (standard + fungi == old fungi_default)
    # ------------------------------------------------------------------

    def test_fungi_preset_min_codons(self):
        cfg = load_config(preset="fungi")
        raw = _fungi_default_raw()
        assert cfg.orf.min_codons == raw["orf"]["min_codons"]

    def test_fungi_preset_max_transcript_length(self):
        cfg = load_config(preset="fungi")
        raw = _fungi_default_raw()
        assert (
            cfg.transcriptomic_filter.max_transcript_length
            == raw["transcriptomic_filter"]["max_transcript_length"]
        )

    def test_fungi_preset_max_intron_length(self):
        cfg = load_config(preset="fungi")
        raw = _fungi_default_raw()
        assert (
            cfg.transcriptomic_filter.max_intron_length
            == raw["transcriptomic_filter"]["max_intron_length"]
        )

    def test_fungi_preset_split_enabled(self):
        cfg = load_config(preset="fungi")
        raw = _fungi_default_raw()
        assert (
            cfg.transcript_splitting.split_enabled == raw["transcript_splitting"]["split_enabled"]
        )

    def test_fungi_preset_backbone_weight(self):
        cfg = load_config(preset="fungi")
        raw = _fungi_default_raw()
        assert cfg.scoring.weights.backbone == raw["scoring"]["weights"]["backbone"]

    def test_fungi_preset_max_isoforms(self):
        cfg = load_config(preset="fungi")
        raw = _fungi_default_raw()
        assert cfg.scoring.max_isoforms_per_locus == raw["scoring"]["max_isoforms_per_locus"]

    def test_fungi_preset_fungal_single_exon_mode(self):
        cfg = load_config(preset="fungi")
        raw = _fungi_default_raw()
        assert cfg.scoring.fungal_single_exon_mode == raw["scoring"]["fungal_single_exon_mode"]

    def test_fungi_preset_utr_caps(self):
        cfg = load_config(preset="fungi")
        raw = _fungi_default_raw()
        assert cfg.utr.max_5p_bp == raw["utr"]["max_5p_bp"]
        assert cfg.utr.max_3p_bp == raw["utr"]["max_3p_bp"]
        assert cfg.utr.max_total_bp == raw["utr"]["max_total_bp"]

    def test_fungi_preset_protein_validation_weights(self):
        cfg = load_config(preset="fungi")
        raw = _fungi_default_raw()
        assert cfg.protein_validation.diamond_weight == raw["protein_validation"]["diamond_weight"]
        assert cfg.protein_validation.psauron_weight == raw["protein_validation"]["psauron_weight"]
        assert cfg.protein_validation.min_score == raw["protein_validation"]["min_score"]
        assert cfg.protein_validation.policy == raw["protein_validation"]["policy"]

    def test_fungi_preset_end_tolerance(self):
        cfg = load_config(preset="fungi")
        raw = _fungi_default_raw()
        assert cfg.utr.end_tolerance_bp == raw["utr"]["end_tolerance_bp"]

    # ------------------------------------------------------------------
    # Part 13d: apicomplexa preset smoke tests
    # ------------------------------------------------------------------

    def test_apicomplexa_preset_loads(self):
        cfg = load_config(preset="apicomplexa")
        assert isinstance(cfg, PipelineConfig)

    def test_apicomplexa_has_tiberius_backbone_label(self):
        cfg = load_config(preset="apicomplexa")
        assert cfg.scoring.backbone_label == "Tiberius"

    def test_apicomplexa_has_larger_max_transcript_length_than_fungi(self):
        cfg_ap = load_config(preset="apicomplexa")
        cfg_fn = load_config(preset="fungi")
        assert cfg_ap.transcriptomic_filter.max_transcript_length > (
            cfg_fn.transcriptomic_filter.max_transcript_length
        )

    def test_apicomplexa_has_lower_backbone_weight_than_fungi(self):
        cfg_ap = load_config(preset="apicomplexa")
        cfg_fn = load_config(preset="fungi")
        assert cfg_ap.scoring.weights.backbone < cfg_fn.scoring.weights.backbone

    def test_apicomplexa_has_genblast_in_skip_orf(self):
        cfg = load_config(preset="apicomplexa")
        assert "GenBlast" in cfg.qc.skip_orf_inference_tracks

    def test_apicomplexa_has_higher_5p_utr_than_fungi(self):
        cfg_ap = load_config(preset="apicomplexa")
        cfg_fn = load_config(preset="fungi")
        assert cfg_ap.utr.max_5p_bp > cfg_fn.utr.max_5p_bp

    def test_apicomplexa_plus_user_override_layers_correctly(self, tmp_path):
        override = tmp_path / "override.yaml"
        override.write_text("orf:\n  min_codons: 25\n")
        cfg = load_config(str(override), preset="apicomplexa")
        assert cfg.orf.min_codons == 25
        assert cfg.scoring.backbone_label == "Tiberius"

    # ------------------------------------------------------------------
    # Part 12: deprecated preset alias
    # ------------------------------------------------------------------

    def test_fungi_default_alias_emits_deprecation_warning(self):
        with pytest.warns(DeprecationWarning, match="fungi_default"):
            cfg = load_config(preset="fungi_default")
        # Values must be identical to the current fungi preset.
        cfg_fungi = load_config(preset="fungi")
        assert cfg.orf.min_codons == cfg_fungi.orf.min_codons
        assert cfg.scoring.fungal_single_exon_mode == cfg_fungi.scoring.fungal_single_exon_mode

    # ------------------------------------------------------------------
    # Part 10: dump_config round-trip
    # ------------------------------------------------------------------

    def test_dump_config_returns_dict(self):
        from gmb.pipeline.config import dump_config

        cfg = load_config()
        d = dump_config(cfg)
        assert isinstance(d, dict)
        assert "orf" in d
        assert "scoring" in d

    def test_dump_config_values_match_cfg(self):
        from gmb.pipeline.config import dump_config

        cfg = load_config(preset="fungi")
        d = dump_config(cfg)
        assert d["orf"]["min_codons"] == cfg.orf.min_codons
        assert d["scoring"]["max_isoforms_per_locus"] == cfg.scoring.max_isoforms_per_locus

    # ------------------------------------------------------------------
    # Part 17: package layout — bundled presets load outside checkout
    # ------------------------------------------------------------------

    def test_fungi_preset_loads_via_importlib(self):
        from importlib.resources import files as pkg_files

        pkg = pkg_files("gmb.configs")
        text = pkg.joinpath("fungi.yaml").read_text(encoding="utf-8")
        assert "fungal_single_exon_mode" in text

    def test_standard_preset_loads_via_importlib(self):
        from importlib.resources import files as pkg_files

        pkg = pkg_files("gmb.configs")
        text = pkg.joinpath("standard.yaml").read_text(encoding="utf-8")
        assert "max_isoforms_per_locus" in text

    def test_apicomplexa_preset_loads_via_importlib(self):
        from importlib.resources import files as pkg_files

        pkg = pkg_files("gmb.configs")
        text = pkg.joinpath("apicomplexa.yaml").read_text(encoding="utf-8")
        assert "Tiberius" in text


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
