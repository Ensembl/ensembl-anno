#!/usr/bin/env python3
"""Production contract tests.

These pin the behaviour a production pipeline is entitled to rely on:

* evidence WEIGHTS resolve by role, never by tool name -- so an arbitrary new
  assembler is scored identically to a known one;
* legacy tool-name config keys still load, with a deprecation warning;
* ``backbone_intron_rescue`` is applicability-gated and defaults to safe;
* long-read evidence is genuinely optional, and a bad long-read track can be
  demoted or rejected;
* preflight classifies splice quality per evidence role;
* every shipped preset resolves to the documented values;
* the new-clade template is loadable and enables nothing.

Tool names used for the generic cases are deliberately invented
(``AssemblerX``, ``LongTool``, ``PredictorQ``) so that a regression to
literal-name matching fails here rather than in production.
"""

from __future__ import annotations

import json
import os
import warnings

import pandas as pd
import pytest

from gmb.pipeline.applicability import (
    ASSEMBLED_CANONICAL_SPLICE_MIN,
    BACKBONE_MULTI_EXON_MAX,
    measure_backbone_resolution,
    normalise_rescue_mode,
    resolve_backbone_intron_rescue,
)
from gmb.pipeline.canonical_evidence import (
    EVIDENCE_CLASS_BACKBONE,
    EVIDENCE_CLASS_LONG_READ,
    EVIDENCE_CLASS_SHORT_READ,
    EvidenceRoles,
)
from gmb.pipeline.config import load_config
from gmb.pipeline.scoring import rescue_enabled, weights_for_role

_CONFIGS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "configs")


def _exon_frame(spec):
    """Build an exon frame: spec = {source: (n_multi_exon, n_single_exon)}."""
    rows = []
    for source, (n_multi, n_single) in spec.items():
        for i in range(n_multi):
            rows += [(source, f"{source}_m{i}"), (source, f"{source}_m{i}")]
        for i in range(n_single):
            rows.append((source, f"{source}_s{i}"))
    return pd.DataFrame(rows, columns=["Source", "transcript_id"])


# ---------------------------------------------------------------------------
# role-based weight resolution
# ---------------------------------------------------------------------------

class TestRoleBasedWeights:
    def test_unknown_tool_name_in_a_known_role_gets_that_roles_weight(self):
        """The whole point of the refactor.

        An assembler nobody has heard of, listed under shortread_labels, must be
        weighted exactly like Scallop -- not dropped to the unknown weight.
        """
        cfg = load_config(preset=None)
        cfg.scoring.shortread_labels = ["AssemblerX", "Scallop"]
        cfg.scoring.longread_label = "LongTool"
        cfg.scoring.backbone_label = "PredictorQ"
        cfg.scoring.weights.backbone = 2.6
        cfg.scoring.weights.short_read = 1.4
        cfg.scoring.weights.long_read = 1.3
        cfg.scoring.weights.unknown = 0.1
        roles = EvidenceRoles.from_config(cfg.scoring)

        assert weights_for_role(cfg.scoring.weights, roles.role_of("AssemblerX")) == 1.4
        assert weights_for_role(cfg.scoring.weights, roles.role_of("Scallop")) == 1.4
        assert weights_for_role(cfg.scoring.weights, roles.role_of("LongTool")) == 1.3
        assert weights_for_role(cfg.scoring.weights, roles.role_of("PredictorQ")) == 2.6

    def test_renaming_every_tool_does_not_change_weights(self):
        """Rename all tools to invented names; weights must be identical."""
        known = load_config(preset=None)
        known.scoring.weights.short_read = 1.7
        known_roles = EvidenceRoles.from_config(known.scoring)
        known_w = [weights_for_role(known.scoring.weights, known_roles.role_of(s))
                   for s in ("Helixer", "Scallop", "StringTie", "Minimap2")]

        renamed = load_config(preset=None)
        renamed.scoring.weights.short_read = 1.7
        renamed.scoring.backbone_label = "PredictorQ"
        renamed.scoring.shortread_labels = ["AssemblerX", "AssemblerY"]
        renamed.scoring.longread_label = "LongTool"
        renamed_roles = EvidenceRoles.from_config(renamed.scoring)
        renamed_w = [weights_for_role(renamed.scoring.weights, renamed_roles.role_of(s))
                     for s in ("PredictorQ", "AssemblerX", "AssemblerY", "LongTool")]

        assert known_w == renamed_w

    def test_source_with_no_role_gets_the_unknown_weight_not_a_literal(self):
        cfg = load_config(preset=None)
        cfg.scoring.weights.unknown = 0.25
        roles = EvidenceRoles.from_config(cfg.scoring)
        assert weights_for_role(cfg.scoring.weights, roles.role_of("SomethingElse")) == 0.25

    def test_legacy_tool_name_weight_keys_still_load(self, tmp_path):
        yaml_file = tmp_path / "legacy.yaml"
        yaml_file.write_text(
            "scoring:\n  weights:\n    helixer: 3.3\n    scallop: 1.4\n"
            "    stringtie: 1.4\n    minimap2: 0.9\n")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            cfg = load_config(str(yaml_file), preset=None)
        assert cfg.scoring.weights.backbone == 3.3
        assert cfg.scoring.weights.short_read == 1.4
        assert cfg.scoring.weights.long_read == 0.9
        assert any(issubclass(w.category, DeprecationWarning) for w in caught)

    def test_conflicting_legacy_shortread_weights_take_the_maximum_and_warn(self, tmp_path):
        """scallop and stringtie are now one role; disagreement must not silently
        down-weight an evidence class."""
        yaml_file = tmp_path / "legacy.yaml"
        yaml_file.write_text("scoring:\n  weights:\n    scallop: 1.0\n    stringtie: 2.5\n")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            cfg = load_config(str(yaml_file), preset=None)
        assert cfg.scoring.weights.short_read == 2.5
        assert any("DIFFERENT" in str(w.message) for w in caught)

    def test_explicit_role_key_beats_legacy_key(self, tmp_path):
        yaml_file = tmp_path / "both.yaml"
        yaml_file.write_text(
            "scoring:\n  weights:\n    short_read: 3.0\n    scallop: 1.0\n")
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            cfg = load_config(str(yaml_file), preset=None)
        assert cfg.scoring.weights.short_read == 3.0


# ---------------------------------------------------------------------------
# backbone intron rescue: off / on / auto
# ---------------------------------------------------------------------------

class TestRescueMode:
    @pytest.mark.parametrize("value,expected", [
        ("off", "off"), ("on", "on"), ("auto", "auto"),
        (False, "off"), (True, "on"), (None, "off"),
        ("OFF", "off"), ("Auto", "auto"), ("true", "on"), ("false", "off"),
    ])
    def test_mode_normalisation(self, value, expected):
        assert normalise_rescue_mode(value) == expected

    def test_typo_is_rejected(self):
        with pytest.raises(ValueError):
            normalise_rescue_mode("atuo")

    def test_off_never_fires(self):
        cfg = load_config(preset=None)
        cfg.scoring.backbone_intron_rescue = "off"
        d = resolve_backbone_intron_rescue(cfg.scoring, _exon_frame(
            {"Helixer": (50, 950), "Scallop": (950, 50)}), 1.0)
        assert d.enabled is False

    def test_on_fires_regardless_of_evidence(self):
        """Expert override: fires even where auto would refuse."""
        cfg = load_config(preset=None)
        cfg.scoring.backbone_intron_rescue = "on"
        d = resolve_backbone_intron_rescue(cfg.scoring, _exon_frame(
            {"Helixer": (950, 50), "Scallop": (950, 50)}), 1.0)
        assert d.enabled is True
        assert "explicitly enabled" in d.reason

    def test_auto_fires_on_a_weak_backbone(self):
        """The measured P. falciparum state: 17% multi-exon backbone, 99% assembled."""
        cfg = load_config(preset=None)
        cfg.scoring.backbone_intron_rescue = "auto"
        d = resolve_backbone_intron_rescue(cfg.scoring, _exon_frame(
            {"Helixer": (173, 827), "Scallop": (987, 13)}), 1.0)
        assert d.enabled is True

    def test_auto_refuses_on_a_strong_backbone(self):
        """The measured Z. tritici state: 72% multi-exon backbone.

        This is the regression that matters most: rescued models were 1.8%
        CDS-exact here against 31.9% for backbone-only.
        """
        cfg = load_config(preset=None)
        cfg.scoring.backbone_intron_rescue = "auto"
        d = resolve_backbone_intron_rescue(cfg.scoring, _exon_frame(
            {"Helixer": (719, 281), "Scallop": (921, 79)}), 1.0)
        assert d.enabled is False
        assert "not under-calling introns" in d.reason

    def test_auto_refuses_when_assembled_splice_quality_is_poor(self):
        """A weak backbone is not enough -- the replacements must be credible."""
        cfg = load_config(preset=None)
        cfg.scoring.backbone_intron_rescue = "auto"
        d = resolve_backbone_intron_rescue(cfg.scoring, _exon_frame(
            {"Helixer": (173, 827), "Scallop": (987, 13)}),
            ASSEMBLED_CANONICAL_SPLICE_MIN - 0.2)
        assert d.enabled is False
        assert "canonically spliced" in d.reason

    def test_auto_refuses_when_splice_quality_was_not_measured(self):
        cfg = load_config(preset=None)
        cfg.scoring.backbone_intron_rescue = "auto"
        d = resolve_backbone_intron_rescue(cfg.scoring, _exon_frame(
            {"Helixer": (173, 827), "Scallop": (987, 13)}), None)
        assert d.enabled is False

    def test_auto_refuses_on_too_little_evidence(self):
        cfg = load_config(preset=None)
        cfg.scoring.backbone_intron_rescue = "auto"
        d = resolve_backbone_intron_rescue(cfg.scoring, _exon_frame(
            {"Helixer": (5, 20), "Scallop": (20, 5)}), 1.0)
        assert d.enabled is False
        assert "too few models" in d.reason

    def test_auto_refuses_with_no_assembled_evidence(self):
        cfg = load_config(preset=None)
        cfg.scoring.backbone_intron_rescue = "auto"
        d = resolve_backbone_intron_rescue(cfg.scoring, _exon_frame(
            {"Helixer": (173, 827)}), 1.0)
        assert d.enabled is False

    def test_gate_is_role_based_not_tool_name_based(self):
        """Same evidence state, invented tool names -> same decision."""
        cfg = load_config(preset=None)
        cfg.scoring.backbone_intron_rescue = "auto"
        cfg.scoring.backbone_label = "PredictorQ"
        cfg.scoring.shortread_labels = ["AssemblerX"]
        d = resolve_backbone_intron_rescue(cfg.scoring, _exon_frame(
            {"PredictorQ": (173, 827), "AssemblerX": (987, 13)}), 1.0)
        assert d.enabled is True

    def test_calibration_thresholds_separate_the_observed_genomes(self):
        """The gate's thresholds must sit between the measured helpful and
        harmful cases, with margin on both sides.

        Measured backbone multi-exon fractions, and whether rescue helped:
            P. falciparum  Tiberius  17.3%   helped (68.6% CDS-exact)
            T. gondii      Tiberius  38.7%   helped (validated)
            Z. tritici     Helixer   71.9%   HARMED  (1.8% CDS-exact)

        A change to BACKBONE_MULTI_EXON_MAX that stops separating these has to
        be a deliberate edit to this test.
        """
        helped = (0.173, 0.387)
        harmed = (0.719,)
        assert max(helped) < BACKBONE_MULTI_EXON_MAX, \
            "threshold must admit every case where rescue helped"
        assert min(harmed) > BACKBONE_MULTI_EXON_MAX, \
            "threshold must exclude every case where rescue harmed"
        # and not by a hair: at least 10 percentage points of margin each side
        assert BACKBONE_MULTI_EXON_MAX - max(helped) > 0.10
        assert min(harmed) - BACKBONE_MULTI_EXON_MAX > 0.10

    def test_assembled_splice_floor_rejects_the_known_bad_track(self):
        """16.4% canonical — the P. falciparum long-read consensus."""
        assert 0.164 < ASSEMBLED_CANONICAL_SPLICE_MIN

    def test_measure_reports_the_numbers_it_decided_on(self):
        cfg = load_config(preset=None)
        stats = measure_backbone_resolution(
            _exon_frame({"Helixer": (200, 800), "Scallop": (900, 100)}), cfg.scoring, 0.99)
        assert stats.backbone_models == 1000
        assert stats.backbone_multi_exon == 200
        assert abs(stats.backbone_multi_exon_fraction - 0.20) < 1e-9
        assert abs(stats.assembled_multi_exon_fraction - 0.90) < 1e-9
        assert abs(stats.assembled_to_backbone_ratio - 4.5) < 1e-9

    def test_unresolved_auto_is_treated_as_off(self):
        """select_isoforms() called without a builder has no run-level evidence."""
        cfg = load_config(preset=None)
        cfg.scoring.backbone_intron_rescue = "auto"
        cfg.scoring.backbone_intron_rescue_resolved = None
        assert rescue_enabled(cfg.scoring) is False

    def test_resolved_decision_wins_over_the_mode(self):
        cfg = load_config(preset=None)
        cfg.scoring.backbone_intron_rescue = "off"
        cfg.scoring.backbone_intron_rescue_resolved = True
        assert rescue_enabled(cfg.scoring) is True


# ---------------------------------------------------------------------------
# long-read handling
# ---------------------------------------------------------------------------

class TestLongReadOptional:
    def test_long_read_absent_is_the_default_and_needs_no_special_case(self):
        for preset in (None, "fungi", "apicomplexa"):
            cfg = load_config(preset=preset)
            assert cfg.scoring.longread_disposition == "primary_structural"

    def test_dispositions_are_the_documented_set(self):
        cfg = load_config(preset=None)
        for value in ("primary_structural", "support_only", "reject"):
            cfg.scoring.longread_disposition = value
            assert cfg.scoring.longread_disposition == value

    def test_guard_cannot_fire_without_a_long_read_role(self):
        """With no long-read source, no structure can be long-read-only."""
        cfg = load_config(preset=None)
        cfg.scoring.longread_structural_guard = True
        roles = EvidenceRoles.from_config(cfg.scoring)
        for source in ("Helixer", "Scallop", "StringTie"):
            assert roles.role_of(source) != EVIDENCE_CLASS_LONG_READ


# ---------------------------------------------------------------------------
# preflight
# ---------------------------------------------------------------------------

class TestPreflightSpliceClassification:
    def _thresholds(self, role):
        from gmb.preflight.checks import _splice_thresholds
        return _splice_thresholds(load_config(preset=None).preflight, role)

    def test_long_read_is_judged_most_strictly(self):
        """A long-read track claims to observe splice structure directly."""
        _, lr_fail = self._thresholds(EVIDENCE_CLASS_LONG_READ)
        _, bb_fail = self._thresholds(EVIDENCE_CLASS_BACKBONE)
        assert lr_fail > bb_fail

    def test_the_known_bad_long_read_fraction_fails(self):
        """The P. falciparum Minimap2 consensus measured 16.4% canonical."""
        _, fail_below = self._thresholds(EVIDENCE_CLASS_LONG_READ)
        assert 0.164 < fail_below

    def test_a_good_assembled_track_passes(self):
        warn_below, _ = self._thresholds(EVIDENCE_CLASS_SHORT_READ)
        assert 1.0 >= warn_below

    def test_protein_alignment_role_is_not_splice_gated(self):
        from gmb.pipeline.canonical_evidence import EVIDENCE_CLASS_PROTEIN_ALIGNMENT
        warn_below, fail_below = self._thresholds(EVIDENCE_CLASS_PROTEIN_ALIGNMENT)
        assert (warn_below, fail_below) == (0.0, 0.0)


# ---------------------------------------------------------------------------
# presets
# ---------------------------------------------------------------------------

class TestShippedPresets:
    def test_standard_is_the_neutral_base(self):
        cfg = load_config(preset="standard")
        s = cfg.scoring
        assert s.structural_corroboration is False
        assert s.longread_structural_guard is False
        assert s.protein_support_mode == "positional"
        assert normalise_rescue_mode(s.backbone_intron_rescue) == "off"
        assert s.weights.backbone == 2.0
        assert s.weights.short_read == 1.0
        assert s.weights.long_read == 1.0

    def test_apicomplexa_resolves_to_the_validated_values(self):
        s = load_config(preset="apicomplexa").scoring
        assert s.backbone_label == "Tiberius"
        assert s.weights.backbone == 2.6
        assert s.weights.long_read == 1.3
        assert s.multi_source_bonus == 1.2
        assert s.max_isoforms_per_locus == 3
        assert s.structural_corroboration is True
        assert s.longread_structural_guard is True
        assert s.protein_support_mode == "cds_span_compatible"
        assert normalise_rescue_mode(s.backbone_intron_rescue) == "auto"

    def test_fungi_resolves_to_the_validated_baseline(self):
        s = load_config(preset="fungi").scoring
        assert s.backbone_label == "Helixer"
        assert s.weights.backbone == 3.1
        assert s.multi_source_bonus == 0.5
        assert s.structural_corroboration is False
        assert s.longread_structural_guard is False
        assert s.protein_support_mode == "positional"
        assert normalise_rescue_mode(s.backbone_intron_rescue) == "off"

    def test_no_preset_silently_loses_its_clade_settings(self):
        """Regression guard for a duplicate top-level YAML key.

        A second `scoring:` block in a preset would make YAML discard the first
        outright, so the preset would load but resolve to neutral defaults.
        """
        assert load_config(preset="apicomplexa").scoring.weights.backbone != \
            load_config(preset="standard").scoring.weights.backbone
        assert load_config(preset="fungi").scoring.weights.backbone != \
            load_config(preset="standard").scoring.weights.backbone

    def test_duplicate_top_level_key_is_rejected(self, tmp_path):
        bad = tmp_path / "dup.yaml"
        bad.write_text("scoring:\n  multi_source_bonus: 9.0\n"
                       "scoring:\n  max_isoforms_per_locus: 1\n")
        with pytest.raises(ValueError, match="Duplicate key"):
            load_config(str(bad), preset=None)


class TestNewCladeTemplate:
    """The blank template must be usable without reading source code."""

    def _path(self):
        return os.path.join(_CONFIGS_DIR, "new_clade_template.yaml")

    def test_template_exists(self):
        assert os.path.exists(self._path()), "new_clade_template.yaml is missing"

    def test_template_loads_on_top_of_the_neutral_preset(self):
        load_config(self._path(), preset="standard")

    def test_template_enables_no_biological_policy(self):
        s = load_config(self._path(), preset="standard").scoring
        assert s.structural_corroboration is False
        assert s.longread_structural_guard is False
        assert normalise_rescue_mode(s.backbone_intron_rescue) == "off"
        assert s.protein_support_mode == "positional"

    def test_template_carries_no_clade_specific_values(self):
        """It must not quietly ship P. falciparum tuning."""
        template = load_config(self._path(), preset="standard").scoring
        neutral = load_config(preset="standard").scoring
        assert template.weights.backbone == neutral.weights.backbone
        assert template.multi_source_bonus == neutral.multi_source_bonus


# ---------------------------------------------------------------------------
# provenance
# ---------------------------------------------------------------------------

class TestRunManifest:
    """The manifest must capture enough to identify and reproduce a run."""

    def _manifest(self, tmp_path, **kw):
        from gmb.provenance import build_run_manifest
        cfg = load_config(preset="apicomplexa")
        genome = tmp_path / "genome.fa"
        genome.write_text(">1\nACGTACGTAC\n")
        backbone = tmp_path / "backbone.gtf"
        backbone.write_text("1\ttest\texon\t1\t9\t.\t+\t.\ttranscript_id \"t1\";\n")
        return build_run_manifest(
            cfg,
            inputs={"genome": str(genome), "Tiberius": str(backbone)},
            output_dir=str(tmp_path),
            **kw,
        )

    def test_records_software_identity(self, tmp_path):
        m = self._manifest(tmp_path)
        assert m.gmb_version
        assert m.python_version
        assert m.hostname
        # git state may be absent outside a checkout, but the key must exist
        assert set(m.git) >= {"commit", "branch", "dirty"}

    def test_hashes_every_input(self, tmp_path):
        m = self._manifest(tmp_path)
        assert len(m.inputs) == 2
        for entry in m.inputs:
            assert entry["sha256"], f"{entry['label']} has no hash"
            assert entry["size_bytes"] is not None

    def test_genome_is_not_reported_as_a_weighted_evidence_track(self, tmp_path):
        """The genome is the coordinate system, not evidence."""
        m = self._manifest(tmp_path)
        genome = [i for i in m.inputs if i["label"] == "genome"][0]
        assert genome["resolved_role"] == "genome_reference_sequence"
        assert genome["resolved_weight"] is None
        assert "genome" not in m.evidence_roles

    def test_records_resolved_roles_and_weights(self, tmp_path):
        m = self._manifest(tmp_path)
        assert m.evidence_roles["Tiberius"] == "backbone"
        assert m.evidence_weights["Tiberius"] == 2.6

    def test_records_the_resolved_rescue_decision_not_just_the_mode(self, tmp_path):
        """`auto` means the outcome is decided at run time; recording the mode
        alone would hide what actually happened."""
        from gmb.pipeline.applicability import resolve_backbone_intron_rescue

        cfg = load_config(preset="apicomplexa")
        decision = resolve_backbone_intron_rescue(
            cfg.scoring, _exon_frame({"Tiberius": (173, 827), "Scallop": (987, 13)}), 1.0)
        m = self._manifest(tmp_path, rescue_decision=decision)
        rescue = m.resolved_policy["backbone_intron_rescue"]
        assert rescue["mode"] == "auto"
        assert rescue["enabled"] is True
        assert "reason" in rescue and rescue["reason"]
        # and the numbers it decided on
        assert rescue["evidence"]["backbone_multi_exon_fraction"] is not None

    def test_writes_both_json_and_tsv(self, tmp_path):
        from gmb.provenance import write_run_manifest

        m = self._manifest(tmp_path)
        json_path, tsv_path = write_run_manifest(m, str(tmp_path))
        assert os.path.exists(json_path) and os.path.exists(tsv_path)
        data = json.load(open(json_path))
        assert data["gmb_version"] == m.gmb_version
        header = open(tsv_path).readline().rstrip("\n").split("\t")
        assert header == ["field", "value"]


# ---------------------------------------------------------------------------
# build -> finalise configuration handoff
# ---------------------------------------------------------------------------

class TestFinaliseConfigHandoff:
    """Finalisation must run under the configuration that produced the build.

    Anything else can silently finalise an annotation under settings that never
    applied to it -- e.g. selecting canonical transcripts under different weights
    than the ones that chose the models.
    """

    def _write_build_config(self, tmp_path, preset="apicomplexa"):
        from gmb.pipeline.config import load_config
        import yaml
        from dataclasses import asdict
        cfg = load_config(preset=preset)
        path = tmp_path / "resolved_config.yaml"
        path.write_text(yaml.safe_dump(asdict(cfg), sort_keys=True))
        return str(path)

    def test_build_config_is_used_by_default(self, tmp_path):
        from gmb.pipeline.finalise import _resolve_finalise_config
        resolved = self._write_build_config(tmp_path, "apicomplexa")
        cfg, source, sha = _resolve_finalise_config(resolved, "standard", None, False)
        assert source == "build_resolved_config"
        assert sha
        assert cfg.scoring.weights.backbone == 2.6          # apicomplexa, not standard

    def test_supplied_config_does_NOT_override_by_default(self, tmp_path):
        """The regression this test exists for.

        Previously finalise used the build config only when --config was absent,
        so passing --config (as every example wrapper does) silently re-resolved
        the configuration independently of the build.
        """
        from gmb.pipeline.finalise import _resolve_finalise_config
        resolved = self._write_build_config(tmp_path, "apicomplexa")
        overlay = tmp_path / "other.yaml"
        overlay.write_text("scoring:\n  weights:\n    backbone: 9.9\n")
        cfg, source, _ = _resolve_finalise_config(resolved, "fungi", [str(overlay)], False)
        assert source == "build_resolved_config"
        assert cfg.scoring.weights.backbone == 2.6, \
            "supplied --config must not silently replace the build configuration"

    def test_explicit_override_is_honoured_and_recorded(self, tmp_path):
        from gmb.pipeline.finalise import _resolve_finalise_config
        resolved = self._write_build_config(tmp_path, "apicomplexa")
        overlay = tmp_path / "other.yaml"
        overlay.write_text("scoring:\n  weights:\n    backbone: 9.9\n")
        cfg, source, _ = _resolve_finalise_config(
            resolved, "standard", [str(overlay)], True)
        assert source == "explicit_override"
        assert cfg.scoring.weights.backbone == 9.9

    def test_falls_back_to_supplied_config_when_build_has_none(self, tmp_path):
        from gmb.pipeline.finalise import _resolve_finalise_config
        missing = str(tmp_path / "does_not_exist.yaml")
        cfg, source, sha = _resolve_finalise_config(missing, "fungi", None, False)
        assert source == "supplied_preset_config"
        assert sha is None
        assert cfg.scoring.weights.backbone == 3.1          # fungi

    def test_divergence_is_reported(self, tmp_path, capsys):
        """A caller passing a different config must be told, not ignored silently."""
        from gmb.pipeline.finalise import _resolve_finalise_config
        resolved = self._write_build_config(tmp_path, "apicomplexa")
        overlay = tmp_path / "other.yaml"
        overlay.write_text("scoring:\n  weights:\n    backbone: 9.9\n")
        _resolve_finalise_config(resolved, "standard", [str(overlay)], False)
        out = capsys.readouterr().out
        assert "WARNING" in out
        assert "override-build-config" in out

    def test_identical_supplied_config_produces_no_warning(self, tmp_path, capsys):
        """The common case -- reproducing the build invocation -- must be quiet."""
        from gmb.pipeline.finalise import _resolve_finalise_config
        resolved = self._write_build_config(tmp_path, "apicomplexa")
        _resolve_finalise_config(resolved, "apicomplexa", None, False)
        assert "WARNING" not in capsys.readouterr().out

    def test_config_differences_finds_nested_keys(self):
        from gmb.pipeline.config import load_config
        from gmb.pipeline.finalise import _config_differences
        a, b = load_config(preset="apicomplexa"), load_config(preset="fungi")
        keys = {k for k, _, _ in _config_differences(a, b)}
        assert "scoring.weights.backbone" in keys
        assert not _config_differences(a, load_config(preset="apicomplexa"))
