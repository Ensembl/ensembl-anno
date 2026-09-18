#!/usr/bin/env python3
"""Selection logic must depend on evidence ROLES, never on literal tool names.

These tests deliberately use invented source names (``PredictorX``,
``AssemblerA`` ...) so that they pass only if the implementation is genuinely
generic. Real tool names appear only in the end-to-end/golden tests.
"""

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from gmb.pipeline.canonical_evidence import (
    EVIDENCE_CLASS_BACKBONE,
    EVIDENCE_CLASS_LONG_READ,
    EVIDENCE_CLASS_OTHER,
    EVIDENCE_CLASS_PROTEIN_ALIGNMENT,
    EVIDENCE_CLASS_SHORT_READ,
    EvidenceRoles,
)
from gmb.pipeline.config import load_config, validate_selection_policy
from gmb.pipeline.applicability import normalise_rescue_mode
from gmb.pipeline.scoring import (
    TIER_BACKBONE_INTRON_RESCUE,
    TIER_BACKBONE_SHORTREAD_AGREEMENT,
    TIER_MULTI_SOURCE_AGREEMENT,
    TIER_SINGLE_SOURCE,
    rank_tier,
    rescue_enabled,
    select_isoforms,
    selection_reason_of,
)

BACKBONE = "PredictorX"
SR_A, SR_B = "AssemblerA", "AssemblerB"
LONGREAD = "LongReaderZ"
PROTEIN = "ProteinAligner1"


@pytest.fixture
def config():
    cfg = load_config()
    s = cfg.scoring
    s.backbone_label = BACKBONE
    s.shortread_labels = [SR_A, SR_B]
    s.longread_label = LONGREAD
    s.protein_alignment_labels = [PROTEIN]
    return cfg


def _model(tid, source, exons):
    return pd.DataFrame({
        "Chromosome": ["1"] * len(exons),
        "Start": [s for s, _e in exons],
        "End": [e for _s, e in exons],
        "Strand": ["+"] * len(exons),
        "Source": [source] * len(exons),
        "transcript_id": [tid] * len(exons),
    })


def _locus(*frames):
    return pd.concat(frames, ignore_index=True)


CHAIN_A = [(100, 300), (500, 800)]
CHAIN_B = [(100, 300), (600, 800)]


class TestEvidenceRoleResolution:
    def test_configured_labels_resolve_to_roles(self, config):
        r = EvidenceRoles.from_config(config.scoring)
        assert r.role_of(BACKBONE) == EVIDENCE_CLASS_BACKBONE
        assert r.role_of(SR_A) == EVIDENCE_CLASS_SHORT_READ
        assert r.role_of(SR_B) == EVIDENCE_CLASS_SHORT_READ
        assert r.role_of(LONGREAD) == EVIDENCE_CLASS_LONG_READ
        assert r.role_of(PROTEIN) == EVIDENCE_CLASS_PROTEIN_ALIGNMENT

    def test_resolution_is_case_insensitive(self, config):
        r = EvidenceRoles.from_config(config.scoring)
        assert r.role_of(BACKBONE.upper()) == EVIDENCE_CLASS_BACKBONE
        assert r.role_of(SR_A.lower()) == EVIDENCE_CLASS_SHORT_READ

    def test_unknown_source_is_other(self, config):
        r = EvidenceRoles.from_config(config.scoring)
        assert r.role_of("SomethingNobodyConfigured") == EVIDENCE_CLASS_OTHER
        assert r.role_of("") == EVIDENCE_CLASS_OTHER
        assert r.role_of(None) == EVIDENCE_CLASS_OTHER

    def test_builtin_fallback_for_unconfigured_known_tools(self):
        """A source the config does not mention still resolves if it is a
        well-known tool, so partial configuration degrades gracefully."""
        r = EvidenceRoles(backbone_label="PredictorX")
        assert r.role_of("Scallop") == EVIDENCE_CLASS_SHORT_READ
        assert r.role_of("PredictorX") == EVIDENCE_CLASS_BACKBONE

    def test_helpers_agree_with_role_of(self, config):
        r = EvidenceRoles.from_config(config.scoring)
        assert r.is_backbone(BACKBONE) and not r.is_backbone(SR_A)
        assert r.is_shortread(SR_A) and not r.is_shortread(LONGREAD)
        assert r.is_longread(LONGREAD) and not r.is_longread(SR_B)
        assert r.is_assembled_transcript(SR_A)
        assert r.is_assembled_transcript(LONGREAD)
        assert not r.is_assembled_transcript(BACKBONE)


class TestSelectionIsRoleDriven:
    def test_backbone_shortread_agreement_uses_configured_labels(self, config):
        """Agreement is detected between arbitrary tool names."""
        config.scoring.structural_corroboration = True
        locus = _locus(_model("b", BACKBONE, CHAIN_A), _model("a", SR_A, CHAIN_A))
        genes = select_isoforms(locus, config, {"b", "a"})
        top = genes[0][0]
        assert top["backbone_shortread_agreement"] is True
        assert top["selection_reason"] == "backbone_shortread_agreement"

    def test_longread_guard_uses_configured_longread_label(self, config):
        config.scoring.longread_structural_guard = True
        locus = _locus(_model("lr", LONGREAD, CHAIN_B), _model("a", SR_A, CHAIN_A))
        genes = select_isoforms(locus, config, {"lr", "a"})
        roles = {m["id"]: m["longread_structural_role"] for g in genes for m in g}
        assert roles["lr"] == "support_only"
        assert roles["a"] == "not_longread"

    def test_backbone_intron_rescue_uses_configured_labels(self, config):
        config.scoring.backbone_intron_rescue = True
        config.scoring.structural_corroboration = True
        backbone_exons = [(1000, 1600)]
        spliced = [(1000, 1300), (1500, 1900)]
        locus = _locus(_model("b", BACKBONE, backbone_exons),
                       _model("a", SR_A, spliced))
        genes = select_isoforms(
            locus, config, {"b", "a"},
            candidate_cds={"b": backbone_exons, "a": spliced},
            canonical_intron_tids={"a"})
        assert genes[0][0]["id"] == "a"
        assert genes[0][0]["selection_reason"] == "backbone_intron_rescue"

    def test_renaming_every_track_does_not_change_the_outcome(self, config):
        """The same topology under different names selects the same structure."""
        config.scoring.structural_corroboration = True
        first = select_isoforms(
            _locus(_model("b", BACKBONE, CHAIN_A), _model("a", SR_A, CHAIN_A),
                   _model("x", SR_B, CHAIN_B)),
            config, {"b", "a", "x"})
        cfg2 = load_config()
        cfg2.scoring.structural_corroboration = True
        cfg2.scoring.backbone_label = "OtherPredictor"
        cfg2.scoring.shortread_labels = ["Asm1", "Asm2"]
        second = select_isoforms(
            _locus(_model("b", "OtherPredictor", CHAIN_A), _model("a", "Asm1", CHAIN_A),
                   _model("x", "Asm2", CHAIN_B)),
            cfg2, {"b", "a", "x"})
        assert [m["id"] for g in first for m in g] == [m["id"] for g in second for m in g]
        assert first[0][0]["selection_reason"] == second[0][0]["selection_reason"]


class TestRankingHierarchy:
    """The hierarchy is defined once and must stay ordered."""

    def test_tier_order(self):
        assert (TIER_BACKBONE_SHORTREAD_AGREEMENT
                > TIER_MULTI_SOURCE_AGREEMENT
                > TIER_BACKBONE_INTRON_RESCUE
                > TIER_SINGLE_SOURCE)

    @staticmethod
    def _struct(**kw):
        base = dict(backbone_shortread_agreement=False, n_structural_support_sources=1,
                    backbone_intron_rescue=False, protein_cds_span=False,
                    longread_demoted=False, is_longread_only=False,
                    rep={"exon_count": 2})
        base.update(kw)
        return base

    def test_tiers_are_assigned_in_order(self):
        assert rank_tier(self._struct(backbone_shortread_agreement=True)) == \
            TIER_BACKBONE_SHORTREAD_AGREEMENT
        assert rank_tier(self._struct(n_structural_support_sources=2)) == \
            TIER_MULTI_SOURCE_AGREEMENT
        assert rank_tier(self._struct(backbone_intron_rescue=True)) == \
            TIER_BACKBONE_INTRON_RESCUE
        assert rank_tier(self._struct()) == TIER_SINGLE_SOURCE

    def test_protein_support_is_not_a_ranking_tier(self):
        """Protein support influences the score and retention, never the tier."""
        assert rank_tier(self._struct(protein_cds_span=True)) == TIER_SINGLE_SOURCE

    def test_reason_matches_tier(self):
        assert selection_reason_of(self._struct(backbone_shortread_agreement=True)) == \
            "backbone_shortread_agreement"
        assert selection_reason_of(self._struct(n_structural_support_sources=2)) == \
            "multi_source_structural_agreement"
        assert selection_reason_of(self._struct(backbone_intron_rescue=True)) == \
            "backbone_intron_rescue"
        assert selection_reason_of(self._struct(protein_cds_span=True)) == \
            "protein_cds_span_support"
        assert selection_reason_of(self._struct()) == "best_single_source"

    def test_demoted_longread_reported_regardless_of_tier(self):
        assert selection_reason_of(
            self._struct(longread_demoted=True, n_structural_support_sources=2)
        ) == "longread_demoted_support_only"


class TestSelectionPolicyValidation:
    def test_clean_config_produces_no_warnings(self, config):
        assert validate_selection_policy(config) == []

    def test_source_in_two_roles_is_fatal(self, config):
        config.scoring.longread_label = SR_A
        with pytest.raises(ValueError, match="two roles"):
            validate_selection_policy(config)

    def test_invalid_protein_support_mode_is_fatal(self, config):
        config.scoring.protein_support_mode = "not_a_mode"
        with pytest.raises(ValueError, match="protein_support_mode"):
            validate_selection_policy(config)

    def test_rescue_without_backbone_warns(self, config):
        config.scoring.backbone_intron_rescue = True
        config.scoring.backbone_label = ""
        assert any("backbone_intron_rescue" in w for w in validate_selection_policy(config))

    def test_guard_without_longread_source_warns(self, config):
        config.scoring.longread_structural_guard = True
        config.scoring.longread_label = ""
        assert any("longread_structural_guard" in w
                   for w in validate_selection_policy(config))

    def test_corroboration_with_one_structural_role_warns(self, config):
        config.scoring.structural_corroboration = True
        config.scoring.shortread_labels = []
        config.scoring.longread_label = ""
        assert any("structural_corroboration" in w
                   for w in validate_selection_policy(config))

    def test_ineffective_config_warns_but_does_not_raise(self, config):
        config.scoring.longread_structural_guard = True
        config.scoring.longread_label = ""
        validate_selection_policy(config)  # must not raise


class TestDefaultsAreBehaviourNeutral:
    def test_all_optional_policies_default_off(self):
        s = load_config().scoring
        assert s.structural_corroboration is False
        assert s.longread_structural_guard is False
        # Now a mode string; the invariant that matters is that it cannot fire.
        assert normalise_rescue_mode(s.backbone_intron_rescue) == "off"
        assert rescue_enabled(s) is False
        assert s.protein_support_mode == "positional"

    def test_neutral_preset_is_explicitly_nameable(self):
        """`--preset standard` must mean the same thing as omitting the preset."""
        assert load_config(preset="standard").scoring == load_config(preset=None).scoring

    def test_fungi_preset_enables_no_experimental_policy(self):
        """The fungal default is the VALIDATED BASELINE behaviour.

        Applying the Apicomplexa policy unchanged to the full Z. tritici genome
        made the annotation worse (CDS exact 4,895 -> 4,611; improvements :
        regressions = 249 : 555). Until a per-flag fungal experiment exists, the
        fungal preset must not enable any of them.
        """
        s = load_config(preset="fungi").scoring
        assert s.structural_corroboration is False
        assert s.longread_structural_guard is False
        assert normalise_rescue_mode(s.backbone_intron_rescue) == "off"
        assert rescue_enabled(s) is False
        assert s.protein_support_mode == "positional"

    def test_apicomplexa_preset_enables_exactly_the_validated_policy(self):
        """Apicomplexa deliberately enables the policy validated on P. falciparum
        and T. gondii -- and nothing beyond it.

        Pinned explicitly so that enabling a further policy, or quietly changing
        one of these, has to be a deliberate edit to this test.
        """
        s = load_config(preset="apicomplexa").scoring
        assert s.structural_corroboration is True
        assert s.longread_structural_guard is True
        assert s.protein_support_mode == "cds_span_compatible"
        # Applicability-gated, NOT unconditionally on: rescue was 68.6% CDS-exact
        # against a Tiberius backbone and 1.8% against a Helixer one.
        assert normalise_rescue_mode(s.backbone_intron_rescue) == "auto"
        # "auto" alone must never fire without run-level evidence to justify it.
        assert rescue_enabled(s) is False

    def test_no_preset_sets_a_long_read_disposition_other_than_default(self):
        """Long-read handling is decided by preflight + operator, not baked into
        a clade preset."""
        for preset in (None, "fungi", "apicomplexa"):
            s = load_config(preset=preset).scoring
            assert s.longread_disposition == "primary_structural"

    def test_removed_options_are_gone(self):
        """Failed experiments must not linger as dead configuration."""
        s = load_config().scoring
        for gone in ("longread_guard_scope",
                     "longread_guard_requires_credible_alternative",
                     "low_confidence_locus_fallback"):
            assert not hasattr(s, gone), f"{gone} should have been removed"
