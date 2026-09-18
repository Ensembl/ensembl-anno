#!/usr/bin/env python3
"""Tests for evidence-aware candidate ranking.

Three independently configurable policies, all defaulting to legacy behaviour:

* ``scoring.structural_corroboration`` -- rank on how many independent sources
  produced the identical intron chain before falling back to the score.
* ``scoring.protein_support_mode`` -- which protein signal may gate retention.
* ``scoring.longread_structural_guard`` -- long-read models do not take primary
  structural priority where a credible multi-exon alternative exists.
"""

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from gmb.pipeline.config import load_config
from gmb.pipeline.scoring import select_isoforms
from gmb.utils.intervals import cds_span_compatible_ids
from gmb.pipeline.applicability import normalise_rescue_mode
from gmb.pipeline.scoring import rescue_enabled


@pytest.fixture
def config():
    cfg = load_config()
    # The default preset's backbone is Helixer; these tests use Tiberius as the
    # loaded backbone track, exactly as the CLI would set it.
    cfg.scoring.backbone_label = "Tiberius"
    return cfg


def _model(tid, source, exons, chrom="1", strand="+"):
    return pd.DataFrame(
        {
            "Chromosome": [chrom] * len(exons),
            "Start": [s for s, _e in exons],
            "End": [e for _s, e in exons],
            "Strand": [strand] * len(exons),
            "Source": [source] * len(exons),
            "transcript_id": [tid] * len(exons),
        }
    )


def _locus(*frames):
    return pd.concat(frames, ignore_index=True)


# Two distinct multi-exon structures at one locus.
CHAIN_A = [(100, 300), (500, 800)]
CHAIN_B = [(100, 300), (600, 800)]


class TestDefaultsAreLegacy:
    def test_new_flags_default_off(self, config):
        s = config.scoring
        assert s.structural_corroboration is False
        assert s.protein_support_mode == "positional"
        assert s.longread_structural_guard is False

    def test_legacy_ranking_unchanged_by_new_code(self, config):
        """Backbone weight still wins under the default policy."""
        locus = _locus(
            _model("tib", "Tiberius", CHAIN_A),
            _model("sc", "Scallop", CHAIN_B),
            _model("st", "StringTie", CHAIN_B),
        )
        genes = select_isoforms(locus, config, {"tib", "sc", "st"})
        top = genes[0][0]
        # Backbone weight (2.6) still beats the Scallop+StringTie structure.
        assert top["structural_support_sources"] == "Tiberius"


class TestStructuralCorroboration:
    def test_backbone_shortread_agreement_outranks_higher_score(self, config):
        """A corroborated structure is preferred over a higher-scoring lone one.

        `lone` is Tiberius-only (backbone weight 2.6); `tib`+`sc` share one
        intron chain, so that structure carries backbone+short-read agreement.
        """
        config.scoring.structural_corroboration = True
        locus = _locus(
            _model("lone", "Tiberius", [(2000, 2300), (2500, 2900)]),
            _model("tib", "Tiberius", CHAIN_A),
            _model("sc", "Scallop", CHAIN_A),
        )
        genes = select_isoforms(locus, config, {"lone", "tib", "sc"})
        top = genes[0][0]
        assert top["backbone_shortread_agreement"] is True
        assert top["selection_reason"] == "backbone_shortread_agreement"

    def test_features_are_recorded(self, config):
        locus = _locus(_model("tib", "Tiberius", CHAIN_A), _model("sc", "Scallop", CHAIN_A))
        genes = select_isoforms(locus, config, {"tib", "sc"})
        m = genes[0][0]
        assert m["n_structural_support_sources"] == 2
        assert set(m["structural_support_sources"].split(",")) == {"Scallop", "Tiberius"}
        assert m["backbone_shortread_agreement"] is True

    def test_single_source_has_no_agreement(self, config):
        locus = _locus(_model("sc", "Scallop", CHAIN_A))
        genes = select_isoforms(locus, config, {"sc"})
        m = genes[0][0]
        assert m["n_structural_support_sources"] == 1
        assert m["backbone_shortread_agreement"] is False


class TestProteinSupportMode:
    def test_positional_mode_uses_positional_support(self, config):
        locus = _locus(_model("sc", "Scallop", CHAIN_A))
        genes = select_isoforms(locus, config, {"sc"}, None, protein_cds_span_tids=set())
        assert genes, "positional support should retain this single-source model"

    def test_cds_span_mode_ignores_bare_positional_support(self, config):
        """Positional-only support no longer satisfies the retention gate."""
        config.scoring.protein_support_mode = "cds_span_compatible"
        config.scoring.require_support_for_single_exon = True
        locus = _locus(_model("sc", "Scallop", [(100, 400)]))
        genes = select_isoforms(locus, config, {"sc"}, None, protein_cds_span_tids=set())
        assert genes == []

    def test_cds_span_mode_accepts_strong_support(self, config):
        config.scoring.protein_support_mode = "cds_span_compatible"
        locus = _locus(_model("sc", "Scallop", [(100, 400)]))
        genes = select_isoforms(locus, config, {"sc"}, None, protein_cds_span_tids={"sc"})
        assert genes
        assert genes[0][0]["protein_cds_span_compatible"] is True

    def test_positional_support_still_recorded_in_cds_span_mode(self, config):
        config.scoring.protein_support_mode = "cds_span_compatible"
        locus = _locus(_model("tib", "Tiberius", CHAIN_A))
        genes = select_isoforms(locus, config, {"tib"}, None, protein_cds_span_tids=set())
        assert genes[0][0]["protein_positional_support"] is True


class TestLongreadStructuralGuard:
    def _mixed_locus(self):
        return _locus(
            _model("mm", "Minimap2", [(100, 900)]),
            _model("sc", "Scallop", CHAIN_A),
        )

    def test_guard_off_leaves_longread_eligible(self, config):
        locus = _locus(_model("mm", "Minimap2", CHAIN_B), _model("sc", "Scallop", CHAIN_A))
        genes = select_isoforms(locus, config, {"mm", "sc"})
        roles = {m["id"]: m["longread_structural_role"] for g in genes for m in g}
        assert roles.get("mm") == "primary"

    def test_guard_demotes_longread_when_multiexon_alternative_exists(self, config):
        config.scoring.longread_structural_guard = True
        locus = _locus(_model("mm", "Minimap2", CHAIN_B), _model("sc", "Scallop", CHAIN_A))
        genes = select_isoforms(locus, config, {"mm", "sc"})
        flat = [m for g in genes for m in g]
        mm = next(m for m in flat if m["id"] == "mm")
        assert mm["longread_structural_role"] == "support_only"
        # the short-read structure must be the primary model
        assert genes[0][0]["id"] == "sc"

    def test_guard_keeps_longread_when_it_is_the_only_candidate(self, config):
        config.scoring.longread_structural_guard = True
        locus = _locus(_model("mm", "Minimap2", CHAIN_A))
        genes = select_isoforms(locus, config, {"mm"})
        assert genes
        assert genes[0][0]["longread_structural_role"] == "primary"

    def test_guard_keeps_longread_when_all_alternatives_single_exon(self, config):
        config.scoring.longread_structural_guard = True
        locus = _locus(
            _model("mm", "Minimap2", CHAIN_A),
            _model("sc", "Scallop", [(100, 400)]),
        )
        genes = select_isoforms(locus, config, {"mm", "sc"})
        flat = [m for g in genes for m in g]
        mm = next(m for m in flat if m["id"] == "mm")
        assert mm["longread_structural_role"] == "primary"

    def test_guard_does_not_delete_longread_models(self, config):
        """Demotion is not removal -- the model stays available."""
        config.scoring.longread_structural_guard = True
        locus = _locus(_model("mm", "Minimap2", CHAIN_B), _model("sc", "Scallop", CHAIN_A))
        genes = select_isoforms(locus, config, {"mm", "sc"})
        assert "mm" in {m["id"] for g in genes for m in g}


class TestCdsSpanCompatibleIds:
    def test_alignment_inside_candidate_span_qualifies(self):
        cands = {"c1": ("1", "+", 100, 1000)}
        prots = {"p1": ("1", "+", 200, 800)}
        assert cds_span_compatible_ids(cands, prots, {"c1"}) == {"c1"}

    def test_alignment_overhanging_candidate_does_not_qualify(self):
        cands = {"c1": ("1", "+", 100, 1000)}
        prots = {"p1": ("1", "+", 900, 1500)}
        assert cds_span_compatible_ids(cands, prots, {"c1"}) == set()

    def test_candidate_without_cds_never_qualifies(self):
        cands = {"c1": ("1", "+", 100, 1000)}
        prots = {"p1": ("1", "+", 200, 800)}
        assert cds_span_compatible_ids(cands, prots, set()) == set()

    def test_opposite_strand_does_not_qualify(self):
        cands = {"c1": ("1", "+", 100, 1000)}
        prots = {"p1": ("1", "-", 200, 800)}
        assert cds_span_compatible_ids(cands, prots, {"c1"}) == set()

    def test_other_chromosome_does_not_qualify(self):
        cands = {"c1": ("1", "+", 100, 1000)}
        prots = {"p1": ("2", "+", 200, 800)}
        assert cds_span_compatible_ids(cands, prots, {"c1"}) == set()

    def test_picks_contained_alignment_among_several(self):
        cands = {"c1": ("1", "+", 100, 1000)}
        prots = {
            "before": ("1", "+", 0, 50),
            "overhang": ("1", "+", 950, 2000),
            "inside": ("1", "+", 300, 700),
        }
        assert cds_span_compatible_ids(cands, prots, {"c1"}) == {"c1"}

    def test_empty_inputs(self):
        assert cds_span_compatible_ids({}, {}, set()) == set()


class TestBackboneIntronRescue:
    """Ab initio backbones under-call introns, emitting one long coding exon
    where the gene is spliced. Where an assembled transcript shows a canonical
    spliced CDS that recovers MORE coding sequence and contains the backbone's
    call, the assembly is the direct observation and takes structural priority.
    """

    # Backbone: one unspliced CDS exon spanning 1000-1600 (600 bp).
    # Assembly:  spliced CDS 1000-1300 + 1500-1900 (700 bp) -- longer, and it
    # covers 500 of the backbone's 600 coding bases.
    BACKBONE_EXONS = [(1000, 1600)]
    ASSEMBLY_EXONS = [(1000, 1300), (1500, 1900)]

    def _locus(self):
        return _locus(
            _model("tib", "Tiberius", self.BACKBONE_EXONS),
            _model("st", "StringTie", self.ASSEMBLY_EXONS),
        )

    def _cds(self):
        return {"tib": self.BACKBONE_EXONS, "st": self.ASSEMBLY_EXONS}

    def _call(self, config, **kw):
        return select_isoforms(
            self._locus(), config, {"tib", "st"}, None,
            candidate_cds=kw.pop("candidate_cds", self._cds()),
            canonical_intron_tids=kw.pop("canonical_intron_tids", {"st"}),
            **kw,
        )

    def test_disabled_by_default(self, config):
        # Now a mode string; assert the rule cannot fire rather than the literal.
        assert normalise_rescue_mode(config.scoring.backbone_intron_rescue) == "off"
        assert rescue_enabled(config.scoring) is False
        genes = self._call(config)
        flat = [m for g in genes for m in g]
        assert all(not m["backbone_intron_rescue"] for m in flat)

    def test_rescue_flags_the_spliced_assembly(self, config):
        config.scoring.backbone_intron_rescue = True
        flat = [m for g in self._call(config) for m in g]
        st = next(m for m in flat if m["id"] == "st")
        assert st["backbone_intron_rescue"] is True
        assert st["selection_reason"] == "backbone_intron_rescue"

    def test_rescued_model_outranks_collapsed_backbone(self, config):
        """Backbone weight (2.6) would otherwise beat a lone StringTie (1.0)."""
        config.scoring.backbone_intron_rescue = True
        config.scoring.structural_corroboration = True
        genes = self._call(config)
        assert genes[0][0]["id"] == "st"

    def test_backbone_wins_without_the_rule(self, config):
        config.scoring.structural_corroboration = True
        genes = self._call(config)
        assert genes[0][0]["id"] == "tib"

    def test_not_triggered_when_assembly_cds_is_shorter(self, config):
        """A fragmentary assembly must not displace a correct single-exon call."""
        config.scoring.backbone_intron_rescue = True
        cds = {"tib": [(1000, 1600)], "st": [(1000, 1100), (1500, 1550)]}  # 150 bp
        flat = [m for g in self._call(config, candidate_cds=cds) for m in g]
        st = next(m for m in flat if m["id"] == "st")
        assert st["backbone_intron_rescue"] is False

    def test_not_triggered_when_backbone_cds_is_spliced(self, config):
        """Only a *collapsed* (unspliced) backbone CDS can be rescued."""
        config.scoring.backbone_intron_rescue = True
        cds = {"tib": [(1000, 1200), (1400, 1600)], "st": self.ASSEMBLY_EXONS}
        flat = [m for g in self._call(config, candidate_cds=cds) for m in g]
        st = next(m for m in flat if m["id"] == "st")
        assert st["backbone_intron_rescue"] is False

    def test_not_triggered_for_noncanonical_introns(self, config):
        config.scoring.backbone_intron_rescue = True
        flat = [m for g in self._call(config, canonical_intron_tids=set()) for m in g]
        st = next(m for m in flat if m["id"] == "st")
        assert st["backbone_intron_rescue"] is False

    def test_not_triggered_when_backbone_cds_not_contained(self, config):
        """A neighbouring spliced gene that merely overlaps must not rescue."""
        config.scoring.backbone_intron_rescue = True
        # Backbone CDS 1000-1600; assembly CDS covers only 1550-1600 of it.
        cds = {"tib": [(1000, 1600)], "st": [(1550, 1700), (1900, 2600)]}
        locus = _locus(
            _model("tib", "Tiberius", [(1000, 1600)]),
            _model("st", "StringTie", [(1550, 1700), (1900, 2600)]),
        )
        genes = select_isoforms(locus, config, {"tib", "st"}, None,
                                candidate_cds=cds, canonical_intron_tids={"st"})
        st = next(m for g in genes for m in g if m["id"] == "st")
        assert st["backbone_intron_rescue"] is False

    def test_rescue_exempts_the_single_source_protein_gate(self, config):
        """The backbone agreeing the locus is coding is corroboration enough."""
        config.scoring.backbone_intron_rescue = True
        config.scoring.protein_support_mode = "cds_span_compatible"
        config.scoring.require_support_for_single_exon = True
        genes = select_isoforms(self._locus(), config, {"tib", "st"}, None,
                                protein_cds_span_tids=set(),
                                candidate_cds=self._cds(),
                                canonical_intron_tids={"st"})
        ids = {m["id"] for g in genes for m in g}
        assert "st" in ids

    def test_rescue_does_not_outrank_corroborated_structure(self, config):
        """A structure two independent sources agree on still wins."""
        config.scoring.backbone_intron_rescue = True
        config.scoring.structural_corroboration = True
        corroborated = [(1000, 1250), (1500, 1950)]
        locus = _locus(
            _model("tib", "Tiberius", self.BACKBONE_EXONS),
            _model("st", "StringTie", self.ASSEMBLY_EXONS),
            _model("sc2", "Scallop", corroborated),
            _model("st2", "StringTie", corroborated),
        )
        cds = {**self._cds(), "sc2": corroborated, "st2": corroborated}
        genes = select_isoforms(locus, config, {"tib", "st", "sc2", "st2"}, None,
                                candidate_cds=cds,
                                canonical_intron_tids={"st", "sc2", "st2"})
        top = genes[0][0]
        assert top["n_structural_support_sources"] > 1
        assert top["selection_reason"] == "multi_source_structural_agreement"

    def test_no_cds_map_disables_the_rule(self, config):
        config.scoring.backbone_intron_rescue = True
        genes = select_isoforms(self._locus(), config, {"tib", "st"}, None,
                                candidate_cds=None, canonical_intron_tids={"st"})
        assert all(not m["backbone_intron_rescue"] for g in genes for m in g)
