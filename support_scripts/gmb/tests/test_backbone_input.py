#!/usr/bin/env python3
"""Any predictor may fill the ``backbone`` evidence role.

Following the convention of ``test_evidence_roles.py``, the generic behaviour is
tested with an invented predictor name (``PredictorZ``) so these tests pass only
if the implementation is genuinely tool-agnostic. A realistic Vipsania-shaped
GFF3 fixture is used only where the point is file-format integration.
"""

import gzip
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
from gmb.pipeline.backbone import (
    DEFAULT_BACKBONE_LABEL,
    BackboneInputError,
    resolve_backbone_input,
    source_label_from_annotation,
)
from gmb.pipeline.canonical_evidence import EVIDENCE_CLASS_BACKBONE, EvidenceRoles
from gmb.pipeline.config import load_config
from gmb.pipeline.scoring import weights_for_role

# A Vipsania-shaped record set: GFF3, gene -> mRNA -> exon/CDS, ID/Parent,
# CDS == exon (no UTR), phase populated, stop codon inside the CDS.
VIPSANIA_GFF3 = """##gff-version 3
1\tPredictorZ\tgene\t29510\t37126\t.\t+\t.\tID=g1;Name=g1
1\tPredictorZ\tmRNA\t29510\t37126\t.\t+\t.\tID=g1_t1;Name=g1_t1;Parent=g1
1\tPredictorZ\tCDS\t29510\t34762\t.\t+\t0\tID=g1_t1.CDS.1;Parent=g1_t1
1\tPredictorZ\texon\t29510\t34762\t.\t+\t0\tID=g1_t1.exon.1;Parent=g1_t1
1\tPredictorZ\tCDS\t35888\t37126\t.\t+\t0\tID=g1_t1.CDS.2;Parent=g1_t1
1\tPredictorZ\texon\t35888\t37126\t.\t+\t0\tID=g1_t1.exon.2;Parent=g1_t1
1\tPredictorZ\tgene\t38982\t40207\t.\t-\t.\tID=g2;Name=g2
1\tPredictorZ\tmRNA\t38982\t40207\t.\t-\t.\tID=g2_t1;Name=g2_t1;Parent=g2
1\tPredictorZ\tCDS\t38982\t39923\t.\t-\t0\tID=g2_t1.CDS.1;Parent=g2_t1
1\tPredictorZ\texon\t38982\t39923\t.\t-\t0\tID=g2_t1.exon.1;Parent=g2_t1
"""

TIBERIUS_GTF = (
    '1\tTiberius\tgene\t29510\t34957\t.\t+\t.\tgene_id "g1";\n'
    '1\tTiberius\ttranscript\t29510\t34957\t.\t+\t.\tgene_id "g1"; transcript_id "g1.t1";\n'
    '1\tTiberius\tCDS\t29510\t34181\t.\t+\t0\tgene_id "g1"; transcript_id "g1.t1";\n'
    '1\tTiberius\texon\t29510\t34181\t.\t+\t0\tgene_id "g1"; transcript_id "g1.t1";\n'
)


@pytest.fixture
def backbone_gff3(tmp_path):
    """A generic predictor's GFF3 backbone."""
    path = tmp_path / "predictorz.gff3"
    path.write_text(VIPSANIA_GFF3, encoding="utf8")
    return str(path)


@pytest.fixture
def tiberius_gtf(tmp_path):
    """A Tiberius GTF backbone, for regression cover."""
    path = tmp_path / "tiberius.gtf"
    path.write_text(TIBERIUS_GTF, encoding="utf8")
    return str(path)


class TestSourceLabelFromAnnotation:
    """The GFF/GTF source column is the natural home of the predictor's name."""

    def test_reads_gff3_source_column(self, backbone_gff3):
        assert source_label_from_annotation(backbone_gff3) == "PredictorZ"

    def test_reads_gtf_source_column(self, tiberius_gtf):
        assert source_label_from_annotation(tiberius_gtf) == "Tiberius"

    def test_reads_gzipped_annotation(self, tmp_path):
        path = tmp_path / "predictorz.gff3.gz"
        with gzip.open(path, "wt", encoding="utf8") as handle:
            handle.write(VIPSANIA_GFF3)
        assert source_label_from_annotation(str(path)) == "PredictorZ"

    def test_mixed_sources_give_no_label(self, tmp_path):
        path = tmp_path / "mixed.gff3"
        path.write_text(
            VIPSANIA_GFF3 + "1\tOtherTool\tgene\t1\t100\t.\t+\t.\tID=g9\n", encoding="utf8"
        )
        assert source_label_from_annotation(str(path)) is None

    def test_placeholder_source_gives_no_label(self, tmp_path):
        path = tmp_path / "dots.gff3"
        path.write_text("1\t.\tgene\t1\t100\t.\t+\t.\tID=g1\n", encoding="utf8")
        assert source_label_from_annotation(str(path)) is None

    def test_missing_file_gives_no_label(self, tmp_path):
        assert source_label_from_annotation(str(tmp_path / "absent.gff3")) is None


class TestResolveBackboneInput:
    """One backbone role; several ways to supply it."""

    def test_generic_backbone_takes_label_from_the_file(self, backbone_gff3):
        path, label = resolve_backbone_input(backbone=backbone_gff3)
        assert (path, label) == (backbone_gff3, "PredictorZ")

    def test_explicit_label_wins_over_the_file(self, backbone_gff3):
        _, label = resolve_backbone_input(backbone=backbone_gff3, backbone_label="Renamed")
        assert label == "Renamed"

    def test_unlabelled_backbone_falls_back_to_default(self, tmp_path):
        path = tmp_path / "nosource.gff3"
        path.write_text("1\t.\tgene\t1\t100\t.\t+\t.\tID=g1\n", encoding="utf8")
        _, label = resolve_backbone_input(backbone=str(path))
        assert label == DEFAULT_BACKBONE_LABEL

    def test_tiberius_flag_unchanged(self, tiberius_gtf):
        """Regression: the historical flag keeps its historical label."""
        assert resolve_backbone_input(tiberius=tiberius_gtf) == (tiberius_gtf, "Tiberius")

    def test_helixer_flag_unchanged(self, backbone_gff3):
        """Regression: --helixer keeps labelling as Helixer regardless of file content."""
        assert resolve_backbone_input(helixer=backbone_gff3) == (backbone_gff3, "Helixer")

    def test_explicit_label_can_rename_a_legacy_flag(self, tiberius_gtf):
        _, label = resolve_backbone_input(tiberius=tiberius_gtf, backbone_label="PredictorZ")
        assert label == "PredictorZ"

    def test_no_backbone_keeps_historical_default_label(self):
        """A run without any backbone must behave exactly as before."""
        assert resolve_backbone_input() == (None, "Helixer")

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"tiberius": "a.gtf", "backbone": "b.gff3"},
            {"helixer": "a.gff3", "backbone": "b.gff3"},
            {"helixer": "a.gff3", "tiberius": "b.gtf"},
            {"helixer": "a.gff3", "tiberius": "b.gtf", "backbone": "c.gff3"},
        ],
    )
    def test_more_than_one_backbone_is_refused(self, kwargs):
        with pytest.raises(BackboneInputError):
            resolve_backbone_input(**kwargs)


class TestBackboneRoleAndWeight:
    """A custom label must resolve to the backbone ROLE and the backbone WEIGHT."""

    def test_custom_label_resolves_to_backbone_role(self):
        roles = EvidenceRoles(backbone_label="PredictorZ")
        assert roles.role_of("PredictorZ") == EVIDENCE_CLASS_BACKBONE
        assert roles.is_backbone("predictorz")  # case-insensitive

    def test_custom_backbone_gets_the_role_keyed_weight(self):
        config = load_config(None, "apicomplexa")
        config.scoring.backbone_label = "PredictorZ"
        roles = EvidenceRoles.from_config(config.scoring)
        weight = weights_for_role(config.scoring.weights, roles.role_of("PredictorZ"))
        assert weight == config.scoring.weights.backbone

    def test_known_backbones_still_resolve_without_a_configured_label(self):
        """The built-in fallback map is untouched."""
        roles = EvidenceRoles()
        assert roles.role_of("Tiberius") == EVIDENCE_CLASS_BACKBONE
        assert roles.role_of("Helixer") == EVIDENCE_CLASS_BACKBONE

    def test_preset_backbone_labels_unchanged(self):
        """Presets must not shift because a new flag exists.

        The apicomplexa preset declares Tiberius as its backbone label; a build
        that passes a backbone flag overrides this afterwards, which is what lets
        another predictor fill the role without editing the preset.
        """
        assert load_config(None, "standard").scoring.backbone_label == "Helixer"
        assert load_config(None, "apicomplexa").scoring.backbone_label == "Tiberius"


class TestEvidenceLoading:
    """A GFF3 backbone must normalise to the same candidate representation."""

    def test_gff3_backbone_parses_with_its_own_source_label(self, backbone_gff3):
        from gmb.pipeline.builder import load_evidence

        path, label = resolve_backbone_input(backbone=backbone_gff3)
        exons, cds = load_evidence(path, label)
        assert exons is not None and cds is not None
        assert set(exons["Source"].unique()) == {"PredictorZ"}
        assert set(cds["Source"].unique()) == {"PredictorZ"}
        # Parent relationships become transcript_id, namespaced by source label
        # (existing behaviour, which also carries the attribution into the IDs).
        assert set(exons["transcript_id"].unique()) == {"PredictorZ_g1_t1", "PredictorZ_g2_t1"}
        # Strand is carried through for both orientations.
        assert set(exons["Strand"].unique()) == {"+", "-"}
        # CDS == exon for this CDS-only predictor.
        assert len(cds) == len(exons) == 3

    def test_gtf_backbone_still_parses(self, tiberius_gtf):
        """Regression: the GTF path is unaffected."""
        from gmb.pipeline.builder import load_evidence

        path, label = resolve_backbone_input(tiberius=tiberius_gtf)
        exons, cds = load_evidence(path, label)
        assert set(exons["Source"].unique()) == {"Tiberius"}
        assert set(exons["transcript_id"].unique()) == {"Tiberius_g1.t1"}


class TestPreflightArgumentMapping:
    """gmb-preflight must present the same labels gmb-build will use."""

    def test_generic_backbone_appears_with_its_resolved_label(self, backbone_gff3):
        from gmb.cli.preflight import _tracks_from_args, build_parser

        args = build_parser().parse_args(
            ["--genome", "g.fa", "--backbone", backbone_gff3]
        )
        tracks = _tracks_from_args(args)
        assert {"label": "PredictorZ", "path": backbone_gff3} in tracks

    def test_tiberius_flag_still_maps_to_tiberius(self, tiberius_gtf):
        from gmb.cli.preflight import _tracks_from_args, build_parser

        args = build_parser().parse_args(["--genome", "g.fa", "--tiberius", tiberius_gtf])
        assert {"label": "Tiberius", "path": tiberius_gtf} in _tracks_from_args(args)


class TestPublicApiBackboneSlot:
    """``gmb.run_gene_model_builder`` must expose the same generic slot."""

    def test_generic_slot_accepted_and_unknown_rejected(self):
        import inspect

        from gmb.api import run_gene_model_builder

        params = inspect.signature(run_gene_model_builder).parameters
        assert params["backbone_kind"].default == "helixer"   # unchanged default
        assert "backbone_label" in params
        with pytest.raises(ValueError):
            run_gene_model_builder("g.fa", backbone_kind="not_a_slot")
