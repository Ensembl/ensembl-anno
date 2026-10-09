"""Regression tests for the evidence-integration changes.

Covers junction-level evidence (gmb.pipeline.junction_support), the opt-in
``scoring.primary_selection: junction_supported`` rule, canonical selection's
``prefer_build_primary``, the now-functional ``scoring.min_cds_bp`` floor and the
configurable protein-validation penalty. Defaults must leave behaviour unchanged.
See docs/release/evidence_integration.md for the Z. tritici measurements.
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from gmb.pipeline.canonical_selection import select_canonical_for_gene
from gmb.pipeline.config import load_config
from gmb.pipeline.junction_support import (
    ProteinIntronIndex,
    observed_introns,
    protein_compatibility,
    protein_compatibility_counts,
    unsupported_intron_counts,
)
from gmb.pipeline.scoring import score_model, select_isoforms

GMB_DIR = Path(__file__).resolve().parent.parent
FIXTURE = GMB_DIR / "tests" / "fixtures" / "z_tritici_region1"


def _exons(tid, source, exons, strand="+", chrom="1"):
    return pd.DataFrame({
        "Chromosome": chrom, "Start": [s for s, _ in exons], "End": [e for _, e in exons],
        "Strand": strand, "transcript_id": tid, "Source": source, "Feature": "exon",
    })


# ---------------------------------------------------------------------------
# junction_support
# ---------------------------------------------------------------------------


class TestJunctionSupport:
    def test_unsupported_introns_counted_against_observed_junctions(self):
        tx = _exons("st1", "StringTie", [(100, 300), (400, 700), (800, 1000)])
        cand = pd.concat([
            _exons("hx", "Helixer", [(100, 300), (420, 700), (800, 1000)]),
            _exons("st1", "StringTie", [(100, 300), (400, 700), (800, 1000)]),
            _exons("one", "Helixer", [(100, 1000)]),
        ])
        counts = unsupported_intron_counts(cand, observed_introns(tx))
        assert counts["hx"] == (2, 1)       # (300,420) unobserved, (700,800) observed
        assert counts["st1"] == (2, 0)
        assert counts["one"] == (0, 0)

    def test_observed_introns_are_strand_specific(self):
        tx = _exons("st1", "StringTie", [(100, 300), (400, 700)], strand="-")
        cand = _exons("hx", "Helixer", [(100, 300), (400, 700)], strand="+")
        assert unsupported_intron_counts(cand, observed_introns(tx))["hx"] == (1, 1)

    def test_protein_compatibility(self):
        prot = pd.concat([
            _exons("same", "OrthoDB", [(150, 300), (400, 600)]),          # agrees
            _exons("other", "OrthoDB", [(150, 320), (400, 600)]),         # other donor
            _exons("unspliced", "OrthoDB", [(420, 650)]),                 # inside one exon
            _exons("graze", "OrthoDB", [(0, 160)]),                       # < 30 bp overlap
            _exons("antisense", "OrthoDB", [(150, 300), (400, 600)], strand="-"),
        ])
        index = ProteinIntronIndex(prot)
        assert protein_compatibility([(150, 300), (400, 700)], "1", "+", index) == (1, 1)
        cand = _exons("c", "StringTie", [(100, 300), (400, 800)])
        assert protein_compatibility_counts(cand, {"c": [(150, 300), (400, 700)]}, prot) == {
            "c": (1, 1)}

    def test_no_cds_or_no_proteins_gives_nothing(self):
        cand = _exons("c", "StringTie", [(100, 300), (400, 800)])
        assert protein_compatibility_counts(cand, {}, _exons("p", "OrthoDB", [(100, 200)])) == {}
        assert protein_compatibility_counts(cand, {"c": [(150, 300)]}, pd.DataFrame()) == {}


# ---------------------------------------------------------------------------
# primary_selection: junction_supported
# ---------------------------------------------------------------------------

HX_EXONS = [(100, 300), (400, 700), (800, 1000)]
SR_EXONS = [(100, 300), (420, 700), (800, 1000)]   # alternative acceptor, shares (700,800)
HX_CDS = [(150, 300), (400, 700), (800, 950)]      # 600 bp
SR_CDS = [(150, 300), (420, 700), (800, 950)]      # 580 bp


def _locus():
    return pd.concat([
        _exons("hx", "Helixer", HX_EXONS),
        _exons("sc", "Scallop", SR_EXONS),
        _exons("st", "StringTie", SR_EXONS),
    ], ignore_index=True)


def _select(mode="junction_supported", junctions=None, compat=None, complete=None,
            cds=None, min_frac=None):
    cfg = load_config(None, "fungi")
    cfg.scoring.primary_selection = mode
    if min_frac is not None:
        cfg.scoring.junction_primary_min_cds_fraction = min_frac
    return select_isoforms(
        _locus(), cfg, set(),
        candidate_cds=cds or {"hx": HX_CDS, "sc": SR_CDS, "st": SR_CDS},
        junction_support=junctions or {"hx": (2, 1), "sc": (2, 0), "st": (2, 0)},
        protein_compatibility=compat or {},
        complete_orf_tids={"hx", "sc", "st"} if complete is None else complete,
    )


class TestJunctionSupportedPrimary:
    def test_default_score_mode_keeps_backbone_primary(self):
        genes = _select(mode="score")
        assert len(genes) == 1
        assert genes[0][0]["id"] == "hx"
        assert genes[0][0]["introns_without_transcript_support"] == 1  # still reported

    def test_promotes_transcript_supported_alternate_and_keeps_backbone(self):
        genes = _select()
        assert len(genes) == 1
        assert genes[0][0]["id"] == "sc"
        assert genes[0][0]["is_primary"] is True
        assert genes[0][0]["selection_reason"] == "transcript_junction_support"
        assert [m["id"] for m in genes[0][1:]] == ["hx"]
        assert genes[0][1]["is_primary"] is False

    def test_no_promotion_when_backbone_introns_all_observed(self):
        genes = _select(junctions={"hx": (2, 0), "sc": (2, 0), "st": (2, 0)})
        assert genes[0][0]["id"] == "hx"

    def test_no_promotion_without_complete_orf(self):
        assert _select(complete={"hx"})[0][0]["id"] == "hx"

    def test_no_promotion_when_protein_evidence_worse(self):
        genes = _select(compat={"hx": (3, 0), "sc": (0, 2), "st": (0, 2)})
        assert genes[0][0]["id"] == "hx"

    def test_no_promotion_when_cds_much_shorter(self):
        short = [(800, 950)]
        assert _select(cds={"hx": HX_CDS, "sc": short, "st": short})[0][0]["id"] == "hx"

    def test_no_promotion_when_cds_is_a_different_coding_locus(self):
        # Regression: a read-through alternate admitted through its span had its
        # CDS beyond the backbone's. Promoting it lost the backbone model.
        far = [(1500, 2100)]
        genes = _select(cds={"hx": HX_CDS, "sc": far, "st": far}, min_frac=0.0)
        assert genes[0][0]["id"] == "hx"


# ---------------------------------------------------------------------------
# scoring.isoform_cds_overlap
# ---------------------------------------------------------------------------

READTHROUGH_EXONS = [(100, 300), (420, 700), (800, 1600)]  # shares intron (700,800)
READTHROUGH_CDS = [(1000, 1500)]                            # no base shared with HX_CDS


def _select_iso(mode, protein_tids=()):
    cfg = load_config(None, "fungi")
    cfg.scoring.isoform_cds_overlap = mode
    locus = pd.concat([
        _exons("hx", "Helixer", HX_EXONS),
        _exons("sc", "Scallop", READTHROUGH_EXONS),
        _exons("st", "StringTie", READTHROUGH_EXONS),
    ], ignore_index=True)
    return select_isoforms(
        locus, cfg, set(protein_tids),
        candidate_cds={"hx": HX_CDS, "sc": READTHROUGH_CDS, "st": READTHROUGH_CDS})


class TestIsoformCdsOverlap:
    def test_default_off_admits_a_different_orf_as_isoform(self):
        assert load_config(None, "fungi").scoring.isoform_cds_overlap == "off"
        genes = _select_iso("off")
        assert [[m["id"] for m in g] for g in genes] == [["hx", "sc"]]

    def test_drop_discards_structure_without_shared_coding_sequence(self):
        assert [[m["id"] for m in g] for g in _select_iso("drop")] == [["hx"]]

    def test_new_gene_gives_it_its_own_gene(self):
        assert sorted(m["id"] for g in _select_iso("new_gene") for m in g) == ["hx", "sc"]
        assert len(_select_iso("new_gene")) == 2

    def test_drop_never_discards_a_backbone_model(self):
        # Protein support lets the read-through outscore the backbone; the
        # backbone's CDS then shares nothing with the primary's.
        genes = _select_iso("drop", protein_tids={"sc", "st"})
        assert sorted(m["id"] for g in genes for m in g) == ["hx", "sc"]
        assert len(genes) == 2

    def test_bare_yaml_off_is_accepted(self, tmp_path):
        p = tmp_path / "off.yaml"
        p.write_text("scoring:\n  isoform_cds_overlap: off\n")
        assert load_config([str(p)], "fungi").scoring.isoform_cds_overlap == "off"


# ---------------------------------------------------------------------------
# canonical_selection.prefer_build_primary
# ---------------------------------------------------------------------------


def _crec(tid, sources, gmb_score):
    return {
        "gene_id": "G1", "transcript_id": tid, "evidence_sources": sources, "exon_count": 2,
        "cds_bp": 900, "transcript_span_bp": 1200, "gmb_score": gmb_score,
        "diamond_hit": None, "diamond_pident": None, "diamond_qcov": None,
        "diamond_scov": None, "diamond_bitscore": None, "diamond_evalue": None,
        "psauron_score": None, "protein_length": 300, "orf_label": None,
        "is_partial_5": False, "is_partial_3": False, "internal_stop_count": 0,
        "protein_coding_score": None,
    }


class TestPreferBuildPrimary:
    RECORDS = [_crec("G1.t1", "Helixer", 5.1), _crec("G1.t2", "Scallop,StringTie", 3.0)]

    def test_default_lets_two_short_read_sources_displace_build_primary(self):
        cfg = load_config(None, "fungi").canonical_selection
        assert cfg.prefer_build_primary is False
        result = select_canonical_for_gene("G1", self.RECORDS, cfg, "Helixer")
        assert result["canonical_transcript_id"] == "G1.t2"

    def test_prefer_build_primary_keeps_t1(self):
        cfg = load_config(None, "fungi").canonical_selection
        cfg.prefer_build_primary = True
        result = select_canonical_for_gene("G1", self.RECORDS, cfg, "Helixer")
        assert result["canonical_transcript_id"] == "G1.t1"

    def test_absorbed_t1_of_another_gene_is_not_this_genes_primary(self):
        cfg = load_config(None, "fungi").canonical_selection
        cfg.prefer_build_primary = True
        records = [_crec("G2.t1", "Helixer", 5.1), _crec("G1.t2", "Scallop,StringTie", 3.0)]
        result = select_canonical_for_gene("G1", records, cfg, "Helixer")
        assert result["canonical_transcript_id"] == "G1.t2"


# ---------------------------------------------------------------------------
# scoring.min_cds_bp and protein_validation.penalty
# ---------------------------------------------------------------------------


class TestMinCdsBp:
    def _run(self, min_cds_bp):
        cfg = load_config(None, "fungi")
        cfg.scoring.min_cds_bp = min_cds_bp
        locus = pd.concat([_exons("sc", "Scallop", SR_EXONS), _exons("st", "StringTie", SR_EXONS)])
        cds = [(150, 270)]  # 120 bp
        return select_isoforms(locus, cfg, set(), candidate_cds={"sc": cds, "st": cds})

    def test_shipped_default_is_off(self):
        assert load_config(None, "fungi").scoring.min_cds_bp == 0
        assert load_config(None, "apicomplexa").scoring.min_cds_bp == 0
        assert len(self._run(0)) == 1

    def test_floor_now_applies(self):
        assert self._run(150) == []


class TestProteinValidationPenalty:
    def test_penalty_is_configurable_and_defaults_to_former_value(self):
        cfg = load_config(None, "fungi")
        assert cfg.protein_validation.penalty == 5.0
        cfg.protein_validation.enabled = True
        cfg.protein_validation.policy = "penalize"
        model = {"id": "t", "source": "Helixer", "combined_evidence": "Helixer",
                 "exon_count": 1, "protein_coding_score": 0.1}
        base = score_model(dict(model, protein_coding_score=0.99), cfg, set())
        assert base - score_model(dict(model), cfg, set()) == pytest.approx(5.0)
        cfg.protein_validation.penalty = 1.5
        assert base - score_model(dict(model), cfg, set()) == pytest.approx(1.5)

    def test_negative_penalty_rejected(self, tmp_path):
        p = tmp_path / "bad.yaml"
        p.write_text("protein_validation:\n  penalty: -1\n")
        with pytest.raises(ValueError):
            load_config([str(p)], "fungi")


def test_invalid_primary_selection_rejected(tmp_path):
    p = tmp_path / "bad.yaml"
    p.write_text("scoring:\n  primary_selection: best\n")
    with pytest.raises(ValueError):
        load_config([str(p)], "fungi")


# ---------------------------------------------------------------------------
# Fixture builds (the bundled 500 kb Z. tritici region)
# ---------------------------------------------------------------------------


def _build(out_dir, overlay_text=None):
    cmd = [sys.executable, "-m", "gmb.cli.build", "--preset", "fungi",
           "--scallop", str(FIXTURE / "scallop_geneset.gtf"),
           "--stringtie", str(FIXTURE / "stringtie_geneset.gtf"),
           "--helixer", str(FIXTURE / "helixer_remapped.gff3"),
           "--orthodb", str(FIXTURE / "orthodb_geneset.gtf"),
           "--uniprot", str(FIXTURE / "uniprot_geneset.gtf"),
           "--genome", str(FIXTURE / "genome.fa"),
           "--output-dir", str(out_dir), "--gene-prefix", "ZT", "--no-log-file"]
    if overlay_text:
        overlay = out_dir.parent / f"{out_dir.name}.yaml"
        overlay.write_text(overlay_text)
        cmd += ["--config", str(overlay)]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(GMB_DIR))
    assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-2000:]
    return pd.read_csv(out_dir / "evidence_attribution.tsv", sep="\t")


def _structures(att):
    """Set of (structure, evidence) for every output transcript, ignoring IDs."""
    return sorted(zip(att.exon_count, att.cds_bp, att.transcript_span_bp, att.evidence_sources))


@pytest.mark.integration
def test_fixture_junction_supported_reorders_but_keeps_transcripts(tmp_path):
    base = _build(tmp_path / "score")
    jp = _build(tmp_path / "jp", "scoring:\n  primary_selection: junction_supported\n")
    # Junction evidence is always reported.
    assert base["introns_without_transcript_support"].notna().all()
    assert (base["protein_alignments_compatible"] >= 0).all()
    # The default never promotes; the opt-in mode does, without losing models.
    assert not (base["selection_reason"] == "transcript_junction_support").any()
    promoted = jp[jp["selection_reason"] == "transcript_junction_support"]
    assert len(promoted) > 0
    assert promoted["transcript_id"].str.endswith(".t1").all()
    assert len(jp) == len(base)
    assert jp["gene_id"].nunique() == base["gene_id"].nunique()
    assert _structures(jp) == _structures(base)
