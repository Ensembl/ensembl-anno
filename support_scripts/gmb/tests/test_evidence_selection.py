"""Unit tests for evidence features and the opt-in evidence-score primary rule (P4d)."""

import pandas as pd
import pytest

from gmb.pipeline.config import dump_config, load_config, validate_evidence_score
from gmb.pipeline.evidence_features import (
    CandidateEvidence,
    EvidenceParams,
    EvidenceScorer,
    EvidenceSelectionError,
    decide_primary,
    require_selection_index,
    same_coding_locus,
)
from gmb.pipeline.junction_support import ProteinIntronIndex
from gmb.pipeline.scoring import select_isoforms
from gmb.pipeline.star_junctions import JunctionIndex

PARAMS = EvidenceParams()

# ----------------------------------------------------------------- synthetic world
# One + strand sequence. CDS of the backbone: (101,200) (300,400) (500,601), i.e.
# introns (200,300) and (400,500), 300 bp in frame. The transcript alternative uses
# acceptor 303 instead of 300 for the first intron (one codon shorter, NAGNAG-like).
HX_EXONS = [(50, 200), (300, 400), (500, 700)]
TX_EXONS = [(50, 200), (303, 400), (500, 700)]
HX_CDS = [(101, 200), (300, 400), (500, 601)]
TX_CDS = [(101, 200), (303, 400), (500, 601)]


def _genome():
    seq = list("GCC" * 1200)
    for s, e in [(200, 300), (200, 303), (400, 500), (700, 1800)]:
        seq[s:s + 2] = "GT"
        seq[e - 2:e] = "AG"
    seq[598:601] = "TAA"  # stop codon closes every CDS above
    return {"1": "".join(seq)}


GENOME = _genome()


def sj_rows(rows):
    return "".join("\t".join(str(x) for x in r) + "\n" for r in rows)


def make_index(tmp_path, junctions, n_groups=5, roles=("construction",), lib_roles=None):
    """*junctions*: {(start1, end): unique reads}; every library observes every junction."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    lines = ["library_id\tpath\trole\tstudy\tstrain\tcondition\tcount_unit\tmax_intron\tmin_intron"
             "\tnoncanonical_filtered"]
    for g in range(n_groups):
        role = (lib_roles or {}).get(g, "construction")
        (tmp_path / f"L{g}.sj.tab").write_text(sj_rows(
            ("1", s, e, 1, 1, 1, reads, 0, 30) for (s, e), reads in sorted(junctions.items())))
        lines.append(f"L{g}\tL{g}.sj.tab\t{role}\tS{g}\tX\tc{g}\tmate\t1000\t21\ttrue")
    m = tmp_path / "manifest.tsv"
    m.write_text("\n".join(lines) + "\n")
    return JunctionIndex.from_manifest(str(m), roles, {"1": len(GENOME["1"])}, checksum=False)


def scorer(index, protein=None, params=PARAMS):
    return EvidenceScorer(index, GENOME, protein, params, {"Helixer"},
                          {"StringTie", "Scallop", "Minimap2"})


def ev(index, cds, sources=("Helixer",), complete=True, protein=None):
    return scorer(index, protein).evaluate("1", "+", cds, set(sources), complete)


# ----------------------------------------------------------------- feature formulas
class TestFeatures:
    def test_supported_breadth_saturates(self, tmp_path):
        idx = make_index(tmp_path, {(201, 300): 50, (401, 500): 50}, n_groups=5)
        e = ev(idx, HX_CDS)
        assert [r["J"] for r in e.introns] == [1.0, 1.0]
        assert e.J == 2.0 and e.n_negative_J == 0 and e.integrity == 0

    def test_breadth_counts_groups_not_reads(self, tmp_path):
        idx = make_index(tmp_path, {(201, 300): 5000, (401, 500): 5000}, n_groups=2)
        assert [r["J"] for r in ev(idx, HX_CDS).introns] == [0.4, 0.4]

    def test_competing_site_contradicts_and_usage_rewards(self, tmp_path):
        # Only the shifted acceptor (201-304) and the second intron are observed.
        idx = make_index(tmp_path, {(201, 303): 100, (401, 500): 100})
        hx, tx = ev(idx, HX_CDS), ev(idx, TX_CDS, ("StringTie", "Scallop"))
        assert hx.introns[0]["status"] == "ALTERNATIVE_SITE_OBSERVED"
        assert hx.introns[0]["J"] == -1.0
        assert tx.introns[0]["J"] == 1.0
        assert tx.evidence - hx.evidence == pytest.approx(2.0)

    def test_minor_site_scores_by_usage(self, tmp_path):
        idx = make_index(tmp_path, {(201, 300): 1, (201, 303): 999, (401, 500): 100})
        r = ev(idx, HX_CDS).introns[0]
        assert r["status"] == "SUPPORTED" and r["site_conflict"] == "minor_site"
        assert r["J"] == pytest.approx(1.0 * 0.001 - 1.0 * 0.999, abs=1e-3)

    def test_unobserved_intron_penalised_only_with_own_expression(self, tmp_path):
        idx = make_index(tmp_path, {(401, 500): 100})          # first intron never seen
        r = ev(idx, HX_CDS).introns
        assert r[0]["status"] == "UNOBSERVED_EXPRESSED_LOCUS" and r[0]["J"] == -0.5
        # Same intron in a model with no supported intron: no expression context, no penalty.
        lone = ev(idx, [(101, 200), (300, 380)])
        assert lone.introns[0]["J"] == 0.0

    def test_missing_evidence_is_neutral(self, tmp_path):
        idx = make_index(tmp_path, {(1201, 1300): 10})          # nothing near the model
        e = ev(idx, HX_CDS)
        assert {r["status"] for r in e.introns} == {"UNOBSERVED_NO_LOCAL_EVIDENCE"}
        assert e.J == 0.0 and e.n_negative_J == 0

    def test_long_intron_unobservable_not_penalised(self, tmp_path):
        idx = make_index(tmp_path, {(201, 300): 50})
        e = ev(idx, [(101, 200), (300, 700), (1800, 1900)])   # (700,1800) = 1,100 bp
        long = [r for r in e.introns if r["length"] > 1000][0]
        assert long["status"] == "UNOBSERVABLE_LONG" and long["J"] == 0.0

    def test_noncanonical_motif_penalised_when_not_contradicted(self, tmp_path):
        idx = make_index(tmp_path, {(1201, 1300): 10})
        e = ev(idx, [(101, 210), (290, 601)])                 # motif not canonical
        assert e.introns[0]["genome_motif"] != "GT/AG" and e.motif == -0.5

    def test_protein_agreement_and_conflict(self, tmp_path):
        idx = make_index(tmp_path, {(201, 300): 50, (401, 500): 50})
        prot = ProteinIntronIndex(pd.DataFrame({
            "transcript_id": ["p1", "p1", "p2", "p2", "p3"], "Chromosome": "1", "Strand": "+",
            "Start": [120, 300, 120, 303, 120], "End": [200, 380, 200, 380, 380]}))
        e = ev(idx, HX_CDS, protein=prot)
        r = e.introns[0]
        # p1 agrees, p2 has a different overlapping intron, p3 is unspliced (uninformative).
        assert (r["protein_agree"], r["protein_conflict"]) == (1, 1)
        assert r["P"] == 0.0
        assert e.introns[1]["P"] == 0.0                         # no informative alignment

    def test_integrity_partial_orf_and_internal_stop(self, tmp_path):
        idx = make_index(tmp_path, {(201, 300): 50})
        assert ev(idx, HX_CDS, complete=False).integrity == -1.0
        stop = ev(idx, [(101, 200), (300, 400), (500, 604)])  # reads through the TAA
        assert stop.internal_stops == 1 and stop.integrity == -5.0
        frame = ev(idx, [(101, 200), (300, 401), (500, 601)])
        assert not frame.frame_ok and frame.integrity <= -5.0

    def test_priors_bounded_and_transcript_class_counted_once(self, tmp_path):
        idx = make_index(tmp_path, {(201, 300): 50})
        a = ev(idx, HX_CDS, ("Helixer", "StringTie", "Scallop", "Minimap2"))
        assert (a.prior_backbone, a.prior_transcript) == (1.0, 0.25)
        assert a.total == pytest.approx(a.evidence + 1.25)

    def test_single_exon_has_no_intron_terms(self, tmp_path):
        idx = make_index(tmp_path, {(201, 300): 50})
        e = ev(idx, [(101, 599)])
        assert e.introns == [] and e.J == 0 and e.evidence == 0


# ----------------------------------------------------------------- the P4d rule
def cand(name, total, evidence, chain=(), complete=True, backbone=False, neg=0, integrity=0.0,
         J=(), stops=0, frame_ok=True):
    """A CandidateEvidence with given scores; *chain* introns get the given J values."""
    cds = [(0, 10)]
    for i, (a, b) in enumerate(chain):
        cds[-1] = (cds[-1][0], a)
        cds.append((b, b + 10))
    c = CandidateEvidence(key=("1", "+", tuple(cds)), chrom="1", strand="+", cds=tuple(cds),
                          sources=frozenset({name}), complete=complete, total=total,
                          evidence=evidence, n_negative_J=neg, integrity=integrity,
                          has_backbone=backbone, internal_stops=stops, frame_ok=frame_ok)
    c.introns = [{"intron_start": a + 1, "intron_end": b, "J": j} for (a, b), j in zip(chain, J)]
    return c


class TestDecision:
    def test_protected_backbone_without_negative_evidence(self):
        b = cand("hx", 1.0, 0.0, backbone=True)
        d = decide_primary(b, [cand("x", 5.0, 5.0)], PARAMS)
        assert d.outcome == "protected_backbone" and d.chosen is b

    def test_positive_evidence_improvement_replaces(self):
        b = cand("hx", 0.0, -1.0, backbone=True, neg=1)
        x = cand("x", 1.25, 1.0)
        assert decide_primary(b, [x], PARAMS).chosen is x

    def test_prior_alone_cannot_promote(self):
        b = cand("st", 0.25, 0.5)
        x = cand("hx", 1.5, 0.5, backbone=True)           # +1.25 total, 0 evidence
        d = decide_primary(b, [x], PARAMS)
        assert d.outcome == "kept_below_margin" and d.chosen is b

    def test_evidence_below_margin_is_kept(self):
        b = cand("hx", 1.0, 0.0, backbone=True, neg=1)
        assert decide_primary(b, [cand("x", 1.6, 0.4)], PARAMS).outcome == "kept_below_margin"

    def test_negative_evidence_change_never_replaces(self):
        b = cand("hx", 2.0, 1.0, backbone=True, neg=1)
        assert not decide_primary(b, [cand("x", 2.6, 0.9)], PARAMS).replaced

    def test_strongly_supported_intron_cannot_be_dropped(self):
        intron = (10, 50)
        b = cand("hx", 0.0, -1.0, chain=[intron], J=[0.8], backbone=True, neg=1)
        x = cand("x", 3.0, 3.0)                               # drops the J=0.8 intron
        d = decide_primary(b, [x], PARAMS)
        assert d.outcome == "kept_no_eligible" and d.n_excluded_supported_intron == 1
        weak = cand("hx", 0.0, -1.0, chain=[intron], J=[0.4], backbone=True, neg=1)
        assert decide_primary(weak, [x], PARAMS).replaced

    def test_partial_orf_and_internal_stop_ineligible(self):
        b = cand("hx", 0.0, -2.0, backbone=True, neg=2)
        d = decide_primary(b, [cand("p", 5, 5, complete=False), cand("s", 5, 5, stops=1),
                               cand("f", 5, 5, frame_ok=False)], PARAMS)
        assert d.outcome == "kept_no_eligible"
        assert (d.n_excluded_partial_orf, d.n_excluded_stop_or_frame) == (1, 2)

    def test_ties_keep_the_earlier_candidate(self):
        b = cand("hx", 0.0, -1.0, backbone=True, neg=1)
        x1, x2 = cand("a", 2.0, 2.0), cand("b", 2.0, 2.0)
        assert decide_primary(b, [x1, x2], PARAMS).chosen is x1
        assert decide_primary(b, [x2, x1], PARAMS).chosen is x2

    def test_best_by_total_then_evidence(self):
        b = cand("hx", 0.0, -1.0, backbone=True, neg=1)
        x1, x2 = cand("a", 2.0, 1.0), cand("b", 2.0, 1.5)
        assert decide_primary(b, [x1, x2], PARAMS).chosen is x2

    def test_same_coding_locus(self):
        assert same_coding_locus([(0, 100)], [(0, 60)], 0.5)
        assert not same_coding_locus([(0, 100)], [(0, 40)], 0.5)


# ----------------------------------------------------------------- role isolation
class TestRoleIsolation:
    def test_construction_index_accepted(self, tmp_path):
        require_selection_index(make_index(tmp_path, {(201, 300): 5}))

    @pytest.mark.parametrize("roles", [("holdout",), ("construction", "holdout"), ("all",)])
    def test_other_roles_refused(self, tmp_path, roles):
        idx = make_index(tmp_path, {(201, 300): 5}, n_groups=2, roles=roles,
                         lib_roles={1: "holdout"})
        with pytest.raises(EvidenceSelectionError):
            require_selection_index(idx)

    def test_failed_construction_library_refused(self, tmp_path):
        idx = make_index(tmp_path, {(201, 300): 5}, n_groups=2)
        (tmp_path / "L1.sj.tab").write_text("")                # zero-byte = failed run
        idx = JunctionIndex.from_manifest(str(tmp_path / "manifest.tsv"), ("construction",),
                                          None, checksum=False)
        with pytest.raises(EvidenceSelectionError, match="failed to load"):
            require_selection_index(idx)

    def test_scorer_refuses_holdout_index(self, tmp_path):
        idx = make_index(tmp_path, {(201, 300): 5}, n_groups=2, roles=("holdout",),
                         lib_roles={0: "holdout", 1: "holdout"})
        with pytest.raises(EvidenceSelectionError):
            scorer(idx)

    def test_builder_refuses_missing_or_wrong_roles(self, tmp_path):
        from types import SimpleNamespace

        from gmb.pipeline.builder import _require_evidence_selection_inputs

        with pytest.raises(SystemExit, match="requires STAR junction evidence"):
            _require_evidence_selection_inputs(SimpleNamespace(star_junctions=None))
        make_index(tmp_path, {(201, 300): 5}, n_groups=1)
        m = str(tmp_path / "manifest.tsv")
        for roles in ("holdout", "construction,holdout", "all", ""):
            with pytest.raises(SystemExit, match="construction"):
                _require_evidence_selection_inputs(
                    SimpleNamespace(star_junctions=m, star_junction_roles=roles))
        _require_evidence_selection_inputs(
            SimpleNamespace(star_junctions=m, star_junction_roles="construction"))

    def test_api_refuses_before_any_stage(self, tmp_path):
        from gmb.api import _require_evidence_score_inputs

        overlay = tmp_path / "es.yaml"
        overlay.write_text("scoring:\n  primary_selection: evidence_score\n")
        with pytest.raises(ValueError, match="requires star_junctions"):
            _require_evidence_score_inputs([str(overlay)], "fungi_latest_settings", None,
                                           "construction")
        with pytest.raises(ValueError, match="construction"):
            _require_evidence_score_inputs([str(overlay)], "fungi_latest_settings", "m.tsv",
                                           "holdout")
        # Any other policy: no check, as before.
        _require_evidence_score_inputs([], "fungi_latest_settings", None, "holdout")


# ----------------------------------------------------------------- configuration
class TestConfig:
    def test_defaults_are_the_frozen_p4d_values(self):
        e = load_config(None, "fungi_latest_settings").scoring.evidence_score
        assert (e.breadth_saturation_groups, e.unobserved_expressed_weight, e.protein_intron_weight,
                e.noncanonical_penalty, e.partial_orf_penalty, e.internal_stop_penalty,
                e.backbone_prior, e.transcript_class_prior, e.replacement_margin,
                e.keep_supported_intron_J, e.same_coding_locus_fraction,
                e.min_candidate_cds_bp) == (5, 0.5, 0.5, 0.5, 1, 5, 1, 0.25, 0.5, 0.6, 0.5, 150)

    @pytest.mark.parametrize("preset", ["standard", "fungi", "fungi_latest_settings", "apicomplexa"])
    def test_no_preset_enables_it_and_resolved_config_unchanged(self, preset):
        cfg = load_config(None, preset)
        assert cfg.scoring.primary_selection != "evidence_score"
        assert "evidence_score" not in dump_config(cfg)["scoring"]

    def test_section_written_when_active(self, tmp_path):
        o = tmp_path / "o.yaml"
        o.write_text("scoring:\n  primary_selection: evidence_score\n"
                     "  evidence_score:\n    replacement_margin: 0.75\n")
        d = dump_config(load_config(str(o), "fungi_latest_settings"))
        assert d["scoring"]["evidence_score"]["replacement_margin"] == 0.75

    @pytest.mark.parametrize("key,value", [
        ("breadth_saturation_groups", 0), ("replacement_margin", 0), ("replacement_margin", -1),
        ("keep_supported_intron_J", 1.5), ("same_coding_locus_fraction", 0),
        ("backbone_prior", -0.1), ("protein_intron_weight", "high"), ("min_candidate_cds_bp", True),
    ])
    def test_invalid_values_rejected(self, tmp_path, key, value):
        o = tmp_path / "o.yaml"
        o.write_text(f"scoring:\n  evidence_score:\n    {key}: {value!r}\n".replace("'", ""))
        with pytest.raises(ValueError, match=key):
            load_config(str(o), "fungi_latest_settings")

    def test_unknown_key_rejected(self, tmp_path):
        o = tmp_path / "o.yaml"
        o.write_text("scoring:\n  evidence_score:\n    skipped_junction_penalty: 1\n")
        with pytest.raises(ValueError, match="Unknown configuration key"):
            load_config(str(o), "fungi_latest_settings")

    def test_validator_accepts_defaults(self):
        validate_evidence_score(load_config(None, "standard").scoring.evidence_score)


# ----------------------------------------------------------------- integration in select_isoforms
def _exons(tid, source, exons):
    return pd.DataFrame({"Chromosome": "1", "Start": [s for s, _ in exons],
                         "End": [e for _, e in exons], "Strand": "+", "transcript_id": tid,
                         "Source": source, "Feature": "exon"})


def _select(tmp_path, mode="evidence_score", max_iso=3, junctions=None, extra=(), cds=None):
    cfg = load_config(None, "fungi_latest_settings")
    cfg.scoring.primary_selection = mode
    cfg.scoring.max_isoforms_per_locus = max_iso
    idx = make_index(tmp_path, junctions or {(201, 303): 100, (401, 500): 100})
    locus = pd.concat([_exons("hx", "Helixer", HX_EXONS), _exons("st", "StringTie", TX_EXONS),
                       *extra], ignore_index=True)
    cds_map = cds or {"hx": HX_CDS, "st": TX_CDS}
    return select_isoforms(
        locus, cfg, {"hx", "st"}, GENOME, candidate_cds=cds_map,
        junction_support={"hx": (2, 0), "st": (2, 0)}, protein_compatibility={},
        complete_orf_tids=set(cds_map), evidence_scorer=scorer(idx))


class TestSelection:
    def test_contradicted_backbone_replaced_and_kept_as_alternate(self, tmp_path):
        genes = _select(tmp_path)
        assert len(genes) == 1
        g = genes[0]
        assert g[0]["id"] == "st" and g[0]["selection_reason"] == "evidence_score"
        assert [m["id"] for m in g[1:]] == ["hx"] and g[1]["is_primary"] is False
        rec = g[0]["evidence_selection"]
        assert rec["outcome"] == "replaced" and rec["old_candidate"] == "hx"
        assert rec["delta_evidence"] >= 0.5 and rec["representative"] == "retained_isoform"

    def test_candidate_outside_retained_isoforms_is_eligible(self, tmp_path):
        genes = _select(tmp_path, max_iso=1)                # GMB keeps only the backbone
        g = genes[0]
        assert [m["id"] for m in g] == ["st"]                # gene size never grows
        rec = g[0]["evidence_selection"]
        assert "new_primary_was_not_a_retained_isoform" in rec["concerns"]
        assert "max_isoforms_per_locus" in rec["isoforms_dropped"]

    def test_supported_backbone_unchanged(self, tmp_path):
        genes = _select(tmp_path, junctions={(201, 300): 100, (401, 500): 100})
        assert genes[0][0]["id"] == "hx"
        assert genes[0][0]["evidence_selection"]["outcome"] == "protected_backbone"

    def test_other_policies_ignore_the_scorer(self, tmp_path):
        for mode in ("score", "junction_supported"):
            genes = _select(tmp_path / mode, mode=mode)
            assert genes[0][0]["id"] == "hx"
            assert all("evidence_selection" not in m for g in genes for m in g)

    def test_missing_scorer_fails_without_fallback(self):
        cfg = load_config(None, "fungi_latest_settings")
        cfg.scoring.primary_selection = "evidence_score"
        with pytest.raises(ValueError, match="never falls back"):
            select_isoforms(_exons("hx", "Helixer", HX_EXONS), cfg, set(), GENOME,
                            candidate_cds={"hx": HX_CDS}, complete_orf_tids={"hx"})

    def test_deterministic(self, tmp_path):
        a = _select(tmp_path / "a")
        b = _select(tmp_path / "b")
        def strip(gs):
            return [[(m["id"], m.get("evidence_selection", {}).get("outcome")) for m in g]
                    for g in gs]

        assert strip(a) == strip(b)


def test_cross_gene_structure_never_duplicated(tmp_path):
    """A CDS already used by another gene of the locus is not offered to this gene."""
    from gmb.pipeline.scoring import _promote_evidence_supported

    idx = make_index(tmp_path, {(201, 303): 100, (401, 500): 100})
    locus = pd.concat([_exons("hx", "Helixer", HX_EXONS), _exons("st", "StringTie", TX_EXONS)],
                      ignore_index=True)

    def model(tid, src, exons, chain):
        df = locus[locus.transcript_id == tid]
        return {"id": tid, "source": src, "chrom": "1", "strand": "+", "df": df, "start": exons[0][0],
                "end": exons[-1][1], "exon_count": len(exons), "intron_chain": chain, "score": 5.0}

    hx = model("hx", "Helixer", HX_EXONS, "200-300,400-500")
    st = model("st", "StringTie", TX_EXONS, "200-303,400-500")
    genes = [[hx], [st]]                                       # st is another gene's primary
    _promote_evidence_supported(
        genes, locus, [hx, st], {}, [], {"hx": HX_CDS, "st": TX_CDS}, {"hx", "st"}, scorer(idx),
        load_config(None, "fungi_latest_settings").scoring, lambda a, b: True, lambda a, b: True,
        "drop", {}, {})
    assert genes[0][0]["id"] == "hx"
    rec = genes[0][0]["evidence_selection"]
    assert rec["n_excluded_other_gene"] == 1 and rec["outcome"] in ("kept_no_eligible",)
