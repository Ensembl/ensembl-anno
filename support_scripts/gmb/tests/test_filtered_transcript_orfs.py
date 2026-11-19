"""Unit tests for the experimental filtered-transcript ORF units
(scoring.evidence_score.filtered_transcript_orfs, GMB Phase 4 proof of concept)."""

import pandas as pd
import pytest

import gmb.pipeline.filtered_transcript_orfs as FT
from gmb.pipeline.annotate_cds_utrs import annotate_all_transcripts
from gmb.pipeline.config import dump_config, load_config
from gmb.pipeline.secondary_orfs import SecondaryOrf, discover

F = "GCC"  # filler: no ATG and no stop codon in any frame
INTRON = "GT" + "C" * 56 + "AG"


def orf(n_codons):
    return "ATG" + F * n_codons + "TAA"


ORF1 = orf(80)  # 246 bp: the removed model's own GMB ORF
ORF2_A, ORF2_B = "ATG" + F * 30, F * 30 + "TAA"  # 186 bp spliced downstream unit
LEFT = F * 20
E1 = F * 10 + ORF1 + F * 5 + ORF2_A
E2 = ORF2_B + F * 10
SEQ = LEFT + E1 + INTRON + E2 + F * 20
EX = [
    (len(LEFT), len(LEFT) + len(E1)),
    (len(LEFT) + len(E1) + len(INTRON), len(LEFT) + len(E1) + len(INTRON) + len(E2)),
]
I0, I1 = EX[0][1], EX[1][0]
ORF1_CDS = ((len(LEFT) + 30, len(LEFT) + 30 + len(ORF1)),)
ORF2_CDS = ((I0 - len(ORF2_A), I0), (I1, I1 + len(ORF2_B)))
GENOME = {"1": SEQ}
SOURCES = {"StringTie", "Scallop", "Minimap2"}


def exons_df(tid, exons, source="StringTie", chrom="1", strand="+"):
    return pd.DataFrame(
        {
            "Chromosome": chrom,
            "Strand": strand,
            "Start": [a for a, _ in exons],
            "End": [b for _, b in exons],
            "transcript_id": tid,
            "Source": source,
        }
    )


def _discover(df, existing=(), max_intron=3000, max_span=20000):
    ann = annotate_all_transcripts(df, GENOME, None, min_codons=33)
    return FT.discover_removed(
        df, ann, GENOME, SOURCES, 33, 150, set(existing), max_intron, max_span, discover
    )


def test_removed_transcripts_are_the_ones_the_filter_dropped():
    all_ = pd.concat([exons_df("a", EX), exons_df("b", EX)])
    kept = exons_df("a", EX)
    out = FT.removed_transcripts(all_, kept)
    assert set(out["transcript_id"]) == {"b"}
    assert FT.removed_transcripts(all_, kept.iloc[0:0])["transcript_id"].nunique() == 2


def test_units_are_the_own_orf_and_the_flanking_orf():
    found, stats = _discover(exons_df("RT", EX), max_span=100)
    by = {c.cds: c for c in found.values()}
    assert set(by) == {ORF1_CDS, ORF2_CDS}
    own, down = by[ORF1_CDS], by[ORF2_CDS]
    assert own.cid == "RT.ftown1" and down.cid == "RT.ftorf1"
    assert [p["relation"] for p in own.parents] == ["removed_own"]
    assert [p["relation"] for p in down.parents] == ["removed_downstream"]
    assert all(p["filter_reason"] == "span" for c in found.values() for p in c.parents)
    assert stats["removed_models"] == 1 and stats["units"] == 2


def test_unit_equal_to_an_existing_structure_is_not_generated():
    found, stats = _discover(exons_df("RT", EX), existing={("1", "+", ORF2_CDS)})
    assert {c.cds for c in found.values()} == {ORF1_CDS}
    assert stats["flanking_equal_to_existing_candidate"] == 1


def test_unit_must_itself_pass_the_filter():
    found, stats = _discover(exons_df("RT", EX), max_intron=10)
    assert found[("1", "+", ORF2_CDS)].status == "exceeds_filter"
    assert found[("1", "+", ORF1_CDS)].status == "discovered"
    assert stats["units_exceeding_filter"] == 1


class _NoJ:
    def find(self, *a):
        return None


class _Cov:
    def __init__(self, cover):
        self.c = cover

    def cover(self, chrom, strand, cds):
        return self.c, 1


def _unit(cds, tid="RT", status="discovered"):
    c = SecondaryOrf(cid=f"{tid}.ftorf1", chrom="1", strand="+", cds=tuple(cds), status=status)
    c.parents.append({"tid": tid, "source": "StringTie", "relation": "removed_downstream"})
    return c


def test_admission_uses_the_frozen_rule_after_locus_assignment():
    near = ((I0 - len(ORF2_A), I0), (I1 - 3, I1 + len(ORF2_B)))  # an existing structure
    structures = {7: {("1", "+", near)}}
    ok = _unit(ORF2_CDS)
    far = _unit(((5000, 5300),), tid="RT2")
    found = {ok.key: ok, far.key: far}
    st = FT.admit_removed(found, structures, index=_NoJ(), coverage=_Cov(0.9), fraction=0.5)
    assert ok.status == "admitted" and ok.loci == {7} and ok.admitted_by == "protein"
    assert far.status == "no_shared_coding_locus"
    assert st["admitted"] == 1 and st["no_shared_coding_locus"] == 1


def test_no_protein_and_no_junctions_is_not_admitted():
    structures = {7: {("1", "+", ORF2_CDS[:1] + ((I1, I1 + len(ORF2_B) + 3),))}}
    c = _unit(ORF2_CDS)
    FT.admit_removed({c.key: c}, structures, index=_NoJ(), coverage=_Cov(0.1), fraction=0.5)
    assert c.status == "not_admitted"


def test_a_unit_sharing_coding_loci_of_two_loci_is_rejected_as_a_fusion():
    structures = {1: {("1", "+", ORF2_CDS)}, 2: {("1", "+", ORF2_CDS[:1] + ((I1, I1 + 90),))}}
    c = _unit(((ORF2_CDS[0][0], ORF2_CDS[0][1]), (I1, I1 + len(ORF2_B) - 3)))
    FT.admit_removed({c.key: c}, structures, index=None, coverage=_Cov(0.9), fraction=0.5)
    assert c.status == "multi_locus" and c.loci == {1, 2}


def test_exceeds_filter_is_never_admitted():
    c = _unit(ORF2_CDS, status="exceeds_filter")
    st = FT.admit_removed({c.key: c}, {7: {("1", "+", ORF2_CDS)}}, None, _Cov(1.0), 0.5)
    assert c.status == "exceeds_filter" and st["exceeds_filter"] == 1


def test_per_locus_bound_keeps_the_best_supported_units(monkeypatch):
    monkeypatch.setattr(FT, "MAX_ADMITTED_PER_LOCUS", 2)
    base = (1000, 1300)
    structures = {3: {("1", "+", (base,))}}
    units = []
    for i in range(4):
        u = SecondaryOrf(cid=f"u{i}", chrom="1", strand="+", cds=((1000, 1300 - 3 * i),))
        for k in range(i + 1):  # u3 has the most parents
            u.parents.append(
                {"tid": f"t{i}_{k}", "source": "StringTie", "relation": "removed_own"}
            )
        units.append(u)
    found = {u.key: u for u in units}
    st = FT.admit_removed(found, structures, None, _Cov(0.9), 0.5)
    assert [u.status for u in units] == ["locus_cap", "locus_cap", "admitted", "admitted"]
    assert st["admitted"] == 2 and st["locus_cap"] == 2 and st["loci_at_cap"] == 1


def _cfg(tmp_path, body):
    p = tmp_path / "c.yaml"
    p.write_text(body)
    return load_config(str(p), preset="fungi_latest_settings")


def test_config_requires_secondary_orfs(tmp_path):
    with pytest.raises(ValueError, match="requires scoring.evidence_score.secondary_orfs"):
        _cfg(
            tmp_path,
            "scoring:\n  primary_selection: evidence_score\n  evidence_score:\n"
            "    filtered_transcript_orfs: true\n",
        )


def test_config_rejects_non_bool(tmp_path):
    with pytest.raises(ValueError, match="filtered_transcript_orfs must be true or false"):
        _cfg(
            tmp_path,
            "scoring:\n  primary_selection: evidence_score\n  evidence_score:\n"
            "    secondary_orfs: true\n    filtered_transcript_orfs: 1\n",
        )


def test_resolved_config_omits_the_switch_when_off(tmp_path):
    on = dump_config(
        _cfg(
            tmp_path,
            "scoring:\n  primary_selection: evidence_score\n"
            "  evidence_score:\n    secondary_orfs: true\n"
            "    filtered_transcript_orfs: true\n",
        )
    )
    off = dump_config(
        _cfg(
            tmp_path,
            "scoring:\n  primary_selection: evidence_score\n"
            "  evidence_score:\n    secondary_orfs: true\n"
            "    filtered_transcript_orfs: false\n",
        )
    )
    plain = dump_config(
        _cfg(
            tmp_path,
            "scoring:\n  primary_selection: evidence_score\n"
            "  evidence_score:\n    secondary_orfs: true\n",
        )
    )
    assert on["scoring"]["evidence_score"]["filtered_transcript_orfs"] is True
    assert off == plain and "filtered_transcript_orfs" not in off["scoring"]["evidence_score"]
