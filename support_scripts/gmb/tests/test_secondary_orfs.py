"""Unit tests for opt-in secondary-ORF candidates (scoring.evidence_score.secondary_orfs)."""

import pandas as pd
import pytest

from gmb.pipeline.annotate_cds_utrs import (
    annotate_all_transcripts,
    find_best_orf,
    find_complete_orfs,
    reverse_complement,
    translate,
)
from gmb.pipeline.config import dump_config, load_config, validate_selection_policy
from gmb.pipeline.evidence_features import EvidenceParams, EvidenceScorer
from gmb.pipeline.scoring import _promote_evidence_supported
from gmb.pipeline.secondary_orfs import (
    ORIGIN,
    ProteinCoverage,
    SecondaryOrf,
    admit,
    attribution_rows,
    build_transcript,
    discover,
    validate_cds,
)
from gmb.pipeline.star_junctions import JunctionIndex

F = "GCC"  # filler: no ATG and no stop codon in any frame
INTRON = "GT" + "C" * 56 + "AG"


def orf(n_codons):
    return "ATG" + F * n_codons + "TAA"


# A read-through transcript: own (longest) ORF, then a second, shorter ORF downstream whose
# body is interrupted by an intron. Spliced cDNA = E1 + E2.
ORF1 = orf(80)  # 246 bp
ORF2_A, ORF2_B = "ATG" + F * 30, F * 30 + "TAA"  # 186 bp, above the 150 bp pool floor
E1 = F * 10 + ORF1 + F * 5 + ORF2_A
E2 = ORF2_B + F * 10
LEFT = F * 20
PLUS = LEFT + E1 + INTRON + E2 + F * 20
P_EX = [
    (len(LEFT), len(LEFT) + len(E1)),
    (len(LEFT) + len(E1) + len(INTRON), len(LEFT) + len(E1) + len(INTRON) + len(E2)),
]
MINUS = reverse_complement(PLUS)
L = len(PLUS)
M_EX = sorted((L - b, L - a) for a, b in P_EX)
GENOME = {"1": PLUS, "2": MINUS}
SOURCES = {"StringTie", "Scallop", "Minimap2"}


def exons_df(tid, source, chrom, strand, exons):
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


def world(extra=()):
    rows = [
        exons_df("st+", "StringTie", "1", "+", P_EX),
        exons_df("st-", "StringTie", "2", "-", M_EX),
        *extra,
    ]
    ex = pd.concat(rows, ignore_index=True)
    ann = annotate_all_transcripts(ex, GENOME, None, min_codons=10)
    txs = [
        (t, s, c, st)
        for t, s, c, st in ex.drop_duplicates("transcript_id")[
            ["transcript_id", "Source", "Chromosome", "Strand"]
        ].itertuples(index=False)
    ]
    return ex, ann, txs


def existing_keys(ex, ann):
    meta = ex.drop_duplicates("transcript_id").set_index("transcript_id")
    return {
        (str(meta.at[t, "Chromosome"]), str(meta.at[t, "Strand"]), tuple(a["cds"]))
        for t, a in ann.items()
        if a["cds"]
    }


# ----------------------------------------------------------------- ORF enumeration
class TestFindCompleteOrfs:
    def test_one_maximal_orf_per_stop(self):
        seq = F * 3 + "ATG" + F * 2 + "ATG" + F * 20 + "TAA" + F * 3
        orfs = find_complete_orfs(seq, min_codons=1)
        assert orfs == [(9, 9 + 3 + 6 + 3 + 60 + 3, 0)]  # the nested ATG is not a second ORF

    def test_length_floors(self):
        seq = orf(5) + F * 2 + orf(40)
        assert len(find_complete_orfs(seq, min_codons=1)) == 2
        assert [o[1] - o[0] for o in find_complete_orfs(seq, min_codons=10)] == [126]
        assert find_complete_orfs(seq, min_codons=1, min_bp=127) == []

    def test_open_orf_without_stop_is_not_complete(self):
        assert find_complete_orfs(orf(20)[:-3] + F * 5, min_codons=1) == []

    def test_best_orf_unchanged(self):
        seq = F * 10 + ORF1 + F * 5 + ORF2_A + ORF2_B
        s, e, p5, p3 = find_best_orf(seq, min_codons=10)
        assert (s, e, p5, p3) == (30, 30 + len(ORF1), False, False)


# ----------------------------------------------------------------- discovery
class TestDiscover:
    def test_downstream_spliced_orf_both_strands(self):
        ex, ann, txs = world()
        found, stats = discover(ann, txs, GENOME, SOURCES, 10, 60, existing_keys(ex, ann))
        assert stats["unique_structures"] == 2 and stats["transcripts_scanned"] == 2
        plus = next(c for c in found.values() if c.chrom == "1")
        minus = next(c for c in found.values() if c.chrom == "2")
        # projection: the ORF spans the intron on the genome, on both strands
        a0 = len(LEFT) + len(F * 10 + ORF1 + F * 5)
        assert plus.cds == ((a0, P_EX[0][1]), (P_EX[1][0], P_EX[1][0] + len(ORF2_B)))
        assert minus.cds == tuple(sorted((L - b, L - a) for a, b in plus.cds))
        for c, seq in ((plus, PLUS), (minus, MINUS)):
            assert validate_cds(seq, c.strand, c.cds, ORF2_A + ORF2_B) == (True, "ok")
            assert c.parents[0]["relation"] == "downstream" and c.cid.endswith(".sorf1")
            assert c.n_introns == 1 and c.cds_bp == len(ORF2_A + ORF2_B)
        assert translate(ORF2_A + ORF2_B) == "M" + "A" * 60 + "*"

    def test_duplicates_merged_with_provenance(self):
        twin = exons_df("sc+", "Scallop", "1", "+", P_EX)
        ex, ann, txs = world([twin])
        found, stats = discover(ann, txs, GENOME, SOURCES, 10, 60, existing_keys(ex, ann))
        plus = next(c for c in found.values() if c.chrom == "1")
        assert [p["tid"] for p in plus.parents] == ["st+", "sc+"]
        assert plus.sources == frozenset({"StringTie", "Scallop"}) and plus.cid == "st+.sorf1"
        assert stats["duplicate_emissions_merged"] == 1

    def test_existing_structure_not_generated(self):
        ex, ann, txs = world()
        first, _ = discover(ann, txs, GENOME, SOURCES, 10, 60, set())
        keys = existing_keys(ex, ann) | set(first)
        found, stats = discover(ann, txs, GENOME, SOURCES, 10, 60, keys)
        assert not found and stats["equal_to_existing_candidate"] == 2

    def test_backbone_and_unknown_sources_not_scanned(self):
        hx = exons_df("hx", "Helixer", "1", "+", P_EX)
        ex, ann, txs = world([hx])
        found, _ = discover(ann, [t for t in txs if t[0] == "hx"], GENOME, SOURCES, 10, 60, set())
        assert not found

    def test_overlapping_other_frame_orf_not_generated(self):
        # +1-frame ORF inside the own ORF: ATG at an offset not divisible by 3
        body = F * 5 + "A" + "ATG" + "GCC" * 15 + "TAG" + F * 40
        own = "ATG" + body + "TAA"
        seq = F * 10 + own + F * 10
        g = {"9": seq}
        ex = exons_df("t", "StringTie", "9", "+", [(0, len(seq))])
        ann = annotate_all_transcripts(ex, g, None, min_codons=10)
        found, stats = discover(ann, [("t", "StringTie", "9", "+")], g, SOURCES, 10, 30, set())
        assert not found and stats["overlapping_other_frame_not_generated"] >= 1

    def test_floor_applies(self):
        ex, ann, txs = world()
        found, _ = discover(
            ann, txs, GENOME, SOURCES, 10, len(ORF2_A + ORF2_B) + 3, existing_keys(ex, ann)
        )
        assert not found

    def test_missing_annotation_or_bad_strand_skipped(self):
        ex, ann, txs = world()
        found, stats = discover({}, txs, GENOME, SOURCES, 10, 60, set())
        assert not found and stats["transcripts_scanned"] == 0
        found, _ = discover(ann, [("st+", "StringTie", "1", ".")], GENOME, SOURCES, 10, 60, set())
        assert not found


class TestValidate:
    def test_rejections(self):
        c = ((len(LEFT) + 30, len(LEFT) + 30 + len(ORF1)),)
        assert validate_cds(PLUS, "+", c, ORF1) == (True, "ok")
        assert validate_cds(PLUS, "+", c, ORF1[:-3] + "TAG")[1] == "projection_mismatch"
        stop = "ATG" + "TAA" + F * 5 + "TAA"
        assert validate_cds(stop, "+", ((0, len(stop)),), stop)[1] == "internal_stop"
        noatg = "GCC" + F * 5 + "TAA"
        assert validate_cds(noatg, "+", ((0, len(noatg)),), noatg)[1] == "no_atg"
        short = orf(5)[:-1]
        assert validate_cds(short, "+", ((0, len(short)),), short)[1] == "frame"


# ----------------------------------------------------------------- admission (A_ev)
def make_index(tmp_path, junctions, n_groups=5):
    tmp_path.mkdir(parents=True, exist_ok=True)
    lines = [
        "library_id\tpath\trole\tstudy\tstrain\tcondition\tcount_unit\tmax_intron\tmin_intron"
        "\tnoncanonical_filtered"
    ]
    for g in range(n_groups):
        (tmp_path / f"L{g}.sj.tab").write_text(
            "2\t11\t40\t1\t1\t1\t1\t0\t30\n"
            + "".join(
                f"1\t{s}\t{e}\t1\t1\t1\t{r}\t0\t30\n" for (s, e), r in sorted(junctions.items())
            )
        )
        lines.append(f"L{g}\tL{g}.sj.tab\tconstruction\tS{g}\tX\tc{g}\tmate\t1000\t21\ttrue")
    (tmp_path / "m.tsv").write_text("\n".join(lines) + "\n")
    return JunctionIndex.from_manifest(
        str(tmp_path / "m.tsv"),
        ("construction",),
        {k: len(v) for k, v in GENOME.items()},
        checksum=False,
    )


def _plus_candidate():
    ex, ann, txs = world()
    found, _ = discover(ann, txs, GENOME, SOURCES, 10, 60, existing_keys(ex, ann))
    return {k: c for k, c in found.items() if c.chrom == "1"}, ex, ann


INTRON_1BASED = (P_EX[0][1] + 1, P_EX[1][0])


class TestAdmit:
    def test_protein_clause(self, tmp_path):
        found, ex, ann = _plus_candidate()
        c = next(iter(found.values()))
        prot = exons_df("p1", "OrthoDB", "1", "+", list(c.cds))
        idx = make_index(tmp_path, {})
        stats = admit(found, {"st+": {0}}, {0: {c.key}}, idx, ProteinCoverage(prot), 0.5)
        assert c.status == "admitted" and c.admitted_by == "protein" and c.protein_cover == 1.0
        assert stats["admitted"] == 1 and c.protein_alignments == 1

    def test_junction_clause_and_not_admitted(self, tmp_path):
        found, ex, ann = _plus_candidate()
        c = next(iter(found.values()))
        admit(
            found,
            {"st+": {0}},
            {0: {c.key}},
            make_index(tmp_path / "j", {INTRON_1BASED: 9}),
            ProteinCoverage(None),
            0.5,
        )
        assert c.status == "admitted" and c.admitted_by == "junctions" and c.introns_supported == 1
        admit(
            found,
            {"st+": {0}},
            {0: {c.key}},
            make_index(tmp_path / "n", {}),
            ProteinCoverage(None),
            0.5,
        )
        assert c.status == "not_admitted"

    def test_half_cover_is_enough_but_not_less(self, tmp_path):
        found, _, _ = _plus_candidate()
        c = next(iter(found.values()))
        a, b = c.cds[0]
        idx = make_index(tmp_path, {})
        half = exons_df("p", "OrthoDB", "1", "+", [(a, a + c.cds_bp // 2)])
        admit(found, {"st+": {0}}, {0: {c.key}}, idx, ProteinCoverage(half), 0.5)
        assert c.status == "admitted"
        less = exons_df("p", "OrthoDB", "1", "+", [(a, a + c.cds_bp // 2 - 3)])
        admit(found, {"st+": {0}}, {0: {c.key}}, idx, ProteinCoverage(less), 0.5)
        assert c.status == "not_admitted"
        other = exons_df("p", "OrthoDB", "1", "-", list(c.cds))  # opposite strand never counts
        admit(found, {"st+": {0}}, {0: {c.key}}, idx, ProteinCoverage(other), 0.5)
        assert c.status == "not_admitted"

    def test_no_shared_coding_locus_is_not_evaluated(self, tmp_path):
        found, ex, ann = _plus_candidate()
        c = next(iter(found.values()))
        own = ("1", "+", tuple(ann["st+"]["cds"]))
        admit(
            found, {"st+": {0}}, {0: {own}}, make_index(tmp_path, {}), ProteinCoverage(None), 0.5
        )
        assert c.status == "no_shared_coding_locus" and c.protein_cover is None


# ----------------------------------------------------------------- transcript representation
class TestTranscript:
    def test_exons_are_the_cds_and_no_utr(self):
        found, _, _ = _plus_candidate()
        c = next(iter(found.values()))
        build_transcript(c, GENOME, 10)
        a = c.annotation
        assert sorted(a["exons"]) == list(c.cds) and sorted(a["cds"]) == list(c.cds)
        assert a["five_prime_utr"] == [] and a["three_prime_utr"] == []
        assert a["protein"] == "M" + "A" * 60 and a["frame_ok"] and not a["is_partial_5"]
        m = c.model
        assert m["candidate_origin"] == ORIGIN and m["combined_evidence"] == "StringTie"
        assert m["secondary_orf_parents"] == "st+" and m["has_backbone"] is False
        assert m["intron_chain"] == f"{P_EX[0][1]}-{P_EX[1][0]}"
        assert m["start"] == c.cds[0][0] and m["end"] == c.cds[-1][1]


# ----------------------------------------------------------------- selection
def _scorer(idx):
    return EvidenceScorer(idx, GENOME, None, EvidenceParams(), {"Helixer"}, SOURCES)


def _gene(tid, src, exons, cds):
    df = exons_df(tid, src, "1", "+", exons)
    chain = (
        ",".join(f"{exons[i][1]}-{exons[i + 1][0]}" for i in range(len(exons) - 1))
        or "single-exon"
    )
    return {
        "id": tid,
        "source": src,
        "chrom": "1",
        "strand": "+",
        "df": df,
        "start": exons[0][0],
        "end": exons[-1][1],
        "exon_count": len(exons),
        "intron_chain": chain,
        "score": 5.0,
        "candidate_cds": list(cds),
        "combined_evidence": src,
    }, df


class TestSelection:
    def _run(self, tmp_path, junctions, helixer_cds, others=()):
        found, _, _ = _plus_candidate()
        c = next(iter(found.values()))
        build_transcript(c, GENOME, 10)
        c.status = "admitted"
        hx, hdf = _gene("hx", "Helixer", [(a, b) for a, b in helixer_cds], helixer_cds)
        genes = [[hx], *[[o[0]] for o in others]]
        locus = pd.concat([hdf, *[o[1] for o in others]], ignore_index=True)
        cds_map = {"hx": helixer_cds, **{o[0]["id"]: o[0]["candidate_cds"] for o in others}}
        _promote_evidence_supported(
            genes,
            locus,
            [g[0] for g in genes],
            {},
            [],
            cds_map,
            set(cds_map),
            _scorer(make_index(tmp_path, junctions)),
            load_config(None, "fungi_latest_settings").scoring,
            lambda a, b: True,
            lambda a, b: True,
            "drop",
            {},
            {},
            secondary_orfs=[c],
        )
        return genes, c

    def test_generated_replaces_contradicted_backbone(self, tmp_path):
        found, _, _ = _plus_candidate()
        c = next(iter(found.values()))
        # Helixer: same ORF but an unobserved alternative acceptor 3 bp upstream (shares the
        # donor of the observed junction, so the junction contradicts it)
        hx_cds = [c.cds[0], (c.cds[1][0] - 3, c.cds[1][1])]
        genes, c = self._run(tmp_path, {INTRON_1BASED: 50}, hx_cds)
        g = genes[0]
        assert g[0]["id"] == c.cid and g[0]["candidate_origin"] == ORIGIN
        rec = g[0]["evidence_selection"]
        assert rec["outcome"] == "replaced" and rec["representative"] == "secondary_orf"
        assert rec["new_candidate_origin"] == ORIGIN
        assert c.pool_entries == [dict(c.pool_entries[0], primary="hx", selected=True)]
        rows = attribution_rows({c.key: c})
        assert rows[0]["selected"] and rows[0]["eligible_genes"] == 1

    def test_not_selected_is_never_an_isoform(self, tmp_path):
        found, _, _ = _plus_candidate()
        c = next(iter(found.values()))
        genes, c = self._run(tmp_path, {}, list(c.cds)[:1] + [(c.cds[1][0], c.cds[1][1])])
        assert all(m["id"] != c.cid for g in genes for m in g)

    def test_other_gene_structure_excluded_and_ineligible_locus(self, tmp_path):
        found, _, _ = _plus_candidate()
        c = next(iter(found.values()))
        far = [
            (P_EX[0][0] + 30, P_EX[0][0] + 30 + len(ORF1))
        ]  # the own ORF: no shared coding locus
        genes, c = self._run(tmp_path, {INTRON_1BASED: 50}, far)
        assert genes[0][0]["id"] == "hx" and c.pool_entries == []


# ----------------------------------------------------------------- configuration
class TestConfig:
    def test_default_off_and_not_dumped(self):
        cfg = load_config(None, "fungi_latest_settings")
        assert cfg.scoring.evidence_score.secondary_orfs is False
        cfg.scoring.primary_selection = "evidence_score"
        assert "secondary_orfs" not in dump_config(cfg)["scoring"]["evidence_score"]
        cfg.scoring.evidence_score.secondary_orfs = True
        assert dump_config(cfg)["scoring"]["evidence_score"]["secondary_orfs"] is True

    def test_requires_evidence_score(self):
        cfg = load_config(None, "fungi_latest_settings")
        cfg.scoring.evidence_score.secondary_orfs = True
        with pytest.raises(ValueError, match="requires scoring.primary_selection: evidence_score"):
            validate_selection_policy(cfg)

    def test_must_be_boolean(self):
        cfg = load_config(None, "fungi_latest_settings")
        cfg.scoring.primary_selection = "evidence_score"
        cfg.scoring.evidence_score.secondary_orfs = "yes"
        with pytest.raises(ValueError, match="true or false"):
            validate_selection_policy(cfg)

    @pytest.mark.parametrize(
        "preset", ["standard", "fungi", "fungi_latest_settings", "apicomplexa"]
    )
    def test_no_preset_enables_it(self, preset):
        assert load_config(None, preset).scoring.evidence_score.secondary_orfs is False


def test_secondary_orf_dataclass_properties():
    c = SecondaryOrf(
        cid="x",
        chrom="1",
        strand="+",
        cds=((0, 30), (40, 100)),
        parents=[{"tid": "a", "source": "Scallop"}],
    )
    assert c.key == ("1", "+", ((0, 30), (40, 100))) and c.cds_bp == 90 and c.n_introns == 1
