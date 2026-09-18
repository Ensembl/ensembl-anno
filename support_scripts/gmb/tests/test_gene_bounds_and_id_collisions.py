"""Regression tests for two generic correctness defects found in cross-species testing.

1. Gene records kept coordinates inherited from an earlier candidate set after
   later stages removed transcripts, so a gene could be arbitrarily wider than
   the transcripts it actually contained.
2. Evidence loading grouped exons by ``transcript_id`` alone. Predictors that
   restart gene numbering on every sequence (Tiberius emits ``g1``, ``g2``, ...
   per contig) therefore had unrelated models on different sequences fused into
   one chimeric transcript, joining unrelated coding fragments in unrelated
   reading frames.

Synthetic coordinates and invented source labels throughout, so the tests pass
only if the behaviour is genuinely generic.
"""

from __future__ import annotations

import pandas as pd
import pytest

from gmb.pipeline.gff3_validate import recompute_gene_bounds, validate_gene


def _gene(gid, start, end, chrom="ctg1", strand="+"):
    return {
        "Chromosome": chrom, "Source": "GMB", "Feature": "gene", "Start": start,
        "End": end, "Score": ".", "Strand": strand, "Frame": ".", "ID": gid, "Parent": "",
    }


def _mrna(tid, gid, start, end, chrom="ctg1", strand="+"):
    return {
        "Chromosome": chrom, "Source": "GMB", "Feature": "mRNA", "Start": start,
        "End": end, "Score": ".", "Strand": strand, "Frame": ".", "ID": tid, "Parent": gid,
    }


def _exon(tid, start, end, n=1, chrom="ctg1", strand="+"):
    return {
        "Chromosome": chrom, "Source": "GMB", "Feature": "exon", "Start": start,
        "End": end, "Score": ".", "Strand": strand, "Frame": ".",
        "ID": f"{tid}.exon{n}", "Parent": tid,
    }


class TestGeneBoundsContract:
    def test_gene_contracts_after_outer_transcript_removed(self):
        """Selection/validation drops the outer isoform: the gene must shrink."""
        rows = [
            _gene("G1", 100, 90_000),          # span still reflects the dropped t2
            _mrna("G1.t1", "G1", 80_000, 90_000),
            _exon("G1.t1", 80_000, 90_000),
            _mrna("G1.t3", "G1", 82_000, 88_000),
            _exon("G1.t3", 82_000, 88_000),
        ]
        out, stats = recompute_gene_bounds(rows)
        gene = next(r for r in out if r["Feature"] == "gene")
        assert (gene["Start"], gene["End"]) == (80_000, 90_000)
        assert stats["genes_adjusted"] == 1
        assert stats["genes_contracted"] == 1
        assert stats["max_contraction_bp"] == 79_900

    def test_gene_widens_when_dedup_reparents_an_outlying_isoform(self):
        """dedup merges another gene's mRNA in as an isoform: the gene must grow."""
        rows = [
            _gene("G1", 5_000, 6_000),
            _mrna("G1.t1", "G1", 5_000, 6_000),
            _exon("G1.t1", 5_000, 6_000),
            _mrna("G2.t1", "G1", 4_000, 7_500),   # absorbed from a merged gene
            _exon("G2.t1", 4_000, 7_500),
        ]
        out, stats = recompute_gene_bounds(rows)
        gene = next(r for r in out if r["Feature"] == "gene")
        assert (gene["Start"], gene["End"]) == (4_000, 7_500)
        assert stats["genes_widened"] == 1

    def test_correct_gene_is_left_untouched(self):
        rows = [
            _gene("G1", 1_000, 2_000),
            _mrna("G1.t1", "G1", 1_000, 2_000),
            _exon("G1.t1", 1_000, 2_000),
        ]
        out, stats = recompute_gene_bounds(rows)
        gene = next(r for r in out if r["Feature"] == "gene")
        assert (gene["Start"], gene["End"]) == (1_000, 2_000)
        assert stats["genes_adjusted"] == 0
        assert len(out) == len(rows)

    def test_gene_with_no_surviving_transcript_is_dropped_with_its_descendants(self):
        rows = [
            _gene("G1", 100, 900),
            _gene("G2", 1_000, 2_000),
            _mrna("G2.t1", "G2", 1_000, 2_000),
            _exon("G2.t1", 1_000, 2_000),
        ]
        out, stats = recompute_gene_bounds(rows)
        assert stats["genes_dropped_no_transcript"] == 1
        assert {r["ID"] for r in out if r["Feature"] == "gene"} == {"G2"}

    def test_multiple_genes_are_each_scoped_to_their_own_children(self):
        rows = [
            _gene("G1", 0, 50_000),
            _mrna("G1.t1", "G1", 1_000, 2_000),
            _exon("G1.t1", 1_000, 2_000),
            _gene("G2", 60_000, 61_000),
            _mrna("G2.t1", "G2", 60_000, 61_000),
            _exon("G2.t1", 60_000, 61_000),
        ]
        out, _ = recompute_gene_bounds(rows)
        spans = {r["ID"]: (r["Start"], r["End"]) for r in out if r["Feature"] == "gene"}
        assert spans == {"G1": (1_000, 2_000), "G2": (60_000, 61_000)}

    def test_invariant_holds_for_every_gene_after_recompute(self):
        rows = [
            _gene("G1", 0, 99_999),
            _mrna("G1.t1", "G1", 10, 20),
            _exon("G1.t1", 10, 20),
            _gene("G2", 500, 600),
            _mrna("G2.t1", "G2", 400, 900),
            _exon("G2.t1", 400, 900),
        ]
        out, _ = recompute_gene_bounds(rows)
        by_parent = {}
        for r in out:
            if r["Feature"] == "mRNA":
                by_parent.setdefault(r["Parent"], []).append(r)
        for g in [r for r in out if r["Feature"] == "gene"]:
            kids = by_parent[g["ID"]]
            assert g["Start"] == min(k["Start"] for k in kids)
            assert g["End"] == max(k["End"] for k in kids)


class TestValidateGeneEquality:
    """validate_gene must reject an over-wide gene, not only an under-covering one."""

    def test_flags_gene_wider_than_children(self):
        gene = _gene("G1", 100, 90_000)
        mrnas = [_mrna("G1.t1", "G1", 80_000, 90_000)]
        assert validate_gene(gene, mrnas)

    def test_flags_gene_narrower_than_children(self):
        gene = _gene("G1", 80_000, 85_000)
        mrnas = [_mrna("G1.t1", "G1", 80_000, 90_000)]
        assert validate_gene(gene, mrnas)

    def test_accepts_exact_match(self):
        gene = _gene("G1", 80_000, 90_000)
        mrnas = [_mrna("G1.t1", "G1", 80_000, 90_000)]
        assert validate_gene(gene, mrnas) == []


class TestCrossSeqidIdCollisions:
    """Reused per-sequence IDs must not fuse unrelated models into one transcript."""

    @staticmethod
    def _write_gtf(tmp_path, rows):
        path = tmp_path / "predictor.gtf"
        with open(path, "w") as fh:
            for chrom, feat, start, end, strand, gid, tid in rows:
                fh.write(
                    f'{chrom}\tPredictorX\t{feat}\t{start}\t{end}\t.\t{strand}\t0\t'
                    f'gene_id "{gid}"; transcript_id "{tid}";\n'
                )
        return str(path)

    def test_same_id_on_two_sequences_stays_two_transcripts(self, tmp_path):
        from gmb.pipeline.builder import load_evidence

        # "g2.t1" names unrelated models on ctgA and ctgB -- the pattern that
        # produced 40-96 kb pseudo-introns and premature stops.
        gtf = self._write_gtf(tmp_path, [
            ("ctgA", "exon", 1_000, 1_500, "+", "g2", "g2.t1"),
            ("ctgB", "exon", 90_000, 92_000, "-", "g2", "g2.t1"),
        ])
        exons, _ = load_evidence(gtf, "PredictorX")
        assert exons["transcript_id"].nunique() == 2
        per_tid_seqs = exons.groupby("transcript_id", observed=True)["Chromosome"].nunique()
        assert (per_tid_seqs == 1).all()

    def test_unique_ids_are_left_unchanged(self, tmp_path):
        from gmb.pipeline.builder import load_evidence

        gtf = self._write_gtf(tmp_path, [
            ("ctgA", "exon", 1_000, 1_500, "+", "g1", "g1.t1"),
            ("ctgB", "exon", 90_000, 92_000, "-", "g9", "g9.t1"),
        ])
        exons, _ = load_evidence(gtf, "PredictorX")
        assert set(exons["transcript_id"]) == {"PredictorX_g1.t1", "PredictorX_g9.t1"}

    def test_multi_exon_model_on_one_sequence_is_not_split(self, tmp_path):
        from gmb.pipeline.builder import load_evidence

        gtf = self._write_gtf(tmp_path, [
            ("ctgA", "exon", 1_000, 1_500, "+", "g1", "g1.t1"),
            ("ctgA", "exon", 2_000, 2_500, "+", "g1", "g1.t1"),
        ])
        exons, _ = load_evidence(gtf, "PredictorX")
        assert exons["transcript_id"].nunique() == 1
        assert len(exons) == 2


class TestTranslationAcrossJunctions:
    """CDS assembly/translation must be frame-correct on both strands.

    The T. gondii internal-stop models were frame-correct per block but joined
    unrelated fragments; these tests pin the surrounding machinery so a real
    phase or ordering regression would be caught here rather than by whole-genome QC.
    """

    def test_phases_are_consistent_with_cumulative_length_plus_strand(self):
        from gmb.pipeline.builder import compute_cds_phases

        cds = [(0, 100), (200, 311), (400, 448)]
        phases = compute_cds_phases(cds, "+")
        cum = 0
        for (s, e), ph in zip(cds, phases):
            assert ph == (3 - (cum % 3)) % 3
            cum += e - s

    def test_phases_are_consistent_with_cumulative_length_minus_strand(self):
        from gmb.pipeline.builder import compute_cds_phases

        cds = [(0, 100), (200, 311), (400, 448)]
        phases = compute_cds_phases(cds, "-")
        # biological order on the minus strand is descending genomic order
        cum = 0
        for (s, e), ph in zip(reversed(cds), reversed(phases)):
            assert ph == (3 - (cum % 3)) % 3
            cum += e - s

    def test_minus_strand_cds_translates_in_biological_order(self):
        from gmb.pipeline.annotate_cds_utrs import build_spliced_seq, translate

        # Reverse complement of ATG AAA TTT TGA laid out across two exons.
        chrom_seq = "AAAA" + "TCAAAATTTCAT" + "AAAA"
        exons = [(4, 10), (10, 16)]
        spliced = build_spliced_seq(exons, "-", chrom_seq)
        prot = translate(spliced)
        assert prot.startswith("MKF")
        assert prot.rstrip("*").count("*") == 0

    def test_translation_detects_a_genuine_internal_stop(self):
        from gmb.pipeline.annotate_cds_utrs import translate

        assert translate("ATGTGAAAATAA").count("*") == 2
