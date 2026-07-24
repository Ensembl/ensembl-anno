"""Tests for gmb.pipeline.longread_consensus.

Covers: seqname-interleaved splitting, chimera/oversized-span rejection,
intron-chain grouping + min-read-support gating, single-exon collapsing,
and the short-read rescue path for otherwise-unsupported single reads.
"""

import os

from gmb.pipeline.longread_consensus import (
    LongreadConsensusConfig,
    build_consensus_for_seqname,
    split_by_seqname,
)


def _write_gtf(path, lines):
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _transcript_block(chrom, tid, strand, exons):
    """Emit a (deliberately-swapped-coords) transcript line + clean exon lines,
    mirroring the real raw file's ~2.93% swapped-coordinate transcript lines."""
    lines = [
        f'{chrom}\tminimap\ttranscript\t{exons[-1][0]}\t{exons[0][1]}\t.\t{strand}\t.\t'
        f'gene_id "{tid}"; transcript_id "{tid}";'
    ]
    for i, (s, e) in enumerate(exons, start=1):
        lines.append(
            f'{chrom}\tminimap\texon\t{s}\t{e}\t.\t{strand}\t.\t'
            f'gene_id "{tid}"; transcript_id "{tid}"; exon_number "{i}";'
        )
    return lines


class TestSplitBySeqname:
    def test_interleaved_seqnames_split_correctly(self, tmp_path):
        lines = []
        lines += _transcript_block("1", "r1", "+", [(100, 200), (300, 400)])
        lines += _transcript_block("2", "r2", "+", [(500, 600)])
        lines += _transcript_block("1", "r3", "+", [(100, 200), (300, 400)])
        lines += _transcript_block("2", "r4", "+", [(500, 600)])
        raw = tmp_path / "raw.gtf"
        _write_gtf(raw, lines)

        split_dir = tmp_path / "split"
        paths = split_by_seqname(str(raw), str(split_dir))

        assert set(paths) == {"1", "2"}
        with open(paths["1"]) as fh:
            content = fh.read()
        assert "r1" in content and "r3" in content and "r2" not in content

    def test_only_seqname_filters(self, tmp_path):
        lines = []
        lines += _transcript_block("1", "r1", "+", [(100, 200)])
        lines += _transcript_block("2", "r2", "+", [(500, 600)])
        raw = tmp_path / "raw.gtf"
        _write_gtf(raw, lines)

        split_dir = tmp_path / "split"
        paths = split_by_seqname(str(raw), str(split_dir), only_seqname="1")

        assert set(paths) == {"1"}


class TestConsensusBuilding:
    def _cfg(self, **overrides):
        cfg = LongreadConsensusConfig()
        for k, v in overrides.items():
            setattr(cfg, k, v)
        return cfg

    def test_multi_exon_consensus_with_sufficient_support(self, tmp_path):
        # 3 reads sharing the same intron chain -> one consensus model.
        lines = []
        lines += _transcript_block("1", "r1", "+", [(100, 200), (300, 400)])
        lines += _transcript_block("1", "r2", "+", [(95, 200), (300, 405)])
        lines += _transcript_block("1", "r3", "+", [(100, 200), (300, 400)])
        split_path = tmp_path / "1.gtf"
        _write_gtf(split_path, lines)

        cfg = self._cfg(min_read_support_multi_exon=2)
        records, stats = build_consensus_for_seqname("1", str(split_path), cfg)

        assert len(records) == 1
        rec = records[0]
        assert rec["read_support"] == 3
        assert rec["n_exons"] == 2
        # Internal junction (200-300) must be exact; only first start / last end
        # collapse to the median across supporting reads.
        assert rec["exons"][0][1] == 200
        assert rec["exons"][1][0] == 300
        assert stats["input_reads"] == 3
        assert stats["consensus_models_output"] == 1

    def test_single_read_multi_exon_dropped_without_support(self, tmp_path):
        lines = _transcript_block("1", "r1", "+", [(100, 200), (300, 400)])
        split_path = tmp_path / "1.gtf"
        _write_gtf(split_path, lines)

        cfg = self._cfg(
            min_read_support_multi_exon=2, require_shortread_support_for_single_read_models=False
        )
        records, stats = build_consensus_for_seqname("1", str(split_path), cfg)

        assert records == []
        assert stats["dropped_single_read_no_rescue"] == 1

    def test_chimeric_intron_rejected(self, tmp_path):
        # Gap of 100,000bp between exons -- far beyond max_intron_length.
        lines = []
        for i in range(3):
            lines += _transcript_block(
                "1", f"r{i}", "+", [(100, 200), (100_300, 100_400)]
            )
        split_path = tmp_path / "1.gtf"
        _write_gtf(split_path, lines)

        cfg = self._cfg(
            max_intron_length=3000, min_read_support_multi_exon=2, max_span_length=200_000
        )
        records, stats = build_consensus_for_seqname("1", str(split_path), cfg)

        assert records == []
        assert stats["dropped_bad_intron_or_overlap"] == 3
        assert stats["surviving_after_basic_filters"] == 0

    def test_oversized_span_rejected(self, tmp_path):
        lines = []
        for i in range(3):
            lines += _transcript_block("1", f"r{i}", "+", [(100, 200_000)])
        split_path = tmp_path / "1.gtf"
        _write_gtf(split_path, lines)

        cfg = self._cfg(max_span_length=45_000)
        records, stats = build_consensus_for_seqname("1", str(split_path), cfg)

        assert records == []
        assert stats["dropped_oversized_span"] == 3

    def test_single_exon_requires_higher_support(self, tmp_path):
        lines = []
        # Only 2 supporting reads, single-exon threshold is 4 -> dropped.
        for i in range(2):
            lines += _transcript_block("1", f"r{i}", "+", [(1000, 1100)])
        split_path = tmp_path / "1.gtf"
        _write_gtf(split_path, lines)

        cfg = self._cfg(min_read_support_single_exon=4)
        records, stats = build_consensus_for_seqname("1", str(split_path), cfg)

        assert records == []
        assert stats["dropped_low_support"] >= 1

    def test_single_read_rescued_by_shortread_support(self, tmp_path):
        lines = _transcript_block("1", "r1", "+", [(1000, 1100)])
        split_path = tmp_path / "1.gtf"
        _write_gtf(split_path, lines)

        scallop_path = tmp_path / "scallop.gtf"
        _write_gtf(
            scallop_path,
            ['1\tScallop\texon\t1000\t1100\t.\t+\t.\tgene_id "s1"; transcript_id "s1";'],
        )

        cfg = self._cfg(
            min_read_support_single_exon=4,
            require_shortread_support_for_single_read_models=True,
            shortread_overlap_fraction=0.5,
        )
        records, stats = build_consensus_for_seqname(
            "1", str(split_path), cfg, shortread_paths=[str(scallop_path)]
        )

        assert len(records) == 1
        assert records[0]["rescued_by_shortread"] is True
        assert stats["rescued_by_shortread"] == 1

    def test_single_read_not_rescued_without_shortread_overlap(self, tmp_path):
        lines = _transcript_block("1", "r1", "+", [(1000, 1100)])
        split_path = tmp_path / "1.gtf"
        _write_gtf(split_path, lines)

        scallop_path = tmp_path / "scallop.gtf"
        _write_gtf(
            scallop_path,
            # Different locus entirely -- no overlap.
            ['1\tScallop\texon\t50000\t50100\t.\t+\t.\tgene_id "s1"; transcript_id "s1";'],
        )

        cfg = self._cfg(
            min_read_support_single_exon=4,
            require_shortread_support_for_single_read_models=True,
        )
        records, stats = build_consensus_for_seqname(
            "1", str(split_path), cfg, shortread_paths=[str(scallop_path)]
        )

        assert records == []
        assert stats["dropped_single_read_no_rescue"] == 1

    def test_swapped_transcript_line_coords_do_not_break_consensus(self, tmp_path):
        # _transcript_block already emits swapped transcript-line coords by
        # construction (mirrors the real file's ~2.93% malformed lines) --
        # this just asserts the consensus builder ignores them correctly.
        lines = []
        for i in range(2):
            lines += _transcript_block("1", f"r{i}", "-", [(1000, 1100), (1300, 1400)])
        split_path = tmp_path / "1.gtf"
        _write_gtf(split_path, lines)

        cfg = self._cfg(min_read_support_multi_exon=2)
        records, stats = build_consensus_for_seqname("1", str(split_path), cfg)

        assert len(records) == 1
        assert records[0]["exons"][0][0] == 1000
        assert records[0]["exons"][-1][1] == 1400

    def test_gap_splits_cluster(self, tmp_path):
        # Two well-separated single-exon loci on the same strand, each with
        # 4 supporting reads -- should yield 2 independent consensus models,
        # not one merged locus.
        lines = []
        for i in range(4):
            lines += _transcript_block("1", f"a{i}", "+", [(1000, 1100)])
        for i in range(4):
            lines += _transcript_block("1", f"b{i}", "+", [(50000, 50100)])
        split_path = tmp_path / "1.gtf"
        _write_gtf(split_path, lines)

        cfg = self._cfg(min_read_support_single_exon=4, min_intergenic_gap=500)
        records, stats = build_consensus_for_seqname("1", str(split_path), cfg)

        assert len(records) == 2
        assert stats["clusters_found"] == 2

    def test_distant_models_get_distinct_gene_ids(self, tmp_path):
        # Regression test: an earlier version assigned gene_id from the raw
        # read locus-cluster id, which at real coverage depth collapses to
        # ~1 cluster per strand for a whole chromosome -- silently merging
        # every unrelated consensus model into "one gene". Two well-
        # separated multi-exon consensus models (same raw read cluster, if
        # clustering were too coarse) must still end up with different
        # gene_ids based on their own final genomic positions.
        lines = []
        for i in range(3):
            lines += _transcript_block("1", f"a{i}", "+", [(1000, 1100), (1300, 1400)])
        for i in range(3):
            lines += _transcript_block("1", f"b{i}", "+", [(500000, 500100), (500300, 500400)])
        split_path = tmp_path / "1.gtf"
        _write_gtf(split_path, lines)

        cfg = self._cfg(min_read_support_multi_exon=2)
        records, stats = build_consensus_for_seqname("1", str(split_path), cfg)

        assert len(records) == 2
        gene_ids = {r["gene_id"] for r in records}
        assert len(gene_ids) == 2, "distant models must not share a gene_id"
        assert stats["genes_output"] == 2
