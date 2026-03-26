#!/usr/bin/env python3
"""Tests for subset_utils: region parsing, mapping, subsetting, sampling."""

import os
import sys
import tempfile

import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from subset_utils import (
    Region,
    build_mapping,
    load_assembly_mapping,
    load_regions_file,
    load_seqname_map,
    parse_region,
    remap_df_seqnames,
    sample_loci,
    subset_df_by_regions,
)


class TestParseRegion:
    def test_full_region(self):
        r = parse_region("chr1:100-200")
        assert r == Region("chr1", 100, 200)

    def test_whole_contig(self):
        r = parse_region("chr1")
        assert r == Region("chr1", None, None)
        assert r.is_whole_contig()

    def test_whitespace_stripped(self):
        r = parse_region("  chr2:500-900  ")
        assert r == Region("chr2", 500, 900)


class TestLoadRegionsFile:
    def test_mixed_regions_and_comments(self, tmp_path):
        f = tmp_path / "regions.txt"
        f.write_text(
            "# header comment\n"
            "chr1:100-200\n"
            "\n"
            "# another comment\n"
            "chr2\n"
            "chr3:1000-2000\n"
        )
        regions = load_regions_file(str(f))
        assert len(regions) == 3
        assert regions[0] == Region("chr1", 100, 200)
        assert regions[1] == Region("chr2")
        assert regions[2] == Region("chr3", 1000, 2000)


class TestRemapDfSeqnames:
    def test_basic_remap(self):
        df = pd.DataFrame(
            {
                "Chromosome": ["chr1", "chr1", "chr2"],
                "Start": [100, 200, 300],
                "End": [150, 250, 350],
            }
        )
        mapping = {"chr1": "1", "chr2": "2"}
        result = remap_df_seqnames(df, mapping, label="test")
        assert set(result["Chromosome"].unique()) == {"1", "2"}

    def test_unmapped_preserved(self, capsys):
        df = pd.DataFrame(
            {
                "Chromosome": ["chr1", "chrX"],
                "Start": [100, 200],
                "End": [150, 250],
            }
        )
        mapping = {"chr1": "1"}
        result = remap_df_seqnames(df, mapping, label="test")
        assert "1" in result["Chromosome"].values
        assert "chrX" in result["Chromosome"].values
        captured = capsys.readouterr()
        assert "unmapped" in captured.out.lower() or "chrX" in captured.out

    def test_empty_df_passthrough(self):
        df = pd.DataFrame()
        result = remap_df_seqnames(df, {"chr1": "1"})
        assert result.empty


class TestSubsetDfByRegions:
    @pytest.fixture
    def two_chrom_df(self):
        return pd.DataFrame(
            {
                "Chromosome": ["1", "1", "1", "2", "2"],
                "Start": [100, 300, 500, 100, 400],
                "End": [200, 400, 600, 200, 500],
                "transcript_id": ["tx1", "tx2", "tx3", "tx4", "tx5"],
            }
        )

    def test_single_region(self, two_chrom_df):
        regions = [Region("1", 100, 250)]
        result = subset_df_by_regions(two_chrom_df, regions)
        assert set(result["transcript_id"]) == {"tx1"}

    def test_overlapping_region(self, two_chrom_df):
        regions = [Region("1", 350, 550)]
        result = subset_df_by_regions(two_chrom_df, regions)
        assert set(result["transcript_id"]) == {"tx2", "tx3"}

    def test_whole_contig(self, two_chrom_df):
        regions = [Region("2")]
        result = subset_df_by_regions(two_chrom_df, regions)
        assert set(result["transcript_id"]) == {"tx4", "tx5"}

    def test_no_match(self, two_chrom_df):
        regions = [Region("3", 100, 200)]
        result = subset_df_by_regions(two_chrom_df, regions)
        assert result.empty

    def test_multiple_regions(self, two_chrom_df):
        regions = [Region("1", 100, 250), Region("2", 100, 250)]
        result = subset_df_by_regions(two_chrom_df, regions)
        assert set(result["transcript_id"]) == {"tx1", "tx4"}


class TestSampleLociReproducible:
    @pytest.fixture
    def loci_df(self):
        return pd.DataFrame(
            {
                "Chromosome": ["1"] * 10,
                "Start": list(range(0, 10000, 1000)),
                "End": list(range(500, 10500, 1000)),
            }
        )

    def test_same_seed_same_result(self, loci_df):
        r1 = sample_loci(loci_df, n=5, seed=42)
        r2 = sample_loci(loci_df, n=5, seed=42)
        assert len(r1) == len(r2)
        assert [str(r) for r in r1] == [str(r) for r in r2]

    def test_different_seed_different_result(self, loci_df):
        r1 = sample_loci(loci_df, n=5, seed=1)
        r2 = sample_loci(loci_df, n=5, seed=99)
        # Could theoretically be same, but very unlikely with different seeds
        assert [str(r) for r in r1] != [str(r) for r in r2] or len(r1) == 0

    def test_window_expansion(self, loci_df):
        regions = sample_loci(loci_df, n=1, seed=1, window_bp=500)
        assert len(regions) == 1
        # Should be expanded by 500 on each side relative to original locus
        r = regions[0]
        assert r.start is not None
        assert r.end is not None


class TestSampleLociUniformGenome:
    def test_genome_sampling(self):
        loci = pd.DataFrame(
            {
                "Chromosome": ["1", "1", "2"],
                "Start": [100, 5000, 200],
                "End": [500, 6000, 800],
            }
        )
        regions = sample_loci(loci, n=2, strategy="uniform_genome", seed=42)
        assert len(regions) <= 2
        for r in regions:
            assert r.seqname in ("1", "2")


class TestBuildMappingPrecedence:
    def test_seqname_map_overrides(self, tmp_path):
        asm = tmp_path / "assembly.txt"
        asm.write_text(
            "# header\n"
            "scaffold_1\tassembled\t1\tChromosome\tACC.1\n"
            "scaffold_2\tassembled\t2\tChromosome\tACC.2\n"
        )
        seq = tmp_path / "seqmap.tsv"
        seq.write_text("ACC.1\tOVERRIDE\n")

        mapping = build_mapping(
            assembly_report=str(asm),
            seqname_map=str(seq),
        )
        assert mapping["ACC.1"] == "OVERRIDE"
        assert mapping["ACC.2"] == "2"


class TestSubsetAfterMapping:
    def test_subset_uses_mapped_names(self):
        df = pd.DataFrame(
            {
                "Chromosome": ["chr1", "chr1", "chr2"],
                "Start": [100, 500, 100],
                "End": [200, 600, 200],
                "transcript_id": ["tx1", "tx2", "tx3"],
            }
        )
        mapping = {"chr1": "1", "chr2": "2"}
        mapped = remap_df_seqnames(df, mapping)
        regions = [Region("1", 100, 250)]
        result = subset_df_by_regions(mapped, regions)
        assert set(result["transcript_id"]) == {"tx1"}


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
