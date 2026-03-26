#!/usr/bin/env python3
"""Tests for mega-transcript splitting logic."""

import os
import sys

import pytest
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from config import load_config, PipelineConfig
from evidence_filter import split_mega_transcripts


def _make_config(**overrides):
    """Return a PipelineConfig with transcript_splitting overrides."""
    cfg = load_config()
    for k, v in overrides.items():
        setattr(cfg.transcript_splitting, k, v)
    return cfg


def _make_tx_df(transcripts):
    """Build exon-level DataFrame from transcript specs.

    transcripts : list of (tid, gene_id, chrom, strand, exons)
    exons : list of (start, end)
    """
    rows = []
    for tid, gene_id, chrom, strand, exons in transcripts:
        for s, e in exons:
            rows.append({
                'Chromosome': chrom,
                'Start': s,
                'End': e,
                'Strand': strand,
                'transcript_id': tid,
                'gene_id': gene_id,
                'Source': 'StringTie',
                'Feature': 'exon',
            })
    return pd.DataFrame(rows)


class TestTwoClusterSplit:
    def test_two_cluster_split(self):
        """Transcript with two exon clusters separated by large gap → two segments."""
        cfg = _make_config(split_enabled=True, split_gap_bp=1000)
        df = _make_tx_df([
            ('tx1', 'gene1', '1', '+',
             [(100, 200), (250, 350),          # cluster 1
              (5000, 5100), (5200, 5300)]),     # cluster 2, gap = 5000-350-1 = 4649
        ])
        stats = {}
        result = split_mega_transcripts(df, cfg, stats)

        tids = sorted(result['transcript_id'].unique())
        assert len(tids) == 2, f"Expected 2 segments, got {tids}"
        assert tids[0].endswith('__seg0001')
        assert tids[1].endswith('__seg0002')
        assert stats['transcripts_split'] == 1
        assert stats['segments_emitted'] == 2

        # seg0001 should be the genomically earlier cluster
        seg1 = result[result['transcript_id'] == tids[0]]
        assert seg1['Start'].min() == 100
        seg2 = result[result['transcript_id'] == tids[1]]
        assert seg2['Start'].min() == 5000


class TestMixedContigDeterministicSplit:
    def test_mixed_contig_always_splits(self):
        """Transcript with exons on two chromosomes → always split, never drop."""
        cfg = _make_config(split_enabled=True, split_gap_bp=1000)
        df = _make_tx_df([
            ('tx1', 'gene1', '1', '+', [(100, 200)]),
        ])
        # Add a second chromosome's exon manually
        df2 = _make_tx_df([
            ('tx1', 'gene1', '2', '+', [(500, 600)]),
        ])
        df = pd.concat([df, df2], ignore_index=True)

        stats = {}
        result = split_mega_transcripts(df, cfg, stats)

        tids = sorted(result['transcript_id'].unique())
        assert len(tids) == 2, f"Expected 2 segments (one per contig), got {tids}"
        assert stats['transcripts_split'] == 1


class TestNoSplitWhenDisabled:
    def test_disabled(self):
        """With split_enabled=False, no splitting occurs."""
        cfg = _make_config(split_enabled=False, split_gap_bp=100)
        df = _make_tx_df([
            ('tx1', 'gene1', '1', '+',
             [(100, 200), (5000, 5100)]),
        ])
        result = split_mega_transcripts(df, cfg)
        tids = result['transcript_id'].unique()
        assert len(tids) == 1
        assert tids[0] == 'tx1'


class TestSingleExonPassthrough:
    def test_single_exon(self):
        """Single-exon transcripts pass through unmodified."""
        cfg = _make_config(split_enabled=True, split_gap_bp=100)
        df = _make_tx_df([
            ('tx1', 'gene1', '1', '+', [(100, 500)]),
        ])
        result = split_mega_transcripts(df, cfg)
        tids = result['transcript_id'].unique()
        assert len(tids) == 1
        assert tids[0] == 'tx1'


class TestMaxSegmentsDrops:
    def test_exceeds_max_segments_is_dropped(self):
        """Transcript exceeding max_segments_per_transcript is dropped, not capped."""
        cfg = _make_config(
            split_enabled=True,
            split_gap_bp=100,
            max_segments_per_transcript=2,
        )
        # Create transcript with 4 well-separated exon clusters
        df = _make_tx_df([
            ('tx1', 'gene1', '1', '+',
             [(100, 200), (1000, 1100), (2000, 2100), (3000, 3100)]),
        ])
        stats = {}
        result = split_mega_transcripts(df, cfg, stats)

        # tx1 should be entirely dropped
        assert result.empty or 'tx1' not in result['transcript_id'].values
        assert stats['transcripts_dropped_max_segments'] == 1


class TestGeneIdUpdateAndProvenance:
    def test_gene_id_and_provenance(self):
        """Segment rows have updated gene_id + provenance columns."""
        cfg = _make_config(split_enabled=True, split_gap_bp=1000)
        df = _make_tx_df([
            ('tx1', 'gene1', '1', '+',
             [(100, 200), (5000, 5100)]),
        ])
        stats = {}
        result = split_mega_transcripts(df, cfg, stats)

        # Both transcript_id and gene_id should be segment-specific
        for _, row in result.iterrows():
            assert '__seg' in row['transcript_id']
            assert '__seg' in row['gene_id']

        # Provenance columns exist and are correct
        assert 'parent_transcript_id' in result.columns
        assert 'parent_gene_id' in result.columns
        assert (result['parent_transcript_id'] == 'tx1').all()
        assert (result['parent_gene_id'] == 'gene1').all()

        # gene_id suffix matches transcript_id suffix
        for _, row in result.iterrows():
            tid_suffix = row['transcript_id'].split('__')[-1]
            gid_suffix = row['gene_id'].split('__')[-1]
            assert tid_suffix == gid_suffix


class TestConfigLoadsSplittingSection:
    def test_defaults(self):
        """Splitting config loads from fungi_default.yaml."""
        cfg = load_config()
        assert cfg.transcript_splitting.split_enabled is True
        assert cfg.transcript_splitting.split_gap_bp == 3000
        assert cfg.transcript_splitting.split_on_contig_change is True
        assert cfg.transcript_splitting.split_on_strand_change is True
        assert cfg.transcript_splitting.split_on_large_exon_bp == 50000
        assert cfg.transcript_splitting.max_segments_per_transcript == 50

    def test_yaml_override(self, tmp_path):
        """Splitting config can be overridden via YAML."""
        yaml_file = tmp_path / 'test.yaml'
        yaml_file.write_text(
            'transcript_splitting:\n'
            '  split_enabled: true\n'
            '  split_gap_bp: 5000\n'
            '  max_segments_per_transcript: 10\n'
        )
        cfg = load_config(str(yaml_file))
        assert cfg.transcript_splitting.split_enabled is True
        assert cfg.transcript_splitting.split_gap_bp == 5000
        assert cfg.transcript_splitting.max_segments_per_transcript == 10


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
