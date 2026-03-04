#!/usr/bin/env python3
"""Tests for scoring module."""

import os
import sys

import pytest
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from config import load_config
from scoring import score_model, select_isoforms


@pytest.fixture
def config():
    return load_config()


def _make_model(tid, source, chrom='1', strand='+',
                exons=None, combined_evidence=None):
    """Create a model dict for testing."""
    if exons is None:
        exons = [(100, 500)]
    df = pd.DataFrame({
        'Start': [s for s, e in exons],
        'End': [e for s, e in exons],
        'Chromosome': [chrom] * len(exons),
        'Strand': [strand] * len(exons),
    })
    return {
        'id': tid,
        'source': source,
        'chrom': chrom,
        'strand': strand,
        'intron_chain': 'single-exon' if len(exons) == 1 else 'multi',
        'protein_support': False,
        'df': df,
        'start': exons[0][0],
        'end': exons[-1][1],
        'exon_count': len(exons),
        'combined_evidence': combined_evidence or source,
    }


def _make_locus_df(models):
    """Combine model dicts into a locus DataFrame."""
    dfs = []
    for m in models:
        df = m['df'].copy()
        df['Source'] = m['source']
        df['transcript_id'] = m['id']
        df['combined_evidence'] = m['combined_evidence']
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


class TestScoreModel:
    def test_helixer_scores_higher(self, config):
        """Helixer model should score higher than single-source Scallop."""
        helixer = _make_model('hx1', 'Helixer',
                              combined_evidence='Helixer')
        scallop = _make_model('sc1', 'Scallop',
                              combined_evidence='Scallop')
        h_score = score_model(helixer, config, set())
        s_score = score_model(scallop, config, set())
        assert h_score > s_score

    def test_protein_support_bonus(self, config):
        """Model with protein support scores higher."""
        m1 = _make_model('tx1', 'Scallop')
        s_no_prot = score_model(m1, config, set())
        s_with_prot = score_model(m1, config, {'tx1'})
        assert s_with_prot > s_no_prot
        assert s_with_prot - s_no_prot == config.scoring.protein_overlap_bonus

    def test_multi_source_bonus(self, config):
        """Multi-source model scores higher than single-source."""
        single = _make_model('tx1', 'Scallop',
                             combined_evidence='Scallop')
        multi = _make_model('tx2', 'Scallop',
                            combined_evidence='Scallop,StringTie')
        s1 = score_model(single, config, set())
        s2 = score_model(multi, config, set())
        assert s2 > s1


class TestSelectIsoforms:
    def test_max_isoforms_limit(self, config):
        """Only top-N isoforms should be selected per locus."""
        # Create models with different structures
        models = []
        for i in range(5):
            models.append(_make_model(
                f'tx{i}', 'Helixer',
                exons=[(100 + i * 10, 500 + i * 10)],
                combined_evidence='Helixer'))
        locus_df = _make_locus_df(models)
        selected = select_isoforms(locus_df, config, set())
        assert len(selected) <= config.scoring.max_isoforms_per_locus

    def test_protein_supported_preferred(self, config):
        """Protein-supported model should be selected as primary."""
        m1 = _make_model('tx1', 'Scallop',
                         combined_evidence='Scallop')
        m2 = _make_model('tx2', 'Scallop',
                         exons=[(100, 500)],
                         combined_evidence='Scallop')
        locus_df = _make_locus_df([m1, m2])
        selected = select_isoforms(locus_df, config, {'tx1'})
        # tx1 should be selected (protein support)
        selected_ids = {m['id'] for m in selected}
        assert 'tx1' in selected_ids

    def test_empty_locus(self, config):
        """Empty locus → empty result."""
        locus_df = pd.DataFrame(columns=[
            'Chromosome', 'Start', 'End', 'Strand',
            'Source', 'transcript_id', 'combined_evidence'])
        selected = select_isoforms(locus_df, config, set())
        assert selected == []


class TestLowMinCodons:
    def test_short_orf_accepted_with_fungal_config(self):
        """ORFs of 40 codons should be accepted with fungal min_codons=33."""
        from annotate_cds_utrs import find_best_orf
        # 40-codon ORF: ATG + 38 coding + stop = 40 codons
        seq = 'ATG' + 'GCT' * 38 + 'TAA'  # 40 codons
        result = find_best_orf(seq, min_codons=33)
        assert result is not None
        assert (result[1] - result[0]) // 3 >= 33

    def test_very_short_orf_still_found(self):
        """Even 20-codon ORFs should be found (as fallback)."""
        from annotate_cds_utrs import find_best_orf
        seq = 'ATG' + 'GCT' * 18 + 'TAA'  # 20 codons
        result = find_best_orf(seq, min_codons=33)
        assert result is not None  # fallback should still return it


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
