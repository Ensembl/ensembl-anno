#!/usr/bin/env python3
"""Isoform scoring and selection for Gene Model Builder.

Replaces the previous hard-coded analyze_locus() with a configurable
scoring function.  Default parameters are tuned for fungal genomes.

All thresholds are driven by PipelineConfig (YAML).
"""

from collections import defaultdict

from annotate_cds_utrs import (
    build_spliced_seq,
    check_splice_sites,
)


# ---------------------------------------------------------------------------
# Intron-chain utility (reused from gene_model_builder)
# ---------------------------------------------------------------------------

def _get_intron_chain(exon_df):
    """Return a string signature of the intron chain for one transcript."""
    exon_df = exon_df.sort_values('Start')
    if len(exon_df) < 2:
        return 'single-exon'
    ends = exon_df['End'].tolist()
    starts = exon_df['Start'].tolist()
    return ','.join(f'{ends[i]}-{starts[i+1]}' for i in range(len(starts) - 1))


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def score_model(model, config, protein_supported_tids, genome=None):
    """Score a single gene model.

    Parameters
    ----------
    model : dict
        Keys: id, source, chrom, strand, intron_chain, df (exon DataFrame),
        start, end, exon_count, combined_evidence.
    config : PipelineConfig
    protein_supported_tids : set of str
    genome : dict or None
        If provided, enables splice-site scoring.

    Returns
    -------
    float : composite score.
    """
    scfg = config.scoring
    score = 0.0

    # Base evidence weight
    sources = set(model.get('combined_evidence', model['source']).split(','))
    weights = scfg.weights
    for s in sources:
        s_lower = s.strip().lower()
        if s_lower == 'helixer':
            score += weights.helixer
        elif s_lower == 'scallop':
            score += weights.scallop
        elif s_lower == 'stringtie':
            score += weights.stringtie
        else:
            score += 1.0  # unknown source gets base weight

    # Multi-source bonus
    if len(sources) > 1:
        score += scfg.multi_source_bonus * (len(sources) - 1)

    # Protein support bonus
    if model['id'] in protein_supported_tids:
        score += scfg.protein_overlap_bonus
        model['protein_support'] = True
    else:
        model['protein_support'] = model.get('protein_support', False)

    # Splice-site penalty (only if genome provided and multi-exon)
    if genome and model['exon_count'] > 1:
        chrom = model['chrom']
        if chrom in genome:
            exons = sorted(
                zip(model['df']['Start'].values, model['df']['End'].values))
            splice = check_splice_sites(exons, model['strand'], genome[chrom])
            n_noncanonical = sum(
                1 for s in splice if s['class'] != 'canonical')
            score -= scfg.noncanonical_splice_penalty * n_noncanonical

    return score


# ---------------------------------------------------------------------------
# Isoform selection
# ---------------------------------------------------------------------------

def select_isoforms(locus_df, config, protein_supported_tids,
                    genome=None):
    """Score and select isoforms for a single locus.

    Parameters
    ----------
    locus_df : DataFrame
        Exon rows for all models in this locus.
    config : PipelineConfig
    protein_supported_tids : set
    genome : dict or None

    Returns
    -------
    list of model dicts (selected isoforms), sorted by start position.
    """
    scfg = config.scoring

    # Build model dicts
    models = []
    for (source, tid), grp in locus_df.groupby(['Source', 'transcript_id']):
        chain = _get_intron_chain(grp)
        models.append({
            'id': tid,
            'source': source,
            'chrom': grp['Chromosome'].iloc[0],
            'strand': grp['Strand'].iloc[0],
            'intron_chain': chain,
            'protein_support': tid in protein_supported_tids,
            'df': grp,
            'start': grp['Start'].min(),
            'end': grp['End'].max(),
            'exon_count': len(grp),
            'combined_evidence': grp['combined_evidence'].iloc[0]
                if 'combined_evidence' in grp.columns else source,
        })

    if not models:
        return []

    # Merge identical structures across sources
    merged = {}
    for m in models:
        if m['intron_chain'] == 'single-exon':
            key = f"{m['chrom']}:{m['strand']}:{m['start']}-{m['end']}"
        else:
            key = f"{m['chrom']}:{m['strand']}:{m['intron_chain']}"

        if key not in merged:
            merged[key] = {
                'sources': set(),
                'protein_support': False,
                'rep': m,
                'score': 0.0,
            }
        s = merged[key]
        s['sources'].add(m['source'])
        if m['protein_support']:
            s['protein_support'] = True
            if not s['rep']['protein_support']:
                s['rep'] = m

    # Score each merged structure
    for key, s in merged.items():
        rep = s['rep']
        rep['combined_evidence'] = ','.join(sorted(s['sources']))
        if s['protein_support']:
            rep['protein_support'] = True
        s['score'] = score_model(
            rep, config, protein_supported_tids, genome)

    # Quality gate: keep models meeting minimum criteria
    candidates = []
    for s in merged.values():
        keep = False
        if s['protein_support']:
            keep = True
        elif 'Helixer' in s['sources']:
            keep = True
        elif len(s['sources']) > 1:
            keep = True
        elif 'StringTie' in s['sources'] and s['rep']['exon_count'] > 1:
            keep = True
        elif (scfg.fungal_single_exon_mode and
              s['rep']['exon_count'] == 1 and
              s['score'] >= scfg.min_alternate_score):
            # Fungal mode: keep well-supported single-exon models
            keep = True

        if keep:
            candidates.append(s)

    if not candidates:
        return []

    # Sort by score descending
    candidates.sort(key=lambda s: s['score'], reverse=True)

    # Select primary + alternates
    selected = []
    primary = candidates[0]
    primary['rep']['is_primary'] = True
    selected.append(primary['rep'])

    # Alternates: different intron chain, above threshold
    for s in candidates[1:]:
        if len(selected) >= scfg.max_isoforms_per_locus:
            break
        if s['score'] >= scfg.min_alternate_score:
            # Must have different structure from primary
            if s['rep']['intron_chain'] != primary['rep']['intron_chain']:
                s['rep']['is_primary'] = False
                selected.append(s['rep'])

    selected.sort(key=lambda m: m['start'])
    return selected
