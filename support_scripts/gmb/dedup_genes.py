#!/usr/bin/env python3
"""Post-processing gene deduplication for Gene Model Builder.

Clusters output genes by high reciprocal overlap on the same strand and
resolves duplicates by merging identical structures or keeping best-scored.
"""


def _reciprocal_overlap(s1, e1, s2, e2):
    """Compute reciprocal overlap fraction between two intervals."""
    overlap = max(0, min(e1, e2) - max(s1, s2))
    len1 = e1 - s1
    len2 = e2 - s2
    if min(len1, len2) == 0:
        return 0.0
    return overlap / min(len1, len2)


def _get_exon_chain(children):
    """Extract sorted exon chain from child rows."""
    exons = [(r['Start'], r['End']) for r in children if r['Feature'] == 'exon']
    return sorted(exons)


def _chains_match(chain1, chain2, tolerance_bp=10):
    """Check if two exon chains are structurally equivalent within tolerance."""
    if len(chain1) != len(chain2):
        return False
    for (s1, e1), (s2, e2) in zip(chain1, chain2):
        if abs(s1 - s2) > tolerance_bp or abs(e1 - e2) > tolerance_bp:
            return False
    return True


def dedup_genes(gff_rows, config):
    """Deduplicate overlapping genes in GFF3 output.

    Parameters
    ----------
    gff_rows : list of dict
        GFF3 rows with Feature, Start, End, Strand, ID, Parent keys.
    config : PipelineConfig
        Must have .dedup attribute.

    Returns
    -------
    list of dict : deduplicated rows
    dict : stats
    """
    dedup_cfg = config.dedup
    if not dedup_cfg.enabled:
        return gff_rows, {'dedup_enabled': False}

    stats = {
        'genes_input': 0,
        'genes_merged': 0,
        'genes_dropped': 0,
        'genes_output': 0,
    }

    # Index rows by parent/gene
    by_parent = {}
    genes = []
    for r in gff_rows:
        if r.get('Feature') == 'gene':
            genes.append(r)
        parent = r.get('Parent', '')
        if parent:
            by_parent.setdefault(parent, []).append(r)

    stats['genes_input'] = len(genes)

    # Build gene metadata
    gene_info = []
    for g in genes:
        gid = g['ID']
        mrnas = [r for r in by_parent.get(gid, []) if r['Feature'] == 'mRNA']
        all_children = by_parent.get(gid, [])
        # Get first mRNA's exon chain for comparison
        chains = []
        for m in mrnas:
            mc = by_parent.get(m['ID'], [])
            chains.append(_get_exon_chain(mc))
        gene_info.append({
            'gene_row': g,
            'mrnas': mrnas,
            'chains': chains,
            'children': all_children,
            'chrom': g['Chromosome'],
            'strand': g['Strand'],
            'start': g['Start'],
            'end': g['End'],
            'merged_into': None,
            'dropped': False,
        })

    # Cluster genes by overlap on same strand
    threshold = dedup_cfg.reciprocal_overlap_threshold
    tolerance = dedup_cfg.same_structure_tolerance_bp

    for i in range(len(gene_info)):
        gi = gene_info[i]
        if gi['dropped'] or gi['merged_into'] is not None:
            continue
        for j in range(i + 1, len(gene_info)):
            gj = gene_info[j]
            if gj['dropped'] or gj['merged_into'] is not None:
                continue
            if gi['chrom'] != gj['chrom']:
                continue

            # Opposite strand: always keep both
            if gi['strand'] != gj['strand']:
                if dedup_cfg.policy == 'keep_both_if_opposite_strand':
                    continue
                else:
                    continue

            ro = _reciprocal_overlap(gi['start'], gi['end'],
                                     gj['start'], gj['end'])
            if ro < threshold:
                continue

            # High overlap on same strand — check structure
            structures_match = False
            if gi['chains'] and gj['chains']:
                structures_match = _chains_match(
                    gi['chains'][0], gj['chains'][0], tolerance)

            if structures_match:
                if dedup_cfg.policy == 'merge_as_isoforms':
                    gj['merged_into'] = i
                    stats['genes_merged'] += 1
                else:
                    gj['dropped'] = True
                    stats['genes_dropped'] += 1
            else:
                if dedup_cfg.policy == 'keep_best_drop_rest':
                    # Keep the one with more evidence sources
                    ev_i = len(gi.get('mrnas', []))
                    ev_j = len(gj.get('mrnas', []))
                    if ev_j > ev_i:
                        gi['dropped'] = True
                        stats['genes_dropped'] += 1
                    else:
                        gj['dropped'] = True
                        stats['genes_dropped'] += 1

    # Rebuild output
    output_rows = []
    for i, gi in enumerate(gene_info):
        if gi['dropped']:
            continue
        if gi['merged_into'] is not None:
            continue

        output_rows.append(gi['gene_row'])

        # Add own mRNAs and their children
        for m in gi['mrnas']:
            output_rows.append(m)
            output_rows.extend(by_parent.get(m['ID'], []))

        # Merge any absorbed gene's mRNAs as isoforms
        for j, gj in enumerate(gene_info):
            if gj.get('merged_into') == i:
                for m in gj['mrnas']:
                    m_copy = dict(m)
                    m_copy['Parent'] = gi['gene_row']['ID']
                    output_rows.append(m_copy)
                    output_rows.extend(by_parent.get(m['ID'], []))

    stats['genes_output'] = sum(1 for r in output_rows if r.get('Feature') == 'gene')
    return output_rows, stats
