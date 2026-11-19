"""Coding units of filter-removed transcript models (EXPERIMENTAL, Phase 4 proof of concept).

Opt-in: ``scoring.evidence_score.filtered_transcript_orfs: true``. It requires
``scoring.evidence_score.secondary_orfs: true`` (and therefore evidence-score primary
selection and construction STAR junctions). Nothing here runs otherwise, and a build with the
option off is byte-identical to one without the key.

Why
    The transcriptomic filter (``evidence_filter.filter_chimeras``) removes transcript models
    with an intron longer than ``transcriptomic_filter.max_intron_length`` or a span longer than
    ``transcriptomic_filter.max_transcript_length``. Most of the removed models are
    read-through assemblies that join adjacent genes. Removing them keeps fused structures out
    of clustering, but it also discards every coding unit they carry; Phase 3B's secondary-ORF
    discovery only scans the models that survived the filter. On four fungal genomes
    (GMB Phase 4, ``PHASE4_CANDIDATE_CEILING.md``) a reproducible class of missed reference
    genes is a complete ORF of such a removed model, sharing a coding locus with an existing
    gene's candidates.

What it does
    Each removed model is annotated exactly as a candidate would be (``annotate_all_transcripts``)
    and scanned with the frozen secondary-ORF discovery (:func:`gmb.pipeline.secondary_orfs.discover`,
    the same transcript sources, ``orf.min_codons`` and CDS floor). Its own GMB ORF is a unit too
    when it is complete (GMB never saw it). A unit is never a transcript: like a secondary ORF it
    becomes a CDS-only candidate whose exons are its CDS segments, so no read-through exon or UTR
    re-enters clustering.

Safeguards (each recorded as a status in ``secondary_orfs.tsv``)
    * complete ATG..stop ORF, re-read from the genome (``secondary_orfs.validate_cds``);
    * the unit itself would pass the transcript filter it came from: no CDS intron longer than
      ``max_intron_length``, CDS span within ``max_transcript_length`` (``exceeds_filter``);
    * not equal to an existing candidate's CDS or to a secondary ORF of a kept model
      (counted, not generated);
    * no new gene: a unit may only join the evidence-score pool of an existing gene, through the
      frozen A_ev admission (:func:`gmb.pipeline.secondary_orfs.admit`, unchanged): it must share a
      coding locus with an existing CDS structure and have >= 50% protein-alignment cover or every
      CDS intron observed in construction junctions;
    * no fusion: a unit that shares a coding locus with structures of two or more loci is rejected
      (``multi_locus``);
    * bounded: at most :data:`MAX_ADMITTED_PER_LOCUS` admitted units per locus; the units carried
      by the most removed models are kept (then the longest CDS, then discovery order), the rest
      are recorded as ``locus_cap``.
    Selection is unchanged: an admitted unit becomes a primary only by the frozen P4d rule, and is
    marked ``CandidateOrigin=filtered_transcript_orf``.
"""

from __future__ import annotations

import bisect
from collections import Counter, defaultdict

from gmb.pipeline.evidence_features import same_coding_locus
from gmb.pipeline.junction_diagnostics import introns_of
from gmb.pipeline.secondary_orfs import SecondaryOrf, admit, own_orf_interval, validate_cds

#: ``CandidateOrigin`` of a generated unit.
ORIGIN = "filtered_transcript_orf"
#: Safety bound on admitted units per locus (not a tuned parameter; reported when reached).
MAX_ADMITTED_PER_LOCUS = 10


def removed_transcripts(tx_exons, tx_kept):
    """Exon rows of the transcript models the transcriptomic filter removed."""
    if tx_exons is None or tx_exons.empty:
        return tx_exons
    kept = (
        set(tx_kept["transcript_id"].unique())
        if tx_kept is not None and not tx_kept.empty
        else set()
    )
    return tx_exons[~tx_exons["transcript_id"].isin(kept)].copy()


def _filter_reason(exons, max_intron, max_span) -> str:
    ex = sorted(exons)
    why = []
    if max_intron is not None and any(
        ex[i + 1][0] - ex[i][1] > max_intron for i in range(len(ex) - 1)
    ):
        why.append("intron")
    if max_span is not None and ex[-1][1] - ex[0][0] > max_span:
        why.append("span")
    return "+".join(why) or "other"


def discover_removed(
    removed_df,
    annotations: dict,
    genome: dict,
    transcript_sources,
    min_codons: int,
    min_bp: int,
    existing_keys,
    max_intron,
    max_span,
    discover_fn,
):
    """Units of removed transcript models: ``({key: SecondaryOrf}, Counter)``.

    *annotations* holds the ORF annotation of the removed models only. *existing_keys* holds
    every existing candidate CDS structure and every secondary ORF of a kept model; a unit equal
    to one of them is counted, not generated. *discover_fn* is
    :func:`gmb.pipeline.secondary_orfs.discover` (passed in so the frozen function is used as is).
    """
    stats = Counter()
    if removed_df is None or removed_df.empty:
        return {}, stats
    order = removed_df.drop_duplicates("transcript_id")[
        ["transcript_id", "Source", "Chromosome", "Strand"]
    ]
    transcripts = [(t, str(s), str(c), str(st)) for t, s, c, st in order.itertuples(index=False)]
    sources = frozenset(transcript_sources)
    stats["removed_models"] = len(transcripts)
    found, dstats = discover_fn(
        annotations, transcripts, genome, sources, min_codons, min_bp, set(existing_keys)
    )
    for k, v in dstats.items():
        stats[f"flanking_{k}"] += v
    # The removed model's own GMB ORF, when complete, is a unit as well.
    per_parent: Counter = Counter()
    for tid, source, chrom, strand in transcripts:
        if source not in sources or chrom not in genome:
            continue
        ann = annotations.get(tid)
        if not ann or not ann.get("cds") or ann.get("is_partial_5") or ann.get("is_partial_3"):
            continue
        own = own_orf_interval(ann, strand)
        cds = tuple(sorted((int(a), int(b)) for a, b in ann["cds"]))
        if own is None or sum(b - a for a, b in cds) < min_bp:
            continue
        ok, why = validate_cds(genome[chrom], strand, cds, ann["cdna"][own[0] : own[1]])
        if not ok:
            stats[f"own_rejected_{why}"] += 1
            continue
        key = (str(chrom), strand, cds)
        if key in existing_keys:
            stats["own_equal_to_existing"] += 1
            continue
        c = found.get(key)
        if c is None:
            per_parent[tid] += 1
            c = found[key] = SecondaryOrf(
                cid=f"{tid}.ftown{per_parent[tid]}", chrom=str(chrom), strand=strand, cds=cds
            )
        c.parents.append(
            {
                "tid": tid,
                "source": str(source),
                "relation": "own",
                "tx_start": own[0],
                "tx_end": own[1],
                "tx_len": len(ann["cdna"]),
                "own_orf_bp": own[1] - own[0],
            }
        )
    # Provenance: every parent is a removed model; units are renamed so they are recognisable.
    reasons = {}
    for tid, grp in removed_df.groupby("transcript_id", sort=False):
        reasons[tid] = _filter_reason(
            list(zip(grp["Start"].astype(int), grp["End"].astype(int))), max_intron, max_span
        )
    renamed = Counter()
    for c in found.values():
        for p in c.parents:
            p["relation"] = f"removed_{p['relation']}"
            p["filter_reason"] = reasons.get(p["tid"], "other")
        if ".sorf" in c.cid:
            base = c.cid.rsplit(".sorf", 1)[0]
            renamed[base] += 1
            c.cid = f"{base}.ftorf{renamed[base]}"
    # The unit must itself pass the filter its parent failed.
    for c in found.values():
        ints = introns_of(c.cds)
        long_intron = max_intron is not None and any(b - a > max_intron for a, b in ints)
        long_span = max_span is not None and c.cds[-1][1] - c.cds[0][0] > max_span
        if long_intron or long_span:
            c.status = "exceeds_filter"
    stats["units"] = len(found)
    stats["units_exceeding_filter"] = sum(c.status == "exceeds_filter" for c in found.values())
    return found, stats


class _StructureIndex:
    """Same-strand lookup of (CDS structure, locus) by CDS span."""

    def __init__(self, locus_structures: dict):
        by = defaultdict(list)
        for li, keys in locus_structures.items():
            for chrom, strand, cds in keys:
                by[(chrom, strand)].append((cds[0][0], cds[-1][1], cds, li))
        self._by, self._starts, self._max = {}, {}, {}
        for k, v in by.items():
            v.sort(key=lambda x: (x[0], x[1], str(x[3])))
            self._by[k] = v
            self._starts[k] = [x[0] for x in v]
            self._max[k] = max(x[1] - x[0] for x in v)

    def loci_sharing(self, chrom, strand, cds, fraction) -> set:
        k = (chrom, strand)
        if k not in self._by:
            return set()
        v, st = self._by[k], self._starts[k]
        lo, hi = cds[0][0], cds[-1][1]
        i = bisect.bisect_left(st, lo - self._max[k])
        out = set()
        while i < len(v) and v[i][0] < hi:
            if v[i][1] > lo and same_coding_locus(cds, v[i][2], fraction):
                out.add(v[i][3])
            i += 1
        return out


def admit_removed(
    found: dict, locus_structures: dict, index, coverage, fraction: float
) -> Counter:
    """Locus assignment by shared coding locus, then the frozen A_ev admission, then the bounds."""
    stats = Counter()
    six = _StructureIndex(locus_structures)
    for c in found.values():
        if c.status == "exceeds_filter":
            stats["exceeds_filter"] += 1
            continue
        loci = six.loci_sharing(c.chrom, c.strand, c.cds, fraction)
        if len(loci) > 1:
            c.loci = loci
            c.status = "multi_locus"
            stats["multi_locus"] += 1
            continue
        # The frozen rule, applied with this unit's parents mapped to the unit's locus.
        tx_loci = {p["tid"]: set(loci) for p in c.parents}
        stats.update(admit({c.key: c}, tx_loci, locus_structures, index, coverage, fraction))
    # Bound: keep the units carried by the most removed models (then the longest CDS, then
    # discovery order); the rest of a locus over the bound are recorded as ``locus_cap``.
    by_locus = defaultdict(list)
    for i, c in enumerate(found.values()):
        if c.status == "admitted":
            by_locus[next(iter(c.loci))].append(
                (-len({p["tid"] for p in c.parents}), -c.cds_bp, i, c)
            )
    per_locus: Counter = Counter()
    for li, items in by_locus.items():
        for *_rank, c in sorted(items, key=lambda x: x[:3]):
            per_locus[li] += 1
            if per_locus[li] > MAX_ADMITTED_PER_LOCUS:
                c.status = "locus_cap"
                stats["admitted"] -= 1
                stats["locus_cap"] += 1
    stats["loci_at_cap"] = sum(1 for n in per_locus.values() if n > MAX_ADMITTED_PER_LOCUS)
    return stats
