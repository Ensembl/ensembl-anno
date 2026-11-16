"""Evidence-filtered secondary-ORF candidates for evidence-score primary selection.

Opt-in: ``scoring.evidence_score.secondary_orfs: true``, only with
``scoring.primary_selection: evidence_score`` and construction-only STAR junctions
(``--star-junctions``). Nothing here runs otherwise. This is the ``A_ev`` generator
frozen in the GMB Phase 3A experiment (Z. tritici development windows) and evaluated on
B. cinerea; see ``docs/evidence_selection.md``.

Discovery (after ORF annotation)
    Transcript-evidence models (the short- and long-read classes, the same sources the
    evidence scorer gives the transcript-class prior) are scanned for every complete,
    maximal ATG..stop ORF (:func:`~gmb.pipeline.annotate_cds_utrs.find_complete_orfs`)
    of at least ``orf.min_codons`` codons and the evidence-score CDS floor
    (``max(min_candidate_cds_bp, scoring.min_cds_bp)``) that lies entirely upstream or
    downstream of the transcript's own GMB ORF on the cDNA. ORFs overlapping the own ORF
    in another frame are not generated. Each ORF is projected with
    :func:`~gmb.pipeline.annotate_cds_utrs.map_cds_to_genomic` and must re-read
    identically from the genome, start with ATG, end with a stop codon, contain no
    internal stop, be a multiple of three and pass the frame-continuity check. CDS
    structures are deduplicated across parent transcripts (every parent is kept as
    provenance); a structure that already is any candidate's CDS is not generated.

Admission (A_ev; before primary selection)
    Only a structure that shares a coding locus (``same_coding_locus_fraction``,
    reciprocal) with at least one existing CDS structure of its parent's locus can ever
    enter an existing gene's evidence-score pool, so only those are evaluated. One is
    admitted when at least half of its CDS is covered by same-strand protein-alignment
    exons, or it is spliced and every CDS intron is observed (unique reads, same strand
    or unstranded) in the construction junction index. Protein alignments are used as
    supplied: the target species' own proteins must be excluded when the protein set is
    prepared (GMB cannot recognise them).

Selection
    Admitted structures never enter clustering, gene grouping or isoform admission. Each
    joins only the evidence-score pool of an existing gene of its parent's locus, under
    that pool's own rules, and becomes a primary only by the unchanged P4d rule. Its
    transcript is the CDS itself (exons = CDS segments): no UTR is carried over from a
    read-through parent into a neighbouring locus. It keeps its parents' transcript
    evidence class (``Evidence``) and is marked ``CandidateOrigin=secondary_orf``.
"""

from __future__ import annotations

import bisect
from collections import Counter
from dataclasses import dataclass, field

import pandas as pd

from gmb.pipeline.annotate_cds_utrs import (
    START_CODON,
    STOP_CODONS,
    annotate_transcript,
    check_frame_continuity,
    find_complete_orfs,
    map_cds_to_genomic,
    reverse_complement,
    translate,
)
from gmb.pipeline.evidence_features import same_coding_locus
from gmb.pipeline.junction_diagnostics import introns_of

#: ``CandidateOrigin`` value written for generated transcripts.
ORIGIN = "secondary_orf"
#: A_ev protein clause: minimum fraction of the CDS covered by protein-alignment exons.
PROTEIN_COVER_MIN = 0.5


@dataclass
class SecondaryOrf:
    """One generated CDS structure and its provenance."""

    cid: str
    chrom: str
    strand: str
    cds: tuple
    parents: list = field(default_factory=list)
    loci: set = field(default_factory=set)
    status: str = "discovered"
    protein_cover: float | None = None
    protein_alignments: int | None = None
    introns_supported: int | None = None
    admitted_by: str = ""
    model: dict | None = None
    annotation: dict | None = None
    pool_entries: list = field(default_factory=list)

    @property
    def key(self) -> tuple:
        return (self.chrom, self.strand, self.cds)

    @property
    def cds_bp(self) -> int:
        return sum(b - a for a, b in self.cds)

    @property
    def n_introns(self) -> int:
        return len(introns_of(self.cds))

    @property
    def sources(self) -> frozenset:
        return frozenset(p["source"] for p in self.parents)


def _tx_offset(exons, strand, g):
    """Transcript (5'->3') offset of genomic base *g*, or None if it is not exonic."""
    order = sorted(exons) if strand == "+" else sorted(exons, reverse=True)
    cum = 0
    for a, b in order:
        if a <= g < b:
            return cum + (g - a if strand == "+" else b - 1 - g)
        cum += b - a
    return None


def own_orf_interval(ann: dict, strand: str):
    """cDNA interval of the transcript's own GMB ORF, or None when it has none."""
    cds = sorted(ann.get("cds") or [])
    if not cds:
        return None
    g5 = cds[0][0] if strand == "+" else cds[-1][1] - 1
    start = _tx_offset(ann["exons"], strand, g5)
    if start is None:
        return None
    return start, start + sum(b - a for a, b in cds)


def _genomic_cds(chrom_seq: str, strand: str, cds) -> str:
    seq = "".join(chrom_seq[a:b] for a, b in sorted(cds)).upper()
    return reverse_complement(seq) if strand == "-" else seq


def validate_cds(chrom_seq: str, strand: str, cds, orf_seq: str):
    """(ok, reason) for a projected complete ORF."""
    seq = _genomic_cds(chrom_seq, strand, cds)
    if seq != orf_seq.upper():
        return False, "projection_mismatch"
    if len(seq) % 3:
        return False, "frame"
    if seq[:3] != START_CODON:
        return False, "no_atg"
    if seq[-3:] not in STOP_CODONS:
        return False, "no_stop"
    if "*" in translate(seq)[:-1]:
        return False, "internal_stop"
    if not check_frame_continuity(list(cds), strand):
        return False, "phase"
    return True, "ok"


def discover(
    annotations: dict,
    transcripts,
    genome: dict,
    transcript_sources,
    min_codons: int,
    min_bp: int,
    existing_keys,
) -> tuple[dict, Counter]:
    """Generate secondary-ORF structures.

    *transcripts* is an iterable of ``(tid, source, chrom, strand)`` in a deterministic
    order (the candidate order); the first parent to emit a structure names it
    (``<parent>.sorf<n>``). Returns ``({key: SecondaryOrf}, counters)``.
    """
    stats = Counter()
    found: dict = {}
    per_parent: Counter = Counter()
    sources = frozenset(transcript_sources)
    for tid, source, chrom, strand in transcripts:
        if source not in sources or strand not in ("+", "-") or chrom not in genome:
            continue
        ann = annotations.get(tid)
        if not ann or not ann.get("cdna") or not ann.get("exons"):
            continue
        stats["transcripts_scanned"] += 1
        exons = sorted(ann["exons"])
        cdna = ann["cdna"]
        own = own_orf_interval(ann, strand)
        for s, e, _frame in find_complete_orfs(cdna, min_codons=min_codons, min_bp=min_bp):
            if own is not None and (s, e) == own:
                continue
            if own is None:
                relation = "no_gmb_orf"
            elif e <= own[0]:
                relation = "upstream"
            elif s >= own[1]:
                relation = "downstream"
            else:
                stats["overlapping_other_frame_not_generated"] += 1
                continue
            stats["orfs_found"] += 1
            cds = tuple((int(a), int(b)) for a, b in map_cds_to_genomic(s, e, exons, strand))
            ok, why = validate_cds(genome[chrom], strand, cds, cdna[s:e])
            if not ok:
                stats[f"rejected_{why}"] += 1
                continue
            key = (str(chrom), strand, cds)
            if key in existing_keys:
                stats["equal_to_existing_candidate"] += 1
                continue
            c = found.get(key)
            if c is None:
                per_parent[tid] += 1
                c = found[key] = SecondaryOrf(
                    cid=f"{tid}.sorf{per_parent[tid]}", chrom=str(chrom), strand=strand, cds=cds
                )
            else:
                stats["duplicate_emissions_merged"] += 1
            c.parents.append(
                {
                    "tid": tid,
                    "source": str(source),
                    "relation": relation,
                    "tx_start": s,
                    "tx_end": e,
                    "tx_len": len(cdna),
                    "own_orf_bp": (own[1] - own[0]) if own else 0,
                }
            )
    stats["unique_structures"] = len(found)
    return found, stats


class ProteinCoverage:
    """Same-strand protein-alignment exon coverage of CDS structures."""

    def __init__(self, prot_exons: pd.DataFrame | None):
        self._iv: dict = {}
        self._raw: dict = {}
        self._max = 0
        if prot_exons is None or prot_exons.empty:
            return
        by: dict = {}
        cols = prot_exons[["transcript_id", "Chromosome", "Strand", "Start", "End"]]
        for t, c, s, a, b in cols.itertuples(index=False):
            by.setdefault((str(c), str(s)), []).append((int(a), int(b), t))
        for k, v in by.items():
            v.sort()
            merged: list = []
            for a, b, _ in v:
                if merged and a <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], b)
                else:
                    merged.append([a, b])
            self._iv[k] = ([m[0] for m in merged], merged)
            self._raw[k] = ([x[0] for x in v], v)
            self._max = max([self._max] + [b - a for a, b, _ in v])

    def cover(self, chrom: str, strand: str, cds) -> tuple[float, int]:
        """(fraction of CDS bases covered, distinct alignments overlapping the CDS)."""
        ent = self._iv.get((chrom, strand))
        if not ent:
            return 0.0, 0
        starts, merged = ent
        bp = 0
        for a, b in cds:
            i = max(0, bisect.bisect_right(starts, a) - 1)
            while i < len(merged) and merged[i][0] < b:
                bp += max(0, min(b, merged[i][1]) - max(a, merged[i][0]))
                i += 1
        rs, raw = self._raw[(chrom, strand)]
        i = bisect.bisect_left(rs, cds[0][0] - self._max)
        alignments = set()
        while i < len(raw) and raw[i][0] < cds[-1][1]:
            a, b, t = raw[i]
            if any(min(b, y) - max(a, x) > 0 for x, y in cds):
                alignments.add(t)
            i += 1
        return bp / sum(b - a for a, b in cds), len(alignments)


def junction_observed(index, chrom: str, strand: str, start0: int, end: int) -> bool:
    """The intron is a construction junction with unique reads on this strand (or unstranded)."""
    j = index.find(chrom, start0, end)
    code = 1 if strand == "+" else 2
    return j is not None and int(index.j_strand[j]) in (0, code) and int(index.j_unique[j]) > 0


def admit(
    found: dict,
    tx_loci: dict,
    locus_structures: dict,
    index,
    coverage: ProteinCoverage,
    fraction: float,
) -> Counter:
    """Apply the A_ev rule to structures that could enter an existing gene's pool."""
    stats = Counter()
    for c in found.values():
        c.loci = set()
        for p in c.parents:
            c.loci |= tx_loci.get(p["tid"], set())
        if not any(
            k[0] == c.chrom and k[1] == c.strand and same_coding_locus(c.cds, k[2], fraction)
            for li in sorted(c.loci)
            for k in locus_structures.get(li, ())
        ):
            c.status = "no_shared_coding_locus"
            stats[c.status] += 1
            continue
        c.protein_cover, c.protein_alignments = coverage.cover(c.chrom, c.strand, c.cds)
        introns = introns_of(c.cds)
        c.introns_supported = sum(
            junction_observed(index, c.chrom, c.strand, a, b) for a, b in introns
        )
        by_protein = c.protein_cover >= PROTEIN_COVER_MIN
        by_junctions = bool(introns) and c.introns_supported == len(introns)
        if by_protein or by_junctions:
            c.status = "admitted"
            c.admitted_by = "+".join(
                x for x, y in (("protein", by_protein), ("junctions", by_junctions)) if y
            )
        else:
            c.status = "not_admitted"
        stats[c.status] += 1
    return stats


def build_transcript(c: SecondaryOrf, genome: dict, min_codons: int) -> None:
    """Attach the GMB model dict and ORF annotation of a generated transcript (exons = CDS)."""
    starts = [a for a, _ in c.cds]
    ends = [b for _, b in c.cds]
    src = sorted(c.sources)
    df = pd.DataFrame(
        {
            "Chromosome": c.chrom,
            "Source": src[0],
            "Strand": c.strand,
            "Start": starts,
            "End": ends,
            "transcript_id": c.cid,
        }
    )
    ann = annotate_transcript(
        df, c.chrom, c.strand, genome, cds_df=df[["Start", "End"]], min_codons=min_codons
    )
    chain = (
        "single-exon"
        if len(c.cds) < 2
        else ",".join(f"{ends[i]}-{starts[i + 1]}" for i in range(len(starts) - 1))
    )
    evidence = ",".join(src)
    c.annotation = ann
    c.model = {
        "id": c.cid,
        "source": src[0],
        "chrom": c.chrom,
        "strand": c.strand,
        "intron_chain": chain,
        "df": df,
        "start": starts[0],
        "end": ends[-1],
        "exon_count": len(c.cds),
        "combined_evidence": evidence,
        "protein_support": False,
        "protein_positional_support": False,
        "protein_cds_span_compatible": False,
        "protein_sources": set(),
        "protein_evidence": "",
        "protein_coding_score": None,
        "score": None,
        "has_backbone": False,
        "candidate_cds": list(c.cds),
        "candidate_cds_bp": c.cds_bp,
        "backbone_intron_rescue": False,
        "structural_support_sources": evidence,
        "n_structural_support_sources": len(src),
        "backbone_shortread_agreement": False,
        "protein_support_strength": "",
        "longread_structural_role": "not_longread",
        "candidate_origin": ORIGIN,
        "secondary_orf_parents": ",".join(p["tid"] for p in c.parents),
    }


def attribution_rows(found: dict) -> list[dict]:
    """One row per generated structure (secondary_orfs.tsv)."""
    rows = []
    for c in found.values():
        sel = [e for e in c.pool_entries if e["selected"]]
        rows.append(
            {
                "candidate_id": c.cid,
                "chrom": c.chrom,
                "strand": c.strand,
                "cds": ",".join(f"{a + 1}-{b}" for a, b in c.cds),
                "cds_bp": c.cds_bp,
                "n_cds_exons": len(c.cds),
                "status": c.status,
                "admitted_by": c.admitted_by,
                "protein_cover": None if c.protein_cover is None else round(c.protein_cover, 4),
                "protein_alignments": c.protein_alignments,
                "n_introns": c.n_introns,
                "introns_junction_supported": c.introns_supported,
                "n_parents": len(c.parents),
                "parents": ",".join(p["tid"] for p in c.parents),
                "parent_sources": ",".join(sorted(c.sources)),
                "relations": ",".join(sorted({p["relation"] for p in c.parents})),
                "loci": ",".join(str(x) for x in sorted(c.loci)),
                "eligible_genes": len(c.pool_entries),
                "eligible_primaries": ",".join(e["primary"] for e in c.pool_entries),
                "selected": bool(sel),
                "selected_gene": ";".join(e.get("gene_id", "") for e in sel),
                "decision_at_eligible_genes": ";".join(
                    e["outcome"] + (":chosen" if e["selected"] else "") for e in c.pool_entries
                ),
            }
        )
    return rows
