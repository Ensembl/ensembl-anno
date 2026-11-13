"""Quantitative evidence features for CDS candidates, and the evidence-score primary rule.

This module is used only when ``scoring.primary_selection: evidence_score`` is
configured; nothing here runs otherwise. It turns the Phase 1 STAR junction index
(:mod:`gmb.pipeline.star_junctions`), spliced protein alignments and the genome into
bounded, interpretable per-intron and per-model features for any CDS structure.
These are GMB candidates today and locally reconstructed CDS candidates later.
It also implements the opt-in replacement rule ("P4d") that chooses a gene's primary
among the existing candidates of its coding locus.

Features (all bounded; see ``docs/evidence_selection.md``)
-----------------------------------------------------------
Per CDS intron, from :class:`gmb.pipeline.junction_diagnostics.JunctionAnnotator`:

* ``J`` junction term. ``f(g) = min(1, g / breadth_saturation_groups)`` over
  independent condition groups (not libraries, so replicates count once):

  - SUPPORTED: ``f(own) * u - f(best competitor) * (1 - u)``, with ``u`` the lower
    of donor and acceptor usage against junctions sharing that site;
  - ALTERNATIVE_SITE_OBSERVED: ``-f(competitor)``; OPPOSITE_STRAND: ``-f(opposite)``;
  - UNOBSERVED_EXPRESSED_LOCUS: ``-unobserved_expressed_weight * E``, where ``E`` is
    the median ``f`` of the candidate's *own* supported introns (0 if it has none);
  - unobservable (outside the aligner's intron-length range, filtered motif),
    unobserved without local expression, multimapper-only: 0.
* ``P`` protein term: ``protein_intron_weight * (2a / (a + c) - 1)`` over spliced
  protein alignments overlapping the intron that contain it exactly (``a``) or
  contain a different overlapping intron (``c``). Unspliced or merely overlapping
  alignments are uninformative (0). The fraction bounds the term, so repeated
  alignments of one protein cannot add more than one alignment's sign.
* motif: ``-noncanonical_penalty`` for a non-canonical genome motif when ``J >= 0``.

Per model: integrity (``-partial_orf_penalty`` for a partial ORF,
``-internal_stop_penalty`` for an internal stop or broken frame), a bounded backbone
prior and a bounded transcript-class prior. StringTie, Scallop and long-read models
are **one** class, counted once; junction support is never added per assembler.
``evidence`` = J + P + motif + integrity; ``total`` = evidence + priors.

Replacement rule (P4d)
----------------------
A gene's primary *B* is replaced by a candidate *X* from the same coding locus only if:

1. *B* is not a protected backbone: a backbone structure with no intron of negative
   ``J`` and intact integrity is never replaced;
2. *X* has a complete ORF, no internal stop and an intact frame;
3. *X* does not drop any intron of *B* with ``J >= keep_supported_intron_J``;
4. ``total(X) - total(B) >= replacement_margin`` **and**
   ``evidence(X) - evidence(B) >= replacement_margin``. A prior can therefore
   protect *B* but can never be the reason to switch.

Among eligible candidates the highest ``total`` (then ``evidence``) is taken;
remaining ties keep the earlier candidate in locus order.

Evidence roles
--------------
The junction index used for selection must contain only ``construction``
libraries. :func:`require_selection_index` enforces this wherever an index reaches
selection; held-out or auxiliary libraries are refused, never silently dropped.
"""

from __future__ import annotations

import statistics
import time
from dataclasses import dataclass, field

from gmb.pipeline.junction_diagnostics import JunctionAnnotator, ModelSpec, introns_of

#: The only junction-manifest role whose libraries may influence selection.
SELECTION_JUNCTION_ROLE = "construction"
#: Version tag recorded with every decision; bump when a formula changes.
POLICY_VERSION = "p4d-1"
CANONICAL_MOTIFS = ("GT/AG", "GC/AG", "AT/AC")

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_CODON: dict = {}
_bases = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _x in enumerate(_bases):
    for _j, _y in enumerate(_bases):
        for _k, _z in enumerate(_bases):
            _CODON[_x + _y + _z] = _aa[16 * _i + 4 * _j + _k]


class EvidenceSelectionError(RuntimeError):
    """Evidence-score selection was requested but cannot run as specified."""


@dataclass(frozen=True)
class EvidenceParams:
    """Parameters of the evidence-score rule (``scoring.evidence_score``)."""

    breadth_saturation_groups: float = 5.0
    unobserved_expressed_weight: float = 0.5
    protein_intron_weight: float = 0.5
    noncanonical_penalty: float = 0.5
    partial_orf_penalty: float = 1.0
    internal_stop_penalty: float = 5.0
    backbone_prior: float = 1.0
    transcript_class_prior: float = 0.25
    replacement_margin: float = 0.5
    keep_supported_intron_J: float = 0.6
    same_coding_locus_fraction: float = 0.5
    min_candidate_cds_bp: int = 150

    @classmethod
    def from_config(cls, ecfg) -> EvidenceParams:
        return cls(**{k: getattr(ecfg, k) for k in cls.__dataclass_fields__})

    def as_dict(self) -> dict:
        return {k: getattr(self, k) for k in self.__dataclass_fields__}


def require_selection_index(index) -> None:
    """Refuse any junction index that is not built from construction libraries only."""
    st = index.stats
    roles = list(st.get("roles") or [])
    if roles != [SELECTION_JUNCTION_ROLE]:
        raise EvidenceSelectionError(
            f"evidence_score selection requires junction roles exactly "
            f"['{SELECTION_JUNCTION_ROLE}'], got {roles}; held-out and auxiliary libraries "
            "must never influence selection")
    bad = [lib.library_id for lib in index.loaded if lib.role != SELECTION_JUNCTION_ROLE]
    if bad:
        raise EvidenceSelectionError(
            f"junction libraries with a role other than '{SELECTION_JUNCTION_ROLE}' were "
            f"loaded for selection: {bad[:5]}")
    if not st.get("libraries_loaded"):
        raise EvidenceSelectionError(
            "evidence_score selection: no construction junction library could be loaded "
            f"({st.get('libraries_by_status')})")
    if st.get("libraries_loaded") != st.get("libraries_selected"):
        failed = [f"{lib.library_id}:{lib.status}" for lib in index.libraries
                  if lib.role == SELECTION_JUNCTION_ROLE and lib.status != "LOADED"]
        raise EvidenceSelectionError(
            "evidence_score selection: construction junction libraries failed to load "
            f"(fix or remove them from the manifest): {failed[:5]}")


def breadth(groups: float, saturation: float) -> float:
    return min(1.0, groups / saturation)


def cds_sequence(genome: dict, chrom: str, strand: str, cds) -> str:
    seq = "".join(genome[chrom][a:b] for a, b in sorted(cds)).upper()
    return seq.translate(_COMP)[::-1] if strand == "-" else seq


def translate(seq: str) -> str:
    return "".join(_CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def shared_cds_bp(a, b) -> int:
    return sum(max(0, min(e1, e2) - max(s1, s2)) for s1, e1 in a for s2, e2 in b)


def same_coding_locus(a, b, fraction: float) -> bool:
    """Both CDSs share at least *fraction* of their own length with the other."""
    shared = shared_cds_bp(a, b)
    return (shared >= fraction * sum(e - s for s, e in a)
            and shared >= fraction * sum(e - s for s, e in b))


@dataclass
class CandidateEvidence:
    """Features of one CDS structure (``key`` = (chrom, strand, CDS segments))."""

    key: tuple
    chrom: str
    strand: str
    cds: tuple
    sources: frozenset
    complete: bool
    introns: list = field(default_factory=list)
    J: float = 0.0
    P: float = 0.0
    motif: float = 0.0
    integrity: float = 0.0
    prior_backbone: float = 0.0
    prior_transcript: float = 0.0
    evidence: float = 0.0
    total: float = 0.0
    n_negative_J: int = 0
    internal_stops: int = 0
    frame_ok: bool = True
    model_class: str = ""
    has_backbone: bool = False

    @property
    def cds_bp(self) -> int:
        return sum(b - a for a, b in self.cds)

    @property
    def chain(self) -> tuple:
        return tuple(introns_of(self.cds))

    def strongly_supported_introns(self, threshold: float) -> set:
        return {(r["intron_start"] - 1, r["intron_end"]) for r in self.introns if r["J"] >= threshold}

    def components(self) -> dict:
        return {"J": self.J, "P": self.P, "motif": self.motif, "integrity": self.integrity,
                "prior_backbone": self.prior_backbone, "prior_transcript": self.prior_transcript,
                "evidence": self.evidence, "total": self.total}


class EvidenceScorer:
    """Computes :class:`CandidateEvidence` for CDS structures (cached per structure).

    *protein_index* is a :class:`gmb.pipeline.junction_support.ProteinIntronIndex`
    (or ``None``: protein terms are then 0). *backbone_sources* and
    *transcript_sources* are the source labels holding the backbone and the
    transcript-evidence roles.
    """

    def __init__(self, index, genome: dict, protein_index, params: EvidenceParams,
                 backbone_sources, transcript_sources):
        require_selection_index(index)
        self.index = index
        self.genome = genome
        self.protein_index = protein_index
        self.params = params
        self.backbone_sources = frozenset(backbone_sources)
        self.transcript_sources = frozenset(transcript_sources)
        self.annotator = JunctionAnnotator(index, genome)
        self._cache: dict = {}
        self.n_scored = 0
        self.seconds = 0.0

    def evaluate(self, chrom: str, strand: str, cds, sources, complete: bool) -> CandidateEvidence:
        cds = tuple(sorted((int(a), int(b)) for a, b in cds))
        sources = frozenset(sources)
        ck = (chrom, strand, cds, sources, bool(complete))
        hit = self._cache.get(ck)
        if hit is not None:
            return hit
        t0 = time.perf_counter()
        p = self.params
        rows, summary = self.annotator.model(ModelSpec(
            model_id="", source=",".join(sorted(sources)), chrom=chrom, strand=strand,
            segments=list(cds)))
        own = [breadth(r["n_groups"], p.breadth_saturation_groups) for r in rows
               if r["status"] == "SUPPORTED"]
        expression = statistics.median(own) if own else 0.0
        ix = self.index
        for r in rows:
            status = r["status"]
            alt = breadth(int(r.get("best_alternative_n_groups") or 0), p.breadth_saturation_groups)
            if status == "SUPPORTED":
                usage = [u for u in (r["donor_usage"], r["acceptor_usage"]) if u is not None]
                u = min(usage) if usage else 1.0
                j = breadth(r["n_groups"], p.breadth_saturation_groups) * u - alt * (1 - u)
            elif status == "ALTERNATIVE_SITE_OBSERVED":
                j = -alt
            elif status == "OPPOSITE_STRAND":
                jid = ix.find(chrom, r["intron_start"] - 1, r["intron_end"])
                j = -breadth(int(ix.j_ngroup[jid]), p.breadth_saturation_groups)
            elif status == "UNOBSERVED_EXPRESSED_LOCUS":
                j = -p.unobserved_expressed_weight * expression
            else:  # MULTIMAP_ONLY, UNOBSERVED_NO_LOCAL_EVIDENCE, UNOBSERVABLE_*
                j = 0.0
            a = c = 0
            if self.protein_index is not None:
                a, c = self.protein_index.intron_agreement(
                    chrom, strand, r["intron_start"] - 1, r["intron_end"])
            prot = p.protein_intron_weight * (2 * a / (a + c) - 1) if a + c else 0.0
            nc = -p.noncanonical_penalty if (r["genome_motif"] not in CANONICAL_MOTIFS and j >= 0) else 0.0
            r.update(J=round(j, 4), P=round(prot, 4), protein_agree=a, protein_conflict=c,
                     motif_penalty=nc)
        seq = cds_sequence(self.genome, chrom, strand, cds) if chrom in self.genome else ""
        aa = translate(seq)
        stops = aa[:-1].count("*") if aa else 0
        frame_ok = len(seq) % 3 == 0
        integrity = 0.0
        if not complete:
            integrity -= p.partial_orf_penalty
        if stops or not frame_ok:
            integrity -= p.internal_stop_penalty
        ev = CandidateEvidence(
            key=(chrom, strand, cds), chrom=chrom, strand=strand, cds=cds, sources=sources,
            complete=bool(complete), introns=rows, internal_stops=stops, frame_ok=frame_ok,
            model_class=summary["model_class"])
        ev.J = round(sum(r["J"] for r in rows), 4)
        ev.P = round(sum(r["P"] for r in rows), 4)
        ev.motif = sum(r["motif_penalty"] for r in rows)
        ev.integrity = integrity
        ev.has_backbone = bool(sources & self.backbone_sources)
        ev.prior_backbone = p.backbone_prior if ev.has_backbone else 0.0
        ev.prior_transcript = p.transcript_class_prior if sources & self.transcript_sources else 0.0
        ev.n_negative_J = sum(1 for r in rows if r["J"] < 0)
        ev.evidence = round(ev.J + ev.P + ev.motif + ev.integrity, 4)
        ev.total = round(ev.evidence + ev.prior_backbone + ev.prior_transcript, 4)
        self._cache[ck] = ev
        self.n_scored += 1
        self.seconds += time.perf_counter() - t0
        return ev


@dataclass
class Decision:
    """Outcome of :func:`decide_primary` for one gene."""

    outcome: str               # replaced | protected_backbone | kept_below_margin | kept_no_eligible
    primary: CandidateEvidence
    chosen: CandidateEvidence
    best_alternative: CandidateEvidence | None = None
    n_pool: int = 0
    n_excluded_stop_or_frame: int = 0
    n_excluded_partial_orf: int = 0
    n_excluded_supported_intron: int = 0
    n_eligible: int = 0

    @property
    def replaced(self) -> bool:
        return self.outcome == "replaced"


def decide_primary(primary: CandidateEvidence, pool: list, params: EvidenceParams) -> Decision:
    """The P4d rule. *pool* holds the other same-coding-locus candidates in locus order."""
    d = Decision(outcome="kept_no_eligible", primary=primary, chosen=primary, n_pool=len(pool))
    intact = [c for c in pool if not c.internal_stops and c.frame_ok]
    d.n_excluded_stop_or_frame = len(pool) - len(intact)
    if primary.has_backbone and primary.n_negative_J == 0 and primary.integrity == 0:
        d.outcome = "protected_backbone"
        return d
    complete = [c for c in intact if c.complete]
    d.n_excluded_partial_orf = len(intact) - len(complete)
    keep = primary.strongly_supported_introns(params.keep_supported_intron_J)
    eligible = [c for c in complete if not (keep - set(c.chain))]
    d.n_excluded_supported_intron = len(complete) - len(eligible)
    d.n_eligible = len(eligible)
    if not eligible:
        return d
    best = eligible[0]
    for c in eligible[1:]:
        if (c.total, c.evidence) > (best.total, best.evidence):
            best = c
    d.best_alternative = best
    m = params.replacement_margin
    if best.total - primary.total >= m and best.evidence - primary.evidence >= m:
        d.outcome, d.chosen = "replaced", best
    else:
        d.outcome = "kept_below_margin"
    return d


def fmt_cds(cds) -> str:
    return ",".join(f"{a + 1}-{b}" for a, b in cds)


def decision_concerns(d: Decision, old_prot_net: int | None, new_prot_net: int | None) -> list:
    """Flags a reviewer should look at; they never change the decision."""
    if not d.replaced:
        return []
    old, new = d.primary, d.chosen
    out = []
    if new.cds_bp < 0.9 * old.cds_bp:
        out.append("cds_shortened_over_10pct")
    if new.cds_bp > 1.1 * old.cds_bp:
        out.append("cds_lengthened_over_10pct")
    o5 = old.cds[0][0] if old.strand == "+" else old.cds[-1][1]
    n5 = new.cds[0][0] if new.strand == "+" else new.cds[-1][1]
    o3 = old.cds[-1][1] if old.strand == "+" else old.cds[0][0]
    n3 = new.cds[-1][1] if new.strand == "+" else new.cds[0][0]
    if o5 != n5:
        out.append("start_codon_changed")
    if o3 != n3:
        out.append("stop_codon_changed")
    if old_prot_net is not None and new_prot_net is not None and new_prot_net < old_prot_net:
        out.append("protein_compatibility_worse")
    return out
