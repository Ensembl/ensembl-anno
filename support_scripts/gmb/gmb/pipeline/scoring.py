#!/usr/bin/env python3
"""Isoform scoring and selection for Gene Model Builder.

Replaces the previous hard-coded analyze_locus() with a configurable
scoring function.  Default parameters are tuned for fungal genomes.

All thresholds are driven by PipelineConfig (YAML).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from gmb.pipeline.config import PipelineConfig

import pandas as pd

from gmb.pipeline.annotate_cds_utrs import check_splice_sites
from gmb.pipeline.canonical_evidence import (
    EVIDENCE_CLASS_BACKBONE,
    EVIDENCE_CLASS_LONG_READ,
    EVIDENCE_CLASS_PROTEIN_ALIGNMENT,
    EVIDENCE_CLASS_SHORT_READ,
    EvidenceRoles,
)

# ---------------------------------------------------------------------------
# Intron-chain utility (reused from gene_model_builder)
# ---------------------------------------------------------------------------


def _get_intron_chain(exon_df: pd.DataFrame) -> str:
    """Return a string signature of the intron chain for one transcript."""
    exon_df = exon_df.sort_values("Start")
    if len(exon_df) < 2:
        return "single-exon"
    ends = exon_df["End"].tolist()
    starts = exon_df["Start"].tolist()
    return ",".join(f"{ends[i]}-{starts[i + 1]}" for i in range(len(starts) - 1))


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


def rescue_enabled(scoring_config) -> bool:
    """Whether backbone intron rescue is active for this run.

    Prefers the value the builder resolved through the applicability gate. Falls
    back to interpreting the configured mode directly when nothing resolved it --
    "auto" is treated as OFF in that case, because the gate needs run-level
    evidence that a bare select_isoforms() call does not have, and OFF is the
    safe interpretation.
    """
    resolved = getattr(scoring_config, "backbone_intron_rescue_resolved", None)
    if resolved is not None:
        return bool(resolved)
    from gmb.pipeline.applicability import normalise_rescue_mode
    return normalise_rescue_mode(
        getattr(scoring_config, "backbone_intron_rescue", "off")) == "on"


def weights_for_role(weights, role: str) -> float:
    """Base weight for one evidence role.

    The single place that maps a resolved role onto a configured numeric weight.
    An unrecognised role falls back to `weights.unknown` rather than to a
    hard-coded literal, so adding a role to `canonical_evidence` cannot silently
    change scoring.
    """
    return {
        EVIDENCE_CLASS_BACKBONE: weights.backbone,
        EVIDENCE_CLASS_SHORT_READ: weights.short_read,
        EVIDENCE_CLASS_LONG_READ: weights.long_read,
        EVIDENCE_CLASS_PROTEIN_ALIGNMENT: weights.protein_alignment,
    }.get(role, weights.unknown)


def score_model(
    model: dict,
    config: PipelineConfig,
    protein_supported_tids: set[str],
    genome: dict[str, str] | None = None,
    roles: "EvidenceRoles | None" = None,
) -> float:
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
    roles : EvidenceRoles or None
        Resolved evidence roles. Built from `config.scoring` when omitted;
        callers in a loop should build it once and pass it in.

    Returns
    -------
    float
        Composite score.
    """
    scfg = config.scoring
    score = 0.0

    # Base evidence weight, resolved through evidence ROLES -- never through a
    # literal tool name. A source listed under scoring.shortread_labels gets the
    # short_read weight whatever it is called, so StringTie2 / IsoQuant / FLAIR /
    # AssemblerX all behave identically to StringTie.
    sources = set(model.get("combined_evidence", model["source"]).split(","))
    weights = scfg.weights
    if roles is None:
        roles = EvidenceRoles.from_config(scfg)
    for s in sources:
        score += weights_for_role(weights, roles.role_of(s))

    # Multi-source bonus — uses raw named-source count, not biological evidence classes.
    # canonical_selection uses evidence classes for breadth ranking; these are different
    # questions (model retention here vs canonical ranking there), so the difference is
    # intentional, not a bug.
    if len(sources) > 1:
        score += scfg.multi_source_bonus * (len(sources) - 1)

    # Protein support bonus
    if model["id"] in protein_supported_tids:
        score += scfg.protein_overlap_bonus
        model["protein_support"] = True
    else:
        model["protein_support"] = model.get("protein_support", False)

    # Protein Validation Score (if run)
    if "protein_coding_score" in model:
        val_cfg = config.protein_validation
        if val_cfg.enabled and val_cfg.policy in ("penalize", "penalise"):
            # E.g. penalty if it falls below min_score
            if model["protein_coding_score"] < val_cfg.min_score:
                score -= getattr(val_cfg, "penalty", 5.0)
        elif val_cfg.enabled and val_cfg.policy == "bonus":
            score += model["protein_coding_score"]

    # Splice-site penalty (only if genome provided and multi-exon)
    if genome and model["exon_count"] > 1:
        chrom = model["chrom"]
        if chrom in genome:
            exons = sorted(zip(model["df"]["Start"].values, model["df"]["End"].values))
            splice = check_splice_sites(exons, model["strand"], genome[chrom])
            n_noncanonical = sum(1 for s in splice if s["class"] != "canonical")
            score -= scfg.noncanonical_splice_penalty * n_noncanonical

    return score


# ---------------------------------------------------------------------------
# Isoform selection
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Ranking hierarchy
# ---------------------------------------------------------------------------
# One ordered definition of "which evidence pattern beats which". Higher tier
# wins; the numeric score from score_model() breaks ties WITHIN a tier.
#
# Protein support is deliberately NOT a tier. It influences the score (via
# scoring.protein_overlap_bonus) and gates retention, but it does not by itself
# outrank a structurally better-supported model. `selection_reason_of` may still
# report it as the most informative descriptor for a winner that has no
# structural corroboration.
TIER_BACKBONE_SHORTREAD_AGREEMENT = 3
TIER_MULTI_SOURCE_AGREEMENT = 2
TIER_BACKBONE_INTRON_RESCUE = 1
TIER_SINGLE_SOURCE = 0


def rank_tier(s) -> int:
    """Ranking tier of a merged structure (higher wins)."""
    if s["backbone_shortread_agreement"]:
        return TIER_BACKBONE_SHORTREAD_AGREEMENT
    if s["n_structural_support_sources"] > 1:
        return TIER_MULTI_SOURCE_AGREEMENT
    if s["backbone_intron_rescue"]:
        return TIER_BACKBONE_INTRON_RESCUE
    return TIER_SINGLE_SOURCE


def selection_reason_of(s) -> str:
    """Why this structure ranked where it did.

    Mirrors `rank_tier`, then falls back to finer descriptors (protein support,
    long-read role) that describe the winner without affecting its tier.
    """
    if s["longread_demoted"]:
        return "longread_demoted_support_only"
    tier = rank_tier(s)
    if tier == TIER_BACKBONE_SHORTREAD_AGREEMENT:
        return "backbone_shortread_agreement"
    if tier == TIER_MULTI_SOURCE_AGREEMENT:
        return "multi_source_structural_agreement"
    if tier == TIER_BACKBONE_INTRON_RESCUE:
        return "backbone_intron_rescue"
    if s["protein_cds_span"]:
        return "protein_cds_span_support"
    if s["is_longread_only"]:
        return ("single_exon_longread" if s["rep"]["exon_count"] == 1
                else "longread_only_locus")
    return "best_single_source"


def select_isoforms(
    locus_df: pd.DataFrame,
    config: PipelineConfig,
    protein_supported_tids: set[str],
    genome: dict[str, str] | None = None,
    protein_support_sources: dict[str, set[str]] | None = None,
    protein_cds_span_tids: set[str] | None = None,
    candidate_cds: dict[str, list] | None = None,
    canonical_intron_tids: set[str] | None = None,
    junction_support: dict[str, tuple[int, int]] | None = None,
    protein_compatibility: dict[str, tuple[int, int]] | None = None,
    complete_orf_tids: set[str] | None = None,
    excluded: list[dict] | None = None,
    evidence_scorer=None,
) -> list[list[dict]]:
    """Score and select isoforms for a single locus.

    Parameters
    ----------
    locus_df : pd.DataFrame
        Exon rows for all models in this locus.
    config : PipelineConfig
    protein_supported_tids : set of str
    genome : dict or None
    protein_support_sources : dict or None
        ``{candidate_transcript_id: {"<protein track>", ...}}`` naming *which*
        protein-alignment tracks support each candidate.  Attribution only: it is
        recorded on the selected model as ``protein_evidence`` and deliberately
        never enters ``sources`` / ``combined_evidence``, since protein
        alignments are supporting evidence rather than candidate models.
    protein_cds_span_tids : set of str or None
        Candidates whose protein support is *CDS-span compatible*: the alignment
        falls inside the candidate's transcript span and the candidate has a CDS.
        Used as the strong protein signal when
        ``scoring.protein_support_mode == "cds_span_compatible"``.
    candidate_cds : dict or None
        ``{candidate_transcript_id: [(start, end), ...]}`` CDS intervals from the
        ORF stage.  Required by ``scoring.backbone_intron_rescue``.
    canonical_intron_tids : set of str or None
        Candidates whose introns are all canonical.  ``None`` means the caller
        supplied no such information and the check is skipped; an empty set means
        nothing qualifies.  Used by ``backbone_intron_rescue`` as a safety guard.
    junction_support : dict or None
        ``{candidate_transcript_id: (n_introns, n_introns_not_observed)}`` from
        :func:`gmb.pipeline.junction_support.unsupported_intron_counts`. Recorded on
        every selected model; selection-affecting only with
        ``scoring.primary_selection: junction_supported``.
    protein_compatibility : dict or None
        ``{candidate_transcript_id: (compatible, incompatible)}`` spliced protein
        alignments over the candidate's CDS. Same use as *junction_support*.
    complete_orf_tids : set of str or None
        Candidates whose ORF has both an ATG start and a stop codon.
    excluded : list or None
        If given, a record is appended for every retained structure discarded by
        ``scoring.isoform_cds_overlap: drop`` (it would have joined a gene without
        sharing coding sequence with that gene's primary).
    evidence_scorer : gmb.pipeline.evidence_features.EvidenceScorer or None
        Required by, and used only with, ``scoring.primary_selection:
        evidence_score``.

    Returns
    -------
    list of list of dict
        Each inner list is a gene sub-cluster of model dicts (selected
        isoforms), sorted by genomic start position.
    """
    scfg = config.scoring

    # Protein evidence plays three separate roles and they are never conflated:
    #   ranking     -- may decide which of several structures wins
    #   retention   -- may keep a structure that would otherwise be dropped
    #   attribution -- reported, never selection-affecting
    #
    # protein_support_mode selects which signal fills ranking+retention:
    #   "positional"          -- any same-strand overlap (legacy)
    #   "cds_span_compatible" -- the alignment must lie inside the candidate's
    #                            transcript span and the candidate must have a CDS
    # Positional support is recorded for attribution in both modes.
    roles = EvidenceRoles.from_config(scfg)
    guard_on = getattr(scfg, "longread_structural_guard", False)
    span_tids = protein_cds_span_tids or set()
    strong_protein_tids = (
        span_tids if getattr(scfg, "protein_support_mode", "positional")
        == "cds_span_compatible" else protein_supported_tids
    )

    # Build model dicts
    models = []
    for (source, tid), grp in locus_df.groupby(["Source", "transcript_id"]):
        chain = _get_intron_chain(grp)
        models.append(
            {
                "id": tid,
                "source": source,
                "chrom": grp["Chromosome"].iloc[0],
                "strand": grp["Strand"].iloc[0],
                "intron_chain": chain,
                "protein_support": tid in strong_protein_tids,
                "protein_positional_support": tid in protein_supported_tids,
                "protein_cds_span_compatible": tid in span_tids,
                "protein_sources": set(
                    (protein_support_sources or {}).get(tid, ())
                ),
                "df": grp,
                "start": grp["Start"].min(),
                "end": grp["End"].max(),
                "exon_count": len(grp),
                "combined_evidence": (
                    grp["combined_evidence"].iloc[0]
                    if "combined_evidence" in grp.columns
                    else source
                ),
                "protein_coding_score": (
                    grp["protein_coding_score"].iloc[0]
                    if "protein_coding_score" in grp.columns
                    else 0.0
                ),
            }
        )

    if not models:
        return []

    # Merge identical structures across sources
    merged = {}
    for m in models:
        if m["intron_chain"] == "single-exon":
            key = f"{m['chrom']}:{m['strand']}:{m['start']}-{m['end']}"
        else:
            key = f"{m['chrom']}:{m['strand']}:{m['intron_chain']}"

        if key not in merged:
            merged[key] = {
                "sources": set(),
                "protein_support": False,
                "protein_sources": set(),
                "protein_positional": False,
                "protein_cds_span": False,
                "rep": m,
                "score": 0.0,
            }
        s = merged[key]
        s["sources"].add(m["source"])
        s["protein_positional"] |= m.get("protein_positional_support", False)
        s["protein_cds_span"] |= m.get("protein_cds_span_compatible", False)
        # Attribution only -- kept out of s["sources"] on purpose (see docstring).
        s["protein_sources"].update(m.get("protein_sources", ()))
        if m["protein_support"]:
            s["protein_support"] = True
            if not s["rep"]["protein_support"]:
                s["rep"] = m
        # Propagate validation score if available
        if "protein_coding_score" in m:
            s["rep"]["protein_coding_score"] = m["protein_coding_score"]

    # Score each merged structure
    for _key, s in merged.items():
        rep = s["rep"]
        rep["combined_evidence"] = ",".join(sorted(s["sources"]))
        # Protein-alignment attribution: union across every model merged into
        # this structure, recorded separately from combined_evidence.
        rep["protein_evidence"] = ",".join(sorted(s["protein_sources"]))
        if s["protein_support"]:
            rep["protein_support"] = True
        s["score"] = score_model(rep, config, strong_protein_tids, genome, roles=roles)

        # ---- structural corroboration (role-based, no literal tool names) ----
        src_roles = {src: roles.role_of(src) for src in s["sources"]}
        s["has_backbone"] = EVIDENCE_CLASS_BACKBONE in src_roles.values()
        s["has_shortread"] = EVIDENCE_CLASS_SHORT_READ in src_roles.values()
        s["is_longread_only"] = set(src_roles.values()) == {EVIDENCE_CLASS_LONG_READ}
        s["backbone_shortread_agreement"] = s["has_backbone"] and s["has_shortread"]
        # Independent structural sources. The long-read role is excluded while
        # the guard is active: that track contributes locus extent rather than
        # splice structure, so counting it would overstate agreement.
        corroborating = {
            src for src, role in src_roles.items()
            if not (guard_on and role == EVIDENCE_CLASS_LONG_READ)
        }
        s["n_structural_support_sources"] = len(corroborating)
        s["structural_support_sources"] = ",".join(sorted(s["sources"]))
        rep["structural_support_sources"] = s["structural_support_sources"]
        rep["n_structural_support_sources"] = s["n_structural_support_sources"]
        rep["backbone_shortread_agreement"] = s["backbone_shortread_agreement"]
        rep["protein_support_strength"] = (
            "strong" if s["protein_cds_span"]
            else "weak" if s["protein_positional"]
            else "none"
        )

    # ---- backbone intron rescue ----
    # An ab initio backbone that calls a single coding exon where an assembled
    # transcript shows a canonical spliced CDS recovering MORE coding sequence
    # has under-called introns. The assembly is the direct observation, so it
    # takes structural priority and is exempt from the single-source protein
    # gate. The "longer CDS" requirement is what separates a genuine collapse
    # (the true CDS extends past the truncated backbone call) from a fragmentary
    # assembly overlapping a correctly-called single-exon gene.
    cds_map = candidate_cds or {}
    # None => caller supplied no canonical-intron information, so skip that
    # check; an empty set => supplied, and nothing qualifies.
    canon_tids = canonical_intron_tids
    for s in merged.values():
        cds = sorted(cds_map.get(s["rep"]["id"], []) or [])
        s["cds"] = cds
        s["cds_bp"] = sum(e - st for st, e in cds)
        s["n_cds_exons"] = len(cds)
        # Kept under its own key: the retention floor reads rep["cds_bp"].
        s["rep"]["candidate_cds_bp"] = s["cds_bp"]
        s["rep"]["candidate_cds"] = cds
        s["rep"]["has_backbone"] = s["has_backbone"]
        s["backbone_intron_rescue"] = False
        s["rep"]["backbone_intron_rescue"] = False

    if rescue_enabled(scfg) and cds_map:
        collapsed = [s for s in merged.values()
                     if s["has_backbone"] and s["n_cds_exons"] <= 1 and s["cds"]]
        if collapsed:
            for s in merged.values():
                if not s["has_shortread"] or s["n_cds_exons"] < 2:
                    continue
                if canon_tids is not None and s["rep"]["id"] not in canon_tids:
                    continue
                cs, ce = s["cds"][0][0], s["cds"][-1][1]
                for b in collapsed:
                    bs, be = b["cds"][0][0], b["cds"][-1][1]
                    if b["rep"]["chrom"] != s["rep"]["chrom"]:
                        continue
                    if b["rep"]["strand"] != s["rep"]["strand"]:
                        continue
                    if min(ce, be) - max(cs, bs) <= 0:
                        continue
                    if s["cds_bp"] <= b["cds_bp"]:
                        continue
                    # The collapsed call must be a *part of* the recovered CDS:
                    # most of the backbone's coding bases should fall inside the
                    # assembly's CDS. This separates a genuine collapse from a
                    # neighbouring spliced gene that merely overlaps.
                    covered = 0
                    for bstart, bend in b["cds"]:
                        for astart, aend in s["cds"]:
                            covered += max(0, min(bend, aend) - max(bstart, astart))
                    if b["cds_bp"] and covered * 2 < b["cds_bp"]:
                        continue
                    s["backbone_intron_rescue"] = True
                    s["rep"]["backbone_intron_rescue"] = True
                    break

    # ---- retention gate ----
    # RETENTION asks only "is this structure supported enough to exist?".
    # It never decides which structure wins -- that is RANKING, below.
    val_cfg = config.protein_validation

    def structurally_valid(s):
        """Biological validity floor every model must clear."""
        if s["rep"]["strand"] not in ("+", "-"):
            return False
        if val_cfg.enabled and val_cfg.policy == "drop":
            score = s["rep"].get("protein_coding_score")
            if score is not None and score < val_cfg.min_score:
                return False
        # The candidate's own CDS (from the ORF stage). Before this read s["cds_bp"]
        # it read rep["cds_bp"], which nothing sets, so the floor never applied.
        cds_bp = s.get("cds_bp", 0)
        if cds_bp > 0 and cds_bp < scfg.min_cds_bp:
            return False
        return True

    def passes_gate(s):
        # A backbone-intron rescue is corroborated by the backbone agreeing the
        # locus is coding, so it is not "single-source unsupported".
        supported = s["protein_support"] or s["backbone_intron_rescue"]
        backbone_protected = s["has_backbone"] and scfg.keep_backbone_without_support
        keep = False
        if supported or backbone_protected:
            keep = True
        elif len(s["sources"]) > 1:
            keep = True
        elif s["has_shortread"] and s["rep"]["exon_count"] > 1:
            keep = not scfg.require_protein_support_for_single_source
        elif (
            scfg.fungal_single_exon_mode
            and s["rep"]["exon_count"] == 1
            and s["score"] >= scfg.min_alternate_score
        ):
            keep = True

        # Single-exon models require support when configured; a backbone
        # single-exon gene is exempt when the operator opted in to keeping the
        # backbone without support.
        if (
            keep
            and scfg.require_support_for_single_exon
            and s["rep"]["exon_count"] == 1
            and not backbone_protected
            and not supported
            and len(s["sources"]) < 2
        ):
            keep = False
        return keep

    candidates = [s for s in merged.values()
                  if structurally_valid(s) and passes_gate(s)]

    if not candidates:
        return []

    # ---- long-read structural guard (RANKING) ----
    # A long-read-only structure does not take primary structural priority where
    # any multi-exon structure from another role exists. It is not removed: it
    # stays eligible as an alternate isoform and keeps its attribution.
    # `longread_disposition` is the unconditional form; the guard is the
    # conditional one. support_only demotes every long-read-only structure;
    # the guard demotes one only where another multi-exon structure survives.
    disposition = str(getattr(scfg, "longread_disposition", "primary_structural"))
    demote_all_longread = disposition == "support_only"
    if guard_on:
        has_other_multiexon = any(
            (not s["is_longread_only"]) and s["rep"]["exon_count"] > 1
            for s in candidates
        )
    else:
        has_other_multiexon = False
    for s in candidates:
        s["longread_demoted"] = bool(
            s["is_longread_only"]
            and (demote_all_longread or (guard_on and has_other_multiexon)))
        s["rep"]["longread_structural_role"] = (
            "support_only" if s["longread_demoted"]
            else "primary" if s["is_longread_only"]
            else "not_longread"
        )

    # ---- ranking ----
    # RANKING decides which surviving structure wins the locus. The hierarchy is
    # a single ordered list of tiers; the numeric score only breaks ties within
    # a tier. Keep TIERS, rank_tier() and selection_reason_of() in step -- they
    # are the same hierarchy expressed for sorting and for reporting.
    #
    #   RANKING HIERARCHY (strongest first) -- see rank_tier()
    #     1. backbone + assembled-transcript structural agreement
    #     2. two or more independent structural sources agree
    #     3. backbone intron rescue (spliced assembly over a collapsed backbone)
    #     4. numeric score (protein support contributes here, not as a tier)
    #   Demoted long-read structures always sort after everything else.
    if getattr(scfg, "structural_corroboration", False):
        candidates.sort(key=lambda s: (s["longread_demoted"], -rank_tier(s), -s["score"]))
    else:
        # Legacy ranking: numeric score alone, with the two optional policies
        # (guard, rescue) still able to move a structure if enabled.
        candidates.sort(key=lambda s: (s["longread_demoted"],
                                       not s["backbone_intron_rescue"], -s["score"]))

    for s in candidates:
        reason = selection_reason_of(s)
        s["rep"]["selection_reason"] = reason


    same_gene_ovlp_thresh = scfg.same_gene_overlap_threshold

    def is_same_gene(m1, m2):
        if m1["chrom"] != m2["chrom"] or m1["strand"] != m2["strand"]:
            return False
        if m1["intron_chain"] != "single-exon" and m2["intron_chain"] != "single-exon":
            i1 = set(m1["intron_chain"].split(","))
            i2 = set(m2["intron_chain"].split(","))
            if i1.intersection(i2):
                return True
        overlap = min(m1["end"], m2["end"]) - max(m1["start"], m2["start"])
        if overlap > 0:
            len1 = m1["end"] - m1["start"]
            len2 = m2["end"] - m2["start"]
            if (overlap / min(len1, len2)) > same_gene_ovlp_thresh:
                return True
        return False

    # scoring.isoform_cds_overlap: whether an isoform must share coding sequence
    # with its gene's primary (Ensembl groups transcripts into a gene by shared
    # coding bases). "off" keeps the structural test alone, under which UTR and
    # read-through overlaps make a different ORF an "isoform" of the gene.
    cds_mode = getattr(scfg, "isoform_cds_overlap", "off")

    def shares_coding_base(m1, m2):
        c1, c2 = m1.get("candidate_cds") or [], m2.get("candidate_cds") or []
        if not c1 or not c2:
            return True  # no ORF on one side: the structural test decides
        return any(min(e1, e2) > max(s1, s2) for s1, e1 in c1 for s2, e2 in c2)

    genes = []

    for s in candidates:
        r = s["rep"]
        r["score"] = s["score"]
        # Try to assign to an existing gene sub-cluster
        found = -1
        conflicting_primary = None
        for i, g_isoforms in enumerate(genes):
            if is_same_gene(r, g_isoforms[0]):
                if cds_mode != "off" and not shares_coding_base(r, g_isoforms[0]):
                    conflicting_primary = conflicting_primary or g_isoforms[0]
                    continue
                found = i
                break

        if found == -1:
            if conflicting_primary is not None and cds_mode == "drop" and not r.get("has_backbone"):
                if excluded is not None:
                    excluded.append({
                        "candidate_transcript_id": r["id"],
                        "evidence_sources": r.get("combined_evidence", r.get("source", "")),
                        "chrom": r["chrom"], "strand": r["strand"],
                        "start": int(r["start"]), "end": int(r["end"]),
                        "cds_bp": r.get("candidate_cds_bp", 0),
                        "gmb_score": s["score"],
                        "would_join_primary": conflicting_primary["id"],
                        "exclusion_reason": "isoform_cds_overlap_no_shared_coding_base",
                    })
                continue
            r["is_primary"] = True
            genes.append([r])
        else:
            g_isoforms = genes[found]
            primary = g_isoforms[0]
            if len(g_isoforms) < scfg.max_isoforms_per_locus:
                if s["score"] >= scfg.min_alternate_score:
                    if r["intron_chain"] != primary["intron_chain"]:
                        r["is_primary"] = False
                        g_isoforms.append(r)

    junctions = junction_support or {}
    compat = protein_compatibility or {}
    for g in genes:
        for m in g:
            m["n_introns"], m["introns_without_transcript_support"] = junctions.get(
                m["id"], (m["exon_count"] - 1, None))
            m["protein_alignments_compatible"], m["protein_alignments_incompatible"] = (
                compat.get(m["id"], (0, 0)))

    primary_mode = getattr(scfg, "primary_selection", "score")
    if primary_mode in ("junction_supported", "evidence_score"):
        complete = complete_orf_tids or set()
        frac = getattr(scfg, "junction_primary_min_cds_fraction", 0.9)
        for g in genes:
            _promote_junction_supported(g, complete, frac)

    if primary_mode == "evidence_score":
        if evidence_scorer is None:
            raise ValueError(
                "scoring.primary_selection: evidence_score needs STAR junction evidence "
                "(gmb-build --star-junctions MANIFEST); it never falls back to another policy")
        _promote_evidence_supported(
            genes, locus_df, models, merged, candidates, cds_map, complete_orf_tids or set(),
            evidence_scorer, scfg, is_same_gene, shares_coding_base, cds_mode, junctions, compat,
            config.protein_validation)

    # Sort genes by start coordinate
    genes.sort(key=lambda g: min(m["start"] for m in g))
    return genes


def _promote_junction_supported(gene: list[dict], complete_orf_tids: set, min_cds_fraction: float):
    """Replace a backbone primary that has transcript-unsupported introns (in place).

    The backbone wins every locus on score, so an assembled structure that RNA-seq
    supports better can only ever be an alternate. Here the primary is swapped for
    the alternate with the fewest unobserved introns when that alternate:
      * is not itself a backbone structure,
      * has strictly fewer introns that no transcriptomic model observes,
      * has a complete ORF (ATG start and stop codon),
      * keeps at least ``min_cds_fraction`` of the backbone's CDS length,
      * shares at least half of its CDS bases with the backbone CDS and covers at
        least half of the backbone CDS (the same coding locus, not a neighbour
        reached through a read-through span), and
      * is contradicted by no more spliced protein alignments, net of agreeing
        ones, than the backbone is.
    Ties go to better protein compatibility, then the longer CDS. The backbone
    model is kept as an alternate, so no structure is lost.
    """
    primary = gene[0]
    p_bad = primary.get("introns_without_transcript_support")
    if not primary.get("has_backbone") or not p_bad:
        return

    def net(m):
        return m["protein_alignments_compatible"] - m["protein_alignments_incompatible"]

    p_cds = primary.get("candidate_cds_bp") or 0
    p_ivs = primary.get("candidate_cds") or []

    def same_coding_locus(m):
        shared = sum(max(0, min(e1, e2) - max(s1, s2))
                     for s1, e1 in p_ivs for s2, e2 in (m.get("candidate_cds") or []))
        m_cds = m.get("candidate_cds_bp") or 0
        return p_cds > 0 and m_cds > 0 and 2 * shared >= p_cds and 2 * shared >= m_cds

    eligible = [
        m for m in gene[1:]
        if not m.get("has_backbone")
        and m.get("introns_without_transcript_support") is not None
        and m["introns_without_transcript_support"] < p_bad
        and m["id"] in complete_orf_tids
        and (m.get("candidate_cds_bp") or 0) >= min_cds_fraction * p_cds
        and net(m) >= net(primary)
        and same_coding_locus(m)
    ]
    if not eligible:
        return
    best = min(eligible, key=lambda m: (m["introns_without_transcript_support"],
                                        -net(m), -(m.get("candidate_cds_bp") or 0)))
    gene.remove(best)
    best["is_primary"] = True
    best["selection_reason"] = "transcript_junction_support"
    primary["is_primary"] = False
    gene.insert(0, best)


def _merged_key(m: dict) -> str:
    """The structure key select_isoforms merges identical models under."""
    if m["intron_chain"] == "single-exon":
        return f"{m['chrom']}:{m['strand']}:{m['start']}-{m['end']}"
    return f"{m['chrom']}:{m['strand']}:{m['intron_chain']}"


def _promote_evidence_supported(genes, locus_df, models, merged, candidates, cds_map,
                                complete_orf_tids, scorer, scfg, is_same_gene,
                                shares_coding_base, cds_mode, junctions, compat, val_cfg=None):
    """``primary_selection: evidence_score``: re-choose each gene's primary (in place).

    The candidate pool of a gene is every CDS structure of the locus (whether or not
    GMB retained it as an isoform) that shares at least ``same_coding_locus_fraction``
    of its CDS with the current primary in both directions, has at least
    ``min_candidate_cds_bp`` of CDS and passes the structural validity floor, and is
    not used by any other gene of the locus. The P4d rule
    (:func:`gmb.pipeline.evidence_features.decide_primary`) decides. A replaced
    primary is kept as an alternate only if it passes the normal isoform admission
    rules against the new primary; gene grouping is otherwise unchanged. Every gene
    gets an ``evidence_selection`` record on its primary.
    """
    from gmb.pipeline.evidence_features import (
        POLICY_VERSION,
        decide_primary,
        decision_concerns,
        fmt_cds,
        same_coding_locus,
    )

    params = scorer.params
    min_bp = max(params.min_candidate_cds_bp, getattr(scfg, "min_cds_bp", 0) or 0)

    # CDS structures of every candidate transcript in the locus, in locus order.
    order = locus_df.drop_duplicates("transcript_id")[
        ["transcript_id", "Source", "Chromosome", "Strand"]]
    structs: dict = {}
    tid_key: dict = {}
    for tid, src, chrom, strand in order.itertuples(index=False):
        cds = cds_map.get(tid)
        if not cds or strand not in ("+", "-"):
            continue
        key = (str(chrom), str(strand), tuple(sorted((int(a), int(b)) for a, b in cds)))
        st = structs.setdefault(key, {"tids": [], "sources": set(), "complete": False})
        st["tids"].append(tid)
        st["sources"].add(str(src))
        st["complete"] |= tid in complete_orf_tids
        tid_key[tid] = key
    models_by_id = {m["id"]: m for m in models}

    def valid(key):
        cds = key[2]
        if sum(b - a for a, b in cds) < min_bp:
            return False
        if val_cfg is not None and val_cfg.enabled and val_cfg.policy == "drop":
            for t in structs[key]["tids"]:
                sc = models_by_id.get(t, {}).get("protein_coding_score")
                if sc is not None and sc < val_cfg.min_score:
                    return False
        return True

    def evaluate(key):
        st = structs[key]
        return scorer.evaluate(key[0], key[1], key[2], st["sources"], st["complete"])

    # The GFF3 validation limits the build applies later (merged exon length and,
    # when enabled, transcript span): a transcript that would be dropped there
    # cannot become a primary.
    max_exon, max_span = getattr(scorer, "structural_limits", None) or (None, None)

    def fits(tid):
        m = models_by_id.get(tid)
        if m is None:
            return False
        if max_exon is None and max_span is None:
            return True
        ivs = sorted(zip(m["df"]["Start"].astype(int), m["df"]["End"].astype(int)))
        union = []
        for a, b in ivs:
            if union and a <= union[-1][1]:
                union[-1][1] = max(union[-1][1], b)
            else:
                union.append([a, b])
        if max_exon is not None and max(b - a for a, b in union) > max_exon:
            return False
        return max_span is None or union[-1][1] - union[0][0] <= max_span

    def representative(key, gene):
        """The transcript that carries a chosen CDS: prefer what GMB already kept."""
        tids = [t for t in structs[key]["tids"] if fits(t)]
        for m in gene:
            if m["id"] in tids:
                return m, "retained_isoform"
        for s in candidates:
            if s["rep"]["id"] in tids:
                return s["rep"], "admitted_structure"
        for s in merged.values():
            if s["rep"]["id"] in tids:
                return s["rep"], "scored_structure"
        m = models_by_id[tids[0]]
        s = merged[_merged_key(m)]
        r = dict(s["rep"])
        for k in ("id", "source", "df", "start", "end", "exon_count", "intron_chain",
                  "protein_coding_score"):
            if k in m:
                r[k] = m[k]
        r["candidate_cds"] = sorted(cds_map.get(m["id"]) or [])
        r["candidate_cds_bp"] = sum(b - a for a, b in r["candidate_cds"])
        r.setdefault("longread_structural_role",
                     "primary" if s.get("is_longread_only") else "not_longread")
        return r, "member_transcript"

    def net(tid):
        c = compat.get(tid)
        return None if c is None else c[0] - c[1]

    used: dict = {}
    for gi, g in enumerate(genes):
        for m in g:
            k = tid_key.get(m["id"])
            if k is not None:
                used.setdefault(k, set()).add(gi)

    for gi, g in enumerate(genes):
        primary = g[0]
        pk = tid_key.get(primary["id"])
        base = {"policy_version": POLICY_VERSION, "chrom": primary["chrom"], "strand": primary["strand"],
                "old_candidate": primary["id"]}
        if pk is None:
            primary["evidence_selection"] = dict(base, outcome="no_cds_primary")
            continue
        P = evaluate(pk)
        pool, n_other_gene, n_invalid_transcript = [], 0, 0
        for k in structs:
            if k == pk or k[0] != pk[0] or k[1] != pk[1]:
                continue
            if not same_coding_locus(P.cds, k[2], params.same_coding_locus_fraction):
                continue
            if used.get(k, set()) - {gi}:
                n_other_gene += 1
                continue
            if not valid(k):
                continue
            if not any(fits(t) for t in structs[k]["tids"]):
                n_invalid_transcript += 1
                continue
            pool.append(evaluate(k))
        d = decide_primary(P, pool, params)
        rec = dict(base, outcome=d.outcome, protected_backbone=d.outcome == "protected_backbone",
                   old_structure=fmt_cds(P.cds), old_sources=",".join(sorted(P.sources)),
                   n_pool=d.n_pool, n_excluded_other_gene=n_other_gene,
                   n_excluded_invalid_transcript=n_invalid_transcript,
                   n_excluded_stop_or_frame=d.n_excluded_stop_or_frame,
                   n_excluded_partial_orf=d.n_excluded_partial_orf,
                   n_excluded_supported_intron=d.n_excluded_supported_intron,
                   n_eligible=d.n_eligible, margin=params.replacement_margin,
                   protected_introns=fmt_cds(sorted(P.strongly_supported_introns(
                       params.keep_supported_intron_J))),
                   **{f"old_{k}": v for k, v in P.components().items()})
        alt = d.best_alternative
        if alt is not None:
            rec.update(best_alternative=fmt_cds(alt.cds),
                       best_alternative_delta_total=round(alt.total - P.total, 4),
                       best_alternative_delta_evidence=round(alt.evidence - P.evidence, 4))
        if not d.replaced:
            primary["evidence_selection"] = rec
            continue
        new_m, tier = representative(d.chosen.key, g)
        old_chain, new_chain = set(P.chain), set(d.chosen.chain)
        kept, dropped = [], []
        for m in [x for x in g if x is not new_m]:
            if m["intron_chain"] == new_m["intron_chain"]:
                reason = "same_intron_chain_as_new_primary"
            elif len(kept) + 1 >= scfg.max_isoforms_per_locus:
                reason = "max_isoforms_per_locus"
            elif not is_same_gene(m, new_m):
                reason = "not_same_gene_as_new_primary"
            elif cds_mode != "off" and not shares_coding_base(m, new_m):
                reason = "no_shared_coding_base"
            elif (m.get("score") or 0.0) < scfg.min_alternate_score:
                reason = "below_min_alternate_score"
            else:
                reason = None
            if reason:
                dropped.append(f"{m['id']}:{reason}")
            else:
                kept.append(m)
        for m in kept:
            m["is_primary"] = False
        new_m["is_primary"] = True
        new_m["selection_reason"] = "evidence_score"
        if "n_introns" not in new_m:
            new_m["n_introns"], new_m["introns_without_transcript_support"] = junctions.get(
                new_m["id"], (new_m["exon_count"] - 1, None))
        if "protein_alignments_compatible" not in new_m:
            new_m["protein_alignments_compatible"], new_m["protein_alignments_incompatible"] = (
                compat.get(new_m["id"], (0, 0)))
        g[:] = [new_m, *kept]
        concerns = decision_concerns(d, net(primary["id"]), net(new_m["id"]))
        if tier != "retained_isoform":
            concerns.append("new_primary_was_not_a_retained_isoform")
        if not any(m is primary for m in kept):
            concerns.append("old_primary_not_kept_as_isoform")
        rec.update(new_candidate=new_m["id"], new_structure=fmt_cds(d.chosen.cds),
                   new_sources=",".join(sorted(d.chosen.sources)), representative=tier,
                   delta_total=round(d.chosen.total - P.total, 4),
                   delta_evidence=round(d.chosen.evidence - P.evidence, 4),
                   gained_introns=fmt_cds(sorted(new_chain - old_chain)),
                   lost_introns=fmt_cds(sorted(old_chain - new_chain)),
                   old_cds_bp=P.cds_bp, new_cds_bp=d.chosen.cds_bp,
                   old_protein_net=net(primary["id"]), new_protein_net=net(new_m["id"]),
                   isoforms_dropped=";".join(dropped), concerns=";".join(concerns),
                   **{f"new_{k}": v for k, v in d.chosen.components().items()})
        new_m["evidence_selection"] = rec
        primary.pop("evidence_selection", None)
        used.setdefault(d.chosen.key, set()).add(gi)
