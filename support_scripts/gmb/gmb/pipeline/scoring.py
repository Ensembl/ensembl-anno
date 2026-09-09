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


def score_model(
    model: dict,
    config: PipelineConfig,
    protein_supported_tids: set[str],
    genome: dict[str, str] | None = None,
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

    Returns
    -------
    float
        Composite score.
    """
    scfg = config.scoring
    score = 0.0

    # Base evidence weight
    sources = set(model.get("combined_evidence", model["source"]).split(","))
    weights = scfg.weights
    backbone_lower = scfg.backbone_label.strip().lower()
    for s in sources:
        s_lower = s.strip().lower()
        if s_lower == backbone_lower:
            score += weights.backbone
        elif s_lower == "scallop":
            score += weights.scallop
        elif s_lower == "stringtie":
            score += weights.stringtie
        elif s_lower == "minimap2":
            score += weights.minimap2
        else:
            score += 1.0  # unknown source gets base weight

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
                score -= 5.0  # Arbitrary high penalty. Can be tied to config later.
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
        s["score"] = score_model(rep, config, strong_protein_tids, genome)

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
        s["backbone_intron_rescue"] = False
        s["rep"]["backbone_intron_rescue"] = False

    if getattr(scfg, "backbone_intron_rescue", False) and cds_map:
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
        cds_bp = s["rep"].get("cds_bp", 0)
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
    if guard_on:
        has_other_multiexon = any(
            (not s["is_longread_only"]) and s["rep"]["exon_count"] > 1
            for s in candidates
        )
    else:
        has_other_multiexon = False
    for s in candidates:
        s["longread_demoted"] = bool(guard_on and has_other_multiexon
                                     and s["is_longread_only"])
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

    genes = []

    for s in candidates:
        r = s["rep"]
        r["score"] = s["score"]
        # Try to assign to an existing gene sub-cluster
        found = -1
        for i, g_isoforms in enumerate(genes):
            if is_same_gene(r, g_isoforms[0]):
                found = i
                break

        if found == -1:
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

    # Sort genes by start coordinate
    genes.sort(key=lambda g: min(m["start"] for m in g))
    return genes
