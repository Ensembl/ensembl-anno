#!/usr/bin/env python3
"""Isoform scoring and selection for Gene Model Builder.

Replaces the previous hard-coded analyze_locus() with a configurable
scoring function.  Default parameters are tuned for fungal genomes.

All thresholds are driven by PipelineConfig (YAML).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from config import PipelineConfig

import pandas as pd

from annotate_cds_utrs import check_splice_sites

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
    return ",".join(f"{ends[i]}-{starts[i+1]}" for i in range(len(starts) - 1))


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
    for s in sources:
        s_lower = s.strip().lower()
        if s_lower == "helixer":
            score += weights.helixer
        elif s_lower == "scallop":
            score += weights.scallop
        elif s_lower == "stringtie":
            score += weights.stringtie
        else:
            score += 1.0  # unknown source gets base weight

    # Multi-source bonus
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
        if val_cfg.enabled and val_cfg.policy == "penalize":
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


def select_isoforms(
    locus_df: pd.DataFrame,
    config: PipelineConfig,
    protein_supported_tids: set[str],
    genome: dict[str, str] | None = None,
) -> list[list[dict]]:
    """Score and select isoforms for a single locus.

    Parameters
    ----------
    locus_df : pd.DataFrame
        Exon rows for all models in this locus.
    config : PipelineConfig
    protein_supported_tids : set of str
    genome : dict or None

    Returns
    -------
    list of list of dict
        Each inner list is a gene sub-cluster of model dicts (selected
        isoforms), sorted by genomic start position.
    """
    scfg = config.scoring

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
                "protein_support": tid in protein_supported_tids,
                "df": grp,
                "start": grp["Start"].min(),
                "end": grp["End"].max(),
                "exon_count": len(grp),
                "combined_evidence": grp["combined_evidence"].iloc[0]
                if "combined_evidence" in grp.columns
                else source,
                "protein_coding_score": grp["protein_coding_score"].iloc[0]
                if "protein_coding_score" in grp.columns
                else 0.0,
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
                "rep": m,
                "score": 0.0,
            }
        s = merged[key]
        s["sources"].add(m["source"])
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
        if s["protein_support"]:
            rep["protein_support"] = True
        s["score"] = score_model(rep, config, protein_supported_tids, genome)

    # Quality gate: keep models meeting minimum criteria
    candidates = []
    val_cfg = config.protein_validation

    for s in merged.values():
        keep = False

        # Check Protein Validation Strict Drop Policy first
        if val_cfg.enabled and val_cfg.policy == "drop":
            if "protein_coding_score" in s["rep"]:
                if s["rep"]["protein_coding_score"] < val_cfg.min_score:
                    continue  # strictly drops the model

        # CDS length gate: skip models with very short CDS
        cds_bp = s["rep"].get("cds_bp", 0)
        if cds_bp > 0 and cds_bp < scfg.min_cds_bp:
            continue

        # Apply configured gating logic
        if s["protein_support"] or ("Helixer" in s["sources"] and scfg.keep_helixer_without_support):
            keep = True
        elif len(s["sources"]) > 1:
            if scfg.require_protein_support_for_single_source and not s["protein_support"]:
                # Example: configuring strict logic where multi-source without protein drops
                # But since it's multi-source, normally we keep it. For this codebase,
                # we just stick to boolean flags to keep multi-source.
                keep = True
            else:
                keep = True
        elif ("StringTie" in s["sources"] and s["rep"]["exon_count"] > 1) or ("Scallop" in s["sources"] and s["rep"]["exon_count"] > 1):
            if not scfg.require_protein_support_for_single_source:
                keep = True
        elif (
            scfg.fungal_single_exon_mode
            and s["rep"]["exon_count"] == 1
            and s["score"] >= scfg.min_alternate_score
        ):
            # Fungal mode: keep well-supported single-exon models
            keep = True

        # Single-exon models require protein support when configured
        if keep and scfg.require_support_for_single_exon and s["rep"]["exon_count"] == 1:
            if not s["protein_support"] and len(s["sources"]) < 2:
                keep = False

        if keep:
            candidates.append(s)

    if not candidates:
        return []

    # Sort by score descending
    candidates.sort(key=lambda s: s["score"], reverse=True)

    def is_same_gene(m1, m2):
        if m1["intron_chain"] != "single-exon" and m2["intron_chain"] != "single-exon":
            i1 = set(m1["intron_chain"].split(","))
            i2 = set(m2["intron_chain"].split(","))
            if i1.intersection(i2):
                return True
        overlap = min(m1["end"], m2["end"]) - max(m1["start"], m2["start"])
        if overlap > 0:
            len1 = m1["end"] - m1["start"]
            len2 = m2["end"] - m2["start"]
            if (overlap / min(len1, len2)) > 0.15:
                return True
        return False

    genes = []

    for s in candidates:
        r = s["rep"]
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
