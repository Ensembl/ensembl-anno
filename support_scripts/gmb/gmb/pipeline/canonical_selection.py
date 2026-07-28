#!/usr/bin/env python3
"""First-pass canonical transcript selection for multi-isoform GMB genes.

Standalone post-build reporting step -- NOT wired into gmb.cli.build.

Why standalone: canonical selection needs to run after transcripts have
been constructed *and* protein-validation results are available (see the
module's own report), and it must never delete alternative isoforms -- it
is an annotation/reporting layer, not a filter. The safest way to guarantee
both of those is to read GMB's own already-written output files
(``consensus.gff3``, ``evidence_attribution.tsv``, optionally
``protein_validation.tsv``) rather than reach into ``gmb.cli.build``'s
internal state, where a bug could risk mutating the real transcript set.
This also makes the tool trivially re-runnable with different scoring
weights against an unchanged build, and keeps protein validation and domain
evidence both genuinely optional: this tool runs (with reduced-confidence
scoring) even if ``protein_validation.tsv`` was never produced.

Output (all written to --output-dir, source files never touched):
    canonical_transcripts.tsv   -- one row per gene: winner + component
                                   scores + rank/runner-up/reason/confidence.
                                   NOTE: `canonical_total_score` (and its
                                   sibling `runner_up_total_score`) is the
                                   continuous weighted score, reported for
                                   interpreting *how much better* the winner
                                   is -- it is NOT what picked the winner.
                                   Selection uses the lexicographic priority
                                   tuple in `_rank_key` (complete ORF, then
                                   protein-validation subtotal, independent
                                   source count, gmb_score, CDS length,
                                   transcript ID), so a lower total_score can
                                   still win outright (e.g. a complete ORF
                                   beats a fractionally-higher-scoring
                                   partial one). See `_rank_key` and
                                   `reason_code` for the actual decision.
    transcript_ranking.tsv      -- one row per transcript (all isoforms,
                                   winners and non-winners alike)
    canonical_selection_summary.json -- run-level counts
    consensus.canonical_annotated.gff3 (optional, --annotate-gff3) -- a
        COPY of the input GFF3 with `canonical=1`/`canonical=0` added to
        each mRNA's attributes. The original consensus.gff3 is never
        modified in place.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict
from typing import TYPE_CHECKING, Optional

import pandas as pd

from gmb.pipeline.config import load_config
from gmb.utils.gff import parse_gff3_hierarchy
from gmb.utils.io import ensure_dir
from gmb.utils.logging import resolve_log_file, setup_logging

if TYPE_CHECKING:
    from gmb.pipeline.config import CanonicalSelectionConfig, PipelineConfig


# ---------------------------------------------------------------------------
# Data assembly: merge consensus.gff3 + evidence_attribution.tsv +
# (optional) protein_validation.tsv into one per-transcript record dict.
# ---------------------------------------------------------------------------


def load_transcript_records(
    consensus_gff3: str,
    evidence_attribution_tsv: str,
    protein_validation_tsv: Optional[str] = None,
) -> dict[str, list[dict]]:
    """Build ``{gene_id: [transcript_record, ...]}`` from GMB's own outputs.

    Every record carries whatever fields were available; protein-validation
    fields are ``None`` throughout (not an error) when
    ``protein_validation_tsv`` is not supplied or does not exist -- matching
    the constraint that protein validation stays optional end-to-end.
    """
    hierarchy = parse_gff3_hierarchy(consensus_gff3)
    transcripts = hierarchy["transcripts"]

    ev_df = pd.read_csv(evidence_attribution_tsv, sep="\t")
    ev_by_tid = ev_df.set_index("transcript_id").to_dict(orient="index")

    pv_by_tid = {}
    if protein_validation_tsv and os.path.exists(protein_validation_tsv):
        pv_df = pd.read_csv(protein_validation_tsv, sep="\t")
        pv_by_tid = pv_df.set_index("transcript_id").to_dict(orient="index")

    genes: dict[str, list[dict]] = defaultdict(list)
    for tid, tx in transcripts.items():
        ev = ev_by_tid.get(tid, {})
        pv = pv_by_tid.get(tid, {})

        protein_len = tx["cds"] and sum(e - s for s, e in tx["cds"]) // 3 - 1
        record = {
            "gene_id": tx["gene_id"],
            "transcript_id": tid,
            "evidence_sources": ev.get("evidence_sources", tx.get("evidence", "")) or "",
            "exon_count": _nan_to_none(ev.get("exon_count")),
            "cds_bp": _nan_to_none(ev.get("cds_bp")),
            "transcript_span_bp": _nan_to_none(ev.get("transcript_span_bp")),
            "gmb_score": _nan_to_none(ev.get("gmb_score")),
            # protein_validation.tsv, when present, has its own (identical)
            # gmb_score/evidence_sources copy -- prefer evidence_attribution's
            # since it is always written, but fall back to protein_validation's
            # in case only that file is supplied to the CLI directly.
            "diamond_hit": pv.get("diamond_hit"),
            "diamond_pident": _nan_to_none(pv.get("diamond_pident")),
            "diamond_qcov": _nan_to_none(pv.get("diamond_qcov")),
            "diamond_scov": _nan_to_none(pv.get("diamond_scov")),
            "diamond_bitscore": _nan_to_none(pv.get("diamond_bitscore")),
            "diamond_evalue": _nan_to_none(pv.get("diamond_evalue")),
            "psauron_score": _nan_to_none(pv.get("psauron_score")),
            "protein_length": _nan_to_none(pv.get("protein_length")) or (protein_len or None),
            "orf_label": pv.get("orf_label"),
            "is_partial_5": _nan_to_bool(pv.get("is_partial_5")),
            "is_partial_3": _nan_to_bool(pv.get("is_partial_3")),
            "internal_stop_count": _nan_to_none(pv.get("internal_stop_count")) or 0,
            "protein_coding_score": _nan_to_none(pv.get("protein_coding_score")),
        }
        if record["gmb_score"] is None:
            record["gmb_score"] = _nan_to_none(pv.get("gmb_score"))
        if not record["evidence_sources"] and pv.get("evidence_sources"):
            record["evidence_sources"] = pv["evidence_sources"]

        genes[tx["gene_id"]].append(record)

    return genes


def _nan_to_none(val):
    if val is None:
        return None
    try:
        if pd.isna(val):
            return None
    except (TypeError, ValueError):
        pass
    return val


def _nan_to_bool(val) -> Optional[bool]:
    val = _nan_to_none(val)
    if val is None:
        return None
    if isinstance(val, bool):
        return val
    if isinstance(val, str):
        return val.strip().lower() in ("true", "1", "yes")
    return bool(val)


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


def _clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))


def score_transcript(
    record: dict,
    cfg: CanonicalSelectionConfig,
    backbone_label: str,
    gene_gmb_score_range: tuple[float, float],
    domain_metrics: Optional[dict] = None,
) -> dict:
    """Compute all component scores + total for one transcript.

    Every component below is individually reported (not just the total) so
    the weights can be reviewed/re-tuned later without re-running GMB or
    protein validation -- see canonical_transcripts.tsv / transcript_ranking.tsv.
    """
    pv_cfg = cfg.protein_validation
    ev_cfg = cfg.evidence
    st_cfg = cfg.structure
    dm_cfg = cfg.domains

    has_diamond_hit = bool(record.get("diamond_hit"))
    psauron = record.get("psauron_score")
    pident = record.get("diamond_pident")
    qcov = record.get("diamond_qcov")
    scov = record.get("diamond_scov")

    psauron_component = psauron if psauron is not None else 0.0
    identity_component = (pident / 100.0) if pident is not None else 0.0
    qcov_component = (qcov / 100.0) if qcov is not None else 0.0
    scov_component = (scov / 100.0) if scov is not None else 0.0
    hit_component = 1.0 if has_diamond_hit else 0.0

    protein_validation_subtotal = (
        pv_cfg.psauron_weight * psauron_component
        + pv_cfg.diamond_identity_weight * identity_component
        + pv_cfg.diamond_query_coverage_weight * qcov_component
        + pv_cfg.diamond_target_coverage_weight * scov_component
        + pv_cfg.diamond_hit_bonus * hit_component
    )
    has_any_protein_support = psauron is not None or has_diamond_hit

    sources = {
        s.strip() for s in str(record.get("evidence_sources") or "").split(",") if s.strip()
    }
    lower_sources = {s.lower() for s in sources}
    n_sources = len(sources)
    source_count_component = min(n_sources / max(ev_cfg.independent_source_cap, 1), 1.0)
    transcriptomic_flag = 1.0 if ({"scallop", "stringtie"} & lower_sources) else 0.0
    protein_alignment_flag = 1.0 if ({"orthodb", "genblast", "uniprot"} & lower_sources) else 0.0
    longread_flag = 1.0 if any("minimap2" in s for s in lower_sources) else 0.0
    backbone_flag = 1.0 if backbone_label.lower() in lower_sources else 0.0

    gmb_score = record.get("gmb_score")
    lo, hi = gene_gmb_score_range
    if gmb_score is None:
        gmb_score_component = 0.5
    elif hi > lo:
        gmb_score_component = _clamp01((gmb_score - lo) / (hi - lo))
    else:
        gmb_score_component = 0.5  # every isoform tied on gmb_score

    evidence_subtotal = (
        ev_cfg.gmb_score_weight * gmb_score_component
        + ev_cfg.independent_source_weight * source_count_component
        + ev_cfg.transcriptomic_support_weight * transcriptomic_flag
        + ev_cfg.protein_alignment_support_weight * protein_alignment_flag
        + ev_cfg.longread_support_weight * longread_flag
        + ev_cfg.backbone_support_weight * backbone_flag
    )

    is_partial_5 = bool(record.get("is_partial_5"))
    is_partial_3 = bool(record.get("is_partial_3"))
    has_complete_orf = (
        record.get("orf_label") is not None and not is_partial_5 and not is_partial_3
    )
    has_internal_stop = (record.get("internal_stop_count") or 0) > 0
    is_partial = is_partial_5 or is_partial_3

    structure_subtotal = (
        (st_cfg.complete_orf_bonus if has_complete_orf else 0.0)
        - (st_cfg.internal_stop_penalty if has_internal_stop else 0.0)
        - (st_cfg.partial_cds_penalty if is_partial else 0.0)
    )

    domain_subtotal = 0.0
    domain_detail = {}
    if dm_cfg.enabled and domain_metrics:
        dm = domain_metrics.get(record["transcript_id"], {})
        domain_detail = dm
        cov = dm.get("protein_coverage_fraction") or 0.0
        domain_subtotal = (
            dm_cfg.domain_support_weight * (1.0 if dm.get("n_significant_domains") else 0.0)
            + dm_cfg.domain_coverage_weight * cov
            + (dm_cfg.complete_domain_bonus if dm.get("has_complete_domain_coverage") else 0.0)
            - (
                dm_cfg.fragmented_domain_penalty
                if dm.get("n_nonoverlapping_domains", 0) > 1 and cov < 0.5
                else 0.0
            )
            - (dm_cfg.suspicious_fusion_penalty if dm.get("has_suspected_fusion") else 0.0)
            + (
                dm_cfg.cross_provider_agreement_bonus
                if dm.get("cross_provider_agreement")
                else 0.0
            )
        )

    total = protein_validation_subtotal + evidence_subtotal + structure_subtotal + domain_subtotal

    return {
        "psauron_component": round(psauron_component, 4),
        "diamond_identity_component": round(identity_component, 4),
        "diamond_qcov_component": round(qcov_component, 4),
        "diamond_scov_component": round(scov_component, 4),
        "diamond_hit_component": hit_component,
        "protein_validation_subtotal": round(protein_validation_subtotal, 4),
        "has_any_protein_support": has_any_protein_support,
        "gmb_score_component": round(gmb_score_component, 4),
        "source_count_component": round(source_count_component, 4),
        "n_independent_sources": n_sources,
        "transcriptomic_support": bool(transcriptomic_flag),
        "protein_alignment_support": bool(protein_alignment_flag),
        "longread_support": bool(longread_flag),
        "backbone_support": bool(backbone_flag),
        "evidence_subtotal": round(evidence_subtotal, 4),
        "has_complete_orf": has_complete_orf,
        "has_internal_stop": has_internal_stop,
        "is_partial": is_partial,
        "structure_subtotal": round(structure_subtotal, 4),
        "domain_subtotal": round(domain_subtotal, 4),
        "domain_detail": domain_detail,
        "total_score": round(total, 4),
    }


# ---------------------------------------------------------------------------
# Per-gene selection: deterministic tie-break order (documented choice)
# ---------------------------------------------------------------------------

# 1. complete valid ORF                                   (has_complete_orf)
# 2. stronger protein-validation score                    (protein_validation_subtotal)
# 3. broader independent evidence support                 (n_independent_sources)
# 4. higher current GMB score                              (raw gmb_score)
# 5. greater CDS length, among structurally valid models   (cds_bp)
# 6. stable transcript ID lexical ordering                 (transcript_id)
#
# Chosen over ranking purely by the continuous weighted `total_score`
# because the task calls for an explicit, auditable priority order (e.g. a
# complete ORF must never lose to a fractionally-higher-scoring partial one
# just because its evidence weights summed higher) -- `total_score` is still
# computed and reported for magnitude/confidence purposes (the "how much
# better" question), but this tuple decides "who wins".
def _rank_key(rec: dict, scored: dict) -> tuple:
    return (
        not scored["has_complete_orf"],  # False (complete) sorts before True
        -scored["protein_validation_subtotal"],
        -scored["n_independent_sources"],
        -(rec.get("gmb_score") if rec.get("gmb_score") is not None else float("-inf")),
        -(rec.get("cds_bp") if rec.get("cds_bp") is not None else -1),
        rec["transcript_id"],
    )


def _reason_code(winner_rec, winner_scored, runner_rec, runner_scored) -> str:
    if not winner_scored["has_any_protein_support"]:
        return "LOW_CONFIDENCE_NO_PROTEIN_SUPPORT"
    if runner_rec is None:
        return "ONLY_ISOFORM"
    if winner_scored["has_complete_orf"] != runner_scored["has_complete_orf"]:
        return "ONLY_COMPLETE_ORF"
    if (
        winner_scored["protein_validation_subtotal"]
        != runner_scored["protein_validation_subtotal"]
    ):
        # Distinguish which protein-validation signal was the larger
        # contributor to the gap, for a more specific reason code.
        psauron_gap = winner_scored["psauron_component"] - runner_scored["psauron_component"]
        cov_gap = (
            winner_scored["diamond_qcov_component"] + winner_scored["diamond_scov_component"]
        ) - (runner_scored["diamond_qcov_component"] + runner_scored["diamond_scov_component"])
        if abs(psauron_gap) >= abs(cov_gap):
            return "BEST_PSAURON"
        return "BEST_DIAMOND_COVERAGE"
    if winner_scored["n_independent_sources"] != runner_scored["n_independent_sources"]:
        return "BEST_MULTI_SOURCE_SUPPORT"
    if (winner_rec.get("gmb_score") or 0) != (runner_rec.get("gmb_score") or 0):
        return "BEST_PROTEIN_AND_EVIDENCE"
    if (winner_rec.get("cds_bp") or 0) != (runner_rec.get("cds_bp") or 0):
        return "LONGEST_VALID_CDS_TIEBREAK"
    return "LEXICAL_TIEBREAK"


def select_canonical_for_gene(
    gene_id: str,
    records: list[dict],
    cfg: CanonicalSelectionConfig,
    backbone_label: str,
    domain_metrics: Optional[dict] = None,
) -> dict:
    """Rank a gene's isoforms and return the full selection result.

    Never mutates or drops ``records`` -- purely read-only scoring/ranking.
    """
    gmb_scores = [r["gmb_score"] for r in records if r.get("gmb_score") is not None]
    gene_range = (min(gmb_scores), max(gmb_scores)) if gmb_scores else (0.0, 0.0)

    scored_by_tid = {
        r["transcript_id"]: score_transcript(r, cfg, backbone_label, gene_range, domain_metrics)
        for r in records
    }
    ranked = sorted(records, key=lambda r: _rank_key(r, scored_by_tid[r["transcript_id"]]))

    winner = ranked[0]
    winner_scored = scored_by_tid[winner["transcript_id"]]
    runner = ranked[1] if len(ranked) > 1 else None
    runner_scored = scored_by_tid[runner["transcript_id"]] if runner else None

    reason = _reason_code(winner, winner_scored, runner, runner_scored)

    score_gap = (
        winner_scored["total_score"] - runner_scored["total_score"] if runner_scored else None
    )
    highest_raw_gmb_tid = max(
        records,
        key=lambda r: (r.get("gmb_score") if r.get("gmb_score") is not None else float("-inf")),
    )["transcript_id"]

    if not winner_scored["has_any_protein_support"]:
        confidence = "low_confidence"
    elif score_gap is None:
        confidence = "high_confidence"  # only one isoform
    elif score_gap < cfg.low_confidence_score_gap:
        confidence = "low_confidence"
    else:
        confidence = "high_confidence"

    return {
        "gene_id": gene_id,
        "n_isoforms": len(records),
        "canonical_transcript_id": winner["transcript_id"],
        "canonical_total_score": winner_scored["total_score"],
        "runner_up_transcript_id": runner["transcript_id"] if runner else None,
        "runner_up_total_score": runner_scored["total_score"] if runner_scored else None,
        "score_gap": round(score_gap, 4) if score_gap is not None else None,
        "selection_reason": reason,
        "confidence": confidence,
        "differs_from_highest_gmb_score": winner["transcript_id"] != highest_raw_gmb_tid,
        "highest_gmb_score_transcript_id": highest_raw_gmb_tid,
        "ranked_transcript_ids": [r["transcript_id"] for r in ranked],
        "_scored_by_tid": scored_by_tid,
        "_ranked_records": ranked,
    }


# ---------------------------------------------------------------------------
# Top-level run
# ---------------------------------------------------------------------------


def run_canonical_selection(
    consensus_gff3: str,
    evidence_attribution_tsv: str,
    output_dir: str,
    config: PipelineConfig,
    protein_validation_tsv: Optional[str] = None,
    domain_sidecar_rows: Optional[list[dict]] = None,
    annotate_gff3: bool = False,
) -> dict:
    ensure_dir(output_dir)
    cfg = config.canonical_selection

    genes = load_transcript_records(
        consensus_gff3, evidence_attribution_tsv, protein_validation_tsv
    )

    domain_metrics = None
    if cfg.domains.enabled and domain_sidecar_rows:
        from gmb.pipeline.domain_evidence import summarize_transcript_domain_metrics

        protein_lengths = {
            r["transcript_id"]: r["protein_length"]
            for recs in genes.values()
            for r in recs
            if r.get("protein_length")
        }
        domain_metrics = summarize_transcript_domain_metrics(domain_sidecar_rows, protein_lengths)

    gene_results = {}
    for gene_id, records in genes.items():
        gene_results[gene_id] = select_canonical_for_gene(
            gene_id, records, cfg, config.scoring.backbone_label, domain_metrics
        )

    canonical_rows = []
    ranking_rows = []
    for gene_id, result in gene_results.items():
        scored_by_tid = result["_scored_by_tid"]
        winner_scored = scored_by_tid[result["canonical_transcript_id"]]
        canonical_rows.append(
            {
                "gene_id": gene_id,
                "canonical_transcript_id": result["canonical_transcript_id"],
                "canonical_total_score": result["canonical_total_score"],
                "protein_validation_subtotal": winner_scored["protein_validation_subtotal"],
                "evidence_subtotal": winner_scored["evidence_subtotal"],
                "structure_subtotal": winner_scored["structure_subtotal"],
                "domain_subtotal": winner_scored["domain_subtotal"],
                "rank": 1,
                "n_isoforms": result["n_isoforms"],
                "selection_reason": result["selection_reason"],
                "confidence": result["confidence"],
                "runner_up_transcript_id": result["runner_up_transcript_id"],
                "runner_up_total_score": result["runner_up_total_score"],
                "score_gap": result["score_gap"],
                "differs_from_highest_gmb_score": result["differs_from_highest_gmb_score"],
                "highest_gmb_score_transcript_id": result["highest_gmb_score_transcript_id"],
            }
        )
        for rank, rec in enumerate(result["_ranked_records"], start=1):
            scored = scored_by_tid[rec["transcript_id"]]
            ranking_rows.append(
                {
                    "gene_id": gene_id,
                    "transcript_id": rec["transcript_id"],
                    "rank": rank,
                    "is_canonical": rank == 1,
                    "total_score": scored["total_score"],
                    "protein_validation_subtotal": scored["protein_validation_subtotal"],
                    "evidence_subtotal": scored["evidence_subtotal"],
                    "structure_subtotal": scored["structure_subtotal"],
                    "domain_subtotal": scored["domain_subtotal"],
                    "psauron_score": rec.get("psauron_score"),
                    "diamond_hit": rec.get("diamond_hit"),
                    "diamond_pident": rec.get("diamond_pident"),
                    "diamond_qcov": rec.get("diamond_qcov"),
                    "diamond_scov": rec.get("diamond_scov"),
                    "has_complete_orf": scored["has_complete_orf"],
                    "has_internal_stop": scored["has_internal_stop"],
                    "gmb_score": rec.get("gmb_score"),
                    "n_independent_sources": scored["n_independent_sources"],
                    "evidence_sources": rec.get("evidence_sources"),
                    "cds_bp": rec.get("cds_bp"),
                }
            )

    can_path = os.path.join(output_dir, "canonical_transcripts.tsv")
    pd.DataFrame(canonical_rows).to_csv(can_path, sep="\t", index=False)

    rank_path = os.path.join(output_dir, "transcript_ranking.tsv")
    pd.DataFrame(ranking_rows).to_csv(rank_path, sep="\t", index=False)

    n_multi_isoform = sum(1 for r in gene_results.values() if r["n_isoforms"] > 1)
    n_differs = sum(1 for r in gene_results.values() if r["differs_from_highest_gmb_score"])
    reason_counts = defaultdict(int)
    confidence_counts = defaultdict(int)
    for r in gene_results.values():
        reason_counts[r["selection_reason"]] += 1
        confidence_counts[r["confidence"]] += 1

    summary = {
        "n_genes": len(gene_results),
        "n_multi_isoform_genes": n_multi_isoform,
        "n_single_isoform_genes": len(gene_results) - n_multi_isoform,
        "n_canonical_differs_from_highest_gmb_score": n_differs,
        "selection_reason_counts": dict(reason_counts),
        "confidence_counts": dict(confidence_counts),
        "config": {
            "protein_validation": vars(cfg.protein_validation),
            "evidence": vars(cfg.evidence),
            "structure": vars(cfg.structure),
            "domains": vars(cfg.domains),
            "low_confidence_score_gap": cfg.low_confidence_score_gap,
        },
    }
    summary_path = os.path.join(output_dir, "canonical_selection_summary.json")
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2, default=str)

    if annotate_gff3:
        canonical_ids = {r["canonical_transcript_id"] for r in gene_results.values()}
        annotated_path = os.path.join(output_dir, "consensus.canonical_annotated.gff3")
        _write_annotated_gff3(consensus_gff3, canonical_ids, annotated_path)
        summary["annotated_gff3"] = annotated_path

    print(f"Canonical selection: {len(gene_results)} genes ({n_multi_isoform} multi-isoform)")
    print(f"  {can_path}")
    print(f"  {rank_path}")
    print(f"  {summary_path}")

    return summary


def _write_annotated_gff3(input_gff3: str, canonical_ids: set[str], output_path: str) -> None:
    """Write a COPY of ``input_gff3`` with `canonical=1`/`canonical=0` added
    to each mRNA row's attributes. The source file is never modified."""
    with open(input_gff3) as fh_in, open(output_path, "w") as fh_out:
        for line in fh_in:
            if line.startswith("#") or not line.strip():
                fh_out.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "mRNA":
                fh_out.write(line)
                continue
            attrs = parts[8]
            tid = None
            for kv in attrs.split(";"):
                if kv.startswith("ID="):
                    tid = kv[3:]
                    break
            flag = "1" if tid in canonical_ids else "0"
            sep = "" if attrs.endswith(";") else ";"
            parts[8] = f"{attrs}{sep}canonical={flag}"
            fh_out.write("\t".join(parts) + "\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "First-pass canonical transcript selection: reads GMB's own "
            "consensus.gff3/evidence_attribution.tsv (and, optionally, "
            "protein_validation.tsv) and writes a ranking/reporting layer. "
            "Never modifies or drops alternative isoforms."
        )
    )
    parser.add_argument("--consensus-gff3", required=True)
    parser.add_argument("--evidence-attribution", required=True)
    parser.add_argument(
        "--protein-validation", default=None, help="Optional protein_validation.tsv"
    )
    parser.add_argument("--config", default=None, help="YAML with a canonical_selection: section")
    parser.add_argument("--preset", default="fungi")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--annotate-gff3",
        action="store_true",
        help="Also write a copy of the GFF3 with canonical=1/0 mRNA attributes",
    )
    parser.add_argument("--log-file", default=None)
    parser.add_argument("--no-log-file", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ensure_dir(args.output_dir)
    log_file = resolve_log_file(args.output_dir, args.log_file, args.no_log_file)
    setup_logging(log_file=log_file, capture_stdio=log_file is not None)
    if log_file:
        print(f"Logging to {log_file}")

    for label, path in [
        ("--consensus-gff3", args.consensus_gff3),
        ("--evidence-attribution", args.evidence_attribution),
    ]:
        if not os.path.exists(path):
            sys.exit(f"ERROR: {label} path not found: {path}")
    if args.protein_validation and not os.path.exists(args.protein_validation):
        print(
            f"WARNING: --protein-validation path not found ({args.protein_validation}); "
            "continuing with protein-validation components unavailable (0.0/None)."
        )
        args.protein_validation = None

    config = load_config(args.config, args.preset)
    if not config.canonical_selection.enabled:
        sys.exit("canonical_selection.enabled is false in the given config; nothing to do.")

    run_canonical_selection(
        consensus_gff3=args.consensus_gff3,
        evidence_attribution_tsv=args.evidence_attribution,
        output_dir=args.output_dir,
        config=config,
        protein_validation_tsv=args.protein_validation,
        annotate_gff3=args.annotate_gff3,
    )
    print("Done.")


if __name__ == "__main__":
    main()
