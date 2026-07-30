#!/usr/bin/env python3
"""Parse completed InterProScan 6 output and review ambiguous canonical choices.

Execution independence
----------------------
This module NEVER runs InterProScan. It consumes files that already exist,
by explicit path, and is therefore indifferent to how they were produced:
Docker locally, Singularity/Apptainer on a cluster, a Slurm executor, or a
shared results directory. No container runtime, work directory, image name
or site-specific path appears anywhere in this module -- that separation is
deliberate and load-bearing (see the execution-profile examples in the
README rather than adding any of it here).

Canonical machine-readable input
--------------------------------
JSONL is the primary input. Verified against a real InterProScan 6.0.1 run
(InterPro 109.0), it is the only emitted format carrying, together:
  * the `representative` flag per location (GFF3 has it, TSV does not);
  * the integrated InterPro `entry` (TSV has accession/description columns
    but they are "-" for unintegrated matches and carry no entry *type*);
  * the member-database library name AND its release version.
GFF3 is accepted as an alternative for representative locations, but it
does not carry integrated-entry type, so JSONL is preferred where both
exist.

Evidence model
--------------
Raw match count is deliberately NOT scored: overlapping member-database
signatures routinely describe the same biological region (the reference
run produced 3864 locations but only 424 representative ones), so counting
matches would reward redundancy rather than architecture. Scoring is based
on representative locations, integrated InterPro entries, and the fraction
of the protein covered by them.

Absence of matches is never itself a penalty -- a lineage-specific protein
with no InterPro signature is a normal, expected outcome for Apicomplexa,
not evidence against a model.
"""

from __future__ import annotations

import json
import os
from collections import defaultdict
from typing import Optional

# ---------------------------------------------------------------------------
# Feature-type weighting
# ---------------------------------------------------------------------------
#
# InterPro entry/signature types are not equally strong evidence that a
# model's coding structure is right. Observed in the reference run:
# Homologous_superfamily, Region, Domain, Family, Repeat, Conserved_site,
# Coiled_coil.

# Architecture-defining: a Family/Domain/Homologous_superfamily hit says
# something substantive about what the protein *is*.
STRONG_TYPES = {"Family", "Domain", "Homologous_superfamily"}
# Real but weaker: informative about parts, not overall architecture.
SUPPORTING_TYPES = {"Conserved_site", "Repeat", "Binding_site", "Active_site", "PTM"}
# Low-complexity / structural predictions. These fire on almost anything
# and must not be read as functional support.
NON_EVIDENCE_TYPES = {"Coiled_coil", "Region", "Disordered"}
# Member databases whose hits are structural/disorder predictors rather
# than functional signatures.
NON_EVIDENCE_LIBRARIES = {"MobiDB-lite", "COILS", "Phobius", "SignalP", "TMHMM"}
# AntiFam exists to identify spurious ORFs -- a hit here is evidence
# AGAINST the model being a real protein.
NEGATIVE_LIBRARIES = {"AntiFam"}


def classify_library(library: Optional[str], feature_type: Optional[str]) -> str:
    """Bucket a match into strong/supporting/non_evidence/negative."""
    if library in NEGATIVE_LIBRARIES:
        return "negative"
    if library in NON_EVIDENCE_LIBRARIES:
        return "non_evidence"
    if feature_type in STRONG_TYPES:
        return "strong"
    if feature_type in SUPPORTING_TYPES:
        return "supporting"
    if feature_type in NON_EVIDENCE_TYPES:
        return "non_evidence"
    return "supporting"


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------


def parse_interproscan_jsonl(path: str) -> dict:
    """Parse InterProScan 6 JSONL into ``{protein_id: [match_record, ...]}``.

    Tolerant by design: a malformed line is skipped rather than aborting a
    whole review, and every optional field (evalue, score, entry, hmm
    bounds) is allowed to be absent -- PROSITE patterns and COILS, for
    example, legitimately emit no e-value.

    ``protein_id`` is taken from ``xref[].id``; when GMB produced the FASTA
    that is the protein SHA-256 checksum, but nothing here requires that.
    """
    by_protein: dict = defaultdict(list)
    if not path or not os.path.exists(path):
        return dict(by_protein)

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                payload = json.loads(line)
            except json.JSONDecodeError:
                continue
            ips_version = payload.get("interproscan-version")
            ipr_version = payload.get("interpro-version")
            for result in payload.get("results") or []:
                for xref in result.get("xref") or [{}]:
                    protein_id = xref.get("id")
                    if not protein_id:
                        continue
                    by_protein[protein_id].extend(
                        _matches_from_result(result, ips_version, ipr_version)
                    )
    return dict(by_protein)


def _matches_from_result(result: dict, ips_version, ipr_version) -> list:
    records = []
    for match in result.get("matches") or []:
        signature = match.get("signature") or {}
        library_release = signature.get("signatureLibraryRelease") or {}
        entry = signature.get("entry")
        library = library_release.get("library")
        feature_type = signature.get("type")
        for location in match.get("locations") or []:
            records.append(
                {
                    "protein_md5": result.get("md5"),
                    "signature_accession": signature.get("accession"),
                    "signature_description": signature.get("description"),
                    "signature_type": feature_type,
                    "member_database": library,
                    "member_database_version": library_release.get("version"),
                    "interpro_accession": (entry or {}).get("accession"),
                    "interpro_description": (entry or {}).get("description"),
                    "interpro_type": (entry or {}).get("type"),
                    "start": location.get("start"),
                    "end": location.get("end"),
                    "representative": bool(location.get("representative")),
                    # Provider-specific: match-level and location-level
                    # values both exist and are NOT comparable across member
                    # databases, so both are carried through unmodified and
                    # no universal significance threshold is applied.
                    "evalue": location.get("evalue", match.get("evalue")),
                    "score": location.get("score", match.get("score")),
                    "evidence_class": classify_library(library, feature_type),
                    "interproscan_version": ips_version,
                    "interpro_version": ipr_version,
                }
            )
    return records


def parse_interproscan_gff3(path: str) -> dict:
    """Parse InterProScan 6 GFF3 into ``{protein_id: [match_record, ...]}``.

    Secondary input: GFF3 carries `representative` and `type` in its
    attributes but not the integrated InterPro entry type, so records
    produced here have ``interpro_accession``/``interpro_type`` set to
    None. Use JSONL when available.
    """
    by_protein: dict = defaultdict(list)
    if not path or not os.path.exists(path):
        return dict(by_protein)

    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = {}
            for item in parts[8].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attrs[key.strip()] = value.strip()
            feature_type = attrs.get("type")
            library = parts[1]
            by_protein[parts[0]].append(
                {
                    "protein_md5": None,
                    "signature_accession": attrs.get("Name"),
                    "signature_description": attrs.get("signature_desc"),
                    "signature_type": feature_type,
                    "member_database": library,
                    "member_database_version": None,
                    "interpro_accession": None,
                    "interpro_description": None,
                    "interpro_type": None,
                    "start": _safe_int(parts[3]),
                    "end": _safe_int(parts[4]),
                    "representative": attrs.get("representative", "").lower() == "true",
                    "evalue": _safe_float(parts[5]),
                    "score": None,
                    "evidence_class": classify_library(library, feature_type),
                    "interproscan_version": None,
                    "interpro_version": None,
                }
            )
    return dict(by_protein)


def _safe_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _safe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Per-protein architecture summary
# ---------------------------------------------------------------------------


def summarise_architecture(matches: list, protein_length: Optional[int] = None) -> dict:
    """Reduce one protein's matches to an interpretable architecture summary.

    Coverage is computed from merged REPRESENTATIVE locations only, so
    overlapping signatures describing one region contribute once.
    """
    representative = [m for m in matches if m.get("representative")]
    strong = [m for m in representative if m["evidence_class"] == "strong"]
    negative = [m for m in matches if m["evidence_class"] == "negative"]

    intervals = sorted(
        (m["start"], m["end"])
        for m in representative
        if m.get("start") is not None and m.get("end") is not None
    )
    merged = []
    for start, end in intervals:
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    covered = sum(end - start + 1 for start, end in merged)
    coverage_fraction = (
        round(covered / protein_length, 4) if protein_length and protein_length > 0 else None
    )

    integrated = sorted(
        {m["interpro_accession"] for m in representative if m.get("interpro_accession")}
    )
    # Multiple member databases independently supporting one integrated
    # entry is stronger than one database asserting it alone.
    dbs_per_entry: dict = defaultdict(set)
    for m in matches:
        if m.get("interpro_accession") and m.get("member_database"):
            dbs_per_entry[m["interpro_accession"]].add(m["member_database"])
    corroborated = sorted(acc for acc, dbs in dbs_per_entry.items() if len(dbs) > 1)

    architecture = [
        f"{m.get('interpro_accession') or m.get('signature_accession')}"
        f"[{m.get('start')}-{m.get('end')}]"
        for m in sorted(strong, key=lambda m: (m.get("start") or 0))
    ]

    # Span of the STRONG (Family/Domain/Homologous_superfamily) architecture
    # only, used to localise a truncation (N- vs C-terminal) when comparing
    # two candidates -- see decide_canonical()/classify_reason_code(). None
    # when there is no strong representative evidence at all.
    strong_coords = [
        (m["start"], m["end"])
        for m in strong
        if m.get("start") is not None and m.get("end") is not None
    ]
    strong_span = (
        (min(s for s, _ in strong_coords), max(e for _, e in strong_coords))
        if strong_coords
        else None
    )

    return {
        "n_matches_raw": len(matches),
        "n_representative_locations": len(representative),
        "n_strong_representative": len(strong),
        "integrated_entries": integrated,
        "n_integrated_entries": len(integrated),
        "corroborated_entries": corroborated,
        "n_corroborated_entries": len(corroborated),
        "representative_coverage_fraction": coverage_fraction,
        "architecture": architecture,
        "strong_span": strong_span,
        # Individual (not merged) strong-representative intervals, used by
        # detect_fusion to find a region entirely disjoint from the other
        # candidate's span -- strong_span alone (just the overall min/max)
        # cannot distinguish "one wide domain" from "two disjoint chunks
        # with a gap between them".
        "strong_intervals": strong_coords,
        "member_databases": sorted(
            {m["member_database"] for m in matches if m.get("member_database")}
        ),
        "has_negative_evidence": bool(negative),
        "negative_signatures": sorted(
            {m["signature_accession"] for m in negative if m.get("signature_accession")}
        ),
        "has_domains": bool(strong),
    }


# ---------------------------------------------------------------------------
# Resolver verdicts
# ---------------------------------------------------------------------------

VERDICT_SUPPORTS_CURRENT = "supports_current"
VERDICT_SUPPORTS_RUNNER_UP = "supports_runner_up"
VERDICT_UNINFORMATIVE = "uninformative"
VERDICT_CONFLICTING = "conflicting"


def compare_candidates(current: dict, runner_up: dict) -> dict:
    """Compare two architecture summaries and return a verdict + explanation.

    Never recommends a change on absence of evidence alone: if neither
    candidate has domains, the verdict is uninformative and the current
    choice stands.
    """
    reasons = []

    # Negative (AntiFam) evidence is decisive in one direction only.
    if current["has_negative_evidence"] and not runner_up["has_negative_evidence"]:
        return {
            "verdict": VERDICT_SUPPORTS_RUNNER_UP,
            "explanation": (
                "current canonical matches a negative-evidence signature "
                f"({','.join(current['negative_signatures'])}); runner-up does not"
            ),
        }
    if runner_up["has_negative_evidence"] and not current["has_negative_evidence"]:
        return {
            "verdict": VERDICT_SUPPORTS_CURRENT,
            "explanation": (
                "runner-up matches a negative-evidence signature "
                f"({','.join(runner_up['negative_signatures'])}); current canonical does not"
            ),
        }

    if not current["has_domains"] and not runner_up["has_domains"]:
        return {
            "verdict": VERDICT_UNINFORMATIVE,
            "explanation": (
                "neither candidate has representative family/domain matches -- "
                "absence of InterPro signal is not evidence against either model"
            ),
        }

    cur_entries = set(current["integrated_entries"])
    run_entries = set(runner_up["integrated_entries"])
    cur_only = sorted(cur_entries - run_entries)
    run_only = sorted(run_entries - cur_entries)

    if cur_only:
        reasons.append(f"current retains integrated entries absent from runner-up: {cur_only}")
    if run_only:
        reasons.append(f"runner-up retains integrated entries absent from current: {run_only}")

    cur_cov = current["representative_coverage_fraction"] or 0.0
    run_cov = runner_up["representative_coverage_fraction"] or 0.0
    if abs(cur_cov - run_cov) >= 0.10:
        better = "current" if cur_cov > run_cov else "runner-up"
        reasons.append(
            f"{better} has materially higher representative-domain coverage "
            f"({cur_cov:.2f} vs {run_cov:.2f})"
        )

    # Both directions have exclusive integrated entries -> genuinely
    # conflicting architectures; a human should look.
    if cur_only and run_only:
        return {
            "verdict": VERDICT_CONFLICTING,
            "explanation": "; ".join(reasons)
            or "each candidate carries integrated entries the other lacks",
        }

    cur_strength = (len(cur_entries), current["n_strong_representative"], cur_cov)
    run_strength = (len(run_entries), runner_up["n_strong_representative"], run_cov)
    if run_strength > cur_strength:
        return {
            "verdict": VERDICT_SUPPORTS_RUNNER_UP,
            "explanation": "; ".join(reasons)
            or "runner-up has a stronger representative domain architecture",
        }
    if cur_strength > run_strength:
        return {
            "verdict": VERDICT_SUPPORTS_CURRENT,
            "explanation": "; ".join(reasons)
            or "current canonical has a stronger representative domain architecture",
        }
    return {
        "verdict": VERDICT_UNINFORMATIVE,
        "explanation": "; ".join(reasons)
        or "candidates have equivalent representative domain architecture",
    }


# ---------------------------------------------------------------------------
# Conservative replacement policy (Part 10/11 of the InterPro-resolver task)
# ---------------------------------------------------------------------------
#
# The four VERDICT_* constants above answer "which architecture looks
# stronger"; the REASON_* constants below answer the more specific "why,
# in terms an operator/reviewer can act on, would we replace the initial
# canonical transcript". A reason code is produced for EVERY reviewed gene,
# whether or not a replacement actually happens (see decide_canonical) --
# INTERPRO_SUPPORTS_INITIAL/INTERPRO_EQUIVALENT/INTERPRO_UNINFORMATIVE/
# INTERPRO_CONFLICTING never replace; the other six are candidate reasons
# TO replace, still gated by check_safeguards() before they take effect.

REASON_COMPLETE_DOMAIN_RESTORED = "INTERPRO_COMPLETE_DOMAIN_RESTORED"
REASON_N_TERMINAL_TRUNCATION_RESOLVED = "INTERPRO_N_TERMINAL_TRUNCATION_RESOLVED"
REASON_C_TERMINAL_TRUNCATION_RESOLVED = "INTERPRO_C_TERMINAL_TRUNCATION_RESOLVED"
REASON_ARCHITECTURE_RESTORED = "INTERPRO_ARCHITECTURE_RESTORED"
REASON_FUSION_AVOIDED = "INTERPRO_FUSION_AVOIDED"
REASON_ANTIFAM_AVOIDED = "INTERPRO_ANTIFAM_AVOIDED"
REASON_SUPPORTS_INITIAL = "INTERPRO_SUPPORTS_INITIAL"
REASON_EQUIVALENT = "INTERPRO_EQUIVALENT"
REASON_UNINFORMATIVE = "INTERPRO_UNINFORMATIVE"
REASON_CONFLICTING = "INTERPRO_CONFLICTING"

# Reason codes that identify a safety concern DIAMOND/psauron cannot
# reliably see on their own (a spurious ORF can still score a plausible
# hit; a fusion/read-through model still contains the real protein's
# sequence and can still score well) -- these bypass the
# max_protein_validation_override_gap safeguard specifically, though never
# the structural safeguards. See check_safeguards().
_SAFETY_OVERRIDE_REASONS = {REASON_ANTIFAM_AVOIDED, REASON_FUSION_AVOIDED}
# Reason codes that never propose a replacement in the first place.
_NO_REPLACEMENT_REASONS = {
    REASON_SUPPORTS_INITIAL,
    REASON_EQUIVALENT,
    REASON_UNINFORMATIVE,
    REASON_CONFLICTING,
}


def detect_fusion(
    current: dict, runner_up: dict, current_length, runner_up_length, length_ratio_flag: float
) -> bool:
    """Heuristic, conservative suspected-fusion detection.

    A read-through/fusion model can have MORE integrated entries than a
    clean single-domain alternative -- which would make
    ``compare_candidates``'s raw strength tuple prefer it, backwards. This
    check looks for the specific pattern that actually suggests a fusion
    rather than genuinely broader support: the current candidate is
    substantially longer than the runner-up AND carries strong-domain
    evidence entirely outside the runner-up's own strong-domain span (i.e.
    an extra region bolted on, not a longer/more-complete version of the
    same region), while the runner-up itself is not domain-free (so this
    is a real comparison, not "runner-up simply has no evidence").

    Deliberately conservative: both the length-ratio AND the
    non-overlapping-extra-region conditions must hold. A merely longer
    isoform with the SAME domain, more fully covered, is exactly the
    N/C-terminal-truncation-resolved case instead, not a fusion.
    """
    if not current.get("has_domains") or not runner_up.get("has_domains"):
        return False
    if not current_length or not runner_up_length or runner_up_length <= 0:
        return False
    ratio = float(current_length) / float(runner_up_length)
    if ratio < length_ratio_flag:
        return False
    run_span = runner_up.get("strong_span")
    if not run_span:
        return False
    run_start, run_end = run_span
    # Does current carry a strong-domain interval entirely outside the
    # runner-up's own strong-domain span (zero overlap) -- the "bolted-on"
    # signature of a fusion, rather than current simply covering more of
    # the SAME region (which is a truncation-resolution case, handled by
    # classify_reason_code's N/C-terminal localisation instead).
    for start, end in current.get("strong_intervals") or []:
        if end < run_start or start > run_end:
            return True
    return False


def classify_reason_code(current: dict, runner_up: dict, comparison: dict, is_fusion: bool) -> str:
    """Map a verdict + the two architecture summaries to a specific reason code.

    Called only after ``detect_fusion`` has had first refusal (a detected
    fusion always reports REASON_FUSION_AVOIDED regardless of the generic
    verdict, since the generic strength tuple is exactly backwards for a
    fusion -- see detect_fusion's docstring).
    """
    if is_fusion:
        return REASON_FUSION_AVOIDED

    verdict = comparison["verdict"]
    if verdict == VERDICT_CONFLICTING:
        return REASON_CONFLICTING
    if verdict == VERDICT_UNINFORMATIVE:
        if not current["has_domains"] and not runner_up["has_domains"]:
            return REASON_UNINFORMATIVE
        return REASON_EQUIVALENT
    if verdict == VERDICT_SUPPORTS_CURRENT:
        return REASON_SUPPORTS_INITIAL

    # verdict == VERDICT_SUPPORTS_RUNNER_UP -- distinguish why.
    if current["has_negative_evidence"] and not runner_up["has_negative_evidence"]:
        return REASON_ANTIFAM_AVOIDED
    if not current["has_domains"] and runner_up["has_domains"]:
        return REASON_COMPLETE_DOMAIN_RESTORED

    cur_span = current.get("strong_span")
    run_span = runner_up.get("strong_span")
    if cur_span and run_span:
        cur_start, cur_end = cur_span
        run_start, run_end = run_span
        # Runner-up's architecture extends further at one end only ->
        # localise the truncation the initial choice was missing.
        extends_n_term = run_start < cur_start
        extends_c_term = run_end > cur_end
        if extends_n_term and not extends_c_term:
            return REASON_N_TERMINAL_TRUNCATION_RESOLVED
        if extends_c_term and not extends_n_term:
            return REASON_C_TERMINAL_TRUNCATION_RESOLVED
    # Either no localisable single-flank difference, or one/both spans
    # unavailable -- still a genuine architecture improvement per the
    # generic verdict, just not cleanly localised to one terminus.
    return REASON_ARCHITECTURE_RESTORED


def check_safeguards(current_row: dict, runner_up_row: dict, reason_code: str, cfg) -> dict:
    """Evaluate every configured safeguard for a proposed replacement.

    Returns ``{"passed": bool, "checks": {name: bool}, "notes": [str, ...]}``.
    A replacement is only ever applied when every configured (True) check
    passes. Structural safeguards apply unconditionally; the
    protein-validation override-gap safeguard is skipped for the two
    safety-critical reason codes in ``_SAFETY_OVERRIDE_REASONS`` (see their
    docstrings for why DIAMOND/psauron alone cannot be trusted to catch
    those specific failure modes).
    """
    checks: dict = {}
    notes: list = []

    if cfg.require_no_worse_structural_class:
        cur_complete = bool(current_row.get("has_complete_orf"))
        run_complete = bool(runner_up_row.get("has_complete_orf"))
        ok = run_complete or not cur_complete
        checks["no_worse_structural_class"] = ok
        if not ok:
            notes.append(
                "runner-up has an incomplete ORF where the initial canonical was complete"
            )

    if cfg.require_no_new_internal_stops:
        cur_stop = bool(current_row.get("has_internal_stop"))
        run_stop = bool(runner_up_row.get("has_internal_stop"))
        ok = not run_stop or cur_stop
        checks["no_new_internal_stops"] = ok
        if not ok:
            notes.append("runner-up introduces an internal stop the initial canonical did not have")

    if reason_code not in _SAFETY_OVERRIDE_REASONS:
        cur_psauron = current_row.get("psauron_score")
        run_psauron = runner_up_row.get("psauron_score")
        cur_cov = current_row.get("diamond_balanced_coverage")
        run_cov = runner_up_row.get("diamond_balanced_coverage")
        gap_exceeded = False
        if cur_psauron is not None and run_psauron is not None:
            if (float(cur_psauron) - float(run_psauron)) > cfg.max_protein_validation_override_gap:
                gap_exceeded = True
        if cur_cov is not None and run_cov is not None:
            if (float(cur_cov) - float(run_cov)) > cfg.max_protein_validation_override_gap:
                gap_exceeded = True
        checks["no_much_stronger_protein_validation_overridden"] = not gap_exceeded
        if gap_exceeded:
            notes.append(
                "initial canonical has materially stronger DIAMOND/psauron evidence than the "
                f"InterPro-recommended candidate (gap > {cfg.max_protein_validation_override_gap}); "
                "InterPro architecture evidence alone must not override it"
            )

    return {"passed": all(checks.values()) if checks else True, "checks": checks, "notes": notes}


def decide_canonical(
    gene_id: str,
    current_row: dict,
    runner_up_row: dict,
    cur_summary: dict,
    run_summary: dict,
    cfg,
) -> dict:
    """Full per-gene decision: verdict, specific reason code, safeguards, and
    the resulting initial/recommended/final canonical transcript IDs.

    ``final_canonical_transcript_id`` equals ``initial_canonical_transcript_id``
    whenever the reason code is not a replacement candidate, a safeguard
    fails, or ``cfg.apply_replacements`` is False (report-only mode) --
    Part 11's "final = initial when disabled/not applied" contract, made
    explicit at the per-gene level rather than only at the master switch.
    """
    initial_id = current_row["transcript_id"]
    recommended_id = runner_up_row["transcript_id"]

    is_fusion = detect_fusion(
        cur_summary,
        run_summary,
        current_row.get("protein_length"),
        runner_up_row.get("protein_length"),
        cfg.length_ratio_flag,
    )
    if is_fusion:
        comparison = {
            "verdict": VERDICT_SUPPORTS_RUNNER_UP,
            "explanation": (
                "initial canonical is substantially longer than the runner-up and carries "
                "strong-domain evidence entirely outside the runner-up's own domain span -- "
                "suspected fusion/read-through model"
            ),
        }
    else:
        comparison = compare_candidates(cur_summary, run_summary)

    reason_code = classify_reason_code(cur_summary, run_summary, comparison, is_fusion)
    is_replacement_candidate = reason_code not in _NO_REPLACEMENT_REASONS

    safeguards = {"passed": True, "checks": {}, "notes": []}
    replaced = False
    final_id = initial_id
    if is_replacement_candidate:
        safeguards = check_safeguards(current_row, runner_up_row, reason_code, cfg)
        if cfg.apply_replacements and safeguards["passed"]:
            replaced = True
            final_id = recommended_id

    return {
        "gene_id": gene_id,
        "initial_canonical_transcript_id": initial_id,
        "interpro_recommended_transcript_id": recommended_id if is_replacement_candidate else initial_id,
        "final_canonical_transcript_id": final_id,
        "replaced": replaced,
        "interpro_verdict": comparison["verdict"],
        "reason_code": reason_code,
        "explanation": comparison["explanation"],
        "safeguards_passed": safeguards["passed"],
        "safeguard_checks": safeguards["checks"],
        "safeguard_notes": "; ".join(safeguards["notes"]),
    }


# ---------------------------------------------------------------------------
# Review report
# ---------------------------------------------------------------------------


def build_resolver_report(
    manifest_rows: list,
    matches_by_protein: dict,
    cfg=None,
) -> list:
    """One report row per reviewed gene.

    ``matches_by_protein`` is keyed on whatever ID the FASTA used -- for a
    GMB-prepared review that is the protein SHA-256, so a protein shared by
    several transcripts resolves once and is reported for each without
    being treated as independent corroboration.

    ``cfg`` is an ``InterProResolverConfig`` (or ``None``, which falls back
    to that dataclass's defaults -- i.e. report-only behaviour is
    impossible to accidentally skip, since ``apply_replacements`` still
    defaults True but replacement additionally requires the caller to have
    decided ``enabled=True`` before this function is even reached; see
    ``main()``). Every row always carries the full initial/recommended/
    final distinction (Part 11) regardless of ``cfg`` -- when replacement
    is not applied, ``final_canonical_transcript_id`` simply equals
    ``initial_canonical_transcript_id``.
    """
    if cfg is None:
        from gmb.pipeline.config import InterProResolverConfig

        cfg = InterProResolverConfig()

    by_gene: dict = defaultdict(list)
    for row in manifest_rows:
        by_gene[row["gene_id"]].append(row)

    report = []
    for gene_id in sorted(by_gene):
        candidates = sorted(by_gene[gene_id], key=lambda r: r.get("current_rank") or 0)
        if len(candidates) < 2:
            continue
        current, runner_up = candidates[0], candidates[1]

        summaries = {}
        for cand in (current, runner_up):
            checksum = cand.get("protein_sha256")
            summaries[cand["transcript_id"]] = summarise_architecture(
                matches_by_protein.get(checksum, []),
                protein_length=cand.get("protein_length"),
            )

        cur_summary = summaries[current["transcript_id"]]
        run_summary = summaries[runner_up["transcript_id"]]

        decision = decide_canonical(gene_id, current, runner_up, cur_summary, run_summary, cfg)

        report.append(
            {
                "gene_id": gene_id,
                "initial_canonical_transcript_id": decision["initial_canonical_transcript_id"],
                "runner_up_transcript_id": runner_up["transcript_id"],
                "initial_confidence": current.get("canonical_confidence"),
                "initial_selection_reason": current.get("current_selection_reason"),
                "review_reason": current.get("review_reason"),
                "current_n_integrated_entries": cur_summary["n_integrated_entries"],
                "runner_up_n_integrated_entries": run_summary["n_integrated_entries"],
                "current_integrated_entries": ",".join(cur_summary["integrated_entries"]),
                "runner_up_integrated_entries": ",".join(run_summary["integrated_entries"]),
                "current_representative_coverage": cur_summary["representative_coverage_fraction"],
                "runner_up_representative_coverage": run_summary[
                    "representative_coverage_fraction"
                ],
                "current_architecture": ";".join(cur_summary["architecture"]),
                "runner_up_architecture": ";".join(run_summary["architecture"]),
                "current_has_negative_evidence": cur_summary["has_negative_evidence"],
                "runner_up_has_negative_evidence": run_summary["has_negative_evidence"],
                "interpro_verdict": decision["interpro_verdict"],
                "reason_code": decision["reason_code"],
                "interpro_recommended_transcript_id": decision["interpro_recommended_transcript_id"],
                "final_canonical_transcript_id": decision["final_canonical_transcript_id"],
                "replaced": decision["replaced"],
                "recommendation_differs_from_initial": (
                    decision["interpro_recommended_transcript_id"]
                    != decision["initial_canonical_transcript_id"]
                ),
                "safeguards_passed": decision["safeguards_passed"],
                "safeguard_notes": decision["safeguard_notes"],
                "explanation": decision["explanation"],
            }
        )
    return report


# ---------------------------------------------------------------------------
# Attribution outputs: canonical_decisions.tsv + final-canonical GFF3
# ---------------------------------------------------------------------------
#
# canonical_selection.py's own canonical_transcripts.tsv/transcript_ranking.tsv
# and consensus.gff3 are NEVER mutated in place, anywhere in this project
# (see canonical_selection._write_annotated_gff3's own long-standing
# documented invariant) -- the functions below only ever WRITE NEW derived
# files, whether or not any InterPro replacement actually happened. This is
# what keeps "final = initial when disabled" true at the file-content level,
# not just conceptually: a build with the resolver disabled or run in
# report-only mode still gets a canonical_decisions.tsv (trivial rows,
# canonical_selection_stage="initial" throughout) and, if requested, an
# annotated GFF3 whose canonical=1 flags exactly match the pre-existing
# initial-selection-only annotated GFF3.


def build_canonical_decisions(canonical_rows: list, resolver_report: list) -> list:
    """Merge canonical_selection's per-gene initial picks with the
    resolver's per-gene decisions into ONE row per gene, covering the WHOLE
    geneset -- not just the genes that were ambiguous enough to be
    reviewed. Part 13's "exactly one final canonical per gene" contract
    needs every gene represented here, with non-reviewed genes trivially
    carrying ``canonical_selection_stage="initial"`` and
    ``final_canonical_transcript_id == initial_canonical_transcript_id``.

    ``canonical_rows`` is canonical_transcripts.tsv's own rows (as dicts,
    e.g. ``pd.read_csv(..., sep="\\t").to_dict(orient="records")``);
    ``resolver_report`` is ``build_resolver_report``'s output. Raises if the
    two disagree about a reviewed gene's initial canonical transcript --
    that would mean the resolver was run against a stale/mismatched
    canonical_transcripts.tsv, which must never be silently tolerated.
    """
    decisions_by_gene = {d["gene_id"]: d for d in resolver_report}
    rows = []
    for canon_row in canonical_rows:
        gene_id = canon_row["gene_id"]
        initial_id = canon_row["canonical_transcript_id"]
        decision = decisions_by_gene.get(gene_id)
        if decision is None:
            rows.append(
                {
                    "gene_id": gene_id,
                    "initial_canonical_transcript_id": initial_id,
                    "initial_selection_reason": canon_row.get("selection_reason"),
                    "canonical_selection_stage": "initial",
                    "interpro_verdict": None,
                    "interpro_recommended_transcript_id": None,
                    "reason_code": canon_row.get("selection_reason"),
                    "replaced": False,
                    "final_canonical_transcript_id": initial_id,
                    "safeguards_passed": None,
                    "safeguard_notes": "",
                    "explanation": "gene not selected for InterPro review",
                }
            )
            continue

        if decision["initial_canonical_transcript_id"] != initial_id:
            raise ValueError(
                f"canonical_decisions: gene {gene_id!r} initial-canonical mismatch -- "
                f"canonical_transcripts.tsv says {initial_id!r}, the resolver report "
                f"(built from a manifest presumably derived from the same file) says "
                f"{decision['initial_canonical_transcript_id']!r}. Refusing to silently "
                "reconcile; the resolver was likely run against a stale or "
                "mismatched canonical_transcripts.tsv/transcript_ranking.tsv pair."
            )
        rows.append(
            {
                "gene_id": gene_id,
                "initial_canonical_transcript_id": initial_id,
                "initial_selection_reason": canon_row.get("selection_reason"),
                "canonical_selection_stage": "interpro",
                "interpro_verdict": decision["interpro_verdict"],
                "interpro_recommended_transcript_id": decision["interpro_recommended_transcript_id"],
                "reason_code": decision["reason_code"],
                "replaced": decision["replaced"],
                "final_canonical_transcript_id": decision["final_canonical_transcript_id"],
                "safeguards_passed": decision["safeguards_passed"],
                "safeguard_notes": decision["safeguard_notes"],
                "explanation": decision["explanation"],
            }
        )
    return rows


def write_canonical_decisions(rows: list, output_dir: str) -> str:
    """Write canonical_decisions.tsv (one row per gene). Returns the path."""
    import pandas as pd

    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, "canonical_decisions.tsv")
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def write_final_annotated_gff3(input_gff3: str, decision_rows: list, output_path: str) -> None:
    """Write a COPY of ``input_gff3`` with FINAL-canonical attributes added
    to every mRNA row: ``canonical=1/0`` (based on the FINAL canonical
    transcript, which equals the initial one unless InterPro replaced it),
    ``canonical_selection_stage=initial|interpro``,
    ``initial_canonical=<ID>``, ``canonical_reason=<code>``. The source
    file is never modified in place -- mirrors
    ``gmb.pipeline.canonical_selection._write_annotated_gff3``'s own
    documented pattern (kept as a separate function/output file rather than
    replacing that one, so the InterPro-unaware initial-only annotation
    keeps working unchanged for a build that never touches this module at
    all).
    """
    decisions_by_gene = {d["gene_id"]: d for d in decision_rows}
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
            gene_id = None
            for kv in attrs.split(";"):
                if kv.startswith("ID="):
                    tid = kv[3:]
                elif kv.startswith("Parent="):
                    gene_id = kv[7:]
            decision = decisions_by_gene.get(gene_id)
            if decision is None:
                # Gene absent from canonical_decisions.tsv entirely (should
                # not happen for a consistent build) -- leave the row
                # untouched rather than guessing at a canonical flag.
                fh_out.write(line)
                continue
            flag = "1" if tid == decision["final_canonical_transcript_id"] else "0"
            extra = (
                f"canonical={flag};"
                f"canonical_selection_stage={decision['canonical_selection_stage']};"
                f"initial_canonical={decision['initial_canonical_transcript_id']};"
                f"canonical_reason={decision['reason_code']}"
            )
            sep = "" if attrs.endswith(";") else ";"
            parts[8] = f"{attrs}{sep}{extra}"
            fh_out.write("\t".join(parts) + "\n")


def main() -> None:
    import argparse

    import pandas as pd

    from gmb.pipeline.config import load_config

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest", required=True, help="interpro_review_manifest.tsv from the review step"
    )
    parser.add_argument(
        "--interpro-jsonl",
        default=None,
        help="InterProScan 6 JSONL output (preferred machine-readable input)",
    )
    parser.add_argument(
        "--interpro-gff3",
        default=None,
        help="InterProScan 6 GFF3 output (used only if --interpro-jsonl is absent)",
    )
    parser.add_argument(
        "--config",
        action="append",
        default=None,
        help="YAML with a canonical_selection.interpro_resolver: section. May be "
        "repeated to layer several files (each later --config overrides earlier ones).",
    )
    parser.add_argument("--preset", default="fungi")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--canonical-transcripts",
        default=None,
        help="canonical_transcripts.tsv from canonical selection -- when given (together "
        "with --consensus-gff3), also writes canonical_decisions.tsv and a final-canonical "
        "annotated GFF3 covering the WHOLE geneset, not just reviewed genes.",
    )
    parser.add_argument(
        "--consensus-gff3",
        default=None,
        help="GMB's consensus.gff3 for the same run -- source file is never modified in "
        "place; a new consensus.final_canonical_annotated.gff3 copy is written instead.",
    )
    args = parser.parse_args()

    if not args.interpro_jsonl and not args.interpro_gff3:
        parser.error("one of --interpro-jsonl or --interpro-gff3 is required")

    config = load_config(args.config, args.preset)
    cfg = config.canonical_selection.interpro_resolver

    if args.interpro_jsonl:
        matches = parse_interproscan_jsonl(args.interpro_jsonl)
        source = args.interpro_jsonl
    else:
        matches = parse_interproscan_gff3(args.interpro_gff3)
        source = args.interpro_gff3

    manifest_rows = pd.read_csv(args.manifest, sep="\t").to_dict(orient="records")
    report = build_resolver_report(manifest_rows, matches, cfg)

    os.makedirs(args.output_dir, exist_ok=True)
    report_path = os.path.join(args.output_dir, "interpro_resolver_report.tsv")
    pd.DataFrame(report).to_csv(report_path, sep="\t", index=False)

    verdicts: dict = defaultdict(int)
    reason_counts: dict = defaultdict(int)
    n_differ = 0
    n_replaced = 0
    for row in report:
        verdicts[row["interpro_verdict"]] += 1
        reason_counts[row["reason_code"]] += 1
        if row["recommendation_differs_from_initial"]:
            n_differ += 1
        if row["replaced"]:
            n_replaced += 1

    # Attribution outputs (Part 12): only produced when the caller supplies
    # both source files -- optional, since the resolver report alone is a
    # complete, valid output on its own (e.g. for evaluating the evidence
    # model before trusting it to touch attribution at all).
    canonical_outputs_modified = False
    decisions_path = None
    annotated_gff3_path = None
    if args.canonical_transcripts and args.consensus_gff3:
        canonical_rows = pd.read_csv(args.canonical_transcripts, sep="\t").to_dict(
            orient="records"
        )
        decision_rows = build_canonical_decisions(canonical_rows, report)
        decisions_path = write_canonical_decisions(decision_rows, args.output_dir)
        annotated_gff3_path = os.path.join(
            args.output_dir, "consensus.final_canonical_annotated.gff3"
        )
        write_final_annotated_gff3(args.consensus_gff3, decision_rows, annotated_gff3_path)
        # True here means "new derived attribution files were written this
        # run" -- consensus.gff3/canonical_transcripts.tsv themselves are
        # NEVER modified in place, by this or any other GMB module.
        canonical_outputs_modified = True

    summary = {
        "interpro_source": source,
        "proteins_with_matches": len(matches),
        "genes_reviewed": len(report),
        "verdicts": dict(verdicts),
        "reason_code_counts": dict(reason_counts),
        "recommendations_differing_from_initial": n_differ,
        "genes_with_canonical_replaced": n_replaced,
        # See the comment above: True only means canonical_decisions.tsv /
        # consensus.final_canonical_annotated.gff3 were written as NEW
        # files this run -- the original consensus.gff3 and
        # canonical_transcripts.tsv are never modified in place.
        "canonical_outputs_modified": canonical_outputs_modified,
        "canonical_decisions_path": decisions_path,
        "final_annotated_gff3_path": annotated_gff3_path,
        "apply_replacements_enabled": bool(cfg.apply_replacements),
        "resolver_enabled": bool(cfg.enabled),
    }
    summary_path = os.path.join(args.output_dir, "interpro_resolver_summary.json")
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)

    print(f"Reviewed {len(report)} gene(s) from {source}")
    for verdict, count in sorted(verdicts.items()):
        print(f"  {verdict}: {count}")
    for reason, count in sorted(reason_counts.items()):
        print(f"  {reason}: {count}")
    print(f"  recommendations differing from initial: {n_differ}")
    if decisions_path:
        print(f"  {decisions_path}")
    if annotated_gff3_path:
        print(f"  {annotated_gff3_path}")
    print(f"  genes with canonical replaced: {n_replaced}")
    print(f"  {report_path}")
    print(f"  {summary_path}")


if __name__ == "__main__":
    main()
