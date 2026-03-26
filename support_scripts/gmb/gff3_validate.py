#!/usr/bin/env python3
"""GFF3 structural validation, auto-fix, and UTR trimming.

Enforces basic GFF3 hierarchy invariants:
  - gene span covers all mRNA children
  - mRNA span covers all exon/CDS/UTR children
  - every CDS and UTR segment is within an exon
  - exons cover the full transcript structure

Policies (configurable):
  - error:            fail-fast with diagnostics
  - fix:              recompute spans, synthesize missing exons
  - drop_transcript:  remove violating transcripts

Also provides UTR trimming against biological length limits.
"""

import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------


def _merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping intervals into a disjoint set."""
    if not intervals:
        return []
    sorted_ivs = sorted(intervals)
    merged = [sorted_ivs[0]]
    for current in sorted_ivs[1:]:
        last = merged[-1]
        if current[0] <= last[1]:
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)
    return merged


def _interval_in_union(
    union_intervals: List[Tuple[int, int]], target_s: int, target_e: int
) -> bool:
    """Check if target interval is fully contained within ANY SINGLE disjoint interval from the union."""
    for cs, ce in union_intervals:
        if cs <= target_s and ce >= target_e:
            return True
    return False


def _intersect_with_union(union_intervals: List[Tuple[int, int]], row: dict) -> List[dict]:
    """Intersect a single row (like CDS/UTR) with the union of intervals.
    Returns a list of split/trimmed rows bounded strictly within the union geometries."""
    target_s, target_e = row["Start"], row["End"]
    results = []
    for cs, ce in union_intervals:
        s = max(cs, target_s)
        e = min(ce, target_e)
        if s < e:
            new_row = dict(row)
            new_row["Start"] = s
            new_row["End"] = e
            results.append(new_row)
    return results


def validate_transcript(mrna_row, exon_rows, cds_rows, utr_rows, config=None, runtime_params=None):
    """Validate structural integrity of a single transcript.

    Returns
    -------
    list of str : violation messages (empty = valid)
    int : max_drift
    float : transcript_span
    """
    violations = []
    mrna_s, mrna_e = mrna_row["Start"], mrna_row["End"]
    exon_intervals = [(r["Start"], r["End"]) for r in exon_rows]
    exon_union = _merge_intervals(exon_intervals)

    exon_min = exon_union[0][0] if exon_union else 0
    exon_max = exon_union[-1][1] if exon_union else 0

    max_exon_len = 0
    for cs, ce in exon_union:
        if (ce - cs) > max_exon_len:
            max_exon_len = ce - cs

    # Max exon length constraint
    eff_max_exon_len = (
        config.max_exon_len_bp if config and hasattr(config, "max_exon_len_bp") else 15000
    )
    eff_max_span = (
        config.max_transcript_span_bp
        if config and hasattr(config, "max_transcript_span_bp")
        else 100000
    )

    if runtime_params:
        eff_max_exon_len = runtime_params.get("effective_max_exon_len_bp", eff_max_exon_len)
        eff_max_span = runtime_params.get("effective_max_transcript_span_bp", eff_max_span)

    if max_exon_len > eff_max_exon_len:
        violations.append(
            f"Merged exon length {max_exon_len} exceeds maximum allowed ({eff_max_exon_len} bp)"
        )

    transcript_span = 0
    if exon_union:
        transcript_span = exon_max - exon_min
        if (
            config
            and hasattr(config, "max_transcript_span_mode")
            and config.max_transcript_span_mode != "off"
        ):
            if transcript_span > eff_max_span:
                violations.append(
                    f"Transcript span {transcript_span} exceeds maximum allowed ({eff_max_span} bp)"
                )

    max_drift = 0
    for r in cds_rows + utr_rows:
        if r["Start"] < exon_min:
            max_drift = max(max_drift, exon_min - r["Start"])
        if r["End"] > exon_max:
            max_drift = max(max_drift, r["End"] - exon_max)

    # All child features
    all_children = exon_rows + cds_rows + utr_rows
    if all_children:
        child_min = min(r["Start"] for r in all_children)
        child_max = max(r["End"] for r in all_children)
        if mrna_s > child_min:
            violations.append(f"mRNA start {mrna_s} > child min {child_min}")
        if mrna_e < child_max:
            violations.append(f"mRNA end {mrna_e} < child max {child_max}")

    # CDS within single merged exon segment
    for r in cds_rows:
        if not _interval_in_union(exon_union, r["Start"], r["End"]):
            violations.append(f"CDS {r['Start']}-{r['End']} not within any single merged exon")

    # UTR within single merged exon segment
    for r in utr_rows:
        if not _interval_in_union(exon_union, r["Start"], r["End"]):
            violations.append(f"UTR {r['Start']}-{r['End']} not within any single merged exon")

    return violations, max_drift, transcript_span


def validate_gene(gene_row, mrna_rows):
    """Validate gene span covers all mRNA children.

    Returns
    -------
    list of str : violation messages (empty = valid)
    """
    violations = []
    if not mrna_rows:
        return violations
    gene_s, gene_e = gene_row["Start"], gene_row["End"]
    mrna_min = min(r["Start"] for r in mrna_rows)
    mrna_max = max(r["End"] for r in mrna_rows)
    if gene_s > mrna_min:
        violations.append(f"gene start {gene_s} > mRNA min {mrna_min}")
    if gene_e < mrna_max:
        violations.append(f"gene end {gene_e} < mRNA max {mrna_max}")
    return violations


# ---------------------------------------------------------------------------
# Fix helpers
# ---------------------------------------------------------------------------


def fix_transcript(mrna_row, exon_rows, cds_rows, utr_rows):
    """Auto-fix a transcript by recomputing spans and synthesizing missing exons.

    Mutates the rows in-place.

    Returns
    -------
    list of dict : any synthesized exon rows to append
    """
    synth_exons = []

    # We never synthesize new exons from CDS. UTR never did.
    # Exons must be strictly provided by candidates.

    # Recompute mRNA span from all children (exons only effectively, because utr/cds should be contained)
    all_children = exon_rows + cds_rows + utr_rows
    if all_children:
        mrna_row["Start"] = min(r["Start"] for r in all_children)
        mrna_row["End"] = max(r["End"] for r in all_children)

    return synth_exons


def fix_gene(gene_row, mrna_rows):
    """Recompute gene span from mRNA children. Mutates in-place."""
    if mrna_rows:
        gene_row["Start"] = min(r["Start"] for r in mrna_rows)
        gene_row["End"] = max(r["End"] for r in mrna_rows)


# ---------------------------------------------------------------------------
# UTR trimming
# ---------------------------------------------------------------------------


def trim_utrs(utr_5p_rows, utr_3p_rows, cds_rows, exon_rows, config):
    """Apply UTR length constraints.

    Parameters
    ----------
    utr_5p_rows : list of dict
    utr_3p_rows : list of dict
    cds_rows : list of dict
    exon_rows : list of dict
    config : UtrConfig

    Returns
    -------
    (trimmed_5p, trimmed_3p, was_trimmed) : (list, list, bool)
    """
    was_trimmed = False

    # Clip strictly within the transcript's exon span
    exon_min = min(r["Start"] for r in exon_rows) if exon_rows else 0
    exon_max = max(r["End"] for r in exon_rows) if exon_rows else 0

    def _clip_to_span(rows):
        clipped = []
        for r in rows:
            if r["End"] <= exon_min or r["Start"] >= exon_max:
                continue
            if r["Start"] < exon_min or r["End"] > exon_max:
                clipped_r = dict(r)
                clipped_r["Start"] = max(r["Start"], exon_min)
                clipped_r["End"] = min(r["End"], exon_max)
                clipped.append(clipped_r)
            else:
                clipped.append(r)
        return clipped

    utr_5p_rows = _clip_to_span(utr_5p_rows)
    utr_3p_rows = _clip_to_span(utr_3p_rows)

    cds_bp = sum(r["End"] - r["Start"] for r in cds_rows) if cds_rows else 0

    def _total_bp(rows):
        return sum(r["End"] - r["Start"] for r in rows)

    def _hard_cap(rows, max_bp):
        """Trim UTR rows to a hard cap on total bp, trimming from the distal end."""
        if max_bp <= 0 or not rows:
            return []
        remaining = max_bp
        result = []
        for r in rows:
            seg_len = r["End"] - r["Start"]
            if remaining <= 0:
                break
            if seg_len <= remaining:
                result.append(r)
                remaining -= seg_len
            else:
                trimmed = dict(r)
                trimmed["End"] = trimmed["Start"] + remaining
                result.append(trimmed)
                remaining = 0
        return result

    utr_5p_bp = _total_bp(utr_5p_rows)
    utr_3p_bp = _total_bp(utr_3p_rows)

    if config.trim_policy == "drop_utrs":
        if (
            utr_5p_bp > config.max_5p_bp
            or utr_3p_bp > config.max_3p_bp
            or utr_5p_bp + utr_3p_bp > config.max_total_bp
            or (cds_bp > 0 and (utr_5p_bp + utr_3p_bp) / cds_bp > config.max_utr_to_cds_ratio)
        ):
            return [], [], True
        return utr_5p_rows, utr_3p_rows, False

    # hard_cap policy
    out_5p = list(utr_5p_rows)
    out_3p = list(utr_3p_rows)

    if _total_bp(out_5p) > config.max_5p_bp:
        out_5p = _hard_cap(out_5p, config.max_5p_bp)
        was_trimmed = True
    if _total_bp(out_3p) > config.max_3p_bp:
        out_3p = _hard_cap(out_3p, config.max_3p_bp)
        was_trimmed = True

    total = _total_bp(out_5p) + _total_bp(out_3p)
    if total > config.max_total_bp:
        # Proportionally reduce
        ratio = config.max_total_bp / total
        out_5p = _hard_cap(out_5p, int(_total_bp(out_5p) * ratio))
        out_3p = _hard_cap(out_3p, int(_total_bp(out_3p) * ratio))
        was_trimmed = True

    if (
        cds_bp > 0
        and (_total_bp(out_5p) + _total_bp(out_3p)) / cds_bp > config.max_utr_to_cds_ratio
    ):
        max_utr = int(cds_bp * config.max_utr_to_cds_ratio)
        ratio = (
            max_utr / (_total_bp(out_5p) + _total_bp(out_3p))
            if (_total_bp(out_5p) + _total_bp(out_3p)) > 0
            else 0
        )
        out_5p = _hard_cap(out_5p, int(_total_bp(out_5p) * ratio))
        out_3p = _hard_cap(out_3p, int(_total_bp(out_3p) * ratio))
        was_trimmed = True

    return out_5p, out_3p, was_trimmed


# ---------------------------------------------------------------------------
# Main validation pass
# ---------------------------------------------------------------------------


def validate_and_fix_gff3(gff_rows, config, runtime_params=None):
    """Validate and optionally fix GFF3 rows.

    Parameters
    ----------
    gff_rows : list of dict
        GFF3 row dicts with Feature, Start, End, ID, Parent keys.
    config : PipelineConfig
        Must have .validation and .utr attributes.
    runtime_params : dict, optional
        Effective runtime guardrails.

    Returns
    -------
    list of dict : validated/fixed rows
    dict : stats about violations found/fixed/dropped
    """
    val_cfg = config.validation
    utr_cfg = config.utr
    mode = val_cfg.mode

    stats = {
        "genes_checked": 0,
        "transcripts_checked": 0,
        "violations_found": 0,
        "transcripts_fixed": 0,
        "transcripts_dropped": 0,
        "exons_synthesized": 0,
        "utrs_trimmed": 0,
        "transcripts_dropped_drift": 0,
        "transcripts_trimmed_utr": 0,
        "max_utr_length_observed": 0,
        "containment_failures": 0,
        "worst_offenders": [],
    }

    # Index rows by parent
    by_parent = {}
    genes = []
    for r in gff_rows:
        feat = r.get("Feature", "")
        if feat == "gene":
            genes.append(r)
        parent = r.get("Parent", "")
        if parent:
            by_parent.setdefault(parent, []).append(r)

    output_rows = []
    drop_tids = set()

    for gene_row in genes:
        gene_id = gene_row["ID"]
        mrna_rows = [r for r in by_parent.get(gene_id, []) if r["Feature"] == "mRNA"]

        valid_mrna_rows = []
        for mrna_row in mrna_rows:
            tid = mrna_row["ID"]
            children = by_parent.get(tid, [])
            exon_rows = [r for r in children if r["Feature"] == "exon"]
            cds_rows = [r for r in children if r["Feature"] == "CDS"]
            utr_5p = [r for r in children if r["Feature"] == "five_prime_UTR"]
            utr_3p = [r for r in children if r["Feature"] == "three_prime_UTR"]
            utr_rows = utr_5p + utr_3p

            stats["transcripts_checked"] += 1

            # Monitor max UTR length
            utr_len = sum(r["End"] - r["Start"] for r in utr_rows)
            if utr_len > stats["max_utr_length_observed"]:
                stats["max_utr_length_observed"] = utr_len

            # UTR trimming (applied regardless of validation mode)
            if cds_rows and (utr_5p or utr_3p):
                utr_5p, utr_3p, was_trimmed = trim_utrs(
                    utr_5p, utr_3p, cds_rows, exon_rows, utr_cfg
                )
                if was_trimmed:
                    stats["utrs_trimmed"] += 1
                    stats["transcripts_trimmed_utr"] += 1
                utr_rows = utr_5p + utr_3p

            exon_intervals = [(r["Start"], r["End"]) for r in exon_rows]
            exon_union = _merge_intervals(exon_intervals)

            # Feature outside exons policy (trim CDS/UTRs to exon boundaries)
            policy = getattr(val_cfg, "feature_outside_exons_policy", "drop")

            if policy == "trim":
                # Trim CDS and UTRs to exon union
                # _intersect_with_union may split crossing features, creating multiple rows
                # or removing them completely if they have no overlap
                new_cds = []
                for r in cds_rows:
                    intersections = _intersect_with_union(exon_union, r)
                    if len(intersections) < 1:
                        stats["containment_failures"] += 1
                    for ii, iv in enumerate(intersections):
                        iv["ID"] = f"{r['ID']}_{ii}" if ii > 0 else r["ID"]
                        new_cds.append(iv)

                new_utr_5p = []
                for r in utr_5p:
                    intersections = _intersect_with_union(exon_union, r)
                    for ii, iv in enumerate(intersections):
                        iv["ID"] = f"{r['ID']}_{ii}" if ii > 0 else r["ID"]
                        new_utr_5p.append(iv)

                new_utr_3p = []
                for r in utr_3p:
                    intersections = _intersect_with_union(exon_union, r)
                    for ii, iv in enumerate(intersections):
                        iv["ID"] = f"{r['ID']}_{ii}" if ii > 0 else r["ID"]
                        new_utr_3p.append(iv)

                cds_rows = new_cds
                utr_5p = new_utr_5p
                utr_3p = new_utr_3p
                utr_rows = utr_5p + utr_3p

            violations, max_drift, _ = validate_transcript(
                mrna_row,
                exon_rows,
                cds_rows,
                utr_rows,
                config=val_cfg,
                runtime_params=runtime_params,
            )

            local_mode = mode
            if (
                hasattr(val_cfg, "max_feature_drift_bp")
                and max_drift > val_cfg.max_feature_drift_bp
            ):
                local_mode = "drop_transcript"
                stats["transcripts_dropped_drift"] += 1
                # Save to worst offenders
                stats["worst_offenders"].append(
                    {
                        "id": tid,
                        "drift": max_drift,
                        "coords": f"{mrna_row['Chromosome']}:{mrna_row['Start']}-{mrna_row['End']}",
                    }
                )

            if violations:
                stats["violations_found"] += len(violations)
                # Any un-fixed violations are containment failures
                if local_mode != "fix" and local_mode != "drop_transcript":
                    stats["containment_failures"] += 1

                if local_mode == "error" or (
                    policy == "error"
                    and any("not within any single merged exon" in v for v in violations)
                ):
                    msg = f"GFF3 validation error in {tid}:\n"
                    msg += "\n".join(f"  - {v}" for v in violations)
                    print(msg, file=sys.stderr)
                    sys.exit(1)

                elif local_mode == "drop_transcript" or policy == "drop":
                    stats["transcripts_dropped"] += 1
                    drop_tids.add(tid)
                    if val_cfg.log_violations:
                        print(f"  Dropping {tid}: {'; '.join(violations)}")
                    continue

                elif local_mode == "fix":
                    synth = fix_transcript(mrna_row, exon_rows, cds_rows, utr_rows)
                    stats["transcripts_fixed"] += 1
                    stats["exons_synthesized"] += len(synth)
                    exon_rows.extend(synth)

                    # Ensure it's now valid
                    verify_violations, _, _ = validate_transcript(
                        mrna_row,
                        exon_rows,
                        cds_rows,
                        utr_rows,
                        config=val_cfg,
                        runtime_params=runtime_params,
                    )
                    if verify_violations:
                        stats["containment_failures"] += 1
                        if policy == "drop":
                            # Drop if we failed to fix it properly and policy dictates dropping uncontained
                            stats["transcripts_dropped"] += 1
                            drop_tids.add(tid)
                            continue
                        elif policy == "error":
                            sys.exit(
                                f"Failed to fix validation violations for {tid}: {'; '.join(verify_violations)}"
                            )

                    if val_cfg.log_violations:
                        print(f"  Fixed {tid}: {'; '.join(violations)}")

            valid_mrna_rows.append(mrna_row)
            # Rebuild children list
            rebuilt_children = exon_rows + cds_rows + utr_5p + utr_3p
            output_rows.append(mrna_row)
            output_rows.extend(rebuilt_children)

        if valid_mrna_rows:
            stats["genes_checked"] += 1
            # Validate and fix gene span
            gene_violations = validate_gene(gene_row, valid_mrna_rows)
            if gene_violations and mode == "fix":
                fix_gene(gene_row, valid_mrna_rows)
            output_rows.insert(output_rows.index(valid_mrna_rows[0]), gene_row)

    # Sort worst offenders descending by drift
    stats["worst_offenders"].sort(key=lambda x: x["drift"], reverse=True)
    stats["worst_offenders"] = stats["worst_offenders"][:20]

    return output_rows, stats
