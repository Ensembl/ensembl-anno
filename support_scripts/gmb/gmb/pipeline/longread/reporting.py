"""Summary, rescue attribution, dropped-record, and run-manifest outputs."""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from typing import Optional

# ---------------------------------------------------------------------------
# Rescue attribution
# ---------------------------------------------------------------------------


def write_rescue_attribution_tsv(records: list[dict], output_path: str) -> int:
    """Write a per-rescued-transcript attribution TSV.

    Only rows where ``rescued_by_shortread`` is True are written.

    Columns:
    - ``transcript_id``, ``seqname``, ``strand``, ``start``, ``end``, ``n_exons``
    - ``rescue_reason`` — policy code (e.g. ``M1_junction_complete``)
    - ``rescue_sr_tids`` — pipe-separated SR transcript IDs (raw_tid when available)
    - ``rescue_sr_sources`` — pipe-separated assembler names (scallop/stringtie)

    Returns the number of rescued rows written.
    """
    n = 0
    with open(output_path, "w") as fh:
        fh.write(
            "transcript_id\tseqname\tstrand\tstart\tend\tn_exons"
            "\trescue_reason\trescue_sr_tids\trescue_sr_sources\n"
        )
        for rec in records:
            if not rec.get("rescued_by_shortread"):
                continue
            sr_tids = rec.get("rescue_sr_tids") or []
            # Use raw_tid when the tid is namespaced (format: source:idx:raw_tid)
            display_tids = []
            display_sources = []
            for tid in sr_tids:
                parts = tid.split(":", 2)
                if len(parts) == 3:
                    display_tids.append(parts[2])
                    display_sources.append(parts[0])
                else:
                    display_tids.append(tid)
                    display_sources.append("unknown")
            fh.write(
                f"{rec['transcript_id']}\t{rec['seqname']}\t{rec['strand']}\t"
                f"{rec['exons'][0][0]}\t{rec['exons'][-1][1]}\t{rec['n_exons']}\t"
                f"{rec.get('rescue_reason') or ''}\t"
                f"{'|'.join(display_tids)}\t{'|'.join(display_sources)}\n"
            )
            n += 1
    return n


# ---------------------------------------------------------------------------
# Dropped-record attribution
# ---------------------------------------------------------------------------


def write_dropped_records_tsv(dropped: list[dict], output_path: str) -> int:
    """Write per-record rejection log.

    Each row corresponds to one rejected source transcript or group member.

    Columns:
    - ``source_tid``, ``seqname``, ``start``, ``end``, ``strand``, ``n_exons``
    - ``rejection_stage`` — FILTER_SPAN, FILTER_INTRON, FILTER_STRAND, INPUT_PARSE,
                            SUPPORT, SUPPRESSION
    - ``reason_code`` — SPAN_EXCEEDS_MAX, INTRON_EXCEEDS_MAX, UNRESOLVED_STRAND,
                        MIXED_STRAND, NO_EXON_ROWS, DROPPED_LOW_SUPPORT,
                        DROPPED_SINGLE_READ_NO_RESCUE, SUPPRESSED_FRAGMENT_OF_MULTIEXON,
                        SUPPRESSED_CONTAINED, LIKELY_FRAGMENT_OF_MULTIEXON_TRANSCRIPT
    - ``metric`` — relevant numeric value (span, intron length, support count)
    - ``rescue_status`` — rescue reason code or NO_RESCUE or N/A

    Returns the number of rows written.
    """
    n = 0
    with open(output_path, "w") as fh:
        fh.write(
            "source_tid\tseqname\tstart\tend\tstrand\tn_exons"
            "\trejection_stage\treason_code\tmetric\trescue_status\n"
        )
        for rec in dropped:
            fh.write(
                f"{rec.get('source_tid', '')}\t"
                f"{rec.get('seqname', '')}\t"
                f"{rec.get('start', '')}\t"
                f"{rec.get('end', '')}\t"
                f"{rec.get('strand', '')}\t"
                f"{rec.get('n_exons', '')}\t"
                f"{rec.get('rejection_stage', '')}\t"
                f"{rec.get('reason_code', '')}\t"
                f"{rec.get('metric', '')}\t"
                f"{rec.get('rescue_status', '')}\n"
            )
            n += 1
    return n


# ---------------------------------------------------------------------------
# Summary outputs
# ---------------------------------------------------------------------------


_SUMMARY_KEYS = [
    "input_reads",
    "malformed_input_rows",
    "dropped_no_exons",
    "dropped_unresolved_strand",
    "dropped_bad_intron_or_overlap",
    "dropped_oversized_span",
    "surviving_after_basic_filters",
    "clusters_found",
    "chain_or_overlap_groups_found",
    "dropped_low_support",
    "dropped_single_read_no_rescue",
    "suppressed_fragment_of_multiexon",
    "sr_transcripts_loaded",
    "rescued_by_shortread",
    "consensus_models_before_containment_suppression",
    "suppressed_contained_in_higher_support_model",
    "consensus_models_output",
    "genes_output",
]


def write_summary(
    all_stats: list[dict],
    all_records: list[dict],
    args_dict: dict,
    cfg,
    output_dir: str,
    start_time: float,
) -> dict:
    """Write JSON and TSV summaries; return the summary dict."""
    import dataclasses

    totals = {k: sum(s.get(k, 0) for s in all_stats) for k in _SUMMARY_KEYS}

    summary = {
        "command": " ".join(sys.argv),
        "run_date": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "elapsed_s": round(time.time() - start_time, 1),
        "args": args_dict,
        "resolved_config": dataclasses.asdict(cfg),
        "per_seqname": all_stats,
        "totals": totals,
    }

    json_path = os.path.join(output_dir, "longread_consensus_summary.json")
    with open(json_path, "w") as fh:
        json.dump(summary, fh, indent=2, default=str)

    tsv_path = os.path.join(output_dir, "longread_consensus_summary.tsv")
    with open(tsv_path, "w") as fh:
        fh.write("Metric\tValue\n")
        for k, v in totals.items():
            fh.write(f"{k}\t{v}\n")

    return summary


# ---------------------------------------------------------------------------
# Run manifest
# ---------------------------------------------------------------------------


def _sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_run_manifest(
    output_dir: str,
    input_path: Optional[str],
    scallop_paths: list[str],
    stringtie_paths: list[str],
    output_files: dict[str, str],
    summary: dict,
    cfg,
    preset: Optional[str],
    start_time: float,
) -> str:
    """Write a provenance manifest JSON for the run.

    Includes command, input checksums, resolved config, output paths, and
    summary counts.
    """
    import dataclasses

    manifest = {
        "tool": "gmb-longread-consensus",
        "command": " ".join(sys.argv),
        "run_date": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "elapsed_s": round(time.time() - start_time, 1),
        "preset": preset,
        "resolved_config": dataclasses.asdict(cfg),
        "inputs": {},
        "outputs": output_files,
        "shortread_rescue_enabled": cfg.shortread_rescue.enabled,
        "short_read_files": {
            "scallop": scallop_paths,
            "stringtie": stringtie_paths,
        },
        "summary_totals": summary.get("totals", {}),
        "status": "complete",
    }

    if input_path and os.path.exists(input_path):
        manifest["inputs"]["raw_gtf"] = {
            "path": input_path,
            "size_bytes": os.path.getsize(input_path),
            "sha256": _sha256_file(input_path),
        }

    manifest_path = os.path.join(output_dir, "longread_consensus_run_manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2, default=str)
    return manifest_path
