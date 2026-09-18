#!/usr/bin/env python3
"""Helper for the gmb-new-clade-config skill.

Five subcommands, run in this order:

    args           manifest -> GMB CLI evidence arguments + slot map
    measure        reference-free evidence measurements + evidence-state flags
    freeze         hash the overlay + the build's resolved config BEFORE any reference is used
    summarise-run  reference-free hard QC + diagnostics for a finished run
    evaluate       compare a FROZEN config's run against a reference (gmb-compare output)

Design rules
------------
* This script never runs GMB itself; the agent runs the documented CLI
  (gmb-preflight / gmb-build / gmb-finalise / gmb-compare).
* It does not re-implement preflight checks or the backbone-intron-rescue gate.
  It reads their outputs. The only things computed here are distributions
  preflight does not report (intron / span / CDS-length percentiles) and
  summaries of documented build outputs.
* ``args``, ``measure`` and ``summarise-run`` accept NO reference annotation.
  Only ``evaluate`` touches reference-derived data, and it refuses to run unless
  both the overlay and the build's resolved config were frozen first and are
  byte-identical to the frozen versions. The reference must be GFF3.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shlex
import statistics
import sys
import time
from collections import Counter, defaultdict

import yaml

# --------------------------------------------------------------------------
# GMB's fixed input slots. The slot decides only the INTERNAL label; the file's
# extension decides how it is parsed; the evidence ROLE decides behaviour.
# --------------------------------------------------------------------------
ROLE_SLOTS = {
    "backbone": ["--helixer", "--tiberius"],          # exactly one backbone
    "short_read_transcriptomic": ["--scallop", "--stringtie"],
    "long_read_transcriptomic": ["--minimap2"],
    "protein_alignment": ["--orthodb", "--uniprot", "--genblast"],
}
SLOT_LABEL = {"--helixer": "Helixer", "--tiberius": "Tiberius", "--scallop": "Scallop",
              "--stringtie": "StringTie", "--minimap2": "Minimap2", "--orthodb": "OrthoDB",
              "--uniprot": "UniProt", "--genblast": "GenBlast"}
STRUCTURAL_ROLES = ("backbone", "short_read_transcriptomic", "long_read_transcriptomic")
REFERENCE_WORDS = ("reference", "evaluation", "truth", "gold")

_COMP = str.maketrans("ACGTacgtnN", "TGCAtgcaNN")

# Keys a machine_overlay may set. Derived from the current GMB config schema: every
# key here says WHERE a tool or database is, or HOW MUCH hardware to use -- none
# chooses which gene model wins. Everything else (weights, policies, thresholds,
# UTR rules, even protein_validation.enabled) is biology and belongs in the clade
# config, where it gets a recorded reason. Version pins for InterProScan
# (workflow / revision / interpro_release) are deliberately excluded: they change
# the domain annotations the resolver scores.
MACHINE_OVERLAY_ALLOWED = frozenset({
    "protein_validation.diamond_path",
    "protein_validation.psauron_path",
    "protein_validation.diamond_db",
    "protein_validation.psauron_use_cpu",
    "qc.parallel",
    "qc.workers",
    "canonical_selection.interpro_resolver.nextflow.nextflow_executable",
    "canonical_selection.interpro_resolver.nextflow.profile",
    "canonical_selection.interpro_resolver.nextflow.config_file",
    "canonical_selection.interpro_resolver.nextflow.data_dir",
    "canonical_selection.interpro_resolver.nextflow.work_dir",
    "canonical_selection.interpro_resolver.nextflow.output_dir",
    "canonical_selection.interpro_resolver.nextflow.cpus",
    "canonical_selection.interpro_resolver.nextflow.max_workers",
})


def flatten(d, prefix=""):
    """{'a': {'b': 1}} -> {'a.b': 1}. Lists and scalars are leaves."""
    out = {}
    for k, v in (d or {}).items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(flatten(v, key + "."))
        else:
            out[key] = v
    return out


def validate_machine_overlay(path: str) -> None:
    """Fail if a machine overlay contains anything but execution/tool-path keys."""
    if not os.path.exists(path):
        die(f"machine_overlay not found: {path}")
    with open(path) as fh:
        data = yaml.safe_load(fh) or {}
    if not isinstance(data, dict):
        die(f"machine_overlay {path} is not a YAML mapping")
    bad = sorted(k for k in flatten(data) if k not in MACHINE_OVERLAY_ALLOWED)
    if bad:
        die(f"machine_overlay {path} contains settings that are not machine/execution "
            f"paths: {bad}. A machine overlay may only say where tools and databases are "
            f"and how much hardware to use; it must not change model selection. Move these "
            f"to the clade config, where each change is recorded with its reason. Allowed "
            f"keys: {sorted(MACHINE_OVERLAY_ALLOWED)}")


def die(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(2)


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def opener(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def pct(values, q):
    if not values:
        return None
    v = sorted(values)
    return v[min(len(v) - 1, int(q / 100 * len(v)))]


# --------------------------------------------------------------------------
# manifest
# --------------------------------------------------------------------------

def sniff_format(path: str) -> str | None:
    """'gtf' or 'gff3' from the attribute column of the first feature line."""
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                return None
            attrs = cols[8]
            if 'transcript_id "' in attrs or 'gene_id "' in attrs:
                return "gtf"
            if "=" in attrs:
                return "gff3"
            return None
    return None


def load_manifest(path: str) -> dict:
    with open(path) as fh:
        m = yaml.safe_load(fh) or {}
    base = os.path.dirname(os.path.abspath(path))

    def resolve(p):
        return p if (p is None or os.path.isabs(p)) else os.path.normpath(os.path.join(base, p))

    m["genome"] = resolve(m.get("genome"))
    for ev in m.get("evidence", []) or []:
        ev["path"] = resolve(ev.get("path"))
    ev_only = m.get("evaluation_only") or {}
    for k in list(ev_only):
        ev_only[k] = resolve(ev_only[k])
    m["evaluation_only"] = ev_only
    if m.get("machine_overlay"):
        m["machine_overlay"] = resolve(m["machine_overlay"])
        validate_machine_overlay(m["machine_overlay"])
    return m


def assign_slots(m: dict) -> list:
    """Validate the manifest and map every evidence file onto a GMB slot."""
    problems, rows = [], []
    if not m.get("clade"):
        problems.append("manifest has no 'clade' name")
    if not m.get("genome") or not os.path.exists(m["genome"]):
        problems.append(f"genome FASTA not found: {m.get('genome')!r}")

    reference_paths = {os.path.abspath(p) for p in (m.get("evaluation_only") or {}).values() if p}
    used = defaultdict(int)
    for ev in m.get("evidence", []) or []:
        tool, role, path = ev.get("tool"), ev.get("role"), ev.get("path")
        if role not in ROLE_SLOTS:
            problems.append(
                f"'{tool}': role {role!r} is not one of {sorted(ROLE_SLOTS)}. A reference "
                f"annotation is never an evidence role -- put it under evaluation_only.")
            continue
        if not path or not os.path.exists(path):
            problems.append(f"'{tool}': file not found: {path!r}")
            continue
        if os.path.abspath(path) in reference_paths:
            problems.append(f"'{tool}': this file is also listed as the evaluation reference. "
                            f"A reference annotation must never be GMB evidence.")
            continue
        if any(w in os.path.basename(path).lower() for w in REFERENCE_WORDS):
            print(f"WARNING: '{os.path.basename(path)}' looks like a reference annotation. "
                  f"Confirm it is genuinely evidence before continuing.", file=sys.stderr)
        fmt = sniff_format(path)
        if fmt is None:
            problems.append(f"'{tool}': cannot tell whether {path} is GTF or GFF3")
            continue
        # gmb-build parses a file as GTF if and only if its name ends in ".gtf".
        if fmt == "gtf" and not path.endswith(".gtf"):
            problems.append(
                f"'{tool}': {os.path.basename(path)} contains GTF but its name does not end in "
                f"'.gtf'; gmb-build would parse it as GFF3. Decompress or rename it to *.gtf.")
            continue
        if fmt == "gff3" and path.endswith(".gtf"):
            problems.append(f"'{tool}': {os.path.basename(path)} contains GFF3 but is named "
                            f"*.gtf; rename it to *.gff3.")
            continue
        slots = ROLE_SLOTS[role]
        if role == "backbone":
            if used[role]:
                problems.append("more than one backbone supplied; GMB takes exactly one")
                continue
            slot = "--tiberius" if fmt == "gtf" else "--helixer"
        else:
            if used[role] >= len(slots):
                problems.append(
                    f"'{tool}': GMB exposes only {len(slots)} {role} slot(s) "
                    f"({', '.join(slots)}). Merge tracks upstream or drop one, and record why.")
                continue
            slot = slots[used[role]]
        used[role] += 1
        rows.append({"tool": tool, "role": role, "path": path, "format": fmt, "slot": slot,
                     "gmb_label": SLOT_LABEL[slot], "rnaseq_source": ev.get("rnaseq_source")})
    if not used["backbone"]:
        problems.append("no backbone supplied (role: backbone)")
    if problems:
        die("manifest problems:\n  - " + "\n  - ".join(problems))
    return rows


def cmd_args(a):
    m = load_manifest(a.manifest)
    rows = assign_slots(m)
    evidence = []
    for r in rows:
        evidence += [r["slot"], r["path"]]
    configs = []
    if m.get("machine_overlay"):
        configs += ["--config", m["machine_overlay"]]
    out = {"clade": m["clade"], "genome": m["genome"], "slots": rows,
           "evidence_args": evidence, "machine_overlay_args": configs}
    if a.out:
        os.makedirs(a.out, exist_ok=True)
        with open(os.path.join(a.out, "slot_map.tsv"), "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["tool", "role", "gmb_slot", "gmb_internal_label", "format",
                        "rnaseq_source", "path"])
            for r in rows:
                w.writerow([r["tool"], r["role"], r["slot"], r["gmb_label"], r["format"],
                            r["rnaseq_source"] or "", r["path"]])
        with open(os.path.join(a.out, "evidence_args.txt"), "w") as fh:
            fh.write(" ".join(shlex.quote(x) for x in configs + ["--genome", m["genome"]]
                              + evidence) + "\n")
    if a.shell:
        print(" ".join(shlex.quote(x) for x in configs + ["--genome", m["genome"]] + evidence))
    else:
        print(f"{'tool':<16}{'role':<28}{'slot':<13}{'GMB label':<11}{'fmt':<6}path")
        for r in rows:
            print(f"{r['tool']:<16}{r['role']:<28}{r['slot']:<13}{r['gmb_label']:<11}"
                  f"{r['format']:<6}{r['path']}")
        shared = Counter(r["rnaseq_source"] for r in rows
                         if r["role"] == "short_read_transcriptomic" and r["rnaseq_source"])
        for src, n in shared.items():
            if n > 1:
                print(f"\nNOTE: {n} short-read tracks derive from the same RNA-seq "
                      f"('{src}'). They are separate named sources to GMB but NOT independent "
                      f"biological evidence -- weigh their agreement accordingly.")
    return out


# --------------------------------------------------------------------------
# measure (reference-free)
# --------------------------------------------------------------------------

def track_distributions(path: str, genome_fmt: str) -> dict:
    """Intron, span and CDS-length distributions -- things preflight does not report."""
    gff3 = genome_fmt == "gff3"
    exons, cds = defaultdict(list), defaultdict(list)
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "CDS"):
                continue
            key = None
            if gff3:
                for kv in f[8].split(";"):
                    if kv.startswith("Parent="):
                        key = kv[7:].split(",")[0]
            else:
                for kv in f[8].split(";"):
                    kv = kv.strip()
                    if kv.startswith("transcript_id "):
                        key = kv.split(" ", 1)[1].strip('"')
            if key is None:
                continue
            (exons if f[2] == "exon" else cds)[(f[0], key)].append((int(f[3]), int(f[4])))
    models = exons or cds
    introns, spans = [], []
    for segs in models.values():
        s = sorted(segs)
        spans.append(s[-1][1] - s[0][0] + 1)
        introns += [s[i + 1][0] - s[i][1] - 1 for i in range(len(s) - 1)]
    cds_codons = [sum(e - b + 1 for b, e in segs) // 3 for segs in cds.values()]
    q = (50, 90, 95, 99, 99.9)
    out = {
        "introns": len(introns),
        "intron_bp": {f"p{p}": pct(introns, p) for p in q} | {"max": max(introns) if introns else None},
        "span_bp": {f"p{p}": pct(spans, p) for p in q} | {"max": max(spans) if spans else None},
    }
    if cds_codons:
        n = len(cds_codons)
        out["cds_models"] = n
        out["cds_codons_median"] = int(statistics.median(cds_codons))
        out["cds_lt_33_codons_frac"] = round(sum(c < 33 for c in cds_codons) / n, 4)
        out["cds_lt_50_codons_frac"] = round(sum(c < 50 for c in cds_codons) / n, 4)
    return out


def gate_category(probe_policy: dict) -> tuple:
    """Interpret -- never re-decide -- the real gate's result from the probe run."""
    fires = bool(probe_policy.get("backbone_intron_rescue_will_fire"))
    st = probe_policy.get("backbone_resolution") or {}
    try:
        from gmb.pipeline.applicability import (
            ASSEMBLED_CANONICAL_SPLICE_MIN as CS_MIN,
            ASSEMBLED_TO_BACKBONE_RATIO_MIN as R_MIN,
            BACKBONE_MULTI_EXON_MAX as B_MAX,
            MIN_MODELS_FOR_DECISION as N_MIN,
        )
    except Exception:
        return ("fires" if fires else "declines_unknown_reason"), None
    thresholds = {"backbone_multi_exon_max": B_MAX, "assembled_to_backbone_ratio_min": R_MIN,
                  "assembled_canonical_splice_min": CS_MIN, "min_models_each_side": N_MIN}
    if fires:
        return "fires", thresholds
    b, r = st.get("backbone_multi_exon_fraction"), st.get("assembled_to_backbone_ratio")
    cs = st.get("assembled_canonical_splice_fraction")
    if not st or (st.get("backbone_models", 0) < N_MIN or st.get("assembled_models", 0) < N_MIN):
        return "not_evaluable_insufficient_evidence", thresholds
    if cs is None or cs < CS_MIN:
        return "declines_assembled_splice_quality", thresholds
    if b is not None and b > B_MAX:
        return "declines_backbone_well_resolved", thresholds
    if r is not None and r < R_MIN:
        return "declines_backbone_not_clearly_under_resolved", thresholds
    return "declines_unknown_reason", thresholds


def cmd_measure(a):
    m = load_manifest(a.manifest)
    slots = assign_slots(m)
    label_to_row = {r["gmb_label"]: r for r in slots}
    pf = json.load(open(os.path.join(a.preflight, "preflight_report.json")))
    probe = json.load(open(os.path.join(a.probe, "preflight_report.json")))

    checks = pf.get("checks", [])
    def verdicts(name):
        return {c["target"]: c["verdict"] for c in checks if c["name"] == name}
    splice = verdicts("splice_quality")

    tracks = []
    for t in pf.get("tracks", []):
        row = label_to_row.get(t["label"], {})
        entry = {
            "tool": row.get("tool", t["label"]), "gmb_label": t["label"], "role": t["role"],
            "weight_under_standard": t.get("weight"), "models": t["models"],
            "multi_exon": t["multi_exon"], "single_exon": t["single_exon"],
            "multi_exon_fraction": t.get("multi_exon_fraction"),
            "canonical_splice_fraction": t.get("canonical_splice_fraction"),
            "introns_examined": t.get("introns_examined"),
            "median_span": t.get("median_transcript_span"), "max_span": t.get("max_transcript_span"),
            "id_collisions_across_seqids": t.get("duplicate_ids_across_seqids"),
            "unstranded_rows": sum(v for k, v in (t.get("strand_counts") or {}).items()
                                   if k not in ("+", "-")),
            "splice_verdict": splice.get(t["label"]),
            "seqid_verdict": verdicts("seqid_compatibility").get(t["label"]),
            "strand_verdict": verdicts("strand_completeness").get(t["label"]),
            "role_verdict": verdicts("role_resolved").get(t["label"]),
        }
        if t["role"] in STRUCTURAL_ROLES and row:
            entry["distributions"] = track_distributions(row["path"], row["format"])
        tracks.append(entry)

    by_role = defaultdict(list)
    for t in tracks:
        by_role[t["role"]].append(t)

    gate, thresholds = gate_category(probe.get("policy", {}))
    probe_fires = bool(probe.get("policy", {}).get("backbone_intron_rescue_will_fire"))

    flags, why = [], {}
    if probe_fires:
        flags.append("B"); why["B"] = "the implemented reference-free rescue gate FIRES"
    if gate in ("declines_backbone_well_resolved", "declines_backbone_not_clearly_under_resolved"):
        flags.append("A"); why["A"] = f"the rescue gate declines: {gate}"
    sr = by_role.get("short_read_transcriptomic", [])
    bad_sr = [t["tool"] for t in sr if t["splice_verdict"] in ("WARN", "FAIL")]
    if not sr or bad_sr:
        flags.append("D")
        why["D"] = ("no short-read transcript evidence" if not sr
                    else f"short-read splice quality WARN/FAIL: {bad_sr}")
    lr = by_role.get("long_read_transcriptomic", [])
    lr_fail = [t["tool"] for t in lr if t["splice_verdict"] == "FAIL"]
    lr_warn = [t["tool"] for t in lr if t["splice_verdict"] == "WARN"]
    if lr_fail:
        flags.append("E"); why["E"] = f"long-read splice quality FAIL: {lr_fail}"
    elif lr_warn:
        flags.append("E?"); why["E?"] = f"long-read splice quality WARN (borderline): {lr_warn}"
    sparse = [c["name"] for c in checks if c["name"] in
              ("backbone_present", "transcript_evidence_present", "protein_evidence_present")
              and c["verdict"] != "PASS"]
    if sparse or gate == "not_evaluable_insufficient_evidence":
        flags.append("F"); why["F"] = f"sparse evidence: {sparse or gate}"

    out = {
        "clade": m["clade"], "genome": m["genome"],
        "preflight_verdict": pf.get("verdict"), "preflight_counts": pf.get("counts"),
        "tracks": tracks,
        "rescue_gate": {"will_fire": probe_fires, "category": gate,
                        "reason": probe.get("policy", {}).get("backbone_intron_rescue_reason"),
                        "evidence": probe.get("policy", {}).get("backbone_resolution"),
                        "thresholds": thresholds},
        "evidence_state_flags": flags, "flag_reasons": why,
        "preflight_recommendations": pf.get("recommendations", []),
    }

    if a.baseline_build:
        out["baseline_build"] = measure_baseline(a.baseline_build, label_to_row)
        sr_tracks = [r for r in slots if r["role"] == "short_read_transcriptomic"]
        if len(sr_tracks) < 2:
            out["baseline_build"]["utr_note"] = (
                "Only one short-read track: utr.end_support_mode=multisource_end_agreement "
                "cannot be satisfied, so UTRs will mostly be dropped. That is the SAFE outcome.")

    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "evidence_measurements.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    write_evidence_summary(out, os.path.join(a.out, "evidence_summary.md"))
    print(open(os.path.join(a.out, "evidence_summary.md")).read())


def measure_baseline(build_dir: str, label_to_row: dict) -> dict:
    """Reference-free measurements from a build made with the NEUTRAL preset."""
    s = json.load(open(os.path.join(build_dir, "summary.json"))).get("summary", {})
    rows = list(csv.DictReader(open(os.path.join(build_dir, "evidence_attribution.tsv")),
                               delimiter="\t"))
    n = len(rows) or 1
    def role(label):
        r = label_to_row.get(label)
        return r["role"] if r else "unknown"
    mix = Counter()
    for r in rows:
        roles = {role(x) for x in r["evidence_sources"].split(",") if x}
        has_bb = "backbone" in roles
        has_sr = "short_read_transcriptomic" in roles
        has_lr = "long_read_transcriptomic" in roles
        mix["backbone + short-read" if has_bb and has_sr else "backbone only" if has_bb and not
            (has_sr or has_lr) else "short-read only" if has_sr and not (has_bb or has_lr) else
            "long-read only" if has_lr and not (has_bb or has_sr) else "other combination"] += 1
    kept5 = [int(r["utr_5p_bp"]) for r in rows if r.get("utr_5p_action") == "kept" and r["utr_5p_bp"]]
    kept3 = [int(r["utr_3p_bp"]) for r in rows if r.get("utr_3p_action") == "kept" and r["utr_3p_bp"]]
    cds = [int(r["cds_bp"]) for r in rows if r.get("cds_bp")]
    ratio = [(int(r["utr_5p_bp"] or 0) + int(r["utr_3p_bp"] or 0)) / int(r["cds_bp"])
             for r in rows if r.get("cds_bp") and int(r["cds_bp"]) > 0]
    return {
        "selected_transcripts": len(rows),
        "backbone_shortread_identical_chain_pct": round(100 * sum(
            r["backbone_shortread_agreement"] == "True" for r in rows) / n, 2),
        "multi_source_structural_pct": round(100 * sum(
            int(r["n_structural_support_sources"] or 0) > 1 for r in rows) / n, 2),
        "selected_source_mix_pct": {k: round(100 * v / n, 1) for k, v in mix.most_common()},
        "selection_reason_pct": {k: round(100 * v / n, 1) for k, v in
                                 Counter(r["selection_reason"] for r in rows).most_common()},
        "protein_support_pct": {k: round(100 * v / n, 1) for k, v in
                                Counter(r["protein_support_strength"] or "none"
                                        for r in rows).most_common()},
        "protein_filter": {k: s.get(k) for k in (
            "protein_input_transcripts", "protein_redundant_collapsed",
            "protein_after_redundancy", "protein_competition_demoted",
            "protein_supported_candidates", "protein_cds_span_compatible_candidates")},
        "chimeras_removed_large_intron": s.get("chimeras_large_intron"),
        "validation_transcripts_dropped": s.get("validation_transcripts_dropped"),
        "utr": {
            "5p_kept": len(kept5), "3p_kept": len(kept3),
            "5p_kept_bp": {f"p{q}": pct(kept5, q) for q in (50, 90, 99)},
            "3p_kept_bp": {f"p{q}": pct(kept3, q) for q in (50, 90, 99)},
            "utr_to_cds_ratio_p99": round(pct(ratio, 99), 2) if ratio else None,
            "utrs_trimmed_by_hard_caps": s.get("validation_utrs_trimmed"),
            "max_utr_observed_before_caps": s.get("validation_max_utr_length_observed"),
        },
        "cds_bp_median": int(statistics.median(cds)) if cds else None,
        "genes_after_dedup": s.get("dedup_genes_output"),
    }


def write_evidence_summary(o: dict, path: str) -> None:
    L = [f"# Evidence measurements — {o['clade']}", "",
         "Reference-free. Nothing in this file used a reference annotation.", "",
         f"Preflight (standard preset): **{o['preflight_verdict']}** {o['preflight_counts']}", "",
         "## Structural and protein tracks", "",
         "| tool | GMB label | role | models | multi-exon | canonical splice | "
         "intron p50/p99/max | span p50/p99 | splice | seqid | strand |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    for t in o["tracks"]:
        d = t.get("distributions") or {}
        ib, sb = d.get("intron_bp", {}), d.get("span_bp", {})
        mf = t["multi_exon_fraction"]; cs = t["canonical_splice_fraction"]
        L.append(
            f"| {t['tool']} | {t['gmb_label']} | {t['role']} | {t['models']:,} | "
            f"{'-' if mf is None else f'{100*mf:.1f}%'} | {'-' if cs is None else f'{100*cs:.1f}%'} | "
            f"{ib.get('p50','-')}/{ib.get('p99','-')}/{ib.get('max','-')} | "
            f"{sb.get('p50','-')}/{sb.get('p99','-')} | {t['splice_verdict'] or '-'} | "
            f"{t['seqid_verdict'] or '-'} | {t['strand_verdict'] or '-'} |")
    g = o["rescue_gate"]
    L += ["", "## Backbone vs assembled transcripts (the implemented rescue gate)", "",
          f"- **will fire:** {g['will_fire']}  (`{g['category']}`)",
          f"- reason: {g['reason']}", f"- measured: `{json.dumps(g['evidence'])}`",
          f"- gate thresholds: `{json.dumps(g['thresholds'])}`", "",
          "## Evidence-state flags", ""]
    for f in o["evidence_state_flags"]:
        L.append(f"- **{f}** — {o['flag_reasons'][f]}")
    if not o["evidence_state_flags"]:
        L.append("- none raised by the measured rules (see SKILL.md phase 3 for C, which is a "
                 "judgement from the agreement numbers below)")
    b = o.get("baseline_build")
    if b:
        L += ["", "## Neutral-preset baseline build (reference-free)", "",
              f"- selected transcripts: {b['selected_transcripts']:,}; genes {b['genes_after_dedup']}",
              f"- backbone ↔ short-read **identical intron chain**: "
              f"**{b['backbone_shortread_identical_chain_pct']}%**  "
              f"(indicative only, not a threshold: earlier builds measured ~1.2–1.3% in two "
              f"intron-collapsed-backbone genomes and 11.4% in a well-resolved-backbone "
              f"genome, using clade presets rather than standard)",
              f"- multi-source structural agreement: {b['multi_source_structural_pct']}%",
              f"- selected-source mix: {b['selected_source_mix_pct']}",
              f"- protein support: {b['protein_support_pct']}",
              f"- protein filter: {b['protein_filter']}",
              f"- chimeric transcripts removed for a too-long intron: {b['chimeras_removed_large_intron']}",
              f"- UTR: {b['utr']}"]
        if b.get("utr_note"):
            L.append(f"- **UTR note:** {b['utr_note']}")
    if o["preflight_recommendations"]:
        L += ["", "## Preflight recommendations", ""] + [f"- {r}" for r in o["preflight_recommendations"]]
    open(path, "w").write("\n".join(L) + "\n")


# --------------------------------------------------------------------------
# freeze / evaluate
# --------------------------------------------------------------------------

def cmd_freeze(a):
    """Freeze BOTH what the developer wrote and what GMB actually ran with.

    The resolved config is the authoritative record of the build: standard +
    clade overlay + machine overlay, as resolved by gmb-build. It is checked
    against the build's run_manifest.json (which records the overlays used and the
    resolved config's hash) and against every value the overlay sets, so a resolved
    config from some other build cannot be frozen by mistake.
    """
    if not os.path.exists(a.config):
        die(f"source overlay not found: {a.config}")
    if not os.path.exists(a.resolved_config):
        die(f"resolved config not found: {a.resolved_config}. Build the candidate first; "
            f"freeze needs <OUT>/candidate/build/resolved_config.yaml.")
    os.makedirs(a.out, exist_ok=True)
    if os.path.exists(os.path.join(a.out, "evaluation")):
        die("an evaluation/ directory already exists here: the reference has already been "
            "consulted, so a config frozen now cannot be called input-derived. Start a new "
            "output directory and a new config version.")

    overlay_path = os.path.abspath(a.config)
    resolved_path = os.path.abspath(a.resolved_config)
    resolved_sha = sha256(resolved_path)
    build_dir = os.path.dirname(resolved_path)

    # 1. the build's own provenance must name this overlay and this resolved config
    manifest = os.path.join(build_dir, "run_manifest.json")
    if os.path.exists(manifest):
        rm = json.load(open(manifest))
        used = [os.path.abspath(x) for x in (rm.get("config_overlays") or [])]
        if overlay_path not in used:
            die(f"the build in {build_dir} was not run with {overlay_path} (its overlays were "
                f"{used}). Freeze the resolved config of the build that used this overlay.")
        if rm.get("resolved_config_sha256") and rm["resolved_config_sha256"] != resolved_sha:
            die(f"{resolved_path} no longer matches the hash its build recorded in "
                f"run_manifest.json -- it was edited after the build. Rebuild.")
    else:
        print(f"WARNING: no run_manifest.json in {build_dir}; cannot confirm from provenance "
              f"that this build used the overlay. Checking values instead.", file=sys.stderr)

    # 2. every value the overlay sets must be what the build resolved
    with open(overlay_path) as fh:
        overlay = flatten(yaml.safe_load(fh) or {})
    with open(resolved_path) as fh:
        resolved = flatten(yaml.safe_load(fh) or {})
    mismatched = [k for k, v in overlay.items() if k not in resolved or resolved[k] != v]
    if mismatched:
        die(f"the resolved config does not contain the overlay's values for {mismatched}. "
            f"Either it came from a different build, or the overlay uses a deprecated key "
            f"name -- use current key names and rebuild.")

    rec = {
        "source_overlay": overlay_path,
        "source_overlay_sha256": sha256(overlay_path),
        "resolved_config": resolved_path,
        "resolved_config_sha256": resolved_sha,
        "build_dir": build_dir,
        "frozen_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "statement": ("Frozen before any reference evaluation, from reference-free "
                      "measurements only. The resolved config is the authoritative record "
                      "of what GMB actually ran with. Any change after this point is a NEW "
                      "configuration hypothesis: new version, new build, new freeze."),
    }
    path = os.path.join(a.out, "config_freeze.json")
    json.dump(rec, open(path, "w"), indent=1)
    print(f"frozen source overlay   {overlay_path}  sha256 {rec['source_overlay_sha256'][:16]}…")
    print(f"frozen resolved config  {resolved_path}  sha256 {resolved_sha[:16]}…")
    print(f"-> {path}")


def _check_frozen(freeze_path: str) -> dict:
    if not os.path.exists(freeze_path):
        die(f"no freeze record at {freeze_path}. Build the candidate, then run "
            f"`gmb_clade.py freeze <overlay> --resolved-config <build>/resolved_config.yaml "
            f"--out <dir>` BEFORE looking at any reference.")
    rec = json.load(open(freeze_path))
    if "resolved_config_sha256" not in rec:
        die(f"{freeze_path} predates resolved-config freezing and records only the overlay. "
            f"It cannot prove what the build ran with. Rebuild and freeze again in a new "
            f"output directory.")
    if not os.path.exists(rec["source_overlay"]):
        die(f"SOURCE OVERLAY missing: {rec['source_overlay']} no longer exists.")
    if sha256(rec["source_overlay"]) != rec["source_overlay_sha256"]:
        die("SOURCE OVERLAY CHANGED since it was frozen. Evaluating now would let the "
            "reference leak into the configuration. Make the change as a new config version, "
            "build it, freeze its resolved config, then evaluate.")
    if not os.path.exists(rec["resolved_config"]):
        die(f"RESOLVED BUILD CONFIG missing: {rec['resolved_config']} no longer exists, so "
            f"the annotation can no longer be tied to the frozen configuration.")
    if sha256(rec["resolved_config"]) != rec["resolved_config_sha256"]:
        die("RESOLVED BUILD CONFIG CHANGED since it was frozen. The annotation being evaluated "
            "is no longer tied to the configuration frozen before the reference was used. "
            "Rebuild, freeze again in a new output directory, then evaluate.")
    return rec


def _ref_cds_exon_counts(ref_path: str) -> dict:
    """Reference gene -> most CDS exons in any one of its transcripts.

    Counted per transcript, then maximised per gene: summing across isoforms would
    turn a gene with two single-exon isoforms into a "multi-exon" gene.
    """
    tx2g, per_tx = {}, Counter()
    with opener(ref_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] in ("mRNA", "transcript"):
                tx2g[at.get("ID", "")] = at.get("Parent", "")
            elif f[2] == "CDS":
                for t in at.get("Parent", "").split(","):
                    per_tx[t] += 1
    per_gene = Counter()
    for t, n in per_tx.items():
        g = tx2g.get(t)
        if g is None:
            continue
        g = g.split(":", 1)[-1] if ":" in g else g     # Ensembl "gene:ID" -> "ID"
        per_gene[g] = max(per_gene[g], n)
    return per_gene


def _comparison_metrics(cmp_dir: str, ref_exons: dict) -> dict:
    summ = {k: v for k, v in (l.rstrip("\n").split("\t") for l in
            open(os.path.join(cmp_dir, "comparison_summary.tsv")) if "\t" in l)}
    rows = list(csv.DictReader(open(os.path.join(cmp_dir, "comparison_details.tsv")),
                               delimiter="\t"))
    ref = [r for r in rows if r["source"] == "reference"]
    cons = [r for r in rows if r["source"] == "consensus"]
    exact = [r for r in ref if r["classification_cds"] == "Exact_Match"]
    single = sum(1 for r in exact if ref_exons.get(r["gene_id"], 0) == 1)
    merges = sum(1 for v in Counter(r["matched_id"] for r in ref
                                    if r["matched_id"] not in ("", "NA")).values() if v > 1)
    splits = sum(1 for v in Counter(r["matched_id"] for r in cons
                                    if r["matched_id"] not in ("", "NA")).values() if v > 1)
    f = lambda k: summ.get(k)
    return {
        "reference_genes": f("total_reference_genes"), "gmb_genes": f("total_consensus_genes"),
        "cds_exact": len(exact), "cds_exact_multi_exon": len(exact) - single,
        "cds_exact_single_exon": single, "cds_intron_chain": f("sens_cds_intron_chain_recovered"),
        "exact_match": f("ref_Exact_Match"), "partial_match": f("ref_Partial_Match"),
        "structural_mismatch": f("ref_Structural_Mismatch"), "strand_mismatch": f("ref_Strand_Mismatch"),
        "locus_detection": f("sens_locus_detection_rate"), "missed": f("sens_missed_count"),
        "novel": f("spec_novel_consensus_count"), "merges": merges, "splits": splits,
        "_ref_class": {r["gene_id"]: r["classification_cds"] for r in ref},
    }


def _validate_reference(path: str) -> None:
    """The evaluation reference must be GFF3. Checked before anything is written."""
    if not path or not os.path.exists(path):
        die(f"evaluation reference not found: {path!r}")
    fmt = sniff_format(path)
    if fmt == "gtf":
        die("Evaluation reference must currently be GFF3. GTF reference parsing is not "
            "implemented safely; convert to GFF3 before evaluation.")
    if fmt != "gff3":
        die(f"cannot identify {path} as GFF3 (no feature line with ID=/Parent= attributes).")


def cmd_evaluate(a):
    # Every check happens before evaluation/ is created: a bad reference or a broken
    # freeze must leave no misleading output behind.
    _validate_reference(a.reference)
    rec = _check_frozen(a.freeze)
    ref_exons = _ref_cds_exon_counts(a.reference)
    if not ref_exons:
        die(f"no gene -> mRNA/transcript -> CDS structure found in {a.reference}; it is not "
            f"a usable protein-coding GFF3 reference.")
    cand = _comparison_metrics(a.comparison, ref_exons)
    base = _comparison_metrics(a.baseline_comparison, ref_exons) if a.baseline_comparison else None
    out_dir = os.path.join(a.out, "evaluation")
    os.makedirs(out_dir, exist_ok=True)
    keys = [k for k in cand if not k.startswith("_")]
    with open(os.path.join(out_dir, "evaluation.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "standard_preset", "candidate_config", "delta"])
        for k in keys:
            bv = base[k] if base else ""
            try:
                d = float(cand[k]) - float(bv) if base else ""
                d = f"{d:+g}" if d != "" else ""
            except (TypeError, ValueError):
                d = ""
            w.writerow([k, bv, cand[k], d])
    result = {"source_overlay": rec["source_overlay"],
              "source_overlay_sha256": rec["source_overlay_sha256"],
              "resolved_config": rec["resolved_config"],
              "resolved_config_sha256": rec["resolved_config_sha256"],
              "candidate": {k: cand[k] for k in keys}}
    if base:
        rank = {"Exact_Match": 4, "Partial_Match": 3, "Structural_Mismatch": 2,
                "Strand_Mismatch": 1, "Missed": 0}
        imp = reg = 0
        for g, c in cand["_ref_class"].items():
            b = base["_ref_class"].get(g)
            if b is None:
                continue
            imp += rank.get(c, -1) > rank.get(b, -1)
            reg += rank.get(c, -1) < rank.get(b, -1)
        result["standard"] = {k: base[k] for k in keys}
        result["per_gene_vs_standard"] = {"improved": imp, "regressed": reg,
                                          "ratio": round(imp / reg, 3) if reg else None}
    json.dump(result, open(os.path.join(out_dir, "evaluation.json"), "w"), indent=1, default=str)
    print(open(os.path.join(out_dir, "evaluation.tsv")).read())
    if base:
        print(f"per reference gene vs standard: {result['per_gene_vs_standard']}")


# --------------------------------------------------------------------------
# summarise-run (reference-free hard QC + diagnostics)
# --------------------------------------------------------------------------

def cmd_summarise_run(a):
    fin, build = a.finalise, a.build
    fq = json.load(open(os.path.join(fin, "fasta_qc_report.json")))
    uq = json.load(open(os.path.join(fin, "utr_qc_report.json")))
    sc = fq.get("sequence_checks", {})
    genes, mrna, cds_n, exons = {}, {}, Counter(), defaultdict(list)
    for line in open(os.path.join(fin, "consensus.gff3")):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "gene":
            genes[at["ID"]] = (f[0], int(f[3]), int(f[4]), f[6])
        elif f[2] in ("mRNA", "transcript"):
            mrna[at["ID"]] = (at.get("Parent", ""), f[0], int(f[3]), int(f[4]), f[6])
        elif f[2] == "CDS":
            for t in at.get("Parent", "").split(","):
                cds_n[t] += 1
        elif f[2] == "exon":
            for t in at.get("Parent", "").split(","):
                exons[t].append((int(f[3]), int(f[4])))
    kids = defaultdict(list)
    for t, (g, c, s, e, st) in mrna.items():
        kids[g].append((c, s, e, st))
    bound_bad = cross = opp = 0
    for g, (c, s, e, st) in genes.items():
        k = kids.get(g, [])
        if not k or s != min(x[1] for x in k) or e != max(x[2] for x in k):
            bound_bad += 1
        cross += sum(x[0] != c for x in k)
        opp += sum(x[3] != st for x in k)
    canon = Counter(r["gene_id"] for r in csv.DictReader(
        open(os.path.join(fin, "canonical", "canonical_transcripts.tsv")), delimiter="\t"))
    one_canon = all(canon.get(g) == 1 for g in genes) and len(canon) == len(genes)

    genome = {}
    name, buf = None, []
    for line in opener(a.genome):
        if line.startswith(">"):
            if name:
                genome[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        genome[name] = "".join(buf)
    canonical = total = 0
    introns = []
    for t, ex in exons.items():
        ex = sorted(ex)
        if len(ex) < 2 or t not in mrna:
            continue
        chrom, strand = mrna[t][1], mrna[t][4]
        seq = genome.get(chrom, "")
        for i in range(len(ex) - 1):
            introns.append(ex[i + 1][0] - ex[i][1] - 1)
            d = seq[ex[i][1]:ex[i][1] + 2].upper()
            ac = seq[ex[i + 1][0] - 3:ex[i + 1][0] - 1].upper()
            site = (d + ac) if strand == "+" else (ac.translate(_COMP)[::-1] + d.translate(_COMP)[::-1])
            total += 1
            canonical += site in ("GTAG", "GCAG", "ATAC")
    glen = [e - s + 1 for (_, s, e, _) in genes.values()]
    ncds = [cds_n[t] for t in mrna if cds_n[t]]
    hard = {
        "cdna_mismatches": sc.get("cdna_mismatches_total"),
        "cds_mismatches": sc.get("cds_mismatches_total"),
        "protein_mismatches": sc.get("protein_mismatches_total"),
        "internal_stops": sc.get("internal_stops_total"),
        "gene_boundary_violations": bound_bad,
        "utr_invariant_violations": uq.get("violations"),
        "one_canonical_per_gene": "PASS" if one_canon else "FAIL",
        "cross_seqid_parent_child": cross, "opposite_strand_parent_child": opp,
    }
    hard_pass = (all(hard[k] == 0 for k in hard if k != "one_canonical_per_gene")
                 and hard["one_canonical_per_gene"] == "PASS" and fq.get("pass") is True)
    diag = {
        "genes": len(genes), "transcripts": len(mrna),
        "single_cds_exon_pct": round(100 * sum(n == 1 for n in ncds) / len(ncds), 2) if ncds else None,
        "canonical_intron_fraction": round(canonical / total, 4) if total else None,
        "introns": total, "intron_bp_p99": pct(introns, 99), "intron_bp_max": max(introns) if introns else None,
        "gene_bp_p99": pct(glen, 99), "gene_bp_max": max(glen) if glen else None,
        "genes_over_100kb": sum(x > 100000 for x in glen),
    }
    out = {"label": a.label, "hard_qc": hard, "hard_qc_pass": hard_pass, "diagnostics": diag}
    rm = os.path.join(build, "run_manifest.json")
    if os.path.exists(rm):
        r = json.load(open(rm))
        out["resolved_policy"] = r.get("resolved_policy")
        out["runtime_seconds"] = r.get("runtime_seconds")
    os.makedirs(a.out, exist_ok=True)
    json.dump(out, open(os.path.join(a.out, f"run_summary_{a.label}.json"), "w"), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))


def main(argv=None):
    p = argparse.ArgumentParser(prog="gmb_clade.py", description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("args", help="manifest -> GMB evidence arguments + slot map")
    s.add_argument("manifest"); s.add_argument("--out")
    s.add_argument("--shell", action="store_true", help="print only the argument string")

    s = sub.add_parser("measure", help="reference-free evidence measurements")
    s.add_argument("manifest")
    s.add_argument("--preflight", required=True, help="gmb-preflight dir, --preset standard")
    s.add_argument("--probe", required=True,
                   help="gmb-preflight dir, --preset standard + probe_rescue_auto.yaml")
    s.add_argument("--baseline-build", help="gmb-build dir made with --preset standard")
    s.add_argument("--out", required=True)

    s = sub.add_parser("freeze", help="freeze the overlay AND the build's resolved config "
                                      "before any reference is used")
    s.add_argument("config", help="the clade overlay, configs/<clade>.yaml")
    s.add_argument("--resolved-config", required=True,
                   help="<OUT>/candidate/build/resolved_config.yaml from the candidate build")
    s.add_argument("--out", required=True)

    s = sub.add_parser("summarise-run", help="reference-free hard QC + diagnostics")
    s.add_argument("--build", required=True); s.add_argument("--finalise", required=True)
    s.add_argument("--genome", required=True); s.add_argument("--label", required=True)
    s.add_argument("--out", required=True)

    s = sub.add_parser("evaluate", help="evaluate a FROZEN config against a reference")
    s.add_argument("--freeze", required=True, help="config_freeze.json from `freeze`")
    s.add_argument("--comparison", required=True, help="gmb-compare dir for the candidate")
    s.add_argument("--baseline-comparison", help="gmb-compare dir for the standard preset")
    s.add_argument("--reference", required=True,
                   help="reference annotation, GFF3 only (evaluation only)")
    s.add_argument("--out", required=True)

    a = p.parse_args(argv)
    {"args": cmd_args, "measure": cmd_measure, "freeze": cmd_freeze,
     "summarise-run": cmd_summarise_run, "evaluate": cmd_evaluate}[a.cmd](a)


if __name__ == "__main__":
    main()
