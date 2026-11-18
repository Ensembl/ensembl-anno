#!/usr/bin/env python3
"""``gmb-protein-qc`` -- independent DIAMOND/PSAURON QC of finished GMB proteins.

Scores the proteins of a completed annotation with the same DIAMOND and PSAURON
code ``gmb-build`` uses for in-build protein validation (one run of each tool on
the unique sequences), without taking part in model selection::

    gmb-protein-qc --proteins final/prot.fa --gff3 final/consensus.gff3 \\
        --build-dir build/ --diamond-db uniprot_fungi.dmnd --psauron-cpu \\
        --output-dir protein_qc/

Writes ``protein_qc.tsv`` (one row per transcript) and ``protein_qc_summary.json``
(input checksums, tool versions, database, distributions). QC only: nothing is
filtered, no annotation is changed, and neither ``gmb-build`` nor ``gmb-finalise``
reads these files. The build may have had protein validation disabled; this tool
enables it on a private copy of the configuration.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import subprocess
import sys
import time
from collections import Counter

import gmb
from gmb.pipeline.config import DEFAULT_PRESET, load_config
from gmb.pipeline.protein_validation import batch_score_proteins, check_dependencies
from gmb.provenance.manifest import _git_state, _gmb_version, _sha256
from gmb.utils.fasta import load_genome as load_fasta

QC_TSV = "protein_qc.tsv"
QC_SUMMARY = "protein_qc_summary.json"

# Annotation context first, then the batch_score_proteins detail fields.
COLUMNS = [
    "gene_id",
    "transcript_id",
    "candidate_origin",
    "canonical",
    "evidence",
    "protein_sha256",
    "protein_length",
    "diamond_hit",
    "diamond_pident",
    "diamond_qcov",
    "diamond_scov",
    "diamond_bitscore",
    "diamond_evalue",
    "psauron_score",
    "psauron_is_protein",
    "protein_coding_score",
    "validation_status",
    "validation_reason",
    "protein_validation_source",
    "protein_validation_reused_from",
    "diamond_version",
    "psauron_version",
    "diamond_db",
]


def read_transcripts(gff3_path: str) -> dict[str, dict]:
    """mRNA ID -> gene, origin (existing|secondary_orf), canonical tag, evidence."""
    out: dict[str, dict] = {}
    any_canonical = False
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "mRNA":
                continue
            attrs = dict(kv.split("=", 1) for kv in parts[8].split(";") if "=" in kv)
            canonical = "Ensembl_canonical" in attrs.get("tag", "").split(",")
            any_canonical |= canonical
            out[attrs["ID"]] = {
                "gene_id": attrs.get("Parent", ""),
                "candidate_origin": attrs.get("CandidateOrigin", "existing"),
                "canonical": canonical,
                "evidence": attrs.get("Evidence", ""),
            }
    if not any_canonical:  # not a canonical-annotated GFF3: leave the column empty
        for rec in out.values():
            rec["canonical"] = None
    return out


def resolve_config(args):
    """The configuration to score under, with validation forced on in a copy."""
    if args.build_dir:
        path = os.path.join(args.build_dir, "resolved_config.yaml")
        if not os.path.exists(path):
            sys.exit(f"ERROR: {path} not found")
        cfg, source = load_config(path, preset=None, warn_inert=False), path
    else:
        cfg = load_config(args.config or None, args.preset)
        source = f"preset={args.preset} config={args.config or []}"
    overrides = {"enabled": True}
    for flag, field in (
        ("diamond_db", "diamond_db"),
        ("diamond_path", "diamond_path"),
        ("psauron_path", "psauron_path"),
        ("psauron_min_length", "psauron_min_length"),
    ):
        if getattr(args, flag) is not None:
            overrides[field] = getattr(args, flag)
    if args.psauron_cpu:
        overrides["psauron_use_cpu"] = True
    pv = dataclasses.replace(cfg.protein_validation, **overrides)
    return dataclasses.replace(cfg, protein_validation=pv), source


def _diamond_dbinfo(pv) -> dict:
    """Best-effort ``diamond dbinfo`` (sequence/letter counts, format version)."""
    db = pv.diamond_db
    path = db if os.path.exists(db) else f"{db}.dmnd"
    info = {"path": os.path.abspath(path), "size_bytes": None, "mtime": None}
    if os.path.exists(path):
        st = os.stat(path)
        info.update(
            size_bytes=st.st_size,
            mtime=time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(st.st_mtime)),
        )
    try:
        r = subprocess.run(
            [pv.diamond_path, "dbinfo", "-d", db], capture_output=True, text=True, timeout=120
        )
        for line in (r.stdout or "").splitlines():
            key, _, val = line.strip().partition("  ")  # right-aligned "Key  value"
            if key.strip() and val.strip():
                info[key.strip().lower().replace(" ", "_")] = val.strip()
    except Exception:
        pass
    return info


def _group_stats(rows: list[dict]) -> dict:
    scores = sorted(r["psauron_score"] for r in rows if r["psauron_score"] is not None)

    def q(f):
        return round(scores[min(len(scores) - 1, int(f * len(scores)))], 4) if scores else None

    calls = Counter(r["psauron_is_protein"] for r in rows)
    with_protein = [r for r in rows if r["protein_length"]]
    return {
        "transcripts": len(rows),
        "with_protein": len(with_protein),
        "diamond_hit_fraction": (
            round(sum(r["diamond_hit"] is not None for r in with_protein) / len(with_protein), 4)
            if with_protein
            else None
        ),
        "psauron_score": {
            "min": q(0),
            "q25": q(0.25),
            "median": q(0.5),
            "q75": q(0.75),
            "max": q(1),
        },
        "psauron_is_protein": {"true": calls[True], "false": calls[False], "missing": calls[None]},
    }


def run_protein_qc(
    proteins_fa: str,
    gff3: str,
    config,
    output_dir: str,
    config_source: str = "",
    overwrite: bool = False,
) -> dict:
    t0 = time.perf_counter()
    tsv_path = os.path.join(output_dir, QC_TSV)
    if os.path.exists(tsv_path) and not overwrite:
        sys.exit(f"ERROR: {tsv_path} exists; pass --overwrite or choose another --output-dir")
    transcripts = read_transcripts(gff3)
    proteins = {k: str(v) for k, v in load_fasta(proteins_fa).items()}
    unknown = sorted(set(proteins) - set(transcripts))
    if unknown:
        sys.exit(
            f"ERROR: {len(unknown)} protein ID(s) in {proteins_fa} are not mRNAs in "
            f"{gff3} (e.g. {unknown[:3]}); pass the GFF3 that produced this FASTA"
        )
    pv = config.protein_validation
    check_dependencies(pv)
    scores, details = batch_score_proteins(
        {tid: proteins.get(tid, "") for tid in transcripts}, config
    )

    rows = []
    for tid, rec in transcripts.items():
        row = {
            **rec,
            "transcript_id": tid,
            **details[tid],
            "protein_coding_score": scores.get(tid),
        }
        rows.append({c: row.get(c) for c in COLUMNS})
    os.makedirs(output_dir, exist_ok=True)
    with open(tsv_path, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join("" if r[c] is None else str(r[c]) for c in COLUMNS) + "\n")

    groups = {"all": rows}
    for origin in sorted({r["candidate_origin"] for r in rows}):
        groups[f"origin:{origin}"] = [r for r in rows if r["candidate_origin"] == origin]
    if any(r["canonical"] for r in rows):
        groups["canonical"] = [r for r in rows if r["canonical"]]
    first = next(iter(details.values()), {})
    summary = {
        "tool": "gmb-protein-qc",
        "gmb_version": _gmb_version(),
        "git": _git_state(os.path.dirname(os.path.dirname(os.path.abspath(gmb.__file__)))),
        "inputs": {
            "proteins": {"path": os.path.abspath(proteins_fa), "sha256": _sha256(proteins_fa)},
            "gff3": {"path": os.path.abspath(gff3), "sha256": _sha256(gff3)},
            "config_source": config_source,
        },
        "diamond": {
            "path": pv.diamond_path,
            "version": first.get("diamond_version"),
            "database": _diamond_dbinfo(pv),
            "min_query_coverage": pv.diamond_min_query_coverage,
            "min_target_coverage": pv.diamond_min_target_coverage,
        },
        "psauron": {
            "path": pv.psauron_path,
            "version": first.get("psauron_version"),
            "min_length": pv.psauron_min_length,
            "use_cpu": pv.psauron_use_cpu,
        },
        "protein_coding_score_weights": {
            "diamond": pv.diamond_weight,
            "psauron": pv.psauron_weight,
        },
        "genes": len({r["gene_id"] for r in rows}),
        "unique_proteins": len({r["protein_sha256"] for r in rows if r["protein_sha256"]}),
        "validation_status": dict(Counter(r["validation_status"] for r in rows)),
        "groups": {name: _group_stats(g) for name, g in groups.items()},
        "seconds": round(time.perf_counter() - t0, 1),
        "notes": [
            "QC only: not read by gmb-build or gmb-finalise; no model was selected, "
            "filtered or changed.",
            "psauron_is_protein is PSAURON's own call at its built-in threshold; "
            "no GMB threshold is applied.",
            "protein_coding_score uses the configured weights for comparison with "
            "in-build validation; it is not a pass/fail call.",
            "DIAMOND hits are independent evidence only if the database excludes the "
            "target species; check the database provenance.",
        ],
    }
    with open(os.path.join(output_dir, QC_SUMMARY), "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
    return summary


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="gmb-protein-qc", description=__doc__.split("\n\n")[0])
    p.add_argument("--proteins", required=True, help="final protein FASTA (prot.fa)")
    p.add_argument(
        "--gff3",
        required=True,
        help="GFF3 the proteins came from (consensus.gff3 or the canonical-"
        "annotated GFF3); gives gene IDs, secondary-ORF origin, canonical tag",
    )
    p.add_argument("--output-dir", required=True)
    p.add_argument(
        "--build-dir", help="read tool/database settings from BUILD_DIR/resolved_config.yaml"
    )
    p.add_argument(
        "--preset",
        default=DEFAULT_PRESET,
        help=f"preset when --build-dir is not given (default {DEFAULT_PRESET})",
    )
    p.add_argument("--config", action="append", help="YAML override(s); repeatable")
    p.add_argument("--diamond-db", help="DIAMOND database (overrides the config)")
    p.add_argument("--diamond-path", help="diamond executable")
    p.add_argument("--psauron-path", help="psauron executable")
    p.add_argument("--psauron-min-length", type=int, help="psauron -m (aa)")
    p.add_argument("--psauron-cpu", action="store_true", help="run psauron with -c (CPU)")
    p.add_argument("--overwrite", action="store_true", help="replace an existing protein_qc.tsv")
    args = p.parse_args(argv)
    config, source = resolve_config(args)
    s = run_protein_qc(
        args.proteins, args.gff3, config, args.output_dir, source, overwrite=args.overwrite
    )
    a = s["groups"]["all"]
    print(
        f"{a['transcripts']} transcripts, {s['unique_proteins']} unique proteins; "
        f"DIAMOND hit fraction {a['diamond_hit_fraction']}, "
        f"PSAURON median {a['psauron_score']['median']}; wrote {args.output_dir}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
