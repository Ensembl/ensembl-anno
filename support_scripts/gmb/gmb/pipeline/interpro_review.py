#!/usr/bin/env python3
"""Select ambiguous canonical choices for InterProScan review, and prepare
a single batched InterProScan input for them.

Scope
-----
This is the *preparation* half of an optional second-stage resolver. It
consumes what GMB has already computed -- the canonical-selection report,
the protein-validation sidecar, and prot.fa -- and never invokes DIAMOND,
psauron, or InterProScan itself. The resolver half
(gmb.pipeline.interpro_resolver) consumes completed InterProScan output.

Explicitly NOT a geneset QC tool: it never scans a whole proteome, never
computes reference-vs-GMB domain coverage, and never reports Pfam coverage
statistics. Only genes whose canonical choice is genuinely ambiguous are
submitted, and only the top few candidates of each.

Why checksum-keyed
------------------
Candidates are deduplicated on the protein SHA-256 recorded in
protein_validation.tsv, so two transcripts translating to the same protein
are scanned once. The resolver maps the single result back to every
transcript that shares the checksum, and never treats the repeated result
as independent supporting evidence.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import warnings
from typing import TYPE_CHECKING, Optional

import pandas as pd

from gmb.pipeline.canonical_evidence import (
    balanced_coverage,
    evidence_classes_for_transcript,
    named_source_evidence_classes,
)
from gmb.pipeline.config import load_config

if TYPE_CHECKING:
    from gmb.pipeline.config import InterProResolverConfig

# ---------------------------------------------------------------------------
# Independent evidence classes
# ---------------------------------------------------------------------------
#
# GMB records *named software sources* (Scallop, StringTie, Tiberius,
# Minimap2, OrthoDB, GenBlast, UniProt). Those are not the same thing as
# independent biological evidence: Scallop and StringTie are two assemblers
# customarily run over the SAME short-read RNA-seq libraries, so counting
# them as two independent sources overstates the support for a transcript.
#
# The class vocabulary and mapping now live in gmb.pipeline.canonical_evidence
# -- the single shared definition used by both this module and
# gmb.pipeline.canonical_selection (they used to each keep their own copy,
# which had already drifted: this module's version excluded
# protein_validation as a class entirely). balanced_coverage() moved there
# for the same reason -- canonical_selection's protein-plausibility scoring
# uses the identical function.


def evidence_classes(evidence_sources: Optional[str], backbone_label: str = "Helixer") -> set:
    """Named-source evidence classes for one transcript (manifest/report use).

    Thin wrapper around the shared
    ``canonical_evidence.evidence_classes_for_transcript`` that does not have
    a protein-validation-support flag to hand (callers here work from
    ranking-table rows, which may or may not carry DIAMOND/psauron columns
    depending on which sidecar was joined) -- so this entry point never adds
    the ``protein_validation`` class. Use
    ``evidence_classes_for_transcript`` directly when that flag is
    available, as ``build_review_set`` below does.
    """
    classes, unknown_names = named_source_evidence_classes(evidence_sources, backbone_label)
    if unknown_names:
        warnings.warn(
            f"Evidence source(s) not recognised by gmb.pipeline.canonical_evidence, "
            f"mapped to the 'other' evidence class: {sorted(unknown_names)}.",
            stacklevel=2,
        )
    return classes


# ---------------------------------------------------------------------------
# Ambiguity detection
# ---------------------------------------------------------------------------


def _fmt(value) -> str:
    return "" if value is None or (isinstance(value, float) and pd.isna(value)) else str(value)


def detect_ambiguity(
    gene_rows: list,
    canonical_row: dict,
    cfg: InterProResolverConfig,
) -> list:
    """Return the list of ambiguity trigger names firing for one gene.

    ``gene_rows`` is the gene's transcripts from transcript_ranking.tsv,
    already rank-ordered; ``canonical_row`` is its canonical_transcripts.tsv
    row. An empty list means the choice is unambiguous and the gene is not
    submitted for review.
    """
    triggers = []
    if len(gene_rows) < cfg.min_isoforms:
        return triggers

    winner = gene_rows[0]
    runner = gene_rows[1]

    # Identical proteins are not an ambiguous *biological* choice -- the
    # two transcripts encode the same product, so a domain scan cannot
    # distinguish them. (Exact structural duplicates are already collapsed
    # upstream; this catches same-protein/different-UTR survivors.)
    if winner.get("protein_sha256") and winner.get("protein_sha256") == runner.get(
        "protein_sha256"
    ):
        return triggers

    confidence = str(canonical_row.get("confidence") or "")
    if cfg.trigger_low_confidence and confidence == "low_confidence":
        triggers.append("low_confidence")

    gap = canonical_row.get("score_gap")
    if (
        cfg.trigger_close_score
        and gap is not None
        and not pd.isna(gap)
        and float(gap) <= cfg.close_call_score_gap
    ):
        triggers.append("close_score_gap")

    # DIAMOND and psauron preferring different candidates is a genuine
    # conflict of fast protein evidence -- exactly what a domain scan can
    # arbitrate.
    if cfg.trigger_diamond_psauron_disagree:
        w_ps, r_ps = winner.get("psauron_score"), runner.get("psauron_score")
        w_cov = balanced_coverage(winner.get("diamond_qcov"), winner.get("diamond_scov"))
        r_cov = balanced_coverage(runner.get("diamond_qcov"), runner.get("diamond_scov"))
        if None not in (w_ps, r_ps, w_cov, r_cov):
            # A tie is indifference, not disagreement: only flag a conflict
            # when BOTH tools express a strict preference and those
            # preferences point at different candidates. Treating a tie as
            # disagreement would send every evenly-matched pair for review.
            psauron_prefers = 0 if w_ps == r_ps else (1 if w_ps > r_ps else -1)
            diamond_prefers = 0 if w_cov == r_cov else (1 if w_cov > r_cov else -1)
            if psauron_prefers and diamond_prefers and psauron_prefers != diamond_prefers:
                triggers.append("diamond_psauron_disagree")

    # A choice made only by CDS length or transcript ID means no biological
    # evidence separated the candidates at all.
    reason = str(canonical_row.get("selection_reason") or "")
    if cfg.trigger_length_or_lexical_only and reason in (
        "LONGEST_VALID_CDS_TIEBREAK",
        "LEXICAL_TIEBREAK",
    ):
        triggers.append("length_or_lexical_only")

    # One candidate substantially longer than the other, without stronger
    # protein support, is a possible fusion / run-on model.
    w_len, r_len = winner.get("protein_length"), runner.get("protein_length")
    if w_len and r_len and r_len > 0:
        ratio = float(w_len) / float(r_len)
        if ratio >= cfg.length_ratio_flag or ratio <= (1.0 / cfg.length_ratio_flag):
            triggers.append("length_disparity")

    # Complete vs internal-stop disagreement between the top two.
    if cfg.trigger_internal_stop_conflict:
        w_stop = bool(winner.get("has_internal_stop"))
        r_stop = bool(runner.get("has_internal_stop"))
        if w_stop != r_stop:
            triggers.append("internal_stop_conflict")

    return triggers


# ---------------------------------------------------------------------------
# Manifest / FASTA preparation
# ---------------------------------------------------------------------------


def _read_fasta(path: str) -> dict:
    seqs, tid, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if tid is not None:
                    seqs[tid] = "".join(buf)
                tid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if tid is not None:
        seqs[tid] = "".join(buf)
    return seqs


def build_review_set(
    canonical_tsv: str,
    ranking_tsv: str,
    protein_validation_tsv: str,
    prot_fa: str,
    cfg: InterProResolverConfig,
    backbone_label: str = "Helixer",
) -> tuple:
    """Select ambiguous genes and build (manifest_rows, checksum_to_seq, stats).

    Returns manifest rows in stable (gene_id, rank) order so repeated runs
    over unchanged input produce byte-identical output.
    """
    canon = pd.read_csv(canonical_tsv, sep="\t")
    ranking = pd.read_csv(ranking_tsv, sep="\t")
    proteins = _read_fasta(prot_fa)

    pv_by_tid = {}
    if protein_validation_tsv and os.path.exists(protein_validation_tsv):
        pv = pd.read_csv(protein_validation_tsv, sep="\t")
        pv_by_tid = pv.set_index("transcript_id").to_dict(orient="index")

    canon_by_gene = canon.set_index("gene_id").to_dict(orient="index")

    rows_by_gene = {}
    for rec in ranking.sort_values(["gene_id", "rank"]).to_dict(orient="records"):
        pv = pv_by_tid.get(rec["transcript_id"], {})
        rec["protein_sha256"] = pv.get("protein_sha256")
        rec["protein_length"] = pv.get("protein_length", rec.get("protein_length"))
        rows_by_gene.setdefault(rec["gene_id"], []).append(rec)

    manifest_rows = []
    checksum_to_seq = {}
    all_unknown_sources: set = set()
    stats = {
        "genes_total": len(rows_by_gene),
        "genes_multi_isoform": 0,
        "genes_selected_for_review": 0,
        "candidates_submitted": 0,
        "unique_protein_sequences": 0,
        "duplicate_candidates_avoided": 0,
        "trigger_counts": {},
    }

    for gene_id in sorted(rows_by_gene):
        gene_rows = rows_by_gene[gene_id]
        if len(gene_rows) > 1:
            stats["genes_multi_isoform"] += 1
        canonical_row = canon_by_gene.get(gene_id, {})
        triggers = detect_ambiguity(gene_rows, canonical_row, cfg)
        if not triggers:
            continue
        stats["genes_selected_for_review"] += 1
        for t in triggers:
            stats["trigger_counts"][t] = stats["trigger_counts"].get(t, 0) + 1

        for rec in gene_rows[: cfg.max_candidates_per_gene]:
            tid = rec["transcript_id"]
            seq = proteins.get(tid)
            if not seq:
                continue
            checksum = rec.get("protein_sha256")
            if not checksum or (isinstance(checksum, float) and pd.isna(checksum)):
                continue
            if checksum in checksum_to_seq:
                stats["duplicate_candidates_avoided"] += 1
            else:
                checksum_to_seq[checksum] = seq
            stats["candidates_submitted"] += 1
            has_protein_validation_support = bool(rec.get("diamond_hit")) or (
                rec.get("psauron_score") is not None
                and not (
                    isinstance(rec.get("psauron_score"), float) and pd.isna(rec["psauron_score"])
                )
            )
            classes, unknown = evidence_classes_for_transcript(
                rec.get("evidence_sources"), backbone_label, has_protein_validation_support
            )
            all_unknown_sources.update(unknown)
            manifest_rows.append(
                {
                    "gene_id": gene_id,
                    "transcript_id": tid,
                    "protein_sha256": checksum,
                    "protein_length": rec.get("protein_length"),
                    "current_rank": rec.get("rank"),
                    "is_current_canonical": bool(rec.get("is_canonical")),
                    "canonical_confidence": canonical_row.get("confidence"),
                    "review_reason": ";".join(triggers),
                    "runner_up_transcript_id": canonical_row.get("runner_up_transcript_id"),
                    "score_gap": canonical_row.get("score_gap"),
                    "diamond_hit": rec.get("diamond_hit"),
                    "diamond_pident": rec.get("diamond_pident"),
                    "diamond_qcov": rec.get("diamond_qcov"),
                    "diamond_scov": rec.get("diamond_scov"),
                    "diamond_balanced_coverage": balanced_coverage(
                        rec.get("diamond_qcov"), rec.get("diamond_scov")
                    ),
                    "psauron_score": rec.get("psauron_score"),
                    "named_evidence_sources": rec.get("evidence_sources"),
                    "evidence_classes": ",".join(sorted(classes)),
                    "n_named_sources": rec.get("n_independent_sources"),
                    "n_evidence_classes": len(classes),
                    "has_complete_orf": rec.get("has_complete_orf"),
                    "has_internal_stop": rec.get("has_internal_stop"),
                    "current_selection_reason": canonical_row.get("selection_reason"),
                }
            )

    stats["unique_protein_sequences"] = len(checksum_to_seq)
    if all_unknown_sources:
        warnings.warn(
            f"Evidence source(s) not recognised by gmb.pipeline.canonical_evidence, "
            f"mapped to the 'other' evidence class for {stats['candidates_submitted']} "
            f"review candidate(s): {sorted(all_unknown_sources)}. Never dropped, but "
            "consider adding them to canonical_evidence._SOURCE_TO_CLASS if they "
            "represent a genuinely new evidence type.",
            stacklevel=2,
        )
        stats["unknown_evidence_sources"] = sorted(all_unknown_sources)
    return manifest_rows, checksum_to_seq, stats


def write_review_inputs(
    manifest_rows: list,
    checksum_to_seq: dict,
    stats: dict,
    output_dir: str,
    provenance: Optional[dict] = None,
) -> dict:
    """Write the batched FASTA, manifest and summary. Returns written paths.

    One FASTA for the whole review set -- never one InterProScan run per
    gene. FASTA headers are the protein checksum (not a transcript ID), so
    the mapping back to transcripts is many-to-one through the manifest.
    """
    os.makedirs(output_dir, exist_ok=True)

    fasta_path = os.path.join(output_dir, "interpro_review_candidates.faa")
    with open(fasta_path, "w") as fh:
        for checksum in sorted(checksum_to_seq):
            fh.write(f">{checksum}\n{checksum_to_seq[checksum]}\n")

    manifest_path = os.path.join(output_dir, "interpro_review_manifest.tsv")
    pd.DataFrame(manifest_rows).to_csv(manifest_path, sep="\t", index=False)

    summary_path = os.path.join(output_dir, "interpro_review_summary.json")
    summary = dict(stats)
    if provenance:
        summary["provenance"] = provenance
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)

    return {"fasta": fasta_path, "manifest": manifest_path, "summary": summary_path}


# ---------------------------------------------------------------------------
# Mode A: GMB launches InterProScan 6 itself, via Nextflow
# ---------------------------------------------------------------------------
#
# Mode B (the default -- run_interproscan: false) is everything above this
# point: GMB only ever prepares a batched FASTA/manifest and later consumes
# whatever InterProScan output already exists, produced any way the
# operator likes. Mode A is this one function: GMB constructs and launches
# a Nextflow command from InterProNextflowExecConfig, fully parameterised,
# with zero hardcoded Docker/Slurm/Singularity/Apptainer image name,
# executor, or site-specific path anywhere in this module -- every such
# value comes from the config the caller supplies (see
# docs/interpro_resolver.md and docs/interproscan6_cluster.example.config
# for how the same function call maps onto a laptop Docker run vs a
# cluster Slurm+Singularity/Apptainer run).


def build_interproscan_nextflow_command(
    nf_cfg,
    input_fasta: str,
    output_dir: str,
) -> list:
    """Build the Nextflow command line for one InterProScan 6 run.

    Pure function (no subprocess call) so the exact command line is
    independently testable without Nextflow installed. ``nf_cfg`` is an
    ``InterProNextflowExecConfig``.
    """
    outdir = nf_cfg.output_dir or os.path.join(output_dir, "interpro_run", "results")
    workdir = nf_cfg.work_dir or os.path.join(output_dir, "interpro_run", "work")
    cmd = [
        nf_cfg.nextflow_executable,
        "run",
        nf_cfg.workflow,
        "-r",
        nf_cfg.revision,
        "-profile",
        nf_cfg.profile,
        "--input",
        input_fasta,
        "--interpro",
        nf_cfg.interpro_release,
        "--outdir",
        outdir,
        "-work-dir",
        workdir,
    ]
    if nf_cfg.config_file:
        cmd += ["-c", nf_cfg.config_file]
    if nf_cfg.data_dir:
        cmd += ["--datadir", nf_cfg.data_dir]
    if nf_cfg.applications:
        cmd += ["--appl", nf_cfg.applications]
    if nf_cfg.cpus:
        cmd += ["--cpus", str(nf_cfg.cpus)]
    if nf_cfg.max_workers:
        cmd += ["--max-workers", str(nf_cfg.max_workers)]
    cmd += list(nf_cfg.extra_args)
    return cmd


def run_interproscan_workflow(nf_cfg, input_fasta: str, output_dir: str) -> dict:
    """Launch InterProScan 6 via Nextflow (Mode A) and return its result paths.

    Raises ``RuntimeError`` on a non-zero Nextflow exit rather than
    returning a silently-empty result -- a failed scan must never be
    mistaken for "ran, found nothing". Requires ``nf_cfg.data_dir`` to be
    set (there is no sane default for a shared InterPro data directory).
    """
    if not nf_cfg.data_dir:
        raise ValueError(
            "run_interproscan_workflow: canonical_selection.interpro_resolver."
            "nextflow.data_dir must be set (the shared InterPro member-database "
            "directory) -- refusing to guess a path."
        )
    outdir = nf_cfg.output_dir or os.path.join(output_dir, "interpro_run", "results")
    cmd = build_interproscan_nextflow_command(nf_cfg, input_fasta, output_dir)
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"InterProScan Nextflow run failed (exit {result.returncode}).\n"
            f"Command: {' '.join(cmd)}\nstderr:\n{result.stderr}"
        )
    return {"command": cmd, "output_dir": outdir, "stdout": result.stdout}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canonical-transcripts", required=True)
    parser.add_argument("--transcript-ranking", required=True)
    parser.add_argument(
        "--protein-validation",
        default=None,
        help="protein_validation.tsv -- required for protein checksums",
    )
    parser.add_argument("--prot-fa", required=True, help="GMB prot.fa for the same run")
    parser.add_argument(
        "--config",
        action="append",
        default=None,
        help="YAML config override(s); may be repeated to layer several files",
    )
    parser.add_argument("--preset", default="fungi")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    config = load_config(args.config, args.preset)
    cfg = config.canonical_selection.interpro_resolver
    if not cfg.enabled:
        print("canonical_selection.interpro_resolver.enabled is false -- nothing to do.")
        return

    manifest_rows, checksum_to_seq, stats = build_review_set(
        args.canonical_transcripts,
        args.transcript_ranking,
        args.protein_validation,
        args.prot_fa,
        cfg,
        backbone_label=config.scoring.backbone_label,
    )
    paths = write_review_inputs(
        manifest_rows,
        checksum_to_seq,
        stats,
        args.output_dir,
        provenance={
            "max_candidates_per_gene": cfg.max_candidates_per_gene,
            "close_call_score_gap": cfg.close_call_score_gap,
            "min_isoforms": cfg.min_isoforms,
            "backbone_label": config.scoring.backbone_label,
        },
    )

    if cfg.run_interproscan and manifest_rows:
        # Mode A: GMB launches InterProScan 6 itself via Nextflow, fully
        # parameterised by canonical_selection.interpro_resolver.nextflow --
        # never the default (run_interproscan defaults False -> Mode B,
        # consume an already-completed run's output instead).
        print(f"run_interproscan is true -- launching InterProScan via Nextflow "
              f"({cfg.nextflow.workflow} -r {cfg.nextflow.revision}, "
              f"profile={cfg.nextflow.profile})")
        result = run_interproscan_workflow(cfg.nextflow, paths["fasta"], args.output_dir)
        print(f"  InterProScan output: {result['output_dir']}")

    print(
        f"Genes: {stats['genes_total']} "
        f"({stats['genes_multi_isoform']} multi-isoform); "
        f"selected for InterPro review: {stats['genes_selected_for_review']}"
    )
    print(
        f"Candidates: {stats['candidates_submitted']} "
        f"({stats['unique_protein_sequences']} unique sequences, "
        f"{stats['duplicate_candidates_avoided']} duplicate submissions avoided)"
    )
    if stats["trigger_counts"]:
        print("Triggers:")
        for name, count in sorted(stats["trigger_counts"].items()):
            print(f"  {name}: {count}")
    for label, path in paths.items():
        print(f"  {label}: {path}")


if __name__ == "__main__":
    main()
