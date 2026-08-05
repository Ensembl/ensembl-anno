"""Post-build finalisation for Gene Model Builder.

Performs: FASTA regeneration from final GFF3, FASTA QC, UTR validation,
canonical selection, canonical GFF3 annotation, and handover manifest.

Usage::

    gmb-finalise \\
        --build-dir /path/to/gmb/output \\
        --genome /path/to/genome.fa \\
        --preset apicomplexa \\
        --output-dir /path/to/final

Or from Python::

    from gmb.pipeline.finalise import run_finalise
    result = run_finalise(build_dir, genome_path, output_dir, preset="apicomplexa")
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import time
from datetime import datetime, timezone


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def _require(path: str, label: str) -> None:
    if not os.path.exists(path):
        sys.exit(f"ERROR: required file missing: {label} ({path})")


def run_finalise(
    build_dir: str,
    genome_path: str,
    output_dir: str,
    preset: str = "standard",
    config_files: list[str] | None = None,
    comparison_gff3: str | None = None,
) -> dict:
    """Run the complete post-build finalisation workflow.

    Returns a summary dict suitable for JSON serialisation.
    """
    from gmb.pipeline.builder import regenerate_final_fasta
    from gmb.pipeline.canonical_selection import run_canonical_selection
    from gmb.pipeline.config import load_config
    from gmb.pipeline.fasta_qc import parse_gff3, validate_fasta
    from gmb.pipeline.utr_validator import summarise_violations, validate_utr_partition
    from gmb.utils.fasta import load_genome

    t_start = time.monotonic()
    os.makedirs(output_dir, exist_ok=True)

    consensus_gff3 = os.path.join(build_dir, "consensus.gff3")
    evidence_tsv = os.path.join(build_dir, "evidence_attribution.tsv")
    prot_val_tsv = os.path.join(build_dir, "protein_validation.tsv")
    resolved_cfg = os.path.join(build_dir, "resolved_config.yaml")

    _require(consensus_gff3, "consensus.gff3")
    _require(evidence_tsv, "evidence_attribution.tsv")
    _require(genome_path, "genome FASTA")

    errors: list[str] = []

    # --- 1. Copy immutable build GFF3 ---
    shutil.copy2(consensus_gff3, os.path.join(output_dir, "consensus.gff3"))

    # --- 2. Parse GFF3 into 0-based rows for FASTA regen ---
    print("Step 1: Regenerating final FASTA from GFF3...")
    gff_data = parse_gff3(consensus_gff3)
    genome = load_genome(genome_path)

    gff_rows_0based: list[dict] = []
    with open(consensus_gff3) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, src, feat, start, end, score, strand, frame, attrs_str = parts
            attrs = {}
            for kv in attrs_str.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k] = v
            gff_rows_0based.append(
                {
                    "Chromosome": chrom,
                    "Source": src,
                    "Feature": feat,
                    "Start": int(start) - 1,
                    "End": int(end),
                    "Score": score,
                    "Strand": strand,
                    "Frame": frame,
                    "ID": attrs.get("ID", ""),
                    "Parent": attrs.get("Parent", ""),
                }
            )

    fasta_stats = regenerate_final_fasta(gff_rows_0based, genome, output_dir)

    # --- 3. FASTA QC ---
    print("Step 2: Running FASTA QC...")
    fasta_qc = validate_fasta(output_dir, genome_path)
    fasta_qc_path = os.path.join(output_dir, "fasta_qc_report.json")
    with open(fasta_qc_path, "w") as fh:
        json.dump(fasta_qc, fh, indent=2)

    if not fasta_qc.get("pass", False):
        failed = fasta_qc.get("failed_checks", [])
        errors.append(f"FASTA QC failed: {', '.join(failed)}")
        print(f"  FAIL: {', '.join(failed)}")
    else:
        print("  PASS")

    # --- 4. UTR validation ---
    print("Step 3: Running UTR validation...")
    utr_violations = validate_utr_partition(gff_rows_0based)
    n_applicable = sum(1 for r in gff_rows_0based if r["Feature"] == "mRNA")
    utr_summary = summarise_violations(utr_violations, n_applicable)
    utr_qc_path = os.path.join(output_dir, "utr_qc_report.json")
    with open(utr_qc_path, "w") as fh:
        json.dump(utr_summary, fh, indent=2)

    if utr_violations:
        utr_tsv_path = os.path.join(output_dir, "utr_qc_violations.tsv")
        with open(utr_tsv_path, "w") as fh:
            fh.write("transcript_id\tseqname\tstrand\tutr5_gap_bp\tutr3_gap_bp\n")
            for v in utr_violations:
                fh.write(
                    f"{v.transcript_id}\t{v.seqname}\t{v.strand}\t"
                    f"{v.utr5_cds_gap_bp}\t{v.utr3_cds_gap_bp}\n"
                )
        errors.append(f"UTR validation: {len(utr_violations)} violations")
        print(f"  FAIL: {len(utr_violations)} violations")
    else:
        print(f"  PASS ({n_applicable} transcripts checked)")

    # --- 5. Canonical selection ---
    print("Step 4: Running canonical selection...")
    if os.path.exists(resolved_cfg) and not config_files:
        config = load_config(resolved_cfg, preset=None)
        print(f"  Using build's resolved_config.yaml")
    else:
        config = load_config(config_files, preset)
    pv_path = prot_val_tsv if os.path.exists(prot_val_tsv) else None
    canonical_dir = os.path.join(output_dir, "canonical")
    canonical_summary = run_canonical_selection(
        consensus_gff3=os.path.join(output_dir, "consensus.gff3"),
        evidence_attribution_tsv=evidence_tsv,
        output_dir=canonical_dir,
        config=config,
        protein_validation_tsv=pv_path,
        annotate_gff3=True,
    )

    # Verify one canonical per gene
    can_tsv = os.path.join(canonical_dir, "canonical_transcripts.tsv")
    annotated_gff3 = os.path.join(canonical_dir, "consensus.canonical_annotated.gff3")
    n_genes = canonical_summary.get("n_genes", 0)
    print(f"  {n_genes} genes, canonical GFF3 written")

    # --- 6. Canonical consistency checks ---
    if os.path.exists(can_tsv):
        import csv

        with open(can_tsv) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            can_rows = list(reader)
        can_tids = {r["canonical_transcript_id"] for r in can_rows}
        can_genes = {r["gene_id"] for r in can_rows}
        gff_genes = set(gff_data["genes"].keys())
        gff_tids = set(gff_data["transcripts"].keys())

        if can_genes != gff_genes:
            diff = can_genes ^ gff_genes
            errors.append(f"Canonical gene set mismatch: {len(diff)} differences")
        if not can_tids.issubset(gff_tids):
            missing = can_tids - gff_tids
            errors.append(f"Canonical transcripts not in GFF3: {len(missing)}")

    # --- 7. Handover manifest ---
    print("Step 5: Writing handover manifest...")
    elapsed = time.monotonic() - t_start
    manifest = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "build_dir": os.path.abspath(build_dir),
        "output_dir": os.path.abspath(output_dir),
        "genome_path": os.path.abspath(genome_path),
        "preset": preset,
        "config_files": config_files or [],
        "gene_count": gff_data["n_genes"],
        "transcript_count": gff_data["n_transcripts"],
        "coding_transcript_count": len(gff_data["cds_transcript_ids"]),
        "canonical_count": n_genes,
        "fasta_qc_pass": fasta_qc.get("pass", False),
        "fasta_qc_failed_checks": fasta_qc.get("failed_checks", []),
        "utr_validation_pass": len(utr_violations) == 0,
        "utr_violations": len(utr_violations),
        "fasta_stats": fasta_stats,
        "resolved_config_present": os.path.exists(resolved_cfg),
        "resolved_config_sha256": _sha256(resolved_cfg) if os.path.exists(resolved_cfg) else None,
        "runtime_seconds": round(elapsed, 1),
        "errors": errors,
        "outputs": {},
    }

    for fname in [
        "consensus.gff3",
        "cdna.fa",
        "cds.fa",
        "prot.fa",
        "fasta_qc_report.json",
        "utr_qc_report.json",
    ]:
        fpath = os.path.join(output_dir, fname)
        if os.path.exists(fpath):
            manifest["outputs"][fname] = {
                "size": os.path.getsize(fpath),
                "sha256": _sha256(fpath),
            }

    for fname in [
        "canonical_transcripts.tsv",
        "transcript_ranking.tsv",
        "canonical_selection_summary.json",
        "consensus.canonical_annotated.gff3",
    ]:
        fpath = os.path.join(canonical_dir, fname)
        if os.path.exists(fpath):
            manifest["outputs"][f"canonical/{fname}"] = {
                "size": os.path.getsize(fpath),
                "sha256": _sha256(fpath),
            }

    manifest_path = os.path.join(output_dir, "handover_manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    manifest_tsv_path = os.path.join(output_dir, "handover_manifest.tsv")
    with open(manifest_tsv_path, "w") as fh:
        fh.write("key\tvalue\n")
        for k, v in manifest.items():
            if isinstance(v, dict):
                for k2, v2 in v.items():
                    fh.write(f"{k}.{k2}\t{v2}\n")
            elif isinstance(v, list):
                fh.write(f"{k}\t{','.join(str(x) for x in v)}\n")
            else:
                fh.write(f"{k}\t{v}\n")

    # --- Final verdict ---
    if errors:
        print(f"\nFINALISATION FAILED: {len(errors)} error(s):")
        for e in errors:
            print(f"  - {e}")
        return manifest

    print("\nFINALISATION PASSED — ready for QC handover")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Post-build finalisation: FASTA regen, QC, canonical selection, manifest"
    )
    parser.add_argument(
        "--build-dir",
        required=True,
        help="Directory containing gmb-build output (consensus.gff3, evidence_attribution.tsv, ...)",
    )
    parser.add_argument("--genome", required=True, help="Genome FASTA")
    parser.add_argument("--output-dir", required=True, help="Finalisation output directory")
    parser.add_argument(
        "--preset",
        default="standard",
        help="Config preset (default: standard). Use 'apicomplexa' for Plasmodium.",
    )
    parser.add_argument(
        "--config",
        action="append",
        default=None,
        help="Additional YAML config file(s). May be repeated.",
    )
    parser.add_argument("--comparison-gff3", default=None, help="Reference GFF3 for comparison")
    args = parser.parse_args()

    manifest = run_finalise(
        build_dir=args.build_dir,
        genome_path=args.genome,
        output_dir=args.output_dir,
        preset=args.preset,
        config_files=args.config,
        comparison_gff3=args.comparison_gff3,
    )

    if manifest.get("errors"):
        sys.exit(1)


if __name__ == "__main__":
    main()
