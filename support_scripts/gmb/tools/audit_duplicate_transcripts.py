#!/usr/bin/env python3
"""Standalone diagnostic: classify duplicate/alternative isoforms in a
completed GMB run and write duplicate_transcript_audit.tsv.

Reads a run's own output files (consensus.gff3, evidence_attribution.tsv,
optional protein_validation.tsv, prot.fa) -- does not require re-running
gmb.cli.build. Uses the same classification engine
(gmb.pipeline.duplicate_transcript_collapse) that the real in-pipeline collapse
uses, so this audit is a preview of what collapse would do, not a
separately-maintained duplicate of that logic.

Usage:
    python tools/audit_duplicate_transcripts.py \
        --consensus-gff3 out/consensus.gff3 \
        --evidence-attribution out/evidence_attribution.tsv \
        --protein-validation out/protein_validation.tsv \
        --output-dir out/
"""

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gmb.pipeline.config import load_config
from gmb.pipeline.duplicate_transcript_collapse import build_duplicate_audit_rows
from gmb.utils.gff import parse_gff3_hierarchy


def _read_prot_fa(path):
    seqs = {}
    if not path or not os.path.exists(path):
        return seqs
    tid = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if tid is not None:
                    seqs[tid] = "".join(buf)
                tid = line[1:]
                buf = []
            else:
                buf.append(line)
        if tid is not None:
            seqs[tid] = "".join(buf)
    return seqs


def load_genes_for_audit(
    consensus_gff3, evidence_attribution_tsv, protein_validation_tsv, prot_fa
):
    hierarchy = parse_gff3_hierarchy(consensus_gff3)
    transcripts = hierarchy["transcripts"]

    ev_by_tid = (
        pd.read_csv(evidence_attribution_tsv, sep="\t")
        .set_index("transcript_id")
        .to_dict(orient="index")
    )
    pv_by_tid = {}
    if protein_validation_tsv and os.path.exists(protein_validation_tsv):
        pv_by_tid = (
            pd.read_csv(protein_validation_tsv, sep="\t")
            .set_index("transcript_id")
            .to_dict(orient="index")
        )
    proteins = _read_prot_fa(prot_fa)

    genes = {}
    for tid, tx in transcripts.items():
        ev = ev_by_tid.get(tid, {})
        pv = pv_by_tid.get(tid, {})
        record = {
            "transcript_id": tid,
            "gene_id": tx["gene_id"],
            "chrom": tx["chrom"],
            "strand": tx["strand"],
            "exons": tx["exons"],
            "cds": tx["cds"],
            "cds_phases": tx.get("cds_phases"),
            "protein": proteins.get(tid),
            "evidence_sources": ev.get("evidence_sources", tx.get("evidence", "")),
            "gmb_score": _nan_to_none(ev.get("gmb_score")),
            "protein_coding_score": _nan_to_none(pv.get("protein_coding_score")),
            "psauron_score": _nan_to_none(pv.get("psauron_score")),
            "diamond_hit": pv.get("diamond_hit"),
            "diamond_qcov": _nan_to_none(pv.get("diamond_qcov")),
            "diamond_scov": _nan_to_none(pv.get("diamond_scov")),
            "protein_length": _nan_to_none(pv.get("protein_length")),
            "orf_label": pv.get("orf_label"),
            "is_partial_5": _nan_to_bool(pv.get("is_partial_5")),
            "is_partial_3": _nan_to_bool(pv.get("is_partial_3")),
        }
        genes.setdefault(tx["gene_id"], []).append(record)
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


def _nan_to_bool(val):
    val = _nan_to_none(val)
    if val is None:
        return None
    if isinstance(val, bool):
        return val
    if isinstance(val, str):
        return val.strip().lower() in ("true", "1", "yes")
    return bool(val)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--consensus-gff3", required=True)
    parser.add_argument("--evidence-attribution", required=True)
    parser.add_argument("--protein-validation", default=None)
    parser.add_argument("--prot-fa", default=None, help="Default: prot.fa next to consensus.gff3")
    parser.add_argument("--config", default=None)
    parser.add_argument("--preset", default="fungi")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    prot_fa = args.prot_fa or os.path.join(os.path.dirname(args.consensus_gff3), "prot.fa")
    config = load_config(args.config, args.preset)

    genes = load_genes_for_audit(
        args.consensus_gff3, args.evidence_attribution, args.protein_validation, prot_fa
    )

    rows = build_duplicate_audit_rows(genes, config.duplicate_transcript_collapse)

    os.makedirs(args.output_dir, exist_ok=True)
    out_path = os.path.join(args.output_dir, "duplicate_transcript_audit.tsv")
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)

    n_genes_multi = sum(1 for t in genes.values() if len(t) > 1)
    print(f"Genes: {len(genes)} ({n_genes_multi} multi-isoform)")
    print(f"Audit rows: {len(rows)}")
    print(f"Wrote {out_path}")

    from collections import Counter

    counts = Counter(r["duplicate_classification"] for r in rows)
    print("\nClassification counts (transcript rows):")
    for cat, n in sorted(counts.items()):
        print(f"  {cat}: {n}")


if __name__ == "__main__":
    main()
