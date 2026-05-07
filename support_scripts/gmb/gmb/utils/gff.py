"""GFF3 parse/write helpers shared across the pipeline."""

from __future__ import annotations

import gzip
from collections import defaultdict
from typing import Any


def parse_gff3_hierarchy(
    gff3_path: str,
) -> dict[str, Any]:
    """Parse consensus GFF3 into structured transcript data.

    Returns
    -------
    dict with keys:
        genes: dict of gene_id -> {chrom, strand, start, end}
        transcripts: dict of tid -> {gene_id, chrom, strand, exons, cds, cds_phases}
        n_genes, n_transcripts, n_cds_transcripts, cds_transcript_ids
    """
    opener = gzip.open if gff3_path.endswith(".gz") else open

    genes: dict[str, dict] = {}
    transcripts: dict[str, dict] = {}
    exons_by_parent: dict[str, list] = defaultdict(list)
    cds_by_parent: dict[str, list] = defaultdict(list)

    with opener(gff3_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, _src, feat, start, end, _score, strand, frame, attrs_str = parts
            start = int(start)
            end = int(end)

            attrs: dict[str, str] = {}
            for kv in attrs_str.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k] = v

            if feat == "gene":
                genes[attrs["ID"]] = {
                    "chrom": chrom,
                    "strand": strand,
                    "start": start,
                    "end": end,
                }
            elif feat == "mRNA":
                tid = attrs["ID"]
                parent = attrs.get("Parent", "")
                transcripts[tid] = {
                    "gene_id": parent,
                    "chrom": chrom,
                    "strand": strand,
                    "exons": [],
                    "cds": [],
                    "evidence": attrs.get("Evidence", ""),
                }
            elif feat == "exon":
                parent = attrs.get("Parent", "")
                exons_by_parent[parent].append(
                    (start - 1, end)
                )
            elif feat == "CDS":
                parent = attrs.get("Parent", "")
                frame_val = frame if frame != "." else "0"
                cds_by_parent[parent].append(
                    (start - 1, end, int(frame_val))
                )

    for tid, tx in transcripts.items():
        tx["exons"] = sorted(exons_by_parent.get(tid, []))
        raw_cds = sorted(cds_by_parent.get(tid, []), key=lambda x: x[0])
        tx["cds"] = [(s, e) for s, e, _p in raw_cds]
        tx["cds_phases"] = [p for _s, _e, p in raw_cds]

    cds_tids = {tid for tid, tx in transcripts.items() if tx["cds"]}

    return {
        "genes": genes,
        "transcripts": transcripts,
        "n_genes": len(genes),
        "n_transcripts": len(transcripts),
        "n_cds_transcripts": len(cds_tids),
        "cds_transcript_ids": cds_tids,
    }
