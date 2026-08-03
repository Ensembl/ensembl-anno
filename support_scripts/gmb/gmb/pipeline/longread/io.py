"""GTF I/O, split-stage handling, short-read transcript loading, and output writing.

Stage 1 — ``split_by_seqname``
    One sequential streaming pass over the raw GTF.  Writes each line to a
    per-seqname temp file (O(1) memory).  Optionally writes a split manifest
    (``split_manifest.tsv``) with record counts, file sizes, and SHA-256
    checksums so that reuse runs can validate the split files without
    re-reading the full raw input.

Stage reuse modes:
    ``--split-only``         run Stage 1 only, write manifest, exit.
    ``--reuse-split-dir``    skip Stage 1; validate and reuse an existing split.
    ``--input-split-gtf``    process exactly one pre-split seqname file.
"""

from __future__ import annotations

import hashlib
import os
import re
import sys
import time
from typing import Optional

from gmb.utils.io import ensure_dir

ATTR_RE = re.compile(r'(\w+) "([^"]*)"')

# ---------------------------------------------------------------------------
# GTF attribute parsing
# ---------------------------------------------------------------------------


def parse_attrs(attrs_str: str) -> dict[str, str]:
    return dict(ATTR_RE.findall(attrs_str))


# ---------------------------------------------------------------------------
# Load raw per-read transcripts from a per-seqname split file
# ---------------------------------------------------------------------------


def load_transcripts_from_split_file(split_path: str) -> tuple[dict[str, dict], list[dict]]:
    """Reconstruct transcripts from a split file's exon rows.

    Ignores the ``transcript`` feature line's own Start/End (empirically
    ~2.93% have swapped coordinates).  Every transcript's span is derived
    purely from its exon rows.

    Returns
    -------
    transcripts : dict[tid -> {strand, gene_id, exons}]
    malformed_records : list of dicts describing rows that could not be parsed
    """
    transcripts: dict[str, dict] = {}
    malformed: list[dict] = []

    with open(split_path) as fh:
        for lineno, line in enumerate(fh, 1):
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                if line.strip():
                    malformed.append(
                        {"lineno": lineno, "reason": "WRONG_FIELD_COUNT", "line": line[:200]}
                    )
                continue
            chrom, _source, feat, start_s, end_s, _score, strand, _frame, attrs = f
            if feat != "exon":
                continue
            try:
                start, end = int(start_s), int(end_s)
            except ValueError:
                malformed.append(
                    {"lineno": lineno, "reason": "INVALID_COORDINATES", "line": line[:200]}
                )
                continue
            if start > end:
                malformed.append(
                    {"lineno": lineno, "reason": "INVERTED_COORDINATES", "line": line[:200]}
                )
                continue
            if start == end:
                malformed.append(
                    {"lineno": lineno, "reason": "ZERO_LENGTH_EXON", "line": line[:200]}
                )
                continue

            attrs_d = parse_attrs(attrs)
            tid = attrs_d.get("transcript_id")
            if not tid:
                malformed.append(
                    {"lineno": lineno, "reason": "MISSING_TRANSCRIPT_ID", "line": line[:200]}
                )
                continue
            gid = attrs_d.get("gene_id", tid)

            if strand not in ("+", "-", "."):
                malformed.append(
                    {"lineno": lineno, "reason": "INVALID_STRAND", "line": line[:200]}
                )
                continue

            rec = transcripts.get(tid)
            if rec is None:
                rec = {"strand": strand, "gene_id": gid, "exons": []}
                transcripts[tid] = rec
            elif rec["strand"] != strand:
                # Mixed strand within one transcript — flag the transcript
                rec["mixed_strand"] = True
            rec["exons"].append((start, end))

    return transcripts, malformed


# ---------------------------------------------------------------------------
# Stage 1: split the raw file by seqname
# ---------------------------------------------------------------------------


def split_by_seqname(
    input_path: str,
    split_dir: str,
    only_seqname: Optional[str] = None,
    write_manifest: bool = True,
) -> dict[str, str]:
    """Streaming pass: write each GTF line to a per-seqname temp file.

    O(1) memory beyond small per-file write buffers.  The raw file is NOT
    assumed to be seqname-sorted (empirically it is not).

    Parameters
    ----------
    input_path : str
        Raw Minimap2 GTF path.
    split_dir : str
        Directory for ``<seqname>.gtf`` files.
    only_seqname : str or None
        When given, only this seqname's lines are written.
    write_manifest : bool
        When True, write ``split_manifest.tsv`` and compute SHA-256 per file.

    Returns
    -------
    dict mapping seqname → filepath for every seqname written.
    """
    ensure_dir(split_dir)
    handles: dict[str, object] = {}
    paths: dict[str, str] = {}
    line_counts: dict[str, int] = {}
    n_lines = 0
    n_skipped = 0
    t0 = time.time()

    with open(input_path) as fh:
        for line in fh:
            n_lines += 1
            if n_lines % 20_000_000 == 0:
                print(
                    f"    ... split: {n_lines:,} lines read ({time.time() - t0:.0f}s)",
                    file=sys.stderr,
                )
            if not line or line.startswith("#"):
                continue
            tab = line.find("\t")
            if tab < 0:
                continue
            chrom = line[:tab]

            if only_seqname is not None and chrom != only_seqname:
                n_skipped += 1
                continue

            handle = handles.get(chrom)
            if handle is None:
                path = os.path.join(split_dir, f"{chrom}.gtf")
                handle = open(path, "w")
                handles[chrom] = handle
                paths[chrom] = path
                line_counts[chrom] = 0
            handle.write(line)
            line_counts[chrom] += 1

    for handle in handles.values():
        handle.close()

    elapsed = time.time() - t0
    print(
        f"    Split complete: {n_lines:,} lines → {len(paths)} seqname file(s) ({elapsed:.0f}s)"
        + (f", {n_skipped:,} lines skipped (not {only_seqname})" if only_seqname else "")
    )

    if write_manifest and paths:
        _write_split_manifest(split_dir, paths, line_counts)

    return paths


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_split_manifest(
    split_dir: str, paths: dict[str, str], line_counts: dict[str, int]
) -> str:
    manifest_path = os.path.join(split_dir, "split_manifest.tsv")
    with open(manifest_path, "w") as fh:
        fh.write("seqname\tfile_path\tline_count\tfile_size_bytes\tsha256\tstatus\n")
        for seqname in sorted(paths):
            p = paths[seqname]
            size = os.path.getsize(p)
            sha = _sha256(p)
            fh.write(f"{seqname}\t{p}\t{line_counts[seqname]}\t{size}\t{sha}\tcomplete\n")
    return manifest_path


def load_split_manifest(split_dir: str) -> dict[str, dict]:
    """Load and return the split manifest as ``{seqname: row_dict}``."""
    manifest_path = os.path.join(split_dir, "split_manifest.tsv")
    if not os.path.exists(manifest_path):
        raise FileNotFoundError(f"No split manifest in {split_dir}: {manifest_path}")
    rows: dict[str, dict] = {}
    with open(manifest_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            row = dict(zip(header, parts))
            rows[row["seqname"]] = row
    return rows


def validate_split_dir(split_dir: str) -> dict[str, str]:
    """Validate split files against the manifest; return ``{seqname: path}``.

    Raises ``ValueError`` if any file is missing or has a checksum mismatch.
    """
    manifest = load_split_manifest(split_dir)
    paths: dict[str, str] = {}
    errors: list[str] = []
    for seqname, row in manifest.items():
        p = row["file_path"]
        if not os.path.exists(p):
            errors.append(f"  Missing split file: {p}")
            continue
        actual_sha = _sha256(p)
        if actual_sha != row["sha256"]:
            errors.append(
                f"  SHA-256 mismatch for {p}: expected {row['sha256']}, got {actual_sha}"
            )
            continue
        paths[seqname] = p
    if errors:
        raise ValueError("Split directory validation failed:\n" + "\n".join(errors))
    return paths


# ---------------------------------------------------------------------------
# Short-read transcript loading (multi-file, named sources)
# ---------------------------------------------------------------------------


def load_shortread_transcripts(
    scallop_paths: list[str],
    stringtie_paths: list[str],
    seqname: str,
) -> dict[str, list[dict]]:
    """Load full short-read transcript structures for one seqname.

    Returns ``{"scallop": [...], "stringtie": [...], "all": [...]}`` where
    each entry is::

        {
            "tid":    str,        # transcript_id with source-namespaced prefix to avoid collisions
            "raw_tid": str,       # original transcript_id from the file
            "strand": str,
            "exons":  list of (start, end) sorted ascending,
            "source": str,        # "scallop" or "stringtie"
            "source_file": str,   # path of the originating GTF
        }

    Transcript IDs are namespaced as ``{source}:{file_idx}:{raw_tid}`` to
    prevent collisions when the same transcript ID appears in multiple input
    files.
    """
    scallop = _load_source(scallop_paths, "scallop", seqname)
    stringtie = _load_source(stringtie_paths, "stringtie", seqname)
    return {"scallop": scallop, "stringtie": stringtie, "all": scallop + stringtie}


def _load_source(paths: list[str], source_name: str, seqname: str) -> list[dict]:
    all_txs: list[dict] = []
    for file_idx, path in enumerate(paths):
        if not path or not os.path.exists(path):
            continue
        by_tid: dict[str, dict] = {}
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) != 9 or f[0] != seqname or f[2] != "exon":
                    continue
                try:
                    start, end = int(f[3]), int(f[4])
                except ValueError:
                    continue
                attrs = parse_attrs(f[8])
                raw_tid = attrs.get("transcript_id")
                if not raw_tid:
                    continue
                # Namespace the tid to prevent cross-file collisions
                ns_tid = f"{source_name}:{file_idx}:{raw_tid}"
                if ns_tid not in by_tid:
                    by_tid[ns_tid] = {
                        "tid": ns_tid,
                        "raw_tid": raw_tid,
                        "strand": f[6],
                        "exons": [],
                        "source": source_name,
                        "source_file": path,
                    }
                by_tid[ns_tid]["exons"].append((start, end))
        for tx in by_tid.values():
            tx["exons"].sort()
        all_txs.extend(by_tid.values())
    return all_txs


# ---------------------------------------------------------------------------
# GTF output
# ---------------------------------------------------------------------------


def write_consensus_gtf(
    records: list[dict],
    output_path: str,
    source_label: str = "Minimap2Consensus",
) -> None:
    """Write consensus records as a GTF with explicit ``gene`` rows.

    Attributes written per transcript:
    - ``gene_id``, ``transcript_id``
    - ``read_support`` — number of reads in the consensus group
    - ``cluster_size`` — reads in the containing locus cluster
    - ``rescue_applied`` — True/False
    - ``rescue_reason`` — policy code or empty string

    Full SR transcript IDs that contributed rescue evidence are omitted from
    the GTF (they belong in the rescue attribution TSV sidecar).
    """
    from collections import defaultdict

    by_gene: dict[str, list[dict]] = defaultdict(list)
    for rec in records:
        by_gene[rec["gene_id"]].append(rec)

    with open(output_path, "w") as fh:
        for gid, gene_records in by_gene.items():
            gene_start = min(r["exons"][0][0] for r in gene_records)
            gene_end = max(r["exons"][-1][1] for r in gene_records)
            strand = gene_records[0]["strand"]
            seqname = gene_records[0]["seqname"]
            fh.write(
                f"{seqname}\t{source_label}\tgene\t{gene_start}\t{gene_end}\t.\t"
                f'{strand}\t.\tgene_id "{gid}";\n'
            )
            for rec in gene_records:
                tid = rec["transcript_id"]
                exons = rec["exons"]
                start = exons[0][0]
                end = exons[-1][1]
                rescue_reason = rec.get("rescue_reason") or ""
                rescue_applied = rec.get("rescued_by_shortread", False)
                common = (
                    f'gene_id "{gid}"; transcript_id "{tid}"; '
                    f'read_support "{rec["read_support"]}"; '
                    f'cluster_size "{rec["cluster_size"]}"; '
                    f'rescue_applied "{rescue_applied}"; '
                    f'rescue_reason "{rescue_reason}";'
                )
                fh.write(
                    f"{seqname}\t{source_label}\ttranscript\t{start}\t{end}\t.\t"
                    f"{strand}\t.\t{common}\n"
                )
                for i, (s, e) in enumerate(exons, start=1):
                    fh.write(
                        f"{seqname}\t{source_label}\texon\t{s}\t{e}\t.\t"
                        f'{strand}\t.\tgene_id "{gid}"; transcript_id "{tid}"; '
                        f'exon_number "{i}"; read_support "{rec["read_support"]}";\n'
                    )
