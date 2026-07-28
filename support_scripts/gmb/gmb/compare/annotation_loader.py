"""Unified annotation loader for multiple GFF3/GTF formats.

Supports:
  - Ensembl GFF3 (ID=gene:...; Parent=transcript:...)
  - Helixer GFF3 (standard ID=/Parent=)
  - ANNEVO GFF3 (standard ID=/Parent=, confidence scores in CDS)
  - Tiberius hybrid (gene/transcript with bare attrs, exon/CDS with GTF-style attrs)

The loader auto-detects format from file content, parses into a normalised
internal representation, and returns genes/transcripts/exons/CDS DataFrames
with consistent columns ready for comparison.

Cross-chromosome ID collisions: some tools (notably Tiberius) restart gene
numbering ("g1", "g2", ...) independently on every chromosome instead of
emitting genome-wide unique IDs. The loader detects gene_id/transcript_id
values that occur on more than one seqid and namespaces them internally
(e.g. "1:g1" vs "2:g1") so distinct per-chromosome genes/transcripts are
never collapsed together. Original IDs are preserved in the
original_gene_id/original_transcript_id columns. This is a no-op for
inputs (e.g. standard Ensembl GFF3/GTF) whose IDs are already
genome-wide unique. See `_namespace_duplicate_ids`.
"""

from __future__ import annotations

import gzip
import re
import warnings
from collections import defaultdict
from typing import Optional

import pandas as pd

# ---------------------------------------------------------------------------
# Format detection
# ---------------------------------------------------------------------------

_GFF3_ATTR = re.compile(r"ID=|Parent=")
_GTF_ATTR = re.compile(r'(gene_id|transcript_id)\s+"')
_TIBERIUS_BARE = re.compile(r"^\S+$")


def detect_format(path: str) -> str:
    """Detect annotation format by inspecting file content.

    Returns one of: "gff3", "gtf", "tiberius"
    """
    opener = gzip.open if path.endswith(".gz") else open

    gff3_score = 0
    gtf_score = 0
    bare_gene_tx = 0
    gtf_child = 0
    lines_checked = 0

    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                if "gff-version" in line and "3" in line:
                    gff3_score += 5
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            feat = parts[2]
            attrs = parts[8].strip()
            lines_checked += 1

            if _GFF3_ATTR.search(attrs):
                gff3_score += 1
            if _GTF_ATTR.search(attrs):
                gtf_score += 1

            if feat in ("gene", "transcript") and _TIBERIUS_BARE.match(attrs):
                bare_gene_tx += 1

            if (
                feat in ("exon", "CDS")
                and _GTF_ATTR.search(attrs)
                and not _GFF3_ATTR.search(attrs)
            ):
                gtf_child += 1

            if lines_checked >= 200:
                break

    if bare_gene_tx >= 2 and gtf_child >= 2:
        return "tiberius"
    if gff3_score > gtf_score:
        return "gff3"
    if gtf_score > 0:
        return "gtf"
    return "gff3"


# ---------------------------------------------------------------------------
# Attribute parsers
# ---------------------------------------------------------------------------


def _parse_gff3_attrs(attrs_str: str) -> dict[str, str]:
    """Parse GFF3 key=value;key=value attributes."""
    attrs: dict[str, str] = {}
    for kv in attrs_str.split(";"):
        kv = kv.strip()
        if "=" in kv:
            k, v = kv.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs


def _parse_gtf_attrs(attrs_str: str) -> dict[str, str]:
    """Parse GTF key "value"; key "value"; attributes."""
    attrs: dict[str, str] = {}
    for m in re.finditer(r'(\w+)\s+"([^"]*)"', attrs_str):
        attrs[m.group(1)] = m.group(2)
    for m in re.finditer(r"(\w+)=(\S+)", attrs_str):
        if m.group(1) not in attrs:
            attrs[m.group(1)] = m.group(2).rstrip(";")
    return attrs


# ---------------------------------------------------------------------------
# Strip Ensembl ID prefixes
# ---------------------------------------------------------------------------

_ENSEMBL_PREFIX = re.compile(r"^(gene|transcript|chromosome|mRNA):")


def _strip_ensembl_prefix(val: str) -> str:
    return _ENSEMBL_PREFIX.sub("", val)


# ---------------------------------------------------------------------------
# Biotype extraction
# ---------------------------------------------------------------------------


def _extract_biotype(row: dict, attrs: dict, feat: str) -> None:
    """Extract gene_biotype and transcript_biotype from parsed attributes."""
    biotype = attrs.get("biotype", "")
    gene_biotype = attrs.get(
        "gene_biotype", biotype if feat in ("gene", "ncRNA_gene", "pseudogene") else ""
    )
    tx_biotype = attrs.get(
        "transcript_biotype",
        biotype if feat not in ("gene", "ncRNA_gene", "pseudogene", "exon", "CDS") else "",
    )
    row["gene_biotype"] = gene_biotype
    row["transcript_biotype"] = tx_biotype
    row["tags"] = attrs.get("tag", "")


# ---------------------------------------------------------------------------
# Line-by-line parser
# ---------------------------------------------------------------------------


def _parse_annotation_lines(path: str, fmt: str, source_label: str) -> list[dict]:
    """Parse annotation file into a list of raw row dicts."""
    opener = gzip.open if path.endswith(".gz") else open
    rows = []

    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            seqname = parts[0]
            source = parts[1]
            feat = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attrs_str = parts[8].strip()

            if feat in (
                "biological_region",
                "chromosome",
                "region",
                "intron",
                "start_codon",
                "stop_codon",
                "five_prime_UTR",
                "three_prime_UTR",
            ):
                continue

            row = {
                "seqname": seqname,
                "source": source,
                "feature": feat,
                "start": start,
                "end": end,
                "strand": strand,
                "gene_id": "",
                "transcript_id": "",
                "ID": "",
                "Parent": "",
                "gene_biotype": "",
                "transcript_biotype": "",
                "tags": "",
            }

            if fmt == "tiberius":
                _parse_tiberius_row(row, feat, attrs_str)
            elif fmt == "gtf":
                attrs = _parse_gtf_attrs(attrs_str)
                row["gene_id"] = attrs.get("gene_id", "")
                row["transcript_id"] = attrs.get("transcript_id", "")
                _extract_biotype(row, attrs, feat)
                if feat == "gene":
                    row["ID"] = row["gene_id"]
                elif feat in ("transcript", "mRNA"):
                    row["ID"] = row["transcript_id"]
                    row["Parent"] = row["gene_id"]
                else:
                    row["Parent"] = row["transcript_id"]
            else:
                attrs = _parse_gff3_attrs(attrs_str)
                raw_id = attrs.get("ID", "")
                raw_parent = attrs.get("Parent", "")
                row["ID"] = _strip_ensembl_prefix(raw_id)
                row["Parent"] = _strip_ensembl_prefix(raw_parent)
                row["gene_id"] = attrs.get("gene_id", "")
                row["transcript_id"] = attrs.get("transcript_id", "")
                _extract_biotype(row, attrs, feat)

            rows.append(row)

    return rows


def _parse_tiberius_row(row: dict, feat: str, attrs_str: str) -> None:
    """Handle Tiberius hybrid format: bare attrs on gene/transcript, GTF on children."""
    attrs_str = attrs_str.strip()

    if feat == "gene":
        gid = attrs_str.strip().rstrip(";")
        row["ID"] = gid
        row["gene_id"] = gid
    elif feat == "transcript":
        tid = attrs_str.strip().rstrip(";")
        row["ID"] = tid
        row["transcript_id"] = tid
        dot_pos = tid.rfind(".")
        if dot_pos > 0:
            row["gene_id"] = tid[:dot_pos]
            row["Parent"] = tid[:dot_pos]
        else:
            row["Parent"] = tid
            row["gene_id"] = tid
    else:
        attrs = _parse_gtf_attrs(attrs_str)
        row["gene_id"] = attrs.get("gene_id", "")
        row["transcript_id"] = attrs.get("transcript_id", "")
        row["Parent"] = row["transcript_id"]


# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------


def _normalise_rows(rows: list[dict], fmt: str) -> list[dict]:
    """Normalise parsed rows into consistent internal representation.

    - transcript -> mRNA
    - Infer missing gene_id / transcript_id from ID / Parent
    - Build parent lookups for GFF3 hierarchies
    """
    id_to_row: dict[str, dict] = {}
    for r in rows:
        if r["ID"]:
            id_to_row[r["ID"]] = r

    for r in rows:
        feat = r["feature"]

        # Normalise transcript -> mRNA
        if feat == "transcript":
            r["feature"] = "mRNA"

        # Normalise biotype transcript features (lnc_RNA, etc.) for comparison
        if feat in (
            "lnc_RNA",
            "ncRNA",
            "rRNA",
            "tRNA",
            "snRNA",
            "snoRNA",
            "miRNA",
            "pre_miRNA",
            "SRP_RNA",
            "RNase_P_RNA",
            "RNase_MRP_RNA",
            "telomerase_RNA",
            "scRNA",
            "processed_transcript",
            "V_gene_segment",
            "D_gene_segment",
            "J_gene_segment",
            "C_gene_segment",
            "pseudogenic_transcript",
        ):
            r["feature"] = "mRNA"

        # Also handle ncRNA_gene, pseudogene as gene
        if feat in ("ncRNA_gene", "pseudogene"):
            r["feature"] = "gene"

        # Fill in gene_id and transcript_id from ID/Parent where missing
        if r["feature"] == "gene" and not r["gene_id"]:
            r["gene_id"] = r["ID"]

        if r["feature"] == "mRNA":
            if not r["transcript_id"]:
                r["transcript_id"] = r["ID"]
            if not r["gene_id"]:
                r["gene_id"] = r["Parent"]

        if r["feature"] in ("exon", "CDS"):
            if not r["transcript_id"] and r["Parent"]:
                r["transcript_id"] = r["Parent"]
            if not r["gene_id"] and r["Parent"]:
                parent_row = id_to_row.get(r["Parent"])
                if parent_row:
                    r["gene_id"] = parent_row.get("gene_id", parent_row.get("Parent", ""))

    # Propagate biotypes down the hierarchy
    gene_biotype_map: dict[str, str] = {}
    tx_biotype_map: dict[str, str] = {}
    for r in rows:
        if r["feature"] == "gene" and r["gene_biotype"]:
            gene_biotype_map[r["gene_id"]] = r["gene_biotype"]
        if r["feature"] == "mRNA" and r["transcript_biotype"]:
            tx_biotype_map[r["transcript_id"]] = r["transcript_biotype"]

    for r in rows:
        if not r["gene_biotype"] and r["gene_id"] in gene_biotype_map:
            r["gene_biotype"] = gene_biotype_map[r["gene_id"]]
        if not r["transcript_biotype"] and r["transcript_id"] in tx_biotype_map:
            r["transcript_biotype"] = tx_biotype_map[r["transcript_id"]]

    return rows


# ---------------------------------------------------------------------------
# Cross-chromosome ID collision handling
# ---------------------------------------------------------------------------
#
# Some gene predictors restart gene/transcript numbering per sequence rather
# than emitting genome-wide unique IDs. Tiberius is the motivating case: its
# GTF output uses "g1", "g2", ... / "g1.t1", "g2.t1", ... independently on
# every chromosome, so gene_id "g1" on chromosome 1 is a completely different
# gene from gene_id "g1" on chromosome 2. Left alone, this collides in every
# downstream step that keys off gene_id/transcript_id (locus classification,
# transcript grouping, canonical-transcript selection), silently merging
# unrelated features from different chromosomes.
#
# This is generic (keyed on "does this ID appear on >1 seqname", not on any
# Tiberius-specific signature), so it is a no-op for normal Ensembl/GFF3/GTF
# inputs where gene_id/transcript_id are already genome-wide unique.


def _namespace_duplicate_ids(rows: list[dict], source_label: str = "") -> list[dict]:
    """Rewrite gene_id/transcript_id (and dependent ID/Parent) so that IDs
    which are reused across more than one seqname become unique internally.

    Runs after `_normalise_rows`, so gene_id/transcript_id are already filled
    in for every feature. IDs that only ever appear on a single seqname are
    left untouched. Renamed IDs are namespaced as "{seqname}:{original_id}";
    the pre-rename value is preserved in `original_gene_id` /
    `original_transcript_id` on every row. The file on disk is never
    modified.
    """
    if not rows:
        return rows

    gene_id_seqnames: dict[str, set] = defaultdict(set)
    tx_id_seqnames: dict[str, set] = defaultdict(set)
    for r in rows:
        if r["gene_id"]:
            gene_id_seqnames[r["gene_id"]].add(r["seqname"])
        if r["transcript_id"]:
            tx_id_seqnames[r["transcript_id"]].add(r["seqname"])

    dup_gene_ids = {gid for gid, seqs in gene_id_seqnames.items() if len(seqs) > 1}
    dup_tx_ids = {tid for tid, seqs in tx_id_seqnames.items() if len(seqs) > 1}

    for r in rows:
        r["original_gene_id"] = r["gene_id"]
        r["original_transcript_id"] = r["transcript_id"]

    if not dup_gene_ids and not dup_tx_ids:
        return rows

    label = f" in {source_label}" if source_label else ""
    print(
        f"  Detected gene_id/transcript_id values reused across multiple "
        f"seqids{label} ({len(dup_gene_ids)} gene_id, {len(dup_tx_ids)} "
        f"transcript_id) — namespacing by seqid to keep per-chromosome "
        f"features distinct (original IDs kept in original_gene_id / "
        f"original_transcript_id)."
    )

    for r in rows:
        seqname = r["seqname"]
        old_gid = r["gene_id"]
        old_tid = r["transcript_id"]

        new_gid = f"{seqname}:{old_gid}" if old_gid in dup_gene_ids else old_gid
        new_tid = f"{seqname}:{old_tid}" if old_tid in dup_tx_ids else old_tid

        if new_gid != old_gid:
            r["gene_id"] = new_gid
            if r["feature"] == "gene" and r["ID"] == old_gid:
                r["ID"] = new_gid
            if r["feature"] == "mRNA" and r["Parent"] == old_gid:
                r["Parent"] = new_gid

        if new_tid != old_tid:
            r["transcript_id"] = new_tid
            if r["feature"] == "mRNA" and r["ID"] == old_tid:
                r["ID"] = new_tid
            if r["feature"] in ("exon", "CDS") and r["Parent"] == old_tid:
                r["Parent"] = new_tid

    return rows


# ---------------------------------------------------------------------------
# DataFrame construction
# ---------------------------------------------------------------------------


def _rows_to_dataframes(
    rows: list[dict], source_label: str
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Convert normalised rows into genes, mRNA, exons, CDS DataFrames.

    All DataFrames share columns: Chromosome, Start, End, Strand, Feature,
    gene_id, transcript_id, ID, Parent, Source_label.
    """
    if not rows:
        empty = pd.DataFrame(
            columns=[
                "Chromosome",
                "Start",
                "End",
                "Strand",
                "Feature",
                "gene_id",
                "transcript_id",
                "ID",
                "Parent",
                "Source_label",
                "gene_biotype",
                "transcript_biotype",
                "tags",
                "original_gene_id",
                "original_transcript_id",
            ]
        )
        return empty.copy(), empty.copy(), empty.copy(), empty.copy()

    df = pd.DataFrame(rows)
    df.rename(
        columns={
            "seqname": "Chromosome",
            "feature": "Feature",
            "start": "Start",
            "end": "End",
            "strand": "Strand",
        },
        inplace=True,
    )
    df["Source_label"] = source_label

    # Fill remaining blanks
    df["transcript_id"] = df["transcript_id"].replace("", pd.NA)
    df["gene_id"] = df["gene_id"].replace("", pd.NA)

    # For exon/CDS rows without transcript_id, try Parent
    child_mask = df["Feature"].isin(["exon", "CDS"])
    df.loc[child_mask, "transcript_id"] = df.loc[child_mask, "transcript_id"].fillna(
        df.loc[child_mask, "Parent"]
    )

    # Last-resort fills
    df["transcript_id"] = df["transcript_id"].fillna(df["ID"]).fillna("unknown")
    df["gene_id"] = df["gene_id"].fillna(df["Parent"]).fillna(df["ID"]).fillna("unknown")

    keep_cols = [
        "Chromosome",
        "Start",
        "End",
        "Strand",
        "Feature",
        "gene_id",
        "transcript_id",
        "ID",
        "Parent",
        "Source_label",
        "gene_biotype",
        "transcript_biotype",
        "tags",
        "original_gene_id",
        "original_transcript_id",
    ]
    for c in keep_cols:
        if c not in df.columns:
            df[c] = ""
    df = df[keep_cols]

    genes = df[df["Feature"] == "gene"].copy()
    mrna = df[df["Feature"] == "mRNA"].copy()
    exons = df[df["Feature"] == "exon"].copy()
    cds = df[df["Feature"] == "CDS"].copy()

    return genes, mrna, exons, cds


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def load_annotation(
    path: str,
    source_label: str = "Annotation",
    format_hint: str = "auto",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load an annotation file and return normalised DataFrames.

    Parameters
    ----------
    path : str
        Path to annotation file (GFF3, GTF, or Tiberius hybrid). May be gzipped.
    source_label : str
        Label for the source (e.g. "Reference", "Helixer").
    format_hint : str
        One of "auto", "gff3", "gtf", "tiberius". Default "auto" detects from content.

    Returns
    -------
    tuple of (genes, mrna, exons, cds) DataFrames.
    Each DataFrame has columns: Chromosome, Start, End, Strand, Feature,
    gene_id, transcript_id, ID, Parent, Source_label.
    """
    if format_hint == "auto":
        fmt = detect_format(path)
    else:
        fmt = format_hint

    print(f"  Loading {source_label} from {path} (format: {fmt})...")

    rows = _parse_annotation_lines(path, fmt, source_label)
    rows = _normalise_rows(rows, fmt)
    rows = _namespace_duplicate_ids(rows, source_label)
    genes, mrna, exons, cds = _rows_to_dataframes(rows, source_label)

    n_genes = len(genes)
    n_tx = mrna["transcript_id"].nunique() if not mrna.empty else 0
    n_exons = len(exons)
    n_cds = len(cds)
    print(
        f"  {source_label}: {n_genes} genes, {n_tx} transcripts, " f"{n_exons} exons, {n_cds} CDS"
    )

    if n_genes == 0:
        warnings.warn(f"No gene features found in {path}. Check format detection.", stacklevel=2)
    if n_exons == 0:
        warnings.warn(f"No exon features found in {path}.", stacklevel=2)

    return genes, mrna, exons, cds


# ---------------------------------------------------------------------------
# Biotype filtering
# ---------------------------------------------------------------------------


def filter_by_biotype(
    genes: pd.DataFrame,
    mrna: pd.DataFrame,
    exons: pd.DataFrame,
    cds: pd.DataFrame,
    gene_biotypes: Optional[list[str]] = None,
    transcript_biotypes: Optional[list[str]] = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Filter annotation DataFrames by gene and/or transcript biotype.

    Returns filtered (genes, mrna, exons, cds).
    """
    if not gene_biotypes and not transcript_biotypes:
        return genes, mrna, exons, cds

    if gene_biotypes:
        keep_genes = genes[genes["gene_biotype"].isin(gene_biotypes)]
        keep_gene_ids = set(keep_genes["gene_id"])
        genes = keep_genes
        mrna = mrna[mrna["gene_id"].isin(keep_gene_ids)]
        exons = exons[exons["gene_id"].isin(keep_gene_ids)]
        cds = cds[cds["gene_id"].isin(keep_gene_ids)]

    if transcript_biotypes:
        keep_tx = mrna[mrna["transcript_biotype"].isin(transcript_biotypes)]
        keep_tids = set(keep_tx["transcript_id"])
        mrna = keep_tx
        exons = exons[exons["transcript_id"].isin(keep_tids)]
        cds = cds[cds["transcript_id"].isin(keep_tids)]
        keep_gene_ids = set(mrna["gene_id"])
        genes = genes[genes["gene_id"].isin(keep_gene_ids)]

    return genes.copy(), mrna.copy(), exons.copy(), cds.copy()


# ---------------------------------------------------------------------------
# Transcript selection
# ---------------------------------------------------------------------------


def select_transcripts(
    genes: pd.DataFrame,
    mrna: pd.DataFrame,
    exons: pd.DataFrame,
    cds: pd.DataFrame,
    mode: str = "all",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Select representative transcripts per gene.

    Modes:
      - all: keep everything (default)
      - canonical: keep Ensembl canonical transcript (tag=Ensembl_canonical),
                   falling back to longest CDS transcript, then longest transcript
      - longest_cds: keep the transcript with the longest total CDS per gene,
                     falling back to longest transcript when no CDS

    Returns filtered (genes, mrna, exons, cds).
    """
    if mode == "all":
        return genes, mrna, exons, cds

    if mrna.empty:
        return genes, mrna, exons, cds

    selected_tids: set[str] = set()

    # Precompute total CDS length per transcript (vectorised, avoids O(n²) loop)
    if not cds.empty:
        cds_lengths = (
            cds.assign(_len=cds["End"] - cds["Start"]).groupby("transcript_id")["_len"].sum()
        )
    else:
        cds_lengths = pd.Series(dtype="int64")

    # Precompute total exon span per transcript for longest-transcript fallback
    if not exons.empty:
        tx_spans = (
            exons.assign(_len=exons["End"] - exons["Start"]).groupby("transcript_id")["_len"].sum()
        )
    else:
        tx_spans = pd.Series(dtype="int64")

    for _gene_id, gene_mrna in mrna.groupby("gene_id"):
        if mode == "canonical":
            canonical = gene_mrna[gene_mrna["tags"].str.contains("Ensembl_canonical", na=False)]
            if not canonical.empty:
                selected_tids.add(canonical.iloc[0]["transcript_id"])
                continue

        # Tier 1: longest CDS transcript
        gene_tids = gene_mrna["transcript_id"].unique()
        best_tid = None
        best_cds_len = -1
        for tid in gene_tids:
            total_cds = cds_lengths.get(tid, 0)
            if total_cds > best_cds_len:
                best_cds_len = total_cds
                best_tid = tid

        # Tier 2: longest transcript by exon span (when no CDS found)
        if best_cds_len <= 0:
            best_span = -1
            for tid in gene_tids:
                span = tx_spans.get(tid, 0)
                if span > best_span:
                    best_span = span
                    best_tid = tid

        # Tier 3: first transcript (should never be reached)
        if best_tid is None:
            best_tid = gene_mrna.iloc[0]["transcript_id"]

        selected_tids.add(best_tid)

    mrna = mrna[mrna["transcript_id"].isin(selected_tids)].copy()
    exons = exons[exons["transcript_id"].isin(selected_tids)].copy()
    cds = cds[cds["transcript_id"].isin(selected_tids)].copy()

    return genes, mrna, exons, cds


# ---------------------------------------------------------------------------
# Evaluation mode presets
# ---------------------------------------------------------------------------


def apply_evaluation_mode(
    mode: str,
    genes: pd.DataFrame,
    mrna: pd.DataFrame,
    exons: pd.DataFrame,
    cds: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Apply a preset evaluation mode to filter reference annotations.

    Modes:
      - all: no filtering
      - protein_coding: only protein-coding genes and transcripts
      - cds_only: only transcripts that have CDS features
      - canonical: one canonical/longest-CDS transcript per gene
    """
    if mode == "all":
        return genes, mrna, exons, cds

    if mode == "protein_coding":
        genes, mrna, exons, cds = filter_by_biotype(
            genes, mrna, exons, cds, gene_biotypes=["protein_coding"]
        )
        return genes, mrna, exons, cds

    if mode == "cds_only":
        tids_with_cds = set(cds["transcript_id"].unique()) if not cds.empty else set()
        mrna = mrna[mrna["transcript_id"].isin(tids_with_cds)].copy()
        exons = exons[exons["transcript_id"].isin(tids_with_cds)].copy()
        cds = cds[cds["transcript_id"].isin(tids_with_cds)].copy()
        keep_gene_ids = set(mrna["gene_id"])
        genes = genes[genes["gene_id"].isin(keep_gene_ids)].copy()
        return genes, mrna, exons, cds

    if mode == "canonical":
        return select_transcripts(genes, mrna, exons, cds, mode="canonical")

    warnings.warn(f"Unknown evaluation mode '{mode}', using 'all'.", stacklevel=2)
    return genes, mrna, exons, cds


# ---------------------------------------------------------------------------
# Filter audit
# ---------------------------------------------------------------------------


def generate_filter_audit(
    genes: pd.DataFrame,
    mrna: pd.DataFrame,
    exons: pd.DataFrame,
    cds: pd.DataFrame,
    output_dir: str,
    evaluation_mode: str = "all",
    transcript_selection: str = "all",
    gene_biotypes_filter: Optional[str] = None,
    transcript_biotypes_filter: Optional[str] = None,
    pre_filter_genes: int = 0,
    pre_filter_transcripts: int = 0,
) -> dict:
    """Generate an audit of the reference filtering for traceability.

    Writes reference_filter_audit.tsv and reference_filter_audit.json.
    Returns the audit dict.
    """
    import csv
    import json
    import os

    gene_biotype_counts = dict(genes["gene_biotype"].value_counts()) if not genes.empty else {}
    tx_biotype_counts = dict(mrna["transcript_biotype"].value_counts()) if not mrna.empty else {}

    tids_with_cds = set(cds["transcript_id"].unique()) if not cds.empty else set()
    cds_tx_count = len(tids_with_cds)

    canonical_count = 0
    if not mrna.empty and "tags" in mrna.columns:
        canonical_count = int(mrna["tags"].str.contains("Ensembl_canonical", na=False).sum())

    post_genes = len(genes)
    post_tx = mrna["transcript_id"].nunique() if not mrna.empty else 0

    audit = {
        "evaluation_mode": evaluation_mode,
        "transcript_selection": transcript_selection,
        "gene_biotypes_filter": gene_biotypes_filter,
        "transcript_biotypes_filter": transcript_biotypes_filter,
        "pre_filter_genes": pre_filter_genes,
        "pre_filter_transcripts": pre_filter_transcripts,
        "post_filter_genes": post_genes,
        "post_filter_transcripts": post_tx,
        "post_filter_exons": len(exons),
        "post_filter_cds": len(cds),
        "cds_containing_transcripts": cds_tx_count,
        "canonical_transcripts": canonical_count,
        "gene_biotype_counts": gene_biotype_counts,
        "transcript_biotype_counts": tx_biotype_counts,
    }

    os.makedirs(output_dir, exist_ok=True)

    json_path = os.path.join(output_dir, "reference_filter_audit.json")
    with open(json_path, "w") as fh:
        json.dump(audit, fh, indent=2, default=str)

    tsv_path = os.path.join(output_dir, "reference_filter_audit.tsv")
    with open(tsv_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["category", "key", "value"])
        w.writerow(["config", "evaluation_mode", evaluation_mode])
        w.writerow(["config", "transcript_selection", transcript_selection])
        w.writerow(["config", "gene_biotypes_filter", gene_biotypes_filter or "none"])
        w.writerow(["config", "transcript_biotypes_filter", transcript_biotypes_filter or "none"])
        w.writerow(["counts", "pre_filter_genes", pre_filter_genes])
        w.writerow(["counts", "pre_filter_transcripts", pre_filter_transcripts])
        w.writerow(["counts", "post_filter_genes", post_genes])
        w.writerow(["counts", "post_filter_transcripts", post_tx])
        w.writerow(["counts", "post_filter_exons", len(exons)])
        w.writerow(["counts", "post_filter_cds", len(cds)])
        w.writerow(["counts", "cds_containing_transcripts", cds_tx_count])
        w.writerow(["counts", "canonical_transcripts", canonical_count])
        for bt, cnt in sorted(gene_biotype_counts.items(), key=lambda x: -x[1]):
            w.writerow(["gene_biotype", bt, cnt])
        for bt, cnt in sorted(tx_biotype_counts.items(), key=lambda x: -x[1]):
            w.writerow(["transcript_biotype", bt, cnt])

    print(f"  Filter audit: {json_path}")
    return audit


# ---------------------------------------------------------------------------
# FASTA validation
# ---------------------------------------------------------------------------


def validate_against_fasta(
    genes: pd.DataFrame,
    exons: pd.DataFrame,
    cds: pd.DataFrame,
    fasta_path: str,
    label: str = "annotation",
) -> list[str]:
    """Validate annotation coordinates against a FASTA file.

    Checks:
    - Every seqname in annotation exists in FASTA
    - Feature end coordinates do not exceed sequence length
    - Detects possible chr-prefix mismatches

    Returns a list of diagnostic messages (empty = all OK).
    """
    from gmb.utils.fasta import load_genome

    genome = load_genome(fasta_path)
    fasta_seqnames = set(genome.keys())
    fasta_lengths = {name: len(seq) for name, seq in genome.items()}

    all_dfs = [df for df in [genes, exons, cds] if not df.empty]
    if not all_dfs:
        return [f"No features to validate for {label}."]

    anno_seqnames = set()
    for df in all_dfs:
        anno_seqnames.update(df["Chromosome"].unique())

    diagnostics: list[str] = []

    # Check seqname presence
    missing = anno_seqnames - fasta_seqnames
    if missing:
        # Try chr-prefix detection
        chr_stripped = {s.replace("chr", "", 1): s for s in missing if s.startswith("chr")}
        chr_added = {f"chr{s}": s for s in missing if not s.startswith("chr")}

        fixable_strip = {orig for bare, orig in chr_stripped.items() if bare in fasta_seqnames}
        fixable_add = {orig for prefixed, orig in chr_added.items() if prefixed in fasta_seqnames}

        if fixable_strip or fixable_add:
            diagnostics.append(
                f"[{label}] Seqname mismatch may be a chr-prefix issue. "
                f"Annotation has {sorted(fixable_strip | fixable_add)}, "
                f"FASTA uses {'stripped' if fixable_strip else 'prefixed'} names. "
                f"Consider --seqname-map."
            )

        truly_missing = missing - fixable_strip - fixable_add
        if truly_missing:
            sample = sorted(truly_missing)[:5]
            diagnostics.append(
                f"[{label}] {len(truly_missing)} seqname(s) not in FASTA: "
                f"{sample}{'...' if len(truly_missing) > 5 else ''}"
            )

    # Check coordinate bounds
    out_of_bounds = 0
    for df in all_dfs:
        for seqname in df["Chromosome"].unique():
            if seqname not in fasta_lengths:
                continue
            seq_len = fasta_lengths[seqname]
            mask = (df["Chromosome"] == seqname) & (df["End"] > seq_len)
            count = mask.sum()
            if count > 0:
                out_of_bounds += count

    if out_of_bounds:
        diagnostics.append(
            f"[{label}] {out_of_bounds} feature(s) have end coordinates "
            f"exceeding FASTA sequence length."
        )

    if not diagnostics:
        print(f"  [{label}] All coordinates validated against FASTA.")
    else:
        for d in diagnostics:
            warnings.warn(d, stacklevel=2)

    return diagnostics
