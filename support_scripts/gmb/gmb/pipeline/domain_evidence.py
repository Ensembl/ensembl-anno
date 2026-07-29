#!/usr/bin/env python3
"""Protein-domain/feature evidence: normalised sidecar table + adapters.

Not required for the first-pass canonical selector (gmb.pipeline.
canonical_selection runs fine with no domain data at all -- see its
``domains.enabled: false`` default). This module exists so that when Pfam/
InterProScan results *are* available, they can be converted into one
normalised table that canonical_selection consumes generically, instead of
each canonical-scoring code path needing to know about every upstream tool's
raw output format.

Normalised sidecar schema (one row per domain/feature hit)
------------------------------------------------------------
    transcript_id     -- must match the GMB consensus.gff3 transcript ID
    provider          -- e.g. "Pfam", "InterPro", "PANTHER", "CDD"
    feature_accession -- e.g. "PF00069"
    feature_name      -- e.g. "Protein kinase domain"
    start             -- 1-based, protein (amino-acid) coordinates
    end               -- 1-based, inclusive
    score             -- provider's own score (bit score, HMM score, ...)
    evalue            -- provider's own E-value, if any
    coverage          -- fraction of the *domain's own* HMM/profile length
                         aligned (0-1), when the source format reports it
    status            -- provider's own call, e.g. "significant" / "partial"
                         (never invented here -- passed through as-is)

Adapters
--------
Each ``parse_*`` function below takes a raw tool output path and returns a
list of dicts matching the schema above -- the only place that needs to
know the raw format. Add a new adapter for any other tool the same way and
nothing in canonical_selection.py needs to change.

Implemented: HMMER domtblout, InterProScan TSV.
Documented but not implemented (same pattern, not needed for this pass):
InterProScan GFF3 -- would follow the same generic-attribute-parsing style
already used in gmb.utils.gff.parse_gff3_hierarchy (split 9 tab-separated
columns, parse ``;``/``=``-delimited attributes for accession/name/score).
"""

from __future__ import annotations

from collections import defaultdict
from typing import Optional

DOMAIN_SIDECAR_COLUMNS = [
    "transcript_id",
    "provider",
    "feature_accession",
    "feature_name",
    "start",
    "end",
    "score",
    "evalue",
    "coverage",
    "status",
]


def parse_hmmer_domtblout(path: str, provider: str = "Pfam") -> list[dict]:
    """Parse HMMER ``hmmscan --domtblout`` output into the sidecar schema.

    domtblout is whitespace-delimited with a fixed 22-column layout; the
    last column (description) may itself contain spaces, so it is not
    parsed here (feature_name comes from column 4, the query/target's own
    short name, which does not).

    Columns used (1-based, per HMMER's documented domtblout format --
    verified against the format spec: target=the HMM/profile being scanned
    against, query=the input protein sequence, in `hmmscan` orientation):
      1  target name        -> feature_accession (fallback: feature_name)
      2  target accession   -> feature_accession (preferred, if not "-")
      3  tlen (target/HMM profile length, NOT the query/protein length --
         used as the coverage denominator below)
      4  query name         -> transcript_id
      12 c-Evalue (conditional E-value; deliberately NOT column 13's
         i-Evalue/independent E-value -- c-Evalue is the standard
         per-domain significance measure when a query may contain
         several domain copies, which independent E-value is not
         designed for; not verified against a real installed HMMER
         version, since none is available in this environment)
      14 this domain's own bit score -> score
      16 hmm coord from, 17 hmm coord to -- used to estimate coverage
         against col 3 (tlen)
      18 ali coord from, 19 ali coord to -- alignment coordinates in the
         QUERY (protein) sequence -> start/end
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 22:
                continue
            target_name = parts[0]
            target_acc = parts[1]
            tlen = _safe_int(parts[2])
            query_name = parts[3]
            hmm_from = _safe_int(parts[15])
            hmm_to = _safe_int(parts[16])
            ali_from = _safe_int(parts[17])
            ali_to = _safe_int(parts[18])
            c_evalue = _safe_float(parts[11])
            dom_score = _safe_float(parts[13])

            coverage = None
            if tlen and hmm_from is not None and hmm_to is not None and tlen > 0:
                coverage = round((hmm_to - hmm_from + 1) / tlen, 4)

            accession = target_acc if target_acc and target_acc != "-" else target_name
            rows.append(
                {
                    "transcript_id": query_name,
                    "provider": provider,
                    "feature_accession": accession,
                    "feature_name": target_name,
                    "start": ali_from,
                    "end": ali_to,
                    "score": dom_score,
                    "evalue": c_evalue,
                    "coverage": coverage,
                    "status": None,
                }
            )
    return rows


def parse_interproscan_tsv(path: str) -> list[dict]:
    """Parse InterProScan's tab-separated output into the sidecar schema.

    InterProScan TSV has no header row. Standard 11-column layout (1-based,
    per InterProScan's documented "TSV output format" -- NOT verified
    against a real installed InterProScan version, since none is available
    in this environment; confirm this mapping against `interproscan.sh
    --formats tsv` output for whatever version is actually used before
    relying on it for anything beyond QC inspection):
      1  Protein accession           -> transcript_id
      4  Analysis (provider), e.g. Pfam/PANTHER/CDD/SUPERFAMILY
      5  Signature accession         -> feature_accession
      6  Signature description       -> feature_name
      7  Start location (protein aa) -> start
      8  Stop location (protein aa)  -> end
      9  Score -- for e-value-based analyses (e.g. Pfam) this column *is*
         the significance value, not a separate bit-score; InterProScan's
         format does not report score and E-value as two distinct fields
         the way HMMER domtblout does. Stored here under ``score`` for
         provenance; ``evalue`` is left ``None`` rather than guessing
         which analyses' column-9 values are safe to treat as E-values
         and which are not (e.g. Coils/TMHMM report non-numeric values
         here) -- summarize_transcript_domain_metrics()'s significance
         filtering therefore currently cannot threshold InterProScan-
         sourced rows on E-value; this is a known limitation to revisit
         once a real installed version's output can be inspected.
      10 Status ("T"/"?") -- passed through under ``status`` rather than
         reinterpreted (previously misread from column 12, which is the
         *InterPro accession*, not the status flag -- fixed).
      11 Date of the run (not used)
      12+ optional InterPro accession/description, then GO/pathway columns
         if --iprlookup/--goterms/--pathways were used -- column count
         varies by version/flags; only the first 10 (always present) are
         relied on here.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            transcript_id = parts[0]
            provider = parts[3]
            accession = parts[4]
            name = parts[5] if len(parts) > 5 else None
            start = _safe_int(parts[6]) if len(parts) > 6 else None
            end = _safe_int(parts[7]) if len(parts) > 7 else None
            score_raw = parts[8] if len(parts) > 8 else None
            score = _safe_float(score_raw) if score_raw not in (None, "-") else None
            status = parts[9] if len(parts) > 9 else None

            rows.append(
                {
                    "transcript_id": transcript_id,
                    "provider": provider,
                    "feature_accession": accession,
                    "feature_name": name,
                    "start": start,
                    "end": end,
                    "score": score,
                    "evalue": None,
                    "coverage": None,
                    "status": status,
                }
            )
    return rows


def _safe_int(x) -> Optional[int]:
    try:
        return int(x)
    except (TypeError, ValueError):
        return None


def _safe_float(x) -> Optional[float]:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Per-transcript summarisation (what canonical_selection.py actually reads)
# ---------------------------------------------------------------------------


def summarize_transcript_domain_metrics(
    sidecar_rows: list[dict],
    protein_lengths: Optional[dict[str, int]] = None,
    significant_evalue: float = 1e-5,
) -> dict[str, dict]:
    """Reduce raw per-domain sidecar rows to one summary dict per transcript.

    canonical_selection.py's scoring code consumes only this summary (never
    the raw per-domain rows), so a new adapter never requires touching the
    scorer -- only this function needs to know how to fold its rows in.

    Parameters
    ----------
    sidecar_rows : list of dict
        Rows matching DOMAIN_SIDECAR_COLUMNS (any mix of providers/adapters).
    protein_lengths : dict or None
        ``{transcript_id: aa_length}``, used to compute coverage fraction
        of the *protein* (distinct from a single domain's own HMM coverage).
    significant_evalue : float
        E-value threshold below which a hit counts as "significant" for
        ``n_significant_domains``. Hits with no evalue are always counted
        (never penalised for a provider that doesn't report one).

    Returns
    -------
    dict
        ``{transcript_id: {n_significant_domains, protein_coverage_fraction,
        n_nonoverlapping_domains, has_complete_domain_coverage,
        has_duplicate_domain, has_suspected_fusion, cross_provider_agreement,
        providers}}``. All future/TBC metrics the module docstring lists
        beyond these are intentionally not computed yet -- add them here
        (not in canonical_selection.py) as real domain data becomes
        available to validate them against.
    """
    by_transcript = defaultdict(list)
    for row in sidecar_rows:
        tid = row.get("transcript_id")
        if tid:
            by_transcript[tid].append(row)

    summaries = {}
    for tid, rows in by_transcript.items():
        n_significant = sum(
            1 for r in rows if r.get("evalue") is None or r["evalue"] <= significant_evalue
        )

        # Non-overlapping domains: greedy interval scheduling by (end asc)
        # among significant hits with known coordinates.
        intervals = sorted(
            (
                (r["start"], r["end"])
                for r in rows
                if r.get("start") is not None
                and r.get("end") is not None
                and (r.get("evalue") is None or r["evalue"] <= significant_evalue)
            ),
            key=lambda iv: iv[1],
        )
        nonoverlapping = []
        last_end = -1
        for s, e in intervals:
            if s > last_end:
                nonoverlapping.append((s, e))
                last_end = e
        n_nonoverlapping = len(nonoverlapping)

        covered_bp = sum(e - s + 1 for s, e in nonoverlapping)
        protein_len = (protein_lengths or {}).get(tid)
        coverage_fraction = (
            round(covered_bp / protein_len, 4) if protein_len and protein_len > 0 else None
        )

        accessions = [r["feature_accession"] for r in rows if r.get("feature_accession")]
        has_duplicate = len(accessions) != len(set(accessions))

        providers = sorted({r["provider"] for r in rows if r.get("provider")})
        # Cross-provider agreement: >=2 distinct providers both called a
        # domain at this transcript at all (not coordinate-matched here --
        # a real implementation would intersect coordinates per accession
        # family; left as a coarse presence check until real multi-provider
        # data is available to validate a finer one against).
        cross_provider_agreement = len(providers) >= 2

        summaries[tid] = {
            "n_significant_domains": n_significant,
            "protein_coverage_fraction": coverage_fraction,
            "n_nonoverlapping_domains": n_nonoverlapping,
            "has_complete_domain_coverage": (
                coverage_fraction is not None and coverage_fraction >= 0.9
            ),
            "has_duplicate_domain": has_duplicate,
            # Suspected fusion: markedly more non-overlapping domains than
            # is typical for a single well-formed protein. Deliberately a
            # crude, disabled-by-default heuristic (see
            # CanonicalDomainConfig.suspicious_fusion_penalty, TBC=0.0) --
            # flagging based on domain *count* alone risks false positives
            # for genuinely multi-domain proteins, so this is reported for
            # inspection, not wired into scoring in v1.
            "has_suspected_fusion": n_nonoverlapping >= 5,
            "cross_provider_agreement": cross_provider_agreement,
            "providers": ",".join(providers),
        }

    return summaries
