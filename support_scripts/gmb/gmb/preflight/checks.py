#!/usr/bin/env python3
"""Preflight checks on a GMB evidence bundle.

See :mod:`gmb.preflight` for why this exists.

Verdicts
--------
``PASS``  nothing to act on.
``WARN``  usable, but the operator should know. Does not stop a build.
``FAIL``  the build should not proceed as configured. ``gmb-preflight`` exits
          non-zero; the operator may override with ``--allow-fail``.

Thresholds live in ``config.preflight`` so they can be tuned per clade without
touching code, and so the values that judged a run are captured in
``resolved_config.yaml`` alongside everything else.
"""

from __future__ import annotations

import gzip
import hashlib
import os
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field

from gmb.pipeline.canonical_evidence import (
    EVIDENCE_CLASS_BACKBONE,
    EVIDENCE_CLASS_LONG_READ,
    EVIDENCE_CLASS_PROTEIN_ALIGNMENT,
    EVIDENCE_CLASS_SHORT_READ,
    EvidenceRoles,
)

PASS, WARN, FAIL = "PASS", "WARN", "FAIL"
_SEVERITY = {PASS: 0, WARN: 1, FAIL: 2}

_CANONICAL_SITES = {"GTAG", "GCAG", "ATAC"}
_COMPLEMENT = str.maketrans("ACGTacgtnN", "TGCAtgcaNN")


def _rc(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _open(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


@dataclass
class CheckResult:
    """One named check against one input (or the bundle as a whole)."""

    name: str
    verdict: str
    message: str
    target: str = ""
    detail: dict = field(default_factory=dict)


@dataclass
class TrackSummary:
    """Structural summary of one evidence track."""

    label: str
    path: str
    role: str
    weight: float | None = None
    models: int = 0
    multi_exon: int = 0
    single_exon: int = 0
    seqids: list = field(default_factory=list)
    strand_counts: dict = field(default_factory=dict)
    feature_counts: dict = field(default_factory=dict)
    duplicate_ids_across_seqids: int = 0
    canonical_splice_fraction: float | None = None
    introns_examined: int = 0
    median_transcript_span: int | None = None
    max_transcript_span: int | None = None
    max_intron: int | None = None
    sha256: str = ""

    @property
    def multi_exon_fraction(self) -> float | None:
        return (self.multi_exon / self.models) if self.models else None

    def to_dict(self) -> dict:
        d = asdict(self)
        d["multi_exon_fraction"] = (
            round(self.multi_exon_fraction, 4)
            if self.multi_exon_fraction is not None else None
        )
        return d


@dataclass
class PreflightReport:
    """Everything preflight learned about one evidence bundle."""

    checks: list = field(default_factory=list)
    tracks: list = field(default_factory=list)
    genome: dict = field(default_factory=dict)
    policy: dict = field(default_factory=dict)
    recommendations: list = field(default_factory=list)

    @property
    def verdict(self) -> str:
        if not self.checks:
            return PASS
        return max((c.verdict for c in self.checks), key=lambda v: _SEVERITY[v])

    def counts(self) -> dict:
        c = Counter(x.verdict for x in self.checks)
        return {PASS: c[PASS], WARN: c[WARN], FAIL: c[FAIL]}

    def to_dict(self) -> dict:
        return {
            "verdict": self.verdict,
            "counts": self.counts(),
            "genome": self.genome,
            "tracks": [t.to_dict() for t in self.tracks],
            "policy": self.policy,
            "recommendations": self.recommendations,
            "checks": [asdict(c) for c in self.checks],
        }

    def to_text(self) -> str:
        lines = [f"GMB PREFLIGHT — {self.verdict}", "=" * 64, ""]
        g = self.genome
        if g:
            lines += [f"Genome: {g.get('path','?')}",
                      f"  {g.get('seqids',0)} sequence(s), {g.get('total_bp',0):,} bp", ""]
        if self.tracks:
            lines.append("Evidence tracks")
            lines.append(f"  {'label':<12}{'role':<28}{'wt':>5}{'models':>9}"
                         f"{'multi%':>9}{'canon%':>9}")
            for t in self.tracks:
                mf = t.multi_exon_fraction
                cs = t.canonical_splice_fraction
                lines.append(
                    f"  {t.label:<12}{t.role:<28}"
                    f"{('-' if t.weight is None else f'{t.weight:g}'):>5}"
                    f"{t.models:>9,}"
                    f"{('-' if mf is None else f'{100*mf:.1f}'):>9}"
                    f"{('-' if cs is None else f'{100*cs:.1f}'):>9}")
            lines.append("")
        if self.policy:
            lines.append("Selection policy as configured")
            for k, v in self.policy.items():
                lines.append(f"  {k:<34}{v}")
            lines.append("")
        by_verdict = defaultdict(list)
        for c in self.checks:
            by_verdict[c.verdict].append(c)
        for verdict in (FAIL, WARN, PASS):
            items = by_verdict.get(verdict, [])
            if not items:
                continue
            lines.append(f"{verdict} ({len(items)})")
            for c in items:
                tgt = f" [{c.target}]" if c.target else ""
                lines.append(f"  - {c.name}{tgt}: {c.message}")
            lines.append("")
        if self.recommendations:
            lines.append("Recommendations")
            for r in self.recommendations:
                lines.append(f"  * {r}")
            lines.append("")
        counts = self.counts()
        lines.append(f"{counts[PASS]} passed, {counts[WARN]} warning(s), "
                     f"{counts[FAIL]} failure(s)")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# parsing
# ---------------------------------------------------------------------------

def _parse_attributes(field9: str, gff3: bool) -> dict:
    out = {}
    if gff3:
        for kv in field9.split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                out[k.strip()] = v.strip()
    else:
        for kv in field9.split(";"):
            kv = kv.strip()
            if not kv:
                continue
            parts = kv.split(" ", 1)
            if len(parts) == 2:
                out[parts[0]] = parts[1].strip().strip('"')
    return out


def _transcript_key(attrs: dict, gff3: bool) -> str | None:
    """Group exons into transcripts.

    GFF3 exons name their transcript in ``Parent``; GTF exons use
    ``transcript_id``. Preferring the wrong one silently turns every exon into
    its own single-exon "transcript", which would make a well-spliced backbone
    look completely unspliced.
    """
    if gff3:
        for k in ("Parent", "transcript_id", "ID"):
            if k in attrs:
                return attrs[k].split(",")[0]
        return None
    for k in ("transcript_id", "Parent", "ID"):
        if k in attrs:
            return attrs[k]
    return None


def summarise_track(label: str, path: str, role: str, genome: dict | None = None,
                    weight: float | None = None) -> TrackSummary:
    """Structural summary of one GTF/GFF3 evidence track."""
    gff3 = ".gff" in os.path.basename(path).lower()
    summary = TrackSummary(label=label, path=path, role=role, weight=weight)
    exons: dict = defaultdict(list)
    seqids_per_tx: dict = defaultdict(set)
    strand_counts: Counter = Counter()
    feature_counts: Counter = Counter()
    seqids = set()

    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            feature_counts[f[2]] += 1
            if f[2] not in ("exon", "CDS"):
                continue
            attrs = _parse_attributes(f[8], gff3)
            tid = _transcript_key(attrs, gff3)
            if tid is None:
                continue
            seqids.add(f[0])
            strand_counts[f[6]] += 1
            seqids_per_tx[tid].add(f[0])
            exons[(tid, f[2])].append((f[0], int(f[3]), int(f[4]), f[6]))

    # Prefer exon rows; protein-alignment tracks often carry CDS only.
    use = "exon" if any(k[1] == "exon" for k in exons) else "CDS"
    models = {tid: segs for (tid, ft), segs in exons.items() if ft == use}

    summary.models = len(models)
    summary.multi_exon = sum(1 for s in models.values() if len(s) > 1)
    summary.single_exon = summary.models - summary.multi_exon
    summary.seqids = sorted(seqids)
    summary.strand_counts = dict(strand_counts)
    summary.feature_counts = dict(feature_counts.most_common(8))
    summary.duplicate_ids_across_seqids = sum(
        1 for v in seqids_per_tx.values() if len(v) > 1)
    summary.sha256 = _sha256(path)

    spans, max_intron = [], 0
    canonical = total_introns = 0
    for segs in models.values():
        ordered = sorted((s, e) for _c, s, e, _st in segs)
        spans.append(ordered[-1][1] - ordered[0][0])
        for i in range(len(ordered) - 1):
            max_intron = max(max_intron, ordered[i + 1][0] - ordered[i][1] - 1)
        if genome and len(ordered) > 1:
            chrom, strand = segs[0][0], segs[0][3]
            seq = genome.get(chrom)
            if seq is None:
                continue
            for i in range(len(ordered) - 1):
                donor = seq[ordered[i][1]:ordered[i][1] + 2].upper()
                acceptor = seq[ordered[i + 1][0] - 3:ordered[i + 1][0] - 1].upper()
                site = (donor + acceptor) if strand == "+" else (_rc(acceptor) + _rc(donor))
                total_introns += 1
                if site in _CANONICAL_SITES:
                    canonical += 1
    if spans:
        spans.sort()
        summary.median_transcript_span = spans[len(spans) // 2]
        summary.max_transcript_span = spans[-1]
    summary.max_intron = max_intron or None
    summary.introns_examined = total_introns
    if total_introns:
        summary.canonical_splice_fraction = canonical / total_introns
    return summary


def load_genome_lengths(path: str) -> dict:
    """Sequence name -> sequence, read once and reused by every splice check."""
    genome: dict = {}
    name, buf = None, []
    with _open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    genome[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        genome[name] = "".join(buf)
    return genome


# ---------------------------------------------------------------------------
# the checks
# ---------------------------------------------------------------------------

def _splice_thresholds(pcfg, role: str) -> tuple[float, float]:
    """(warn_below, fail_below) canonical-splice fractions for one role.

    Roles are judged differently on purpose. A long-read track is *claiming* to
    observe splice structure directly, so a poor canonical fraction there is a
    hard failure. A protein-alignment track is support-only and never supplies a
    candidate structure, so the same number is merely informational.
    """
    table = {
        EVIDENCE_CLASS_BACKBONE: (pcfg.backbone_splice_warn, pcfg.backbone_splice_fail),
        EVIDENCE_CLASS_SHORT_READ: (pcfg.shortread_splice_warn, pcfg.shortread_splice_fail),
        EVIDENCE_CLASS_LONG_READ: (pcfg.longread_splice_warn, pcfg.longread_splice_fail),
    }
    return table.get(role, (0.0, 0.0))


def run_preflight(inputs: dict, config, genome_seqs: dict | None = None) -> PreflightReport:
    """Check an evidence bundle before building.

    Parameters
    ----------
    inputs : dict
        ``{"genome": path, "tracks": [{"label":…, "path":…}, …]}``. Roles are
        resolved from ``config.scoring``, not from the dict.
    config : PipelineConfig
    genome_seqs : dict or None
        Pre-loaded genome. Loaded from ``inputs["genome"]`` when omitted.
    """
    from gmb.pipeline.applicability import (
        normalise_rescue_mode,
        resolve_backbone_intron_rescue,
    )
    from gmb.pipeline.scoring import weights_for_role

    pcfg = config.preflight
    scfg = config.scoring
    roles = EvidenceRoles.from_config(scfg)
    report = PreflightReport()

    def add(name, verdict, message, target="", **detail):
        report.checks.append(CheckResult(name, verdict, message, target, detail))

    # ---- genome -----------------------------------------------------------
    genome_path = inputs.get("genome")
    if not genome_path or not os.path.exists(genome_path):
        add("genome_readable", FAIL, f"genome FASTA not found: {genome_path!r}")
        return report
    if genome_seqs is None:
        genome_seqs = load_genome_lengths(genome_path)
    if not genome_seqs:
        add("genome_readable", FAIL, "genome FASTA contains no sequences", genome_path)
        return report
    report.genome = {
        "path": genome_path,
        "seqids": len(genome_seqs),
        "total_bp": sum(len(s) for s in genome_seqs.values()),
        "sha256": _sha256(genome_path),
    }
    add("genome_readable", PASS,
        f"{len(genome_seqs)} sequence(s), "
        f"{report.genome['total_bp']:,} bp", genome_path)
    genome_ids = set(genome_seqs)

    # ---- per-track --------------------------------------------------------
    for spec in inputs.get("tracks", []):
        label, path = spec["label"], spec.get("path")
        if not path:
            continue
        if not os.path.exists(path):
            add("file_readable", FAIL, f"file not found: {path}", label)
            continue
        role = roles.role_of(label)
        weight = weights_for_role(scfg.weights, role)
        track = summarise_track(label, path, role, genome_seqs, weight)
        report.tracks.append(track)
        structural = role in (EVIDENCE_CLASS_BACKBONE, EVIDENCE_CLASS_SHORT_READ,
                              EVIDENCE_CLASS_LONG_READ)

        add("file_readable", PASS, f"{track.models:,} model(s)", label)

        # role resolution
        if role == "other" or role is None:
            add("role_resolved", FAIL,
                f"source '{label}' resolves to no configured evidence role, so it "
                f"would receive the generic unknown-source weight "
                f"({scfg.weights.unknown}). Add it to scoring.backbone_label / "
                f"shortread_labels / longread_label / protein_alignment_labels.",
                label, role=role)
        else:
            add("role_resolved", PASS, f"role={role}, weight={weight:g}", label,
                role=role, weight=weight)

        # seqid compatibility
        unknown = [s for s in track.seqids if s not in genome_ids]
        if unknown:
            frac = len(unknown) / max(len(track.seqids), 1)
            verdict = FAIL if frac > pcfg.max_unknown_seqid_fraction else WARN
            add("seqid_compatibility", verdict,
                f"{len(unknown)}/{len(track.seqids)} sequence name(s) are absent "
                f"from the genome FASTA (e.g. {unknown[:3]}). Remap upstream — "
                f"models on these sequences cannot be used.",
                label, unknown_seqids=unknown[:20])
        else:
            add("seqid_compatibility", PASS,
                f"all {len(track.seqids)} sequence name(s) present in the genome",
                label)

        # cross-seqid ID collisions
        if track.duplicate_ids_across_seqids:
            add("cross_seqid_id_collisions", WARN,
                f"{track.duplicate_ids_across_seqids} transcript ID(s) are reused "
                f"on more than one sequence. GMB namespaces these automatically, "
                f"but it means the upstream tool numbers models per sequence.",
                label, count=track.duplicate_ids_across_seqids)
        else:
            add("cross_seqid_id_collisions", PASS, "no reused transcript IDs", label)

        # strand completeness
        unstranded = sum(v for k, v in track.strand_counts.items() if k not in ("+", "-"))
        total_rows = sum(track.strand_counts.values()) or 1
        if unstranded:
            frac = unstranded / total_rows
            verdict = WARN if frac <= pcfg.max_unstranded_fraction else FAIL
            add("strand_completeness", verdict,
                f"{unstranded:,}/{total_rows:,} feature row(s) ({frac:.1%}) have no "
                f"strand; these are excluded from selection.", label)
        else:
            add("strand_completeness", PASS, "every feature is stranded", label)

        # splice quality — the check this module exists for
        if structural:
            cs = track.canonical_splice_fraction
            warn_below, fail_below = _splice_thresholds(pcfg, role)
            if cs is None or track.introns_examined < pcfg.min_introns_for_splice_check:
                add("splice_quality", WARN,
                    f"only {track.introns_examined} intron(s) available — too few to "
                    f"judge splice quality. A track with no introns contributes no "
                    f"splice structure.", label,
                    introns=track.introns_examined)
            elif cs < fail_below:
                add("splice_quality", FAIL,
                    f"only {cs:.1%} of {track.introns_examined:,} intron(s) use a "
                    f"canonical splice site (GT-AG/GC-AG/AT-AC), below the "
                    f"{fail_below:.0%} floor for role '{role}'. This track carries "
                    f"little usable splice information — check the upstream "
                    f"alignment/assembly before building on it.",
                    label, canonical_splice_fraction=round(cs, 4))
            elif cs < warn_below:
                add("splice_quality", WARN,
                    f"{cs:.1%} of {track.introns_examined:,} intron(s) are canonical, "
                    f"below the {warn_below:.0%} expectation for role '{role}'.",
                    label, canonical_splice_fraction=round(cs, 4))
            else:
                add("splice_quality", PASS,
                    f"{cs:.1%} canonical over {track.introns_examined:,} intron(s)",
                    label, canonical_splice_fraction=round(cs, 4))

        # implausible spans
        if track.max_transcript_span and track.max_transcript_span > pcfg.max_transcript_span_warn_bp:
            add("transcript_spans", WARN,
                f"longest transcript spans {track.max_transcript_span:,} bp "
                f"(median {track.median_transcript_span:,} bp). Very long models are "
                f"often chimeras; GMB's transcript_splitting settings apply.", label)
        if track.max_intron and track.max_intron > pcfg.max_intron_warn_bp:
            add("intron_lengths", WARN,
                f"longest intron is {track.max_intron:,} bp; check "
                f"transcriptomic_filter.max_intron_length for this clade.", label)

    # ---- bundle-level -----------------------------------------------------
    by_role = defaultdict(list)
    for t in report.tracks:
        by_role[t.role].append(t)

    if not by_role.get(EVIDENCE_CLASS_BACKBONE):
        add("backbone_present", WARN,
            "no ab initio backbone supplied. GMB will run, but the standard "
            "production workflow expects one.")
    elif len(by_role[EVIDENCE_CLASS_BACKBONE]) > 1:
        add("backbone_present", FAIL,
            f"{len(by_role[EVIDENCE_CLASS_BACKBONE])} tracks resolve to the "
            f"backbone role; exactly one is supported.")
    else:
        add("backbone_present", PASS,
            f"one backbone track ({by_role[EVIDENCE_CLASS_BACKBONE][0].label})")

    if not by_role.get(EVIDENCE_CLASS_SHORT_READ):
        add("transcript_evidence_present", WARN,
            "no assembled transcript evidence. Selection will be driven almost "
            "entirely by the backbone.")
    else:
        add("transcript_evidence_present", PASS,
            f"{len(by_role[EVIDENCE_CLASS_SHORT_READ])} assembled transcript track(s)")

    if not by_role.get(EVIDENCE_CLASS_PROTEIN_ALIGNMENT):
        add("protein_evidence_present", WARN,
            "no protein alignments. Protein support gates retention in several "
            "policies; without it those gates cannot fire.")
    else:
        add("protein_evidence_present", PASS,
            f"{len(by_role[EVIDENCE_CLASS_PROTEIN_ALIGNMENT])} protein alignment track(s)")

    # ---- resolved policy --------------------------------------------------
    lr_tracks = by_role.get(EVIDENCE_CLASS_LONG_READ, [])
    lr_disposition = str(getattr(scfg, "longread_disposition", "primary_structural"))
    rescue_mode = normalise_rescue_mode(getattr(scfg, "backbone_intron_rescue", "off"))

    shortread_cs = [t.canonical_splice_fraction for t in by_role.get(EVIDENCE_CLASS_SHORT_READ, [])
                    if t.canonical_splice_fraction is not None]
    assembled_cs = (sum(shortread_cs) / len(shortread_cs)) if shortread_cs else None

    decision = resolve_backbone_intron_rescue(
        scfg, _exon_frame_from_tracks(report.tracks), assembled_cs)

    report.policy = {
        "structural_corroboration": bool(getattr(scfg, "structural_corroboration", False)),
        "protein_support_mode": getattr(scfg, "protein_support_mode", "positional"),
        "longread_structural_guard": bool(getattr(scfg, "longread_structural_guard", False)),
        "longread_disposition": lr_disposition,
        "backbone_intron_rescue_mode": rescue_mode,
        "backbone_intron_rescue_will_fire": decision.enabled,
        "backbone_intron_rescue_reason": decision.reason,
    }
    if decision.stats is not None:
        report.policy["backbone_resolution"] = decision.stats.as_dict()

    add("backbone_intron_rescue", PASS,
        f"mode '{rescue_mode}' -> {'WILL FIRE' if decision.enabled else 'will not fire'}: "
        f"{decision.reason}")

    if not lr_tracks:
        add("longread_policy", PASS,
            "no long-read track supplied; the long-read guard cannot fire and no "
            "long-read special case applies.")
    else:
        lr_fail = [c for c in report.checks
                   if c.name == "splice_quality" and c.verdict == FAIL
                   and c.target in {t.label for t in lr_tracks}]
        if lr_fail and lr_disposition == "primary_structural":
            add("longread_policy", FAIL,
                f"long-read track failed its splice-quality check but "
                f"scoring.longread_disposition is '{lr_disposition}', so it would "
                f"supply primary splice structure. Set it to 'support_only' or "
                f"'reject', or fix the upstream alignment.")
            report.recommendations.append(
                "Set scoring.longread_disposition: support_only (or reject) — the "
                "long-read track's splice structure is not trustworthy. Whichever "
                "you choose is recorded in the run manifest.")
        elif lr_fail:
            add("longread_policy", PASS,
                f"long-read track failed its splice-quality check and "
                f"longread_disposition is already '{lr_disposition}'.")
        else:
            add("longread_policy", PASS,
                f"{len(lr_tracks)} long-read track(s), disposition '{lr_disposition}'")

    # ---- recommendations --------------------------------------------------
    if rescue_mode == "off" and decision.mode == "off":
        stats = _rescue_stats_preview(report.tracks, scfg, assembled_cs)
        if stats is not None and stats.get("would_fire"):
            report.recommendations.append(
                "scoring.backbone_intron_rescue is 'off', but the evidence state "
                "matches the one where it has been validated to help "
                f"(backbone {stats['backbone_multi_exon_fraction']:.1%} multi-exon vs "
                f"assembled {stats['assembled_multi_exon_fraction']:.1%}). Consider "
                "'auto'.")
    for t in report.tracks:
        if t.role in ("other", None):
            report.recommendations.append(
                f"Give '{t.label}' an evidence role in the config; it currently "
                f"falls back to the unknown-source weight.")
    return report


def _exon_frame_from_tracks(tracks):
    """Build the minimal frame the applicability gate needs from track summaries.

    Avoids re-reading every file: the gate only needs per-source multi-exon
    counts, which `summarise_track` already computed.
    """
    import pandas as pd

    rows = []
    for t in tracks:
        for i in range(t.multi_exon):
            rows.append((t.label, f"{t.label}_m{i}", 0))
            rows.append((t.label, f"{t.label}_m{i}", 1))
        for i in range(t.single_exon):
            rows.append((t.label, f"{t.label}_s{i}", 0))
    if not rows:
        return None
    return pd.DataFrame(rows, columns=["Source", "transcript_id", "_exon"])


def _rescue_stats_preview(tracks, scoring_config, assembled_cs):
    """What `auto` would decide, for a run configured with rescue off."""
    from gmb.pipeline.applicability import measure_backbone_resolution

    frame = _exon_frame_from_tracks(tracks)
    if frame is None:
        return None
    stats = measure_backbone_resolution(frame, scoring_config, assembled_cs)
    import copy

    probe = copy.copy(scoring_config)
    probe.backbone_intron_rescue = "auto"
    from gmb.pipeline.applicability import resolve_backbone_intron_rescue

    decision = resolve_backbone_intron_rescue(probe, frame, assembled_cs)
    d = stats.as_dict()
    d["would_fire"] = decision.enabled
    return d
