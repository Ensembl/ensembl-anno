#!/usr/bin/env python3
"""Reference-free applicability gates for biological selection policies.

Some selection policies are only correct in a particular *evidence state*, not
for a particular clade. ``backbone_intron_rescue`` is the clearest case: it
exists to recover introns from an ab initio backbone that systematically
under-calls them, so it helps exactly when the backbone is structurally
under-resolving genes relative to the assembled transcript evidence, and harms
when the backbone is already good.

This module measures that evidence state **from the input evidence alone**. It
never reads a reference annotation, and it never looks at the organism name.

Everything here runs once per run, before the locus loop, and its decision is
recorded in the run manifest.

Calibration
-----------
The thresholds below were measured on the three genomes GMB has been validated
against, before any reference was consulted for the decision:

===============  ==========  ==================  ===================  =======  ==============
genome           backbone    backbone multi-exon  assembled multi-exon  ratio    rescue outcome
===============  ==========  ==================  ===================  =======  ==============
P. falciparum    Tiberius              17.3%               98.7%       5.71x   68.6% CDS-exact (helps)
T. gondii        Tiberius              38.7%               95.4%       2.47x   validated, helps
Z. tritici       Helixer               71.9%               92.1%       1.27x    1.8% CDS-exact (harms)
===============  ==========  ==================  ===================  =======  ==============

Both criteria separate the two helpful cases from the harmful one with a clear
margin, and they agree with each other on all three genomes -- the decision does
not rest on a single knife-edge number.

This is a three-genome calibration, not a survey. ``auto`` is deliberately
conservative: when the evidence is ambiguous or too thin to judge, it returns
OFF, because OFF is the safe state (it reproduces baseline behaviour).
"""

from __future__ import annotations

from dataclasses import dataclass, field

from gmb.pipeline.canonical_evidence import (
    EVIDENCE_CLASS_BACKBONE,
    EVIDENCE_CLASS_SHORT_READ,
    EvidenceRoles,
)

# --- calibrated thresholds -------------------------------------------------
# A backbone this spliced is not collapsing introns and needs no rescue.
# Observed: 0.173 / 0.387 (rescue helps) vs 0.719 (rescue harms).
BACKBONE_MULTI_EXON_MAX = 0.55
# Assembled transcripts must be this many times more spliced than the backbone.
# Observed: 5.71 / 2.47 (rescue helps) vs 1.27 (rescue harms).
ASSEMBLED_TO_BACKBONE_RATIO_MIN = 1.5
# An assembled track whose introns are mostly non-canonical is not a credible
# source of replacement splice structures. This is the Minimap2 lesson applied
# as a precondition: that track was 16.4% canonical.
ASSEMBLED_CANONICAL_SPLICE_MIN = 0.90
# Below this many models on either side the fractions are too noisy to act on.
MIN_MODELS_FOR_DECISION = 200

RESCUE_MODES = ("off", "on", "auto")


@dataclass
class BackboneResolutionStats:
    """Structural resolution of the backbone vs the assembled transcript tracks."""

    backbone_models: int = 0
    backbone_multi_exon: int = 0
    assembled_models: int = 0
    assembled_multi_exon: int = 0
    assembled_canonical_splice_fraction: float | None = None
    backbone_sources: list = field(default_factory=list)
    assembled_sources: list = field(default_factory=list)

    @property
    def backbone_multi_exon_fraction(self) -> float | None:
        if not self.backbone_models:
            return None
        return self.backbone_multi_exon / self.backbone_models

    @property
    def assembled_multi_exon_fraction(self) -> float | None:
        if not self.assembled_models:
            return None
        return self.assembled_multi_exon / self.assembled_models

    @property
    def assembled_to_backbone_ratio(self) -> float | None:
        b, a = self.backbone_multi_exon_fraction, self.assembled_multi_exon_fraction
        if not b or a is None:
            return None
        return a / b

    def as_dict(self) -> dict:
        def r(v, n=4):
            return None if v is None else round(v, n)
        return {
            "backbone_sources": list(self.backbone_sources),
            "backbone_models": self.backbone_models,
            "backbone_multi_exon": self.backbone_multi_exon,
            "backbone_multi_exon_fraction": r(self.backbone_multi_exon_fraction),
            "assembled_sources": list(self.assembled_sources),
            "assembled_models": self.assembled_models,
            "assembled_multi_exon": self.assembled_multi_exon,
            "assembled_multi_exon_fraction": r(self.assembled_multi_exon_fraction),
            "assembled_to_backbone_ratio": r(self.assembled_to_backbone_ratio),
            "assembled_canonical_splice_fraction": r(
                self.assembled_canonical_splice_fraction),
        }


@dataclass
class RescueDecision:
    """Resolved `backbone_intron_rescue` state for one run."""

    mode: str                 # what the config asked for: off | on | auto
    enabled: bool             # what will actually happen
    reason: str               # human-readable justification
    stats: BackboneResolutionStats | None = None

    def as_dict(self) -> dict:
        d = {"mode": self.mode, "enabled": self.enabled, "reason": self.reason}
        if self.stats is not None:
            d["evidence"] = self.stats.as_dict()
        return d


def normalise_rescue_mode(value) -> str:
    """Map a configured value onto one of RESCUE_MODES.

    Accepts the current strings and the legacy booleans:
      True  -> "on"   (what `backbone_intron_rescue: true` always meant)
      False -> "off"
    """
    if isinstance(value, bool):
        return "on" if value else "off"
    if value is None:
        return "off"
    text = str(value).strip().lower()
    if text in ("true", "yes"):
        return "on"
    if text in ("false", "no", "none", ""):
        return "off"
    if text not in RESCUE_MODES:
        raise ValueError(
            f"scoring.backbone_intron_rescue must be one of {list(RESCUE_MODES)} "
            f"(or a legacy boolean), got {value!r}."
        )
    return text


def measure_backbone_resolution(
    exon_df,
    scoring_config,
    canonical_splice_fraction: float | None = None,
) -> BackboneResolutionStats:
    """Measure backbone vs assembled-transcript structural resolution.

    Parameters
    ----------
    exon_df : pandas.DataFrame
        Candidate exon rows, with at least ``Source`` and ``transcript_id``.
        One row per exon, so exon count per transcript is a simple group size.
    scoring_config : ScoringConfig
        Supplies the evidence-role labels.
    canonical_splice_fraction : float or None
        Canonical GT-AG fraction of the assembled tracks, if already measured
        (preflight computes it). ``None`` means "not measured", which makes the
        auto gate refuse to fire.
    """
    roles = EvidenceRoles.from_config(scoring_config)
    stats = BackboneResolutionStats(
        assembled_canonical_splice_fraction=canonical_splice_fraction)
    if exon_df is None or len(exon_df) == 0:
        return stats

    exons_per_tx = exon_df.groupby(["Source", "transcript_id"], observed=True).size()
    for (source, _tid), n_exons in exons_per_tx.items():
        role = roles.role_of(source)
        if role == EVIDENCE_CLASS_BACKBONE:
            stats.backbone_models += 1
            if n_exons > 1:
                stats.backbone_multi_exon += 1
            if source not in stats.backbone_sources:
                stats.backbone_sources.append(source)
        elif role == EVIDENCE_CLASS_SHORT_READ:
            stats.assembled_models += 1
            if n_exons > 1:
                stats.assembled_multi_exon += 1
            if source not in stats.assembled_sources:
                stats.assembled_sources.append(source)
    return stats


def resolve_backbone_intron_rescue(
    scoring_config,
    exon_df=None,
    canonical_splice_fraction: float | None = None,
) -> RescueDecision:
    """Decide whether backbone intron rescue runs, and say why.

    ``off``  -- never fires. The safe default.
    ``on``   -- expert override; fires regardless of the measured evidence state.
    ``auto`` -- fires only when the measured evidence shows a backbone that is
                under-resolving introns relative to credible assembled
                transcripts. Refuses (returns OFF) whenever the evidence is
                missing, too thin, or ambiguous.
    """
    mode = normalise_rescue_mode(getattr(scoring_config, "backbone_intron_rescue", "off"))

    if mode == "off":
        return RescueDecision(mode, False, "disabled by configuration")
    if mode == "on":
        return RescueDecision(
            mode, True,
            "explicitly enabled by configuration (applicability gate bypassed)")

    stats = measure_backbone_resolution(
        exon_df, scoring_config, canonical_splice_fraction)
    b = stats.backbone_multi_exon_fraction
    a = stats.assembled_multi_exon_fraction
    ratio = stats.assembled_to_backbone_ratio

    if b is None or a is None:
        return RescueDecision(
            mode, False,
            "auto: no backbone and/or assembled-transcript evidence to measure; "
            "defaulting to off", stats)
    if (stats.backbone_models < MIN_MODELS_FOR_DECISION
            or stats.assembled_models < MIN_MODELS_FOR_DECISION):
        return RescueDecision(
            mode, False,
            f"auto: too few models to judge "
            f"(backbone {stats.backbone_models}, assembled {stats.assembled_models}; "
            f"need >= {MIN_MODELS_FOR_DECISION} of each); defaulting to off", stats)
    if canonical_splice_fraction is None:
        return RescueDecision(
            mode, False,
            "auto: assembled-transcript splice quality was not measured, so the "
            "replacement structures cannot be trusted; defaulting to off", stats)
    if canonical_splice_fraction < ASSEMBLED_CANONICAL_SPLICE_MIN:
        return RescueDecision(
            mode, False,
            f"auto: assembled transcripts are only "
            f"{canonical_splice_fraction:.1%} canonically spliced "
            f"(need >= {ASSEMBLED_CANONICAL_SPLICE_MIN:.0%}), so they are not a "
            f"credible source of replacement introns; defaulting to off", stats)
    if b > BACKBONE_MULTI_EXON_MAX:
        return RescueDecision(
            mode, False,
            f"auto: backbone is {b:.1%} multi-exon "
            f"(> {BACKBONE_MULTI_EXON_MAX:.0%}), so it is not under-calling "
            f"introns and needs no rescue", stats)
    if ratio is None or ratio < ASSEMBLED_TO_BACKBONE_RATIO_MIN:
        return RescueDecision(
            mode, False,
            f"auto: assembled transcripts are only {ratio:.2f}x more spliced than "
            f"the backbone (need >= {ASSEMBLED_TO_BACKBONE_RATIO_MIN}x); the "
            f"backbone is not clearly under-resolving genes", stats)
    return RescueDecision(
        mode, True,
        f"auto: backbone is {b:.1%} multi-exon (<= {BACKBONE_MULTI_EXON_MAX:.0%}) "
        f"while assembled transcripts are {a:.1%} multi-exon ({ratio:.2f}x more) "
        f"at {canonical_splice_fraction:.1%} canonical splice sites -- the backbone "
        f"is under-resolving introns and credible replacements exist", stats)
