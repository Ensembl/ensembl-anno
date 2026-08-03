"""Configuration dataclasses, YAML loading, preset resolution, and validation.

All thresholds in ``LongreadConsensusConfig`` were derived from GCA_000002765.3
(P. falciparum) and are **not** organism-neutral defaults.  Use a named preset
(``configs/longread_consensus/``) or supply ``--config`` to override them.
"""

from __future__ import annotations

import dataclasses
import os
from dataclasses import dataclass, field
from typing import Optional

import yaml

# ---------------------------------------------------------------------------
# Legacy keys that were removed — fail loudly if present in YAML
# ---------------------------------------------------------------------------

_REMOVED_KEYS: dict[str, str] = {
    "require_shortread_support_for_single_read_models": (
        "Removed: unsafe overlap rescue caused multi-gene read bridges. "
        "Use shortread_rescue.enabled with policy M1/S1 instead."
    ),
    "shortread_overlap_fraction": (
        "Removed: belongs to the old overlap rescue path (see "
        "require_shortread_support_for_single_read_models)."
    ),
    "rescue_multi_exon_policy": ("Moved to shortread_rescue.multi_exon.policy"),
    "rescue_single_exon_policy": ("Moved to shortread_rescue.single_exon.policy"),
    "rescue_junction_tolerance_bp": ("Moved to shortread_rescue.multi_exon.junction_tolerance_bp"),
    "rescue_reciprocal_overlap_min": (
        "Moved to shortread_rescue.single_exon.reciprocal_overlap_min"
    ),
    "rescue_boundary_tolerance_bp": (
        "Moved to shortread_rescue.single_exon.boundary_tolerance_bp"
    ),
}

# ---------------------------------------------------------------------------
# Nested rescue configuration
# ---------------------------------------------------------------------------

_VALID_MULTI_EXON_POLICIES = frozenset({"M1", "M2", "M3"})
_VALID_SINGLE_EXON_POLICIES = frozenset({"S1", "S2", "S3"})


@dataclass
class MultiExonRescueConfig:
    """Configuration for rescuing otherwise-unsupported multi-exon single reads.

    M1: every intron junction in the candidate is confirmed by at least one
        short-read transcript (combined support across transcripts is allowed).
    M2: the candidate's complete ordered intron chain matches at least one
        single short-read transcript (same number of introns, within tolerance).
    M3: the candidate's junctions are independently confirmed by both Scallop
        AND StringTie (M1-level from each assembler separately).
    """

    enabled: bool = True
    policy: str = "M1"
    junction_tolerance_bp: int = 10
    require_both_assemblers: bool = False


@dataclass
class SingleExonRescueConfig:
    """Configuration for rescuing otherwise-unsupported single-exon single reads.

    S1: a single-exon short-read transcript reciprocally overlaps the candidate
        by at least ``reciprocal_overlap_min`` AND both boundaries are within
        ``boundary_tolerance_bp``.  Supporting SR transcript must also be
        single-exon (``require_single_exon_sr_model``).
    S2: both Scallop AND StringTie independently provide an S1-level match.
    S3: at least two independent single-exon short-read transcripts each satisfy S1.
    """

    enabled: bool = False
    policy: str = "S1"
    reciprocal_overlap_min: float = 0.80
    boundary_tolerance_bp: int = 20
    require_single_exon_sr_model: bool = True
    require_both_assemblers: bool = False


@dataclass
class ShortreadRescueConfig:
    """Top-level toggle for short-read-assisted rescue.

    When ``enabled`` is ``False`` (the default), no short-read files are
    loaded, no rescue is attempted, and the output is identical to running
    in pure long-read mode.

    When ``enabled`` is ``True``, at least one of ``multi_exon.enabled``
    or ``single_exon.enabled`` must also be ``True``, and at least one
    ``--scallop`` or ``--stringtie`` file must be supplied.
    """

    enabled: bool = False
    multi_exon: MultiExonRescueConfig = field(default_factory=MultiExonRescueConfig)
    single_exon: SingleExonRescueConfig = field(default_factory=SingleExonRescueConfig)


# ---------------------------------------------------------------------------
# Top-level configuration
# ---------------------------------------------------------------------------


@dataclass
class LongreadConsensusConfig:
    """Thresholds for raw-read → consensus collapsing.

    Defaults are grounded in GCA_000002765.3 (P. falciparum) biology.
    They are not appropriate for mammals, plants, fungi, or other protists
    without re-derivation.  Use a named preset or ``--config`` to override.
    """

    # Per-read rejection (chimera / genomic-DNA / mis-alignment screening)
    max_span_length: int = 45_000
    max_intron_length: int = 3_000

    # Locus clustering
    min_intergenic_gap: int = 500

    # Consensus grouping.  100 bp matches established long-read collapsing
    # practice for alignment-boundary noise.
    collapse_terminal_variation_bp: int = 100

    # Tolerance for snapping intron boundaries before grouping reads by
    # intron-chain signature.  Used for READ GROUPING only (not for rescue
    # junction matching, which uses explicit coordinate comparison).
    splice_site_tolerance_bp: int = 10

    # Minimum independent-read support to keep a consensus model
    min_read_support_multi_exon: int = 2
    min_read_support_single_exon: int = 4

    # Short-read-assisted rescue (disabled by default)
    shortread_rescue: ShortreadRescueConfig = field(default_factory=ShortreadRescueConfig)

    # Post-consensus containment suppression
    suppress_contained_models: bool = True


# ---------------------------------------------------------------------------
# Config loading helpers
# ---------------------------------------------------------------------------


def _load_multi_exon_rescue(raw: dict) -> MultiExonRescueConfig:
    known = {f.name for f in dataclasses.fields(MultiExonRescueConfig)}
    for key in raw:
        if key not in known:
            raise ValueError(f"Unknown shortread_rescue.multi_exon key: '{key}'")
    return MultiExonRescueConfig(**raw)


def _load_single_exon_rescue(raw: dict) -> SingleExonRescueConfig:
    known = {f.name for f in dataclasses.fields(SingleExonRescueConfig)}
    for key in raw:
        if key not in known:
            raise ValueError(f"Unknown shortread_rescue.single_exon key: '{key}'")
    return SingleExonRescueConfig(**raw)


def _load_shortread_rescue(raw: dict) -> ShortreadRescueConfig:
    known = {f.name for f in dataclasses.fields(ShortreadRescueConfig)}
    multi_raw = raw.pop("multi_exon", {}) or {}
    single_raw = raw.pop("single_exon", {}) or {}
    for key in raw:
        if key not in known:
            raise ValueError(f"Unknown shortread_rescue key: '{key}'")
    return ShortreadRescueConfig(
        enabled=raw.get("enabled", False),
        multi_exon=_load_multi_exon_rescue(multi_raw),
        single_exon=_load_single_exon_rescue(single_raw),
    )


def load_longread_consensus_config(
    path: Optional[str] = None,
    preset: Optional[str] = None,
) -> LongreadConsensusConfig:
    """Load a ``LongreadConsensusConfig`` from a YAML file and/or a preset.

    If both ``path`` and ``preset`` are given, the preset is loaded first and
    the explicit ``path`` is merged on top.

    Parameters
    ----------
    path : str or None
        Path to a YAML config file.  Unknown top-level keys raise ValueError.
    preset : str or None
        Name of a preset from ``configs/longread_consensus/``, e.g.
        ``"pfalciparum_pure"`` or ``"pfalciparum_assisted"``.
    """
    raw: dict = {}

    if preset is not None:
        preset_path = _resolve_preset(preset)
        with open(preset_path) as fh:
            raw = yaml.safe_load(fh) or {}

    if path is not None:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Config file not found: {path}")
        with open(path) as fh:
            override = yaml.safe_load(fh) or {}
        _deep_merge(raw, override)

    # Check for removed legacy keys
    for key in list(raw):
        if key in _REMOVED_KEYS:
            raise ValueError(
                f"Config key '{key}' has been removed.\n"
                f"  Reason: {_REMOVED_KEYS[key]}\n"
                f"  See docs/longread_consensus.md for the current config schema."
            )

    # Extract nested rescue block before building main config
    rescue_raw = raw.pop("shortread_rescue", {}) or {}

    # Validate top-level keys
    known_main = {f.name for f in dataclasses.fields(LongreadConsensusConfig)} - {
        "shortread_rescue"
    }
    for key in raw:
        if key not in known_main:
            raise ValueError(
                f"Unknown longread_consensus config key: '{key}'.  "
                f"Known keys: {sorted(known_main)}"
            )

    rescue_cfg = _load_shortread_rescue(dict(rescue_raw))
    return LongreadConsensusConfig(shortread_rescue=rescue_cfg, **raw)


def validate_config(
    cfg: LongreadConsensusConfig,
    scallop_paths: list[str],
    stringtie_paths: list[str],
) -> None:
    """Raise ``ValueError`` for invalid or inconsistent config.

    Called before any large-file processing so bad config fails fast.
    """
    _validate_thresholds(cfg)
    _validate_rescue(cfg, scallop_paths, stringtie_paths)


def _validate_thresholds(cfg: LongreadConsensusConfig) -> None:
    if cfg.max_span_length <= 0:
        raise ValueError(f"max_span_length must be > 0, got {cfg.max_span_length}")
    if cfg.max_intron_length <= 0:
        raise ValueError(f"max_intron_length must be > 0, got {cfg.max_intron_length}")
    if cfg.min_intergenic_gap < 0:
        raise ValueError(f"min_intergenic_gap must be >= 0, got {cfg.min_intergenic_gap}")
    if cfg.collapse_terminal_variation_bp < 0:
        raise ValueError(
            f"collapse_terminal_variation_bp must be >= 0, got "
            f"{cfg.collapse_terminal_variation_bp}"
        )
    if cfg.splice_site_tolerance_bp < 0:
        raise ValueError(
            f"splice_site_tolerance_bp must be >= 0, got {cfg.splice_site_tolerance_bp}"
        )
    if cfg.min_read_support_multi_exon < 1:
        raise ValueError(
            f"min_read_support_multi_exon must be >= 1, got " f"{cfg.min_read_support_multi_exon}"
        )
    if cfg.min_read_support_single_exon < 1:
        raise ValueError(
            f"min_read_support_single_exon must be >= 1, got "
            f"{cfg.min_read_support_single_exon}"
        )


def _validate_rescue(
    cfg: LongreadConsensusConfig,
    scallop_paths: list[str],
    stringtie_paths: list[str],
) -> None:
    rc = cfg.shortread_rescue
    if not rc.enabled:
        return

    me = rc.multi_exon
    se = rc.single_exon

    if not me.enabled and not se.enabled:
        raise ValueError(
            "shortread_rescue.enabled is true but neither multi_exon.enabled "
            "nor single_exon.enabled is true.  Enable at least one subtype."
        )

    if not scallop_paths and not stringtie_paths:
        raise ValueError(
            "shortread_rescue.enabled is true but no short-read files were "
            "supplied.  Provide at least one --scallop or --stringtie file."
        )

    if me.enabled:
        if me.policy not in _VALID_MULTI_EXON_POLICIES:
            raise ValueError(
                f"Unknown multi_exon rescue policy '{me.policy}'.  "
                f"Valid: {sorted(_VALID_MULTI_EXON_POLICIES)}"
            )
        if me.junction_tolerance_bp < 0:
            raise ValueError(
                f"shortread_rescue.multi_exon.junction_tolerance_bp must be "
                f">= 0, got {me.junction_tolerance_bp}"
            )
        if me.policy == "M3" or me.require_both_assemblers:
            if not scallop_paths or not stringtie_paths:
                raise ValueError(
                    f"Multi-exon policy '{me.policy}' (or require_both_assemblers) "
                    "requires both --scallop and --stringtie files."
                )

    if se.enabled:
        if se.policy not in _VALID_SINGLE_EXON_POLICIES:
            raise ValueError(
                f"Unknown single_exon rescue policy '{se.policy}'.  "
                f"Valid: {sorted(_VALID_SINGLE_EXON_POLICIES)}"
            )
        if not (0.0 < se.reciprocal_overlap_min <= 1.0):
            raise ValueError(
                f"shortread_rescue.single_exon.reciprocal_overlap_min must be "
                f"in (0, 1], got {se.reciprocal_overlap_min}"
            )
        if se.boundary_tolerance_bp < 0:
            raise ValueError(
                f"shortread_rescue.single_exon.boundary_tolerance_bp must be "
                f">= 0, got {se.boundary_tolerance_bp}"
            )
        if se.policy == "S2" or se.require_both_assemblers:
            if not scallop_paths or not stringtie_paths:
                raise ValueError(
                    f"Single-exon policy '{se.policy}' (or require_both_assemblers) "
                    "requires both --scallop and --stringtie files."
                )


# ---------------------------------------------------------------------------
# Preset resolution
# ---------------------------------------------------------------------------

_PRESETS_DIR = None  # resolved lazily — see _resolve_preset


def _get_presets_dir() -> str:
    """Return the longread_consensus presets directory, wheel-safe.

    Presets are packaged under ``gmb/configs/longread_consensus/``.  During
    uninstalled development they live in the source-tree
    ``configs/longread_consensus/`` directory three levels above this file.
    The installed location is checked first so that a wheel install is never
    shadowed by a stale source checkout.
    """
    # Installed location (works in editable installs and wheel installs)
    pkg_location = os.path.normpath(
        os.path.join(os.path.dirname(__file__), "..", "..", "configs", "longread_consensus")
    )
    if os.path.isdir(pkg_location):
        return pkg_location
    # Uninstalled source-tree fallback (three levels up from gmb/pipeline/longread/)
    src_location = os.path.normpath(
        os.path.join(os.path.dirname(__file__), "..", "..", "..", "configs", "longread_consensus")
    )
    return src_location


def _resolve_preset(name: str) -> str:
    """Resolve a preset name to an absolute YAML path."""
    presets_dir = _get_presets_dir()
    candidates = [
        os.path.join(presets_dir, f"{name}.yaml"),
        os.path.join(presets_dir, f"{name}"),
    ]
    for p in candidates:
        p = os.path.normpath(p)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        f"Preset '{name}' not found.  Looked in {os.path.normpath(presets_dir)}.  "
        f"Available: {_list_presets()}"
    )


def _list_presets() -> list[str]:
    presets_dir = os.path.normpath(_get_presets_dir())
    if not os.path.isdir(presets_dir):
        return []
    return sorted(os.path.splitext(f)[0] for f in os.listdir(presets_dir) if f.endswith(".yaml"))


def _deep_merge(base: dict, override: dict) -> None:
    """Merge ``override`` into ``base`` in-place.  Dicts deep-merge; other types replace."""
    for key, val in override.items():
        if key in base and isinstance(base[key], dict) and isinstance(val, dict):
            _deep_merge(base[key], val)
        else:
            base[key] = val
