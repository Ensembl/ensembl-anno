#!/usr/bin/env python3
"""Configuration system for Gene Model Builder.

Loads pipeline parameters from a YAML file with fungal defaults.
All clade- or species-specific behaviour is controlled via config,
never hard-coded in pipeline logic.

Usage:
    from config import load_config
    cfg = load_config()                   # fungal defaults
    cfg = load_config("my_config.yaml")   # custom overrides
"""

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import yaml

# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass
class OrfConfig:
    min_codons: int = 33
    allow_partial_5: bool = True
    allow_partial_3: bool = True
    allow_non_atg_start: bool = False
    stop_codon_char: str = "*"
    partial_prefix: str = ""


@dataclass
class ProteinFilterConfig:
    min_protein_aa: int = 30
    min_exon_count_for_short: int = 1
    redundancy_overlap: float = 0.80
    top_n_per_locus: int = 3
    max_span_bp: int = 50_000
    keep_secondary: bool = True
    min_alignment_coverage: float = 0.8
    min_percent_identity: float = 60.0
    min_bitscore: float = 50.0


@dataclass
class TranscriptomicFilterConfig:
    max_transcript_length: int = 20_000
    max_intron_length: int = 3_000
    min_intergenic_gap: int = 500
    allow_single_exon: bool = True
    strand_consistency_check: bool = False


@dataclass
class TranscriptSplittingConfig:
    split_enabled: bool = False
    split_gap_bp: int = 3_000
    split_on_contig_change: bool = True
    split_on_strand_change: bool = True
    split_on_large_exon_bp: Optional[int] = None
    max_segments_per_transcript: int = 50


@dataclass
class HelixerFilterConfig:
    min_cds_bp: int = 90
    max_exons: int = 50
    enabled: bool = True


@dataclass
class ScoringWeights:
    helixer: float = 2.0
    scallop: float = 1.0
    stringtie: float = 1.0


@dataclass
class ScoringConfig:
    weights: ScoringWeights = field(default_factory=ScoringWeights)
    protein_overlap_bonus: float = 2.0
    multi_source_bonus: float = 1.0
    noncanonical_splice_penalty: float = 0.5
    max_isoforms_per_locus: int = 2
    min_alternate_score: float = 3.0
    fungal_single_exon_mode: bool = True
    keep_helixer_without_support: bool = True
    require_protein_support_for_single_source: bool = False
    min_cds_bp: int = 150
    require_support_for_single_exon: bool = True


@dataclass
class ProteinValidationConfig:
    enabled: bool = False
    diamond_path: str = "diamond"
    psauron_path: str = "psauron"
    diamond_db: str = "swissprot.dmnd"
    psauron_model: str = "default"
    diamond_weight: float = 0.5
    psauron_weight: float = 0.5
    min_score: float = 0.5
    policy: str = "drop"  # 'drop' or 'penalize'


@dataclass
class QcConfig:
    max_transcripts_per_track: int = 5
    skip_orf_inference_tracks: List[str] = field(default_factory=lambda: ["OrthoDB", "UniProt"])
    parallel: bool = False
    workers: int = 4


@dataclass
class ValidationConfig:
    mode: str = "drop"  # "error" | "fix" | "drop_transcript"
    log_violations: bool = True
    max_feature_drift_bp: int = 1500
    feature_outside_exons_policy: str = "drop"  # "drop" | "trim" | "error"
    max_exon_len_bp: int = 15000

    # Percentile-based guardrails vs mega-exons / massive transcripts
    max_exon_len_mode: str = "fixed"  # "fixed" | "percentile"
    max_exon_len_percentile: float = 99.5
    max_exon_len_factor: float = 1.5
    max_exon_len_reference: str = "candidates_supported"  # "candidates_supported" | "reference"

    max_transcript_span_mode: str = "off"  # "off" | "fixed" | "percentile"
    max_transcript_span_bp: int = 100_000
    max_transcript_span_percentile: float = 99.5
    max_transcript_span_factor: float = 1.5


@dataclass
class UtrConfig:
    max_5p_bp: int = 2000
    max_3p_bp: int = 3000
    max_total_bp: int = 5000
    max_utr_to_cds_ratio: float = 5.0
    trim_policy: str = "hard_cap"  # "hard_cap" | "drop_utrs"

    # End support logic
    require_end_support: bool = True
    end_support_mode: str = "multisource_end_agreement"
    end_support_sources: List[str] = field(default_factory=lambda: ["Scallop", "StringTie"])
    end_tolerance_bp: int = 50
    require_multisource_for_utr_5p: bool = True
    require_multisource_for_utr_3p: bool = True
    fallback_policy_when_unsupported: str = "drop_utr"

    # Optional extensions
    min_protein_coding_score_for_utr: Optional[float] = None
    max_end_extension_bp: Optional[int] = None

    def __post_init__(self):
        valid_modes = {"multisource_end_agreement", "protein_validated", "either", "off"}
        if self.end_support_mode not in valid_modes:
            raise ValueError(f"Invalid end_support_mode: {self.end_support_mode}")

        valid_fallbacks = {"drop_utr", "hard_cap", "drop_transcript"}
        if self.fallback_policy_when_unsupported not in valid_fallbacks:
            raise ValueError(
                f"Invalid fallback_policy_when_unsupported: {self.fallback_policy_when_unsupported}"
            )


@dataclass
class DedupConfig:
    enabled: bool = True
    reciprocal_overlap_threshold: float = 0.80
    same_structure_tolerance_bp: int = 10
    policy: str = "merge_as_isoforms"  # | "keep_best_drop_rest" | "keep_both_if_opposite_strand"


@dataclass
class ExportConfig:
    write_cdna: bool = True
    write_protein: bool = True
    write_cds: bool = True
    include_partial: bool = True


@dataclass
class ReportingConfig:
    formats: List[str] = field(default_factory=lambda: ["json", "tsv"])


@dataclass
class PipelineConfig:
    preset: str = "fungi"
    orf: OrfConfig = field(default_factory=OrfConfig)
    protein_filter: ProteinFilterConfig = field(default_factory=ProteinFilterConfig)
    transcriptomic_filter: TranscriptomicFilterConfig = field(
        default_factory=TranscriptomicFilterConfig
    )
    transcript_splitting: TranscriptSplittingConfig = field(
        default_factory=TranscriptSplittingConfig
    )
    helixer_filter: HelixerFilterConfig = field(default_factory=HelixerFilterConfig)
    scoring: ScoringConfig = field(default_factory=ScoringConfig)
    protein_validation: ProteinValidationConfig = field(default_factory=ProteinValidationConfig)
    validation: ValidationConfig = field(default_factory=ValidationConfig)
    utr: UtrConfig = field(default_factory=UtrConfig)
    dedup: DedupConfig = field(default_factory=DedupConfig)
    qc: QcConfig = field(default_factory=QcConfig)
    export: ExportConfig = field(default_factory=ExportConfig)
    reporting: ReportingConfig = field(default_factory=ReportingConfig)


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------


def _update_dataclass(dc, d: dict, path_prefix: str = ""):
    """Recursively update a dataclass instance from a dict with strict merge rules.

    Rules:
      - dict -> deep merge
      - list -> replace entirely
      - scalar -> override
      - unknown keys -> raise ValueError
    """
    if d is None:
        return dc
    for key, val in d.items():
        if not hasattr(dc, key):
            raise ValueError(f"Unknown configuration key: '{path_prefix}{key}'")

        current = getattr(dc, key)
        if hasattr(current, "__dataclass_fields__") and isinstance(val, dict):
            _update_dataclass(current, val, path_prefix=f"{path_prefix}{key}.")
        else:
            setattr(dc, key, val)
    return dc


def _validate_dataclass(dc):
    """Recursively run __post_init__ validation on all dataclasses after manual updates."""
    if hasattr(dc, "__post_init__"):
        dc.__post_init__()
    if hasattr(dc, "__dataclass_fields__"):
        for field_name in dc.__dataclass_fields__:
            val = getattr(dc, field_name)
            if hasattr(val, "__dataclass_fields__"):
                _validate_dataclass(val)


def load_config(path: Optional[str] = None, preset: str = "fungi") -> PipelineConfig:
    """Load pipeline configuration.

    Parameters
    ----------
    path : str or None
        Path to a YAML config file.  If None, uses built-in defaults.
    preset : str
        Preset name ('fungi' uses configs/fungi_default.yaml).

    Returns
    -------
    PipelineConfig
    """
    cfg = PipelineConfig(preset=preset)

    # Base configuration based on preset
    if preset == "fungi":
        default_yaml = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "configs", "fungi_default.yaml"
        )
        if os.path.exists(default_yaml):
            with open(default_yaml) as fh:
                data = yaml.safe_load(fh) or {}
            _update_dataclass(cfg, data)
        else:
            raise FileNotFoundError(f"Missing default config preset: {default_yaml}")

    # User overrides
    if path is not None and os.path.exists(path):
        with open(path) as fh:
            data = yaml.safe_load(fh) or {}
        _update_dataclass(cfg, data)

    _validate_dataclass(cfg)
    return cfg
