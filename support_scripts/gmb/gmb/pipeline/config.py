#!/usr/bin/env python3
"""Configuration system for Gene Model Builder.

Loads pipeline parameters in layers:

    standard.yaml (always-on organism-neutral base)
        → optional clade/workflow preset  (--preset fungi | apicomplexa | ...)
            → user-supplied --config files (in order, each later file wins)

Usage:
    from gmb.pipeline.config import load_config, list_build_presets
    cfg = load_config()                         # standard + fungi (default)
    cfg = load_config(preset="apicomplexa")     # standard + apicomplexa
    cfg = load_config(preset=None)              # standard only
    cfg = load_config("my.yaml")                # standard + fungi + my.yaml
    cfg = load_config(["a.yaml", "b.yaml"])     # layered overrides (b wins)
    list_build_presets()                        # ["apicomplexa", "fungi", ...]
"""

import os
import warnings
from dataclasses import dataclass, field
from importlib.resources import files as _pkg_files
from typing import Optional, Union

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
class BackboneFilterConfig:
    """Plausibility filter for the ab initio backbone track (Helixer or
    Tiberius, whichever was loaded -- see PipelineConfig.scoring.backbone_label).
    Named generically since this section is never Helixer-specific despite
    its legacy `helixer_filter` YAML key (still accepted, see
    _DEPRECATED_KEY_ALIASES below).
    """

    min_cds_bp: int = 90
    max_exons: int = 50
    enabled: bool = True


# Deprecated alias for the pre-rename class name -- see _DEPRECATED_KEY_ALIASES
# for the corresponding YAML key alias (helixer_filter -> backbone_filter).
# Remove once no code/config plausibly still imports this name directly.
HelixerFilterConfig = BackboneFilterConfig


@dataclass
class ScoringWeights:
    # backbone: weight for whichever ab initio backbone is actually loaded
    # (Helixer or Tiberius -- matched against scoring.backbone_label, not a
    # fixed source name). Legacy YAML key `helixer` is still accepted (see
    # _DEPRECATED_KEY_ALIASES).
    backbone: float = 2.0
    scallop: float = 1.0
    stringtie: float = 1.0
    minimap2: float = 1.0


@dataclass
class ScoringConfig:
    weights: ScoringWeights = field(default_factory=ScoringWeights)
    protein_overlap_bonus: float = 2.0
    multi_source_bonus: float = 1.0
    noncanonical_splice_penalty: float = 0.5
    max_isoforms_per_locus: int = 2
    min_alternate_score: float = 3.0
    fungal_single_exon_mode: bool = True
    # Legacy YAML key `keep_helixer_without_support` is still accepted (see
    # _DEPRECATED_KEY_ALIASES) -- this has never been Helixer-specific: it
    # protects whichever backbone (Helixer or Tiberius) is actually loaded.
    keep_backbone_without_support: bool = True
    require_protein_support_for_single_source: bool = False
    min_cds_bp: int = 150
    require_support_for_single_exon: bool = True
    same_gene_overlap_threshold: float = 0.15
    # Source label of the ab initio backbone track actually loaded (set
    # automatically by the CLI based on --helixer vs --tiberius). Determines
    # which string in `weights` and `keep_backbone_without_support` applies —
    # not necessarily "Helixer".
    backbone_label: str = "Helixer"


@dataclass
class ProteinValidationConfig:
    enabled: bool = False
    diamond_path: str = "diamond"
    psauron_path: str = "psauron"
    diamond_db: Optional[str] = None
    # Psauron uses one bundled model checkpoint and exposes no model-
    # selection option; -m specifies minimum protein length, not a model.
    psauron_min_length: int = 5  # psauron -m/--minimum-length (aa), psauron's own default
    psauron_use_cpu: bool = False  # psauron -c/--use-cpu; False lets it auto-detect a GPU
    diamond_weight: float = 0.5
    psauron_weight: float = 0.5
    min_score: float = 0.5
    policy: str = "drop"  # 'drop' or 'penalize' (also accepts 'penalise')
    diamond_min_query_coverage: float = 0.0  # 0-100; additional gate before counting a hit
    diamond_min_target_coverage: float = 0.0  # 0-100; additional gate before counting a hit

    def __post_init__(self) -> None:
        if self.enabled and self.diamond_weight > 0 and self.diamond_db is None:
            raise ValueError(
                "protein_validation.diamond_db must be set when enabled=true and "
                "diamond_weight > 0.  Set diamond_db to the path of your DIAMOND "
                "database (e.g. swissprot.dmnd)."
            )


@dataclass
class QcConfig:
    max_transcripts_per_track: int = 5
    skip_orf_inference_tracks: list[str] = field(default_factory=lambda: ["OrthoDB", "UniProt"])
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
    end_support_sources: list[str] = field(default_factory=lambda: ["Scallop", "StringTie"])
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
class DuplicateTranscriptCollapseConfig:
    """Config for gmb.pipeline.duplicate_transcript_collapse.

    Distinct from DedupConfig: dedup_genes.py merges whole *genes* by
    overlap + a single (first-mRNA, tolerance-bp) structural check -- it is
    what brings exact-duplicate fragments together as isoforms of one gene
    in the first place (see that module's docstring for the root cause:
    exon-level PyRanges clustering can split one transcript's own exons
    across separate clusters when no other evidence bridges the intron
    gap, and each fragment is independently re-admitted by the
    keep_backbone_without_support single-exon exception). This section
    controls a stricter, *exact*-equality-only pass that runs after that,
    within each gene's assembled isoform set, so it only ever removes
    transcripts proven structurally and translationally identical -- never
    genuinely distinct isoforms.
    """

    collapse_exact_duplicates: bool = True
    # When True (default), two transcripts with identical CDS but different
    # UTR extents are kept as separate isoforms rather than collapsed --
    # genuinely different supported UTRs are real evidence, not redundancy.
    preserve_distinct_utrs: bool = True
    # When True (default), two transcripts are only ever collapsed if their
    # translated protein sequences also match exactly (belt-and-braces on
    # top of CDS-coordinate equality -- catches the case where identical
    # CDS coordinates on inconsistent reference sequence data could still
    # translate differently, though that should not happen in practice).
    preserve_distinct_proteins: bool = True
    # CDS phase/frame is compared when both records have it recorded;
    # missing phase information on both sides is not itself a mismatch
    # (see module docstring for why this is the safest available rule
    # given the current data model).
    require_matching_cds_phase: bool = True


@dataclass
class ExportConfig:
    write_cdna: bool = True
    write_protein: bool = True
    write_cds: bool = True
    include_partial: bool = True


@dataclass
class ReportingConfig:
    formats: list[str] = field(default_factory=lambda: ["json", "tsv"])


@dataclass
class CanonicalProteinValidationWeights:
    """Protein-plausibility component weights (gmb.pipeline.canonical_selection).

    psauron/diamond_identity/diamond_balanced_coverage are each already
    bounded 0-1 (psauron score; DIAMOND percentages divided by 100 before
    weighting) -- deliberately NOT including raw DIAMOND bitscore, which
    scales with protein length and would let long proteins dominate
    regardless of match quality. diamond_hit_bonus is a flat, bounded bonus
    for "a significant hit exists at all", independent of its strength.

    diamond_balanced_coverage_weight replaced two separate
    diamond_query_coverage_weight/diamond_target_coverage_weight terms: the
    coverage component now uses gmb.pipeline.canonical_evidence.balanced_coverage
    (min(qcov, scov)/100), not the two sides scored and weighted
    independently. A hit covering 100% of a short query but 10% of a long
    target is a fragment match; weighting each side separately let the
    query-coverage term flatter it even though the target side already told
    the real story. The default (0.30) equals the sum of the two weights it
    replaces, so the total achievable protein-validation_subtotal is
    unchanged from before this change.
    """

    psauron_weight: float = 0.35
    diamond_identity_weight: float = 0.15
    diamond_balanced_coverage_weight: float = 0.30
    diamond_hit_bonus: float = 0.2


@dataclass
class CanonicalEvidenceWeights:
    """Annotation-evidence component weights.

    gmb_score is min-max normalised across the gene's own isoforms before
    weighting (it is open-ended/additive in gmb.pipeline.scoring, not
    naturally bounded); the *_support flags are 0/1 (does the transcript's
    evidence-class set contain a class of that kind -- see
    gmb.pipeline.canonical_evidence). independent_source_weight scales with
    the transcript's *evidence-class* breadth (not a raw named-source
    count -- Scallop+StringTie is one short_read_transcriptomic class, not
    two independent sources), capped at independent_source_cap classes, so
    a transcript backed by many redundant tracks of the same underlying
    evidence type cannot out-rank one backed by several genuinely different
    types. The field name is kept as `independent_source_weight` for
    backward-compatible config keys even though it now scores class breadth
    -- see `n_evidence_classes` in the transcript-ranking output for the
    raw count this weight is derived from, and `n_independent_sources` for
    the (still separately reported, no longer scored here) raw named-source
    count.
    """

    gmb_score_weight: float = 0.30
    independent_source_weight: float = 0.15
    independent_source_cap: int = 3
    transcriptomic_support_weight: float = 0.10
    protein_alignment_support_weight: float = 0.10
    longread_support_weight: float = 0.10
    backbone_support_weight: float = 0.10


@dataclass
class CanonicalStructureWeights:
    complete_orf_bonus: float = 0.30
    internal_stop_penalty: float = 0.50
    partial_cds_penalty: float = 0.20


@dataclass
class CanonicalDomainConfig:
    """Protein-domain/feature evidence (Pfam, InterPro, ...).

    Disabled by default: no domain adapter has run for this first pass (see
    gmb.pipeline.domain_evidence). All weights are TBC placeholders --
    plumbed through end-to-end (config -> scorer -> report) but inert at 0
    until real domain data and a considered weighting are available. Domain
    evidence must stay *supportive*, never a requirement: a transcript with
    zero domain hits (e.g. a genuinely lineage-specific or poorly-annotated
    protein) must not be penalised just for having none -- only
    fragmented/suspicious *positive* domain signals are ever penalised here.
    """

    enabled: bool = False
    domain_support_weight: float = 0.0  # TBC
    domain_coverage_weight: float = 0.0  # TBC
    complete_domain_bonus: float = 0.0  # TBC
    fragmented_domain_penalty: float = 0.0  # TBC
    suspicious_fusion_penalty: float = 0.0  # TBC
    cross_provider_agreement_bonus: float = 0.0  # TBC


@dataclass
class InterProNextflowExecConfig:
    """Fully parameterised, portable settings for running InterProScan 6 via
    Nextflow (Mode A -- ``run_interproscan: true``, GMB launches the scan
    itself). Every field is a placeholder the operator supplies for their
    own environment; NONE of this is hardcoded anywhere in GMB Python code
    (no Docker/Slurm/Singularity/Apptainer image name, executor, or
    site-specific path is ever assumed) -- see docs/interpro_resolver.md
    and docs/interproscan6_cluster.example.config for how these map onto a
    laptop Docker run vs a cluster Slurm+Singularity/Apptainer run.

    Irrelevant, and never read, when ``run_interproscan: false`` (Mode B --
    GMB only consumes an already-completed InterProScan run's output files,
    produced any way the operator likes).
    """

    nextflow_executable: str = "nextflow"
    workflow: str = "ebi-pf-team/interproscan6"
    revision: str = "6.0.1"
    # e.g. "docker" (laptop) or "slurm,singularity" / "slurm,apptainer" (cluster).
    profile: str = "docker"
    # -c <file>: an external Nextflow config supplying site-specific values
    # (queue/partition, container cache dir, bind mounts, ...). See
    # docs/interproscan6_cluster.example.config for a fully-placeholdered
    # example; None runs with no extra config (e.g. plain local Docker).
    config_file: Optional[str] = None
    # --datadir: shared InterPro member-database data directory. Required
    # for a real run; left None here since its location is entirely
    # site-specific and must never be assumed.
    data_dir: Optional[str] = None
    # --interpro: the InterPro data release to scan against (e.g. "109.0"),
    # or "latest".
    interpro_release: str = "latest"
    # --appl: comma-separated application/member-database list; None scans
    # with the workflow's own default application set.
    applications: Optional[str] = None
    work_dir: Optional[str] = None  # -work-dir
    output_dir: Optional[str] = None  # --outdir; None defaults to
    # <canonical-selection output-dir>/interpro_run/results (see
    # gmb.pipeline.interpro_review.run_interproscan_workflow).
    cpus: Optional[int] = None
    max_workers: Optional[int] = None
    # Any additional raw Nextflow CLI arguments the operator's environment
    # requires, appended verbatim and unvalidated -- an intentional escape
    # hatch so this dataclass never needs to anticipate every possible
    # execution profile.
    extra_args: list[str] = field(default_factory=list)


@dataclass
class InterProResolverConfig:
    """Optional second-stage InterProScan resolver for ambiguous canonical
    choices, nested under CanonicalSelectionConfig.interpro_resolver.

    InterProScan is a SECOND-STAGE RESOLVER only: it is never part of a
    normal GMB run, never required for completion, and never triggers a
    second DIAMOND/psauron pass (see gmb.pipeline.protein_validation's
    compute-once contract, which this section never touches). This section
    only decides which small subset of genes is worth the cost of a domain
    scan -- the initial canonical selection itself still happens without it
    (see gmb.pipeline.canonical_selection), and `final = initial` whenever
    this resolver is disabled or a review does not lead to a safeguarded
    replacement.

    ``enabled`` is the master switch and defaults to **False**: InterPro
    involvement is opt-in. ``run_interproscan``/``apply_replacements`` are
    independent sub-flags:
      * ``run_interproscan=False`` (default): Mode B -- GMB only consumes
        an already-completed InterProScan run's output files (any
        execution profile). This is the safer default for a first pass on
        a new dataset: review candidates and manifest are prepared, but
        nothing is launched automatically.
      * ``run_interproscan=True``: Mode A -- GMB launches InterProScan 6
        itself via Nextflow, fully parameterised by ``nextflow`` below.
      * ``apply_replacements=True`` (default, but only takes effect when
        ``enabled=True``): the resolver may replace the initial canonical
        transcript when its conservative replacement policy and safeguards
        (see gmb.pipeline.interpro_resolver.decide_canonical) are met.
        ``apply_replacements=False`` keeps the resolver a report-only
        layer (verdicts/recommendations recorded, canonical outputs
        untouched) even while ``enabled=True`` -- useful for evaluating the
        evidence model on a new dataset before trusting it to replace
        anything.
    """

    enabled: bool = False
    run_interproscan: bool = False
    apply_replacements: bool = True
    # Never review a gene with fewer than this many isoforms -- a
    # single-transcript gene has no choice to disambiguate.
    min_isoforms: int = 2
    # Only the top-N ranked candidates per reviewed gene are submitted;
    # sending every isoform would defeat the point of a targeted resolver.
    max_candidates_per_gene: int = 2
    # Winner/runner-up total_score gap at or below which the choice counts
    # as "close" (an ambiguity trigger). Deliberately looser than
    # CanonicalSelectionConfig.low_confidence_score_gap, since this decides
    # what is worth *reviewing*, not what is reported as low confidence.
    close_call_score_gap: float = 0.10
    # Individual ambiguity triggers, each independently switchable so a run
    # can be scoped (e.g. review only genuinely conflicting evidence).
    trigger_low_confidence: bool = True
    trigger_close_score: bool = True
    trigger_diamond_psauron_disagree: bool = True
    trigger_length_or_lexical_only: bool = True
    trigger_internal_stop_conflict: bool = True
    # A candidate this many times longer than its rival, without stronger
    # protein support, suggests possible fusion or a run-on model.
    length_ratio_flag: float = 1.5
    # --- Replacement-policy safeguards (only relevant when
    # apply_replacements=True; see interpro_resolver.decide_canonical) ---
    # Representative-coverage delta below which a coverage-based
    # improvement is NOT considered material enough to replace the initial
    # canonical on its own (categorical cases -- AntiFam avoidance, a
    # complete domain appearing where there was none -- bypass this
    # threshold, since they are not coverage-magnitude judgements).
    min_coverage_delta_for_replacement: float = 0.10
    # Never replace with a candidate that is structurally worse (e.g. an
    # incomplete ORF where the initial choice was complete).
    require_no_worse_structural_class: bool = True
    # Never replace with a candidate that introduces an internal stop the
    # initial choice did not have.
    require_no_new_internal_stops: bool = True
    # Never replace when the initial choice's own DIAMOND/psauron evidence
    # is materially stronger than the replacement candidate's -- InterPro
    # architecture evidence must not silently override much stronger
    # protein-validation evidence already on record. Expressed as a
    # balanced-DIAMOND-coverage-or-psauron-score gap threshold.
    max_protein_validation_override_gap: float = 0.15
    nextflow: InterProNextflowExecConfig = field(default_factory=InterProNextflowExecConfig)


@dataclass
class CanonicalSelectionConfig:
    """Config for gmb.pipeline.canonical_selection (a standalone post-build
    reporting step -- see that module's docstring for why it reads
    consensus.gff3/evidence_attribution.tsv/protein_validation.tsv rather
    than running inside gmb.cli.build itself).
    """

    enabled: bool = True
    protein_validation: CanonicalProteinValidationWeights = field(
        default_factory=CanonicalProteinValidationWeights
    )
    evidence: CanonicalEvidenceWeights = field(default_factory=CanonicalEvidenceWeights)
    structure: CanonicalStructureWeights = field(default_factory=CanonicalStructureWeights)
    domains: CanonicalDomainConfig = field(default_factory=CanonicalDomainConfig)
    # Winner/runner-up total-score gap below which a selection is flagged
    # "low_confidence" in the report (in addition to the separate
    # LOW_CONFIDENCE_NO_PROTEIN_SUPPORT reason code for missing evidence).
    low_confidence_score_gap: float = 0.05
    interpro_resolver: InterProResolverConfig = field(default_factory=InterProResolverConfig)


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
    # Legacy YAML section name `helixer_filter` is still accepted (see
    # _DEPRECATED_KEY_ALIASES).
    backbone_filter: BackboneFilterConfig = field(default_factory=BackboneFilterConfig)
    scoring: ScoringConfig = field(default_factory=ScoringConfig)
    protein_validation: ProteinValidationConfig = field(default_factory=ProteinValidationConfig)
    validation: ValidationConfig = field(default_factory=ValidationConfig)
    utr: UtrConfig = field(default_factory=UtrConfig)
    dedup: DedupConfig = field(default_factory=DedupConfig)
    duplicate_transcript_collapse: DuplicateTranscriptCollapseConfig = field(
        default_factory=DuplicateTranscriptCollapseConfig
    )
    qc: QcConfig = field(default_factory=QcConfig)
    export: ExportConfig = field(default_factory=ExportConfig)
    reporting: ReportingConfig = field(default_factory=ReportingConfig)
    # InterPro second-stage resolver config lives at
    # canonical_selection.interpro_resolver (see InterProResolverConfig) --
    # nested there, not top-level, since it is conceptually part of
    # canonical selection's optional second stage, not a sibling pipeline
    # stage of its own.
    canonical_selection: CanonicalSelectionConfig = field(default_factory=CanonicalSelectionConfig)


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

# Preset names that have been renamed; load_config emits DeprecationWarning
# and remaps automatically so existing --preset arguments keep working.
_DEPRECATED_PRESET_NAMES: dict[str, str] = {
    "fungi_default": "fungi",
}


def list_build_presets() -> list[str]:
    """Return available clade preset names for gmb-build.

    Returns every .yaml file in the ``gmb.configs`` package directory
    whose name is not ``standard`` (which is the internal always-on base,
    not a user-selectable clade preset).  Returns ``[]`` if the package
    data is not installed.
    """
    try:
        pkg = _pkg_files("gmb.configs")
        return sorted(
            p.name.removesuffix(".yaml")
            for p in pkg.iterdir()
            if p.name.endswith(".yaml") and p.name != "standard.yaml"
        )
    except (FileNotFoundError, NotADirectoryError):
        return []


def _load_bundled_config_yaml(name: str) -> dict:
    """Load a bundled YAML from gmb.configs by base name (no extension).

    Raises FileNotFoundError when the file is missing from the installed
    package -- callers translate this into a user-visible error.
    """
    try:
        text = _pkg_files("gmb.configs").joinpath(f"{name}.yaml").read_text(encoding="utf-8")
        return yaml.safe_load(text) or {}
    except FileNotFoundError:
        raise FileNotFoundError(
            f"Bundled config '{name}.yaml' not found in gmb.configs package.  "
            "Run 'pip install -e .' to ensure package data is installed."
        ) from None


# Deprecated config-key aliases, scoped by the *containing* dataclass's type
# name (so an old name is only ever remapped within the section it actually
# belonged to). Each entry exists because a "generalise backbone naming"
# pass renamed a Helixer-specific-sounding key to its accurate generic name
# without breaking configs already using the old one. Remove an entry (and
# its DeprecationWarning) once no config plausibly still sets the old key --
# not before, and not silently.
_DEPRECATED_KEY_ALIASES: dict = {
    "PipelineConfig": {"helixer_filter": "backbone_filter"},
    "ScoringConfig": {"keep_helixer_without_support": "keep_backbone_without_support"},
    "ScoringWeights": {"helixer": "backbone"},
}


def _apply_deprecated_key_aliases(dc, d: dict, path_prefix: str) -> dict:
    """Return `d` with any legacy keys (scoped to `dc`'s type) renamed to
    their current name, each emitting one DeprecationWarning.

    If both the legacy key and its current name are set in the same dict,
    the current name always wins -- the legacy value is dropped (never
    silently combined) with a warning explaining why.
    """
    aliases = _DEPRECATED_KEY_ALIASES.get(type(dc).__name__)
    if not aliases or not any(old in d for old in aliases):
        return d
    resolved = dict(d)
    for old_key, new_key in aliases.items():
        if old_key not in resolved:
            continue
        old_val = resolved.pop(old_key)
        if new_key in resolved:
            warnings.warn(
                f"Config key '{path_prefix}{old_key}' is deprecated and was "
                f"ignored because '{path_prefix}{new_key}' was also set -- "
                "the current key always takes precedence over the legacy one.",
                DeprecationWarning,
                stacklevel=5,
            )
        else:
            warnings.warn(
                f"Config key '{path_prefix}{old_key}' is deprecated; use "
                f"'{path_prefix}{new_key}' instead. This alias will be "
                "removed in a future release.",
                DeprecationWarning,
                stacklevel=5,
            )
            resolved[new_key] = old_val
    return resolved


def _update_dataclass(dc, d: dict, path_prefix: str = ""):
    """Recursively update a dataclass instance from a dict with strict merge rules.

    Rules:
      - dict -> deep merge
      - list -> replace entirely
      - scalar -> override
      - unknown keys -> raise ValueError, unless a documented deprecated
        alias exists for this dataclass (see _DEPRECATED_KEY_ALIASES)
    """
    if d is None:
        return dc
    d = _apply_deprecated_key_aliases(dc, d, path_prefix)
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


def load_config(
    path: Optional[Union[str, list]] = None,
    preset: Optional[str] = "fungi",
) -> "PipelineConfig":
    """Load pipeline configuration in layers.

    The effective configuration is assembled as:

        standard.yaml  (always-on organism-neutral base)
            → clade/workflow preset YAML  (when preset is not None)
                → user-supplied override files  (in order, last wins)

    Parameters
    ----------
    path : str, list of str, or None
        One or more YAML override files applied on top of the preset.
        Each file deep-merges its dicts and replaces list-valued keys.
        A single string is treated as a one-element list for backward
        compatibility.  Missing files raise ``FileNotFoundError``.
    preset : str or None
        Clade preset name.  Available presets are returned by
        ``list_build_presets()``.  ``None`` (or ``"none"``) loads
        ``standard.yaml`` only, with no clade overlay.  Defaults to
        ``"fungi"`` for backward compatibility.

        The deprecated name ``"fungi_default"`` is accepted with a
        ``DeprecationWarning`` and silently remapped to ``"fungi"``.

    Returns
    -------
    PipelineConfig
    """
    # Normalise preset: None / "none" / "" all mean "standard only".
    if preset in (None, "none", ""):
        preset = None

    # Remap deprecated preset names.
    if preset in _DEPRECATED_PRESET_NAMES:
        new_name = _DEPRECATED_PRESET_NAMES[preset]
        warnings.warn(
            f"Preset '{preset}' is deprecated; use '{new_name}' instead.  "
            "This alias will be removed in a future release.",
            DeprecationWarning,
            stacklevel=2,
        )
        preset = new_name

    cfg = PipelineConfig(preset=preset or "none")

    # Layer 1 — standard.yaml: organism-neutral base (always loaded).
    _update_dataclass(cfg, _load_bundled_config_yaml("standard"))

    # Layer 2 — clade preset: organism-specific overrides (optional).
    if preset is not None:
        available = list_build_presets()
        if preset not in available:
            raise FileNotFoundError(
                f"Unknown preset '{preset}'.  "
                f"Available presets: {available or ['none installed']}.\n"
                "Use --list-presets to see what is installed."
            )
        _update_dataclass(cfg, _load_bundled_config_yaml(preset))

    # Layer 3 — user overrides: applied in order (later files win on shared keys).
    paths = [path] if isinstance(path, str) else (path or [])
    for override_path in paths:
        if override_path is None:
            continue
        if not os.path.exists(override_path):
            raise FileNotFoundError(f"Config file not found: {override_path}")
        with open(override_path) as fh:
            data = yaml.safe_load(fh) or {}
        _update_dataclass(cfg, data)

    _validate_dataclass(cfg)
    return cfg


def dump_config(cfg: "PipelineConfig") -> dict:
    """Serialise a PipelineConfig to a plain dict (for YAML output).

    Uses ``dataclasses.asdict`` so nested dataclasses are recursively
    converted.  The result can be round-tripped through ``yaml.safe_dump``
    and then ``load_config`` (as a user override file), provided the
    ``preset`` top-level key is stripped first.
    """
    from dataclasses import asdict

    return asdict(cfg)
