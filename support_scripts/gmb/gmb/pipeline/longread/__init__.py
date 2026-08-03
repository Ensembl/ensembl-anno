"""Long-read consensus transcript collapsing pipeline.

Public API (stable across refactors):

    from gmb.pipeline.longread import (
        LongreadConsensusConfig,
        ShortreadRescueConfig,
        MultiExonRescueConfig,
        SingleExonRescueConfig,
        load_longread_consensus_config,
        validate_config,
        build_consensus_for_seqname,
        split_by_seqname,
        write_consensus_gtf,
        write_rescue_attribution_tsv,
        write_dropped_records_tsv,
    )

The backwards-compatible import path ``gmb.pipeline.longread_consensus`` is
retained as a shim.
"""

from .config import (
    LongreadConsensusConfig,
    MultiExonRescueConfig,
    ShortreadRescueConfig,
    SingleExonRescueConfig,
    list_presets,
    load_longread_consensus_config,
    validate_config,
)
from .consensus import build_consensus_for_seqname
from .io import (
    load_shortread_transcripts,
    load_split_manifest,
    split_by_seqname,
    validate_split_dir,
    write_consensus_gtf,
)
from .reporting import (
    write_dropped_records_tsv,
    write_rescue_attribution_tsv,
    write_run_manifest,
    write_summary,
)
from .rescue import (
    decide_keep_or_rescue,
    extract_junctions,
    junction_matches,
    m1_junction_complete,
    m2_full_chain,
    m3_dual_assembler,
    s1_reciprocal_boundary,
    s2_dual_assembler_boundary,
    s3_independent_evidence,
    try_rescue_multi_exon,
    try_rescue_single_exon,
)

__all__ = [
    # config
    "LongreadConsensusConfig",
    "ShortreadRescueConfig",
    "MultiExonRescueConfig",
    "SingleExonRescueConfig",
    "list_presets",
    "load_longread_consensus_config",
    "validate_config",
    # consensus
    "build_consensus_for_seqname",
    # io
    "split_by_seqname",
    "load_shortread_transcripts",
    "load_split_manifest",
    "validate_split_dir",
    "write_consensus_gtf",
    # rescue
    "extract_junctions",
    "junction_matches",
    "m1_junction_complete",
    "m2_full_chain",
    "m3_dual_assembler",
    "s1_reciprocal_boundary",
    "s2_dual_assembler_boundary",
    "s3_independent_evidence",
    "try_rescue_multi_exon",
    "try_rescue_single_exon",
    "decide_keep_or_rescue",
    # reporting
    "write_rescue_attribution_tsv",
    "write_dropped_records_tsv",
    "write_run_manifest",
    "write_summary",
]
