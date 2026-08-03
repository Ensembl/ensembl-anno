"""Backwards-compatible shim — re-exports the stable public API from the new package.

The implementation has moved to ``gmb.pipeline.longread.*``.  This module
exists only so that code that imports from ``gmb.pipeline.longread_consensus``
continues to work without changes.

Do not add new logic here.  Add it in the appropriate submodule instead.
"""

from gmb.pipeline.longread import (  # noqa: F401
    LongreadConsensusConfig,
    MultiExonRescueConfig,
    ShortreadRescueConfig,
    SingleExonRescueConfig,
    build_consensus_for_seqname,
    decide_keep_or_rescue,
    extract_junctions,
    junction_matches,
    load_longread_consensus_config,
    load_shortread_transcripts,
    load_split_manifest,
    m1_junction_complete,
    m2_full_chain,
    m3_dual_assembler,
    s1_reciprocal_boundary,
    s2_dual_assembler_boundary,
    s3_independent_evidence,
    split_by_seqname,
    try_rescue_multi_exon,
    try_rescue_single_exon,
    validate_config,
    validate_split_dir,
    write_consensus_gtf,
    write_dropped_records_tsv,
    write_rescue_attribution_tsv,
    write_run_manifest,
    write_summary,
)

# Aliases kept for any code that imported the old private-function names.
# These are the new public equivalents:
#   _m1_junction_complete  → m1_junction_complete
#   _m2_full_chain         → m2_full_chain
#   _m3_dual_assembler_chain → m3_dual_assembler
#   _s1_reciprocal_boundary  → s1_reciprocal_boundary
#   _s2_dual_assembler_boundary → s2_dual_assembler_boundary
#   _s3_independent_evidence    → s3_independent_evidence
#   _rescue_multi_exon_structured / _rescue_single_exon_structured → try_rescue_*
#   _keep_or_rescue             → decide_keep_or_rescue
#   _load_shortread_transcripts_for_seqname → load_shortread_transcripts
#   _junctions_from_exons (REMOVED — used coordinate snapping; replaced by extract_junctions)

_m1_junction_complete = m1_junction_complete
_m2_full_chain = m2_full_chain
_m3_dual_assembler_chain = m3_dual_assembler
_s1_reciprocal_boundary = s1_reciprocal_boundary
_s2_dual_assembler_boundary = s2_dual_assembler_boundary
_s3_independent_evidence = s3_independent_evidence
_rescue_multi_exon_structured = try_rescue_multi_exon
_rescue_single_exon_structured = try_rescue_single_exon
_keep_or_rescue = decide_keep_or_rescue
_load_shortread_transcripts_for_seqname = load_shortread_transcripts
