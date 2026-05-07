#!/usr/bin/env python3
"""CDS & UTR annotation — legacy wrapper. All logic now in :mod:`gmb.pipeline.annotate_cds_utrs`."""
import os, sys  # noqa: E401
_pkg_root = os.path.dirname(os.path.abspath(__file__))
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)
from gmb.pipeline.annotate_cds_utrs import (  # noqa: E402, F401
    CODON_TABLE,
    COMPLEMENT,
    START_CODON,
    STOP_CODONS,
    annotate_all_transcripts,
    annotate_transcript,
    annotations_to_gff_rows,
    build_spliced_seq,
    check_frame_continuity,
    check_splice_sites,
    derive_utrs,
    find_best_orf,
    get_start_stop_positions,
    load_genome,
    main,
    map_cds_to_genomic,
    reverse_complement,
    translate,
    write_augmented_gff3,
    _empty_result,
    _make_orf_label,
)

if __name__ == "__main__":
    main()
