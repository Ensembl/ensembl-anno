#!/usr/bin/env python3
"""FASTA QC — legacy wrapper. All logic now in :mod:`gmb.pipeline.fasta_qc`."""
import os, sys  # noqa: E401
_pkg_root = os.path.dirname(os.path.abspath(__file__))
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)
from gmb.pipeline.fasta_qc import (  # noqa: E402, F401
    load_genome,
    main,
    parse_fasta_ids,
    parse_fasta_records,
    parse_gff3,
    print_report,
    run_coverage_checks,
    run_sequence_checks,
    validate_fasta,
)

if __name__ == "__main__":
    main()
