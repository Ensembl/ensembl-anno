#!/usr/bin/env python3
"""GFF3 validation — legacy wrapper. All logic now in :mod:`gmb.pipeline.gff3_validate`."""
import os, sys  # noqa: E401
_pkg_root = os.path.dirname(os.path.abspath(__file__))
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)
from gmb.pipeline.gff3_validate import (  # noqa: E402, F401
    fix_gene,
    fix_transcript,
    trim_utrs,
    validate_and_fix_gff3,
    validate_gene,
    validate_transcript,
)
