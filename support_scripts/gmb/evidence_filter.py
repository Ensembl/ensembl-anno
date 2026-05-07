#!/usr/bin/env python3
"""Evidence filter — legacy wrapper. All logic now in :mod:`gmb.pipeline.evidence_filter`."""
import os, sys  # noqa: E401
_pkg_root = os.path.dirname(os.path.abspath(__file__))
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)
from gmb.pipeline.evidence_filter import (  # noqa: E402, F401
    filter_chimeras,
    filter_helixer_models,
    filter_protein_evidence,
    split_mega_transcripts,
)
