#!/usr/bin/env python3
"""Protein validation — legacy wrapper. All logic now in :mod:`gmb.pipeline.protein_validation`."""
import os, sys  # noqa: E401
_pkg_root = os.path.dirname(os.path.abspath(__file__))
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)
from gmb.pipeline.protein_validation import (  # noqa: E402, F401
    batch_score_proteins,
    check_dependencies,
)
