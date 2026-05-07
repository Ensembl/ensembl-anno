#!/usr/bin/env python3
"""Scoring — legacy wrapper. All logic now in :mod:`gmb.pipeline.scoring`."""
import os, sys  # noqa: E401
_pkg_root = os.path.dirname(os.path.abspath(__file__))
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)
from gmb.pipeline.scoring import score_model, select_isoforms  # noqa: E402, F401
