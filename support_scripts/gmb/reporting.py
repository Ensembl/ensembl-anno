#!/usr/bin/env python3
"""Reporting — legacy wrapper. All logic now in :mod:`gmb.pipeline.reporting`."""
import os, sys  # noqa: E401
_pkg_root = os.path.dirname(os.path.abspath(__file__))
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)
from gmb.pipeline.reporting import generate_report  # noqa: E402, F401
