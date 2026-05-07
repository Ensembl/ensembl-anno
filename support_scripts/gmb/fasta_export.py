#!/usr/bin/env python3
"""FASTA export — legacy wrapper. All logic now in :mod:`gmb.pipeline.fasta_export`."""
import os, sys  # noqa: E401
_pkg_root = os.path.dirname(os.path.abspath(__file__))
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)
from gmb.pipeline.fasta_export import export_fasta  # noqa: E402, F401
