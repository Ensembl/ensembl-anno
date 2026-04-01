"""Pytest configuration: ensure the gmb package directory is importable.

All modules (scoring, config, evidence_filter, etc.) live in the parent
directory of this tests/ folder.  Adding it to sys.path here means every
test file can ``import scoring`` without a separate ``sys.path`` hack.
"""

import os
import sys

# Parent of tests/ == support_scripts/gmb/  (the "package root")
_GMB_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _GMB_DIR not in sys.path:
    sys.path.insert(0, _GMB_DIR)
