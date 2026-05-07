#!/usr/bin/env python3
"""Legacy wrapper — delegates to gmb.compare.compare_annotations.main().

Preserves backward compatibility:
    python compare_annotations.py ...
"""
import os
import sys

# Ensure the package root (support_scripts/gmb/) is importable
_pkg_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)

from gmb.compare.compare_annotations import main  # noqa: E402

if __name__ == "__main__":
    main()
