#!/usr/bin/env python3
"""Legacy wrapper — delegates to gmb.compare.visualize_disagreements.main().

Preserves backward compatibility:
    python visualize_disagreements.py ...
"""
import os
import sys

# Ensure the package root (support_scripts/gmb/) is importable
_pkg_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)

from gmb.compare.visualize_disagreements import main  # noqa: E402

if __name__ == "__main__":
    main()
