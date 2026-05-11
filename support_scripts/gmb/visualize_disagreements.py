#!/usr/bin/env python3
"""Backward-compatible entry point — delegates to gmb.compare.visualize_disagreements.main().

Canonical usage:
    python -m gmb.cli.visualize ...
"""
import os
import sys

_pkg_root = os.path.dirname(os.path.abspath(__file__))
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)

from gmb.compare.visualize_disagreements import main  # noqa: E402

if __name__ == "__main__":
    main()
