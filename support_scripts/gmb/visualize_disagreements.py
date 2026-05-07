#!/usr/bin/env python3
"""Visualize disagreements between annotation tracks.

This is a legacy wrapper. All logic now lives in
:mod:`gmb.compare.visualize_disagreements`.
"""

from __future__ import annotations

import os
import sys

# Ensure the package root (support_scripts/gmb/) is importable
_pkg_root = os.path.dirname(os.path.abspath(__file__))
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)

# Re-export public API for backwards compatibility
from gmb.compare.visualize_disagreements import (  # noqa: E402, F401
    find_disagreements,
    load_data,
    main,
    plot_locus,
)


if __name__ == "__main__":
    main()
