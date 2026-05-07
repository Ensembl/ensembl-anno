#!/usr/bin/env python3
"""
Compare Annotations
===================
Locus-level comparison of consensus annotation against a reference.

This is a legacy wrapper. All logic now lives in
:mod:`gmb.compare.compare_annotations`.
"""

from __future__ import annotations

import os
import sys

# Ensure the package root (support_scripts/gmb/) is importable
_pkg_root = os.path.dirname(os.path.abspath(__file__))
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)

# Re-export public API for backwards compatibility
from gmb.compare.compare_annotations import (  # noqa: E402, F401
    classify_locus_pairs,
    extract_genbank_proteins,
    generate_comparison_plots,
    load_consensus_genes,
    load_gff,
    main,
    plot_comparison_locus,
    write_summary,
)


if __name__ == "__main__":
    main()
