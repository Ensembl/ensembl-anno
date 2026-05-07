#!/usr/bin/env python3
"""Gene Model Builder orchestrator.

Loads evidence, applies filters, performs optional protein validation,
and generates consensus models.

This is a legacy wrapper. All logic now lives in :mod:`gmb.pipeline.builder`.

DataFrame Schema Contract
-------------------------
Evidence DataFrames passed through the pipeline share these columns:

* ``Chromosome`` : str -- sequence/contig name
* ``Start`` : int -- 0-based start (half-open, pyranges convention)
* ``End`` : int -- 1-based end (half-open, pyranges convention)
* ``Strand`` : str -- ``"+"`` or ``"-"``
* ``transcript_id`` : str -- unique transcript identifier (source-prefixed)
* ``gene_id`` : str -- gene-level grouping identifier (source-prefixed)
* ``Source`` : str -- evidence origin, e.g. ``"Scallop"``, ``"Helixer"``
* ``Feature`` : str -- GFF3/GTF feature type (``"exon"``, ``"CDS"``, etc.)

Additional columns may be present depending on the source format
(``Score``, ``Coverage``, ``Identity``, ``combined_evidence``, etc.).
"""

from __future__ import annotations

import os
import sys

# Ensure the package root (support_scripts/gmb/) is importable
_pkg_root = os.path.dirname(os.path.abspath(__file__))
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)

# Re-export public API for backwards compatibility
from gmb.pipeline.builder import (  # noqa: E402, F401
    compute_percentile_guardrails,
    compute_utr_end_support,
    load_evidence,
    main,
    parse_args,
)


if __name__ == "__main__":
    main()
