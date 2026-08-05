"""gmb-finalise CLI entry point.

Post-build finalisation: regenerates FASTA from the final GFF3, runs
FASTA QC, UTR validation, canonical selection, and produces a handover
manifest.  Designed to be run after ``gmb-build`` completes.
"""

from __future__ import annotations

from gmb.pipeline.finalise import main

__all__ = ["main"]
