"""gmb-finalise CLI entry point.

Post-build finalisation: regenerates FASTA from the final GFF3, runs
FASTA QC, UTR validation, canonical selection, and produces a handover
manifest.  Designed to be run after ``gmb-build`` completes.
"""

from __future__ import annotations

from gmb.pipeline.finalise import main

__all__ = ["main"]


if __name__ == "__main__":
    # `python -m gmb.cli.finalise` used to import main and exit 0 without running.
    raise SystemExit(main())
