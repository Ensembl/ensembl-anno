#!/usr/bin/env python3
"""CLI entry point for canonical transcript selection.

Usage:
    python -m gmb.cli.canonical_selection \
        --consensus-gff3 out/consensus.gff3 \
        --evidence-attribution out/evidence_attribution.tsv \
        --protein-validation out/protein_validation.tsv \
        --output-dir canonical_out/
"""

from gmb.pipeline.canonical_selection import main

if __name__ == "__main__":
    main()
