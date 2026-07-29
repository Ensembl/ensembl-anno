#!/usr/bin/env python3
"""CLI entry point for raw long-read -> consensus preprocessing.

Usage:
    python -m gmb.cli.longread_consensus --input raw_minimap2.gtf --output-dir out/
"""

from gmb.pipeline.longread_consensus import main

if __name__ == "__main__":
    main()
