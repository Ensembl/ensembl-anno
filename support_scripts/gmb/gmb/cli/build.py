#!/usr/bin/env python3
"""CLI entry point for the Gene Model Builder pipeline.

Usage:
    python -m gmb.cli.build --genome genome.fa --output-dir out/ ...
"""

from gmb.pipeline.builder import main

if __name__ == "__main__":
    main()
