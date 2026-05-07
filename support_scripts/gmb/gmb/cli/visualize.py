#!/usr/bin/env python3
"""CLI entry point for disagreement visualization.

Usage:
    python -m gmb.cli.visualize --output-dir qc/ ...
"""

from gmb.compare.visualize_disagreements import main

if __name__ == "__main__":
    main()
