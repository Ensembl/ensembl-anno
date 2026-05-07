#!/usr/bin/env python3
"""CLI entry point for annotation comparison.

Usage:
    python -m gmb.cli.compare --consensus out/consensus.gff3 --reference ref.gff3 ...
"""

from gmb.compare.compare_annotations import main

if __name__ == "__main__":
    main()
