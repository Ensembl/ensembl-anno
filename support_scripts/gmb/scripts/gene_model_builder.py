#!/usr/bin/env python3
"""Legacy wrapper — delegates to gmb.pipeline.builder.main().

Preserves backward compatibility:
    python gene_model_builder.py ...
"""
import os
import sys

_pkg_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)

from gmb.pipeline.builder import main  # noqa: E402

if __name__ == "__main__":
    main()
