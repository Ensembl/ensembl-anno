"""General I/O helpers for the GMB pipeline."""

from __future__ import annotations

import os


def ensure_dir(path: str) -> str:
    """Create directory (and parents) if it doesn't exist. Returns the path."""
    os.makedirs(path, exist_ok=True)
    return path
