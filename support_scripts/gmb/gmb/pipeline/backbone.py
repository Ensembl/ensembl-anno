"""Resolving which annotation file supplies the ``backbone`` evidence role.

GMB has always treated the ab initio backbone as "a path plus a source label":
the label is what :class:`~gmb.pipeline.canonical_evidence.EvidenceRoles`
resolves to the ``backbone`` role, and the numeric weight is looked up by role,
not by tool name. Only the *command line* was tool-specific, offering exactly
two flags (``--helixer`` / ``--tiberius``) and deriving the label from whichever
was used.

This module adds a generic ``--backbone`` input alongside them, so any predictor
can fill the backbone role while keeping its own name in the output attribution.
It introduces no tool-specific behaviour: the resolved label is fed through the
same ``scoring.backbone_label`` field the existing flags already set.
"""

from __future__ import annotations

import gzip
import os
from typing import Optional, Tuple

__all__ = [
    "DEFAULT_BACKBONE_LABEL",
    "FLAG_LABELS",
    "BackboneInputError",
    "source_label_from_annotation",
    "resolve_backbone_input",
]

#: Label used when a generic backbone is supplied and no label can be determined.
DEFAULT_BACKBONE_LABEL = "Backbone"

#: Labels implied by the two historical backbone flags.
FLAG_LABELS = {"helixer": "Helixer", "tiberius": "Tiberius"}

#: Label assumed when no backbone track is supplied at all. Preserves the
#: historical default so runs without a backbone are unchanged.
_NO_BACKBONE_LABEL = "Helixer"


class BackboneInputError(ValueError):
    """Raised when more than one backbone track is supplied."""


def _open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf8", errors="replace")
    return open(path, "r", encoding="utf8", errors="replace")


def source_label_from_annotation(path: str, max_rows: int = 5000) -> Optional[str]:
    """Return the GFF3/GTF ``source`` column value when the file uses just one.

    Column 2 of a GFF3/GTF record names the program that produced it, which is
    exactly the attribution wanted here (``Vipsania``, ``Tiberius``, ...).
    Returns None when the file is unreadable, empty, or mixes several sources,
    so the caller can fall back to an explicit label.
    """
    if not path or not os.path.exists(path):
        return None
    seen = set()
    try:
        with _open_text(path) as handle:
            for i, line in enumerate(handle):
                if i >= max_rows:
                    break
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) < 3:
                    continue
                source = fields[1].strip()
                if source and source != ".":
                    seen.add(source)
                if len(seen) > 1:
                    return None
    except OSError:
        return None
    if len(seen) != 1:
        return None
    return seen.pop()


def resolve_backbone_input(
    helixer: Optional[str] = None,
    tiberius: Optional[str] = None,
    backbone: Optional[str] = None,
    backbone_label: Optional[str] = None,
) -> Tuple[Optional[str], str]:
    """Resolve the single ab initio backbone track to ``(path, label)``.

    Exactly one of *helixer*, *tiberius* and *backbone* may be supplied.
    ``--helixer`` and ``--tiberius`` keep their historical labels. For a generic
    ``--backbone`` the label is, in order of preference: an explicit
    *backbone_label*, the annotation's own source column, then
    :data:`DEFAULT_BACKBONE_LABEL`.

    An explicit *backbone_label* overrides the label for whichever flag was used,
    so an operator can rename a track without changing how it is supplied.

    Raises:
        BackboneInputError: if more than one backbone track is given.
    """
    supplied = [(name, value) for name, value in
                (("--helixer", helixer), ("--tiberius", tiberius), ("--backbone", backbone))
                if value]
    if len(supplied) > 1:
        names = ", ".join(name for name, _ in supplied)
        raise BackboneInputError(
            f"pass only one ab initio backbone track; got {names}. "
            "The backbone is a single evidence role."
        )

    if not supplied:
        return None, (backbone_label or _NO_BACKBONE_LABEL)

    flag, path = supplied[0]
    if backbone_label:
        return path, backbone_label
    if flag == "--backbone":
        return path, (source_label_from_annotation(path) or DEFAULT_BACKBONE_LABEL)
    return path, FLAG_LABELS[flag.lstrip("-")]
