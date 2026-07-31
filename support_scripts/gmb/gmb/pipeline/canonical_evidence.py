#!/usr/bin/env python3
"""Shared evidence-class vocabulary and balanced-coverage helper.

Single source of truth for two things that MUST be computed identically
wherever they are used, rather than duplicated ad hoc:

1. The mapping from named software sources (Scallop, StringTie, Tiberius,
   Helixer, Minimap2, OrthoDB, GenBlast, UniProt, ...) to independent
   *biological* evidence classes. Before this module existed,
   ``gmb.pipeline.canonical_selection`` and ``gmb.pipeline.interpro_review``
   each computed their own version of this, and they had drifted apart
   (``interpro_review``'s ``evidence_classes()`` deliberately excluded
   protein_validation as a class; ``canonical_selection`` counted raw named
   sources instead of classes at all). Both now import from here.

2. ``balanced_coverage()``, the ``min(qcov, scov)``-based DIAMOND coverage
   measure used by both canonical selection's protein-plausibility scoring
   and the InterPro-review ambiguity/manifest code.

Why classes, not named-source counts
-------------------------------------
Scallop and StringTie are two assemblers customarily run over the SAME
short-read RNA-seq libraries -- counting them as two independent pieces of
support overstates evidence breadth. Grouping named sources into biological
evidence *classes* avoids that, while still keeping named-source provenance
available separately (every consumer can report ``named_evidence_sources``
alongside the class set -- named sources are never hidden or deleted, only
additionally summarised).

Long-read and short-read transcriptomic evidence are deliberately kept as
two distinct classes even though both are "transcriptomic": a long-read
consensus model and a short-read-only model are not interchangeable
evidence, and collapsing them would hide exactly the distinction a reviewer
most wants to see.

``protein_validation`` (whether DIAMOND/psauron support the translated
protein) is a class here for evidence-*breadth* purposes (Part 3 of the
task this module was built for requires it as one of five classes), but it
is computed from the transcript's own validation result
(``has_protein_validation_support``), never from ``evidence_sources`` --
DIAMOND/psauron are not annotation evidence tracks, they are a property of
the translated sequence. Canonical selection's continuous score already has
a dedicated, more granular ``protein_validation_subtotal`` for the
*strength* of that support; the class here only marks "was there any",
which is what breadth needs.

Unknown-source handling
------------------------
A named source that is neither the configured backbone label nor in
``_SOURCE_TO_CLASS`` is mapped to ``EVIDENCE_CLASS_OTHER`` -- never
silently dropped -- and its raw name is returned alongside the class set so
a caller can warn once per run (not once per transcript) about sources this
module does not yet recognise. See ``named_source_evidence_classes``.
"""

from __future__ import annotations

from typing import Optional

# ---------------------------------------------------------------------------
# Evidence class vocabulary (5 required classes + the unknown-source bucket)
# ---------------------------------------------------------------------------

EVIDENCE_CLASS_BACKBONE = "backbone"
EVIDENCE_CLASS_SHORT_READ = "short_read_transcriptomic"
EVIDENCE_CLASS_LONG_READ = "long_read_transcriptomic"
EVIDENCE_CLASS_PROTEIN_ALIGNMENT = "protein_alignment"
EVIDENCE_CLASS_PROTEIN_VALIDATION = "protein_validation"
# Bucket for a named source this module does not recognise and that is not
# the configured backbone label. Deliberately one shared bucket (not
# per-source-name, unlike an earlier "other:<name>" scheme) -- the raw name
# is still preserved and returned separately for warning/audit purposes, but
# the class itself should not fragment into an unbounded number of ad hoc
# classes just because a new named source appears.
EVIDENCE_CLASS_OTHER = "other"

# The five classes Part 3 of the canonical-selection evidence-class task
# requires -- used by tests/tools that want to assert "exactly these five
# (plus backbone) are ever produced for a known source set".
NAMED_SOURCE_EVIDENCE_CLASSES = (
    EVIDENCE_CLASS_BACKBONE,
    EVIDENCE_CLASS_SHORT_READ,
    EVIDENCE_CLASS_LONG_READ,
    EVIDENCE_CLASS_PROTEIN_ALIGNMENT,
)
ALL_EVIDENCE_CLASSES = NAMED_SOURCE_EVIDENCE_CLASSES + (EVIDENCE_CLASS_PROTEIN_VALIDATION,)

# Named source (lower-cased) -> evidence class. The two known backbone
# sources (Helixer, Tiberius) are listed here directly so that the mapping
# works correctly in standalone canonical-selection runs where the caller
# may not have propagated backbone_label from the build step. The
# backbone_label check in named_source_evidence_classes still fires first
# (so a custom backbone label also works), and since both paths produce
# EVIDENCE_CLASS_BACKBONE the result is identical when they agree.
_SOURCE_TO_CLASS = {
    # Known backbone sources
    "helixer": EVIDENCE_CLASS_BACKBONE,
    "tiberius": EVIDENCE_CLASS_BACKBONE,
    # Short-read transcriptomic assemblers
    "scallop": EVIDENCE_CLASS_SHORT_READ,
    "stringtie": EVIDENCE_CLASS_SHORT_READ,
    # Long-read transcriptomic
    "minimap2": EVIDENCE_CLASS_LONG_READ,
    "minimap2consensus": EVIDENCE_CLASS_LONG_READ,
    # Protein alignments
    "orthodb": EVIDENCE_CLASS_PROTEIN_ALIGNMENT,
    "genblast": EVIDENCE_CLASS_PROTEIN_ALIGNMENT,
    "uniprot": EVIDENCE_CLASS_PROTEIN_ALIGNMENT,
}


def named_source_evidence_classes(
    evidence_sources: Optional[str], backbone_label: str = "Helixer"
) -> tuple[set, set]:
    """Map a comma-separated named-source string to evidence classes.

    Returns ``(classes, unknown_names)``: ``classes`` never includes
    ``protein_validation`` (see module docstring -- that is added by
    ``evidence_classes_for_transcript`` from the transcript's own
    protein-validation result, not from this string). ``unknown_names`` is
    the set of raw (lower-cased) source names that were not recognised and
    were mapped into ``EVIDENCE_CLASS_OTHER`` -- callers should aggregate
    this across a whole run and emit ONE warning naming them, rather than
    warning per transcript.
    """
    classes = set()
    unknown_names = set()
    backbone_lower = backbone_label.strip().lower()
    for raw in str(evidence_sources or "").split(","):
        name = raw.strip().lower()
        if not name:
            continue
        if name == backbone_lower:
            classes.add(EVIDENCE_CLASS_BACKBONE)
        elif name in _SOURCE_TO_CLASS:
            classes.add(_SOURCE_TO_CLASS[name])
        else:
            classes.add(EVIDENCE_CLASS_OTHER)
            unknown_names.add(name)
    return classes, unknown_names


def evidence_classes_for_transcript(
    evidence_sources: Optional[str],
    backbone_label: str,
    has_protein_validation_support: bool,
) -> tuple[set, set]:
    """Full per-transcript evidence-class set, including ``protein_validation``.

    ``has_protein_validation_support`` should be True iff the transcript has
    a DIAMOND hit and/or a psauron score (i.e. the same test as
    ``score_transcript``'s ``has_any_protein_support`` in
    ``canonical_selection``) -- this function does not recompute that from
    raw DIAMOND/psauron fields itself, to avoid a second definition of
    "has support" drifting from the scorer's own.

    Returns ``(classes, unknown_names)`` -- see ``named_source_evidence_classes``.
    """
    classes, unknown_names = named_source_evidence_classes(evidence_sources, backbone_label)
    if has_protein_validation_support:
        classes.add(EVIDENCE_CLASS_PROTEIN_VALIDATION)
    return classes, unknown_names


# ---------------------------------------------------------------------------
# Balanced DIAMOND coverage
# ---------------------------------------------------------------------------


def balanced_coverage(qcov: Optional[float], scov: Optional[float]) -> Optional[float]:
    """Combine DIAMOND query and target coverage into one bounded 0-1 value.

    Uses the MINIMUM of the two, not the mean: a hit covering 100% of a
    short query but 10% of a long target is a fragment match, and averaging
    (0.55) would flatter it. The minimum reports the weaker side, which is
    the side that actually limits how much the alignment demonstrates. Raw
    bitscore is deliberately not used for this -- it scales with protein
    length and so systematically favours longer isoforms.

    Both inputs are DIAMOND's native 0-100 percentages; the result is 0-1.
    """
    if qcov is None or scov is None:
        return None
    return round(min(qcov, scov) / 100.0, 4)
