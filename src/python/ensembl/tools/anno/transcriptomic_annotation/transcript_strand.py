# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Transcript strand resolution for spliced long-read alignments.

A spliced aligner reports two independent pieces of strand information:

* the **alignment orientation** (SAM FLAG 0x10), i.e. whether the read sequence
  had to be reverse-complemented to match the reference. For an unstranded cDNA
  library this is essentially a coin flip and says nothing about the transcript.
* the **transcript strand**, which Minimap2 infers from the splice-site motifs it
  chose and reports in the ``ts:A:`` tag.

Minimap2 writes ``ts`` relative to the read as stored in the record, so the
transcript strand on the genome is ``ts`` combined with the alignment
orientation::

    genomic transcript strand = ts XOR alignment_is_reverse

Using the SAM FLAG on its own as the transcript strand randomises the annotated
strand of every model built from an unstranded library.
"""

__all__ = [
    "STRAND_FORWARD",
    "STRAND_REVERSE",
    "STRAND_UNSTRANDED",
    "MissingTranscriptStrandError",
    "MissingTsPolicy",
    "transcript_strand",
]

import logging
from typing import Optional

logger = logging.getLogger("__name__." + __name__)

STRAND_FORWARD = "+"
STRAND_REVERSE = "-"
#: BED and GTF both spell "strand not known" as a full stop.
STRAND_UNSTRANDED = "."

_VALID_TS = (STRAND_FORWARD, STRAND_REVERSE)


class MissingTranscriptStrandError(ValueError):
    """Raised when ``ts`` is absent or unusable and the policy is ``error``."""


class MissingTsPolicy:
    """What to do when ``ts`` is absent or not one of ``+`` / ``-``.

    ``UNSTRANDED``
        Emit ``.``. This is the default. Minimap2 only reports ``ts`` when it has
        splice-site evidence, so an unspliced alignment legitimately has no
        transcript strand and saying so is honest.
    ``ERROR``
        Raise :class:`MissingTranscriptStrandError`. For callers that require
        every model to be stranded.
    ``FLAG``
        Fall back to the alignment orientation. This reproduces the historical
        behaviour and is **wrong** for unstranded libraries; it exists only so a
        caller that genuinely depends on the old output can opt in explicitly.
    """

    UNSTRANDED = "unstranded"
    ERROR = "error"
    FLAG = "flag"

    ALL = (UNSTRANDED, ERROR, FLAG)


def transcript_strand(
    ts_tag: Optional[str],
    is_reverse: bool,
    on_missing: str = MissingTsPolicy.UNSTRANDED,
    read_name: Optional[str] = None,
) -> str:
    """Resolve the genomic transcript strand of a spliced alignment.

    Args:
        ts_tag: Value of the Minimap2 ``ts:A:`` tag (``+`` or ``-``), or None when
            the tag is absent. Any other value is treated as missing.
        is_reverse: True when SAM FLAG 0x10 is set, i.e. the read aligned to the
            reverse strand of the reference.
        on_missing: One of :class:`MissingTsPolicy`.
        read_name: Optional read name, used only to make warnings actionable.

    Returns:
        ``+``, ``-`` or ``.``.

    Raises:
        ValueError: If ``on_missing`` is not a recognised policy.
        MissingTranscriptStrandError: If ``ts`` is unusable under the ``error`` policy.
    """
    if on_missing not in MissingTsPolicy.ALL:
        raise ValueError(f"unknown on_missing policy {on_missing!r}; expected one of {MissingTsPolicy.ALL}")

    if ts_tag in _VALID_TS:
        # ts is relative to the stored read, so undo the alignment orientation.
        forward = (ts_tag == STRAND_FORWARD) != bool(is_reverse)
        return STRAND_FORWARD if forward else STRAND_REVERSE

    where = f" for read {read_name}" if read_name else ""
    if on_missing == MissingTsPolicy.ERROR:
        raise MissingTranscriptStrandError(
            f"no usable Minimap2 ts tag{where} (got {ts_tag!r}); cannot determine transcript strand"
        )
    if on_missing == MissingTsPolicy.FLAG:
        return STRAND_REVERSE if is_reverse else STRAND_FORWARD
    return STRAND_UNSTRANDED
