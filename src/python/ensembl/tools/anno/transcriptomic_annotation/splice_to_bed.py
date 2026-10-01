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
Convert spliced SAM/BAM alignments to BED12.

This replaces ``paftools.js splice2bed`` in the long-read path for one reason:
``splice2bed`` writes the SAM FLAG orientation into the BED strand column and
never reads Minimap2's ``ts:A:`` tag, so every model built from an unstranded
cDNA library gets a randomly-signed strand. See
:mod:`.transcript_strand` for the strand rule.

Everything else about the BED12 record is kept byte-compatible with
``splice2bed`` so that swapping the implementation changes the strand and
nothing else:

* one BED record per retained alignment, named after the read (``/1``, ``/2``
  appended for paired records, as ``splice2bed`` does);
* score 1000; thickStart/thickEnd equal to chromStart/chromEnd;
* ``itemRgb`` ``0,128,255`` for a read with a single primary alignment,
  ``255,0,0`` when a read has more than one, ``0,192,0`` for secondary records;
* blocks split on CIGAR ``N`` only.

Coordinates: SAM POS is 1-based inclusive, BED is 0-based half-open.
``chromStart = POS - 1`` and ``chromEnd = chromStart + reference span``.

No third-party dependency is required: SAM text, gzipped SAM and BGZF/BAM are
all read natively.
"""

__all__ = [
    "AlignmentRecord",
    "ConversionStats",
    "iter_alignments",
    "splice_to_bed",
]

import gzip
import logging
import struct
from dataclasses import dataclass, field
from pathlib import Path
from typing import BinaryIO, Dict, Iterator, List, Optional, Sequence, Tuple

from .transcript_strand import MissingTsPolicy, transcript_strand

logger = logging.getLogger("__name__." + __name__)

FLAG_PAIRED = 0x1
FLAG_UNMAPPED = 0x4
FLAG_REVERSE = 0x10
FLAG_FIRST_IN_PAIR = 0x40
FLAG_SECONDARY = 0x100
FLAG_SUPPLEMENTARY = 0x800

#: CIGAR operations in BAM numeric order.
_CIGAR_OPS = "MIDNSHP=X"
#: Operations that consume reference bases inside an aligned block. Unlike
#: ``splice2bed`` (which only counts M and D) ``=`` and ``X`` are handled too;
#: Minimap2 emits M by default, so this cannot change its output.
_REF_CONSUMING = frozenset("MD=X")

_COLOUR_SINGLE_PRIMARY = "0,128,255"
_COLOUR_MULTI_PRIMARY = "255,0,0"
_COLOUR_SECONDARY = "0,192,0"


@dataclass
class AlignmentRecord:
    """The handful of SAM fields this conversion needs."""

    query_name: str
    flag: int
    reference_name: str
    position: int  # 1-based, as in SAM
    cigar: Sequence[Tuple[str, int]]
    ts_tag: Optional[str] = None

    @property
    def is_unmapped(self) -> bool:
        """True when FLAG 0x4 is set."""
        return bool(self.flag & FLAG_UNMAPPED)

    @property
    def is_reverse(self) -> bool:
        """True when FLAG 0x10 is set."""
        return bool(self.flag & FLAG_REVERSE)

    @property
    def is_secondary(self) -> bool:
        """True when FLAG 0x100 is set."""
        return bool(self.flag & FLAG_SECONDARY)

    @property
    def is_supplementary(self) -> bool:
        """True when FLAG 0x800 is set."""
        return bool(self.flag & FLAG_SUPPLEMENTARY)

    @property
    def bed_name(self) -> str:
        """Read name, with the mate suffix ``splice2bed`` appends for pairs."""
        if self.flag & FLAG_PAIRED:
            return f"{self.query_name}/{(self.flag >> 6) & 3}"
        return self.query_name


@dataclass
class ConversionStats:  # pylint:disable=too-many-instance-attributes
    """Counts describing one conversion, for logging and for the caller."""

    alignments_read: int = 0
    records_written: int = 0
    skipped_unmapped: int = 0
    skipped_secondary: int = 0
    skipped_no_cigar: int = 0
    strand_from_ts: int = 0
    strand_unstranded: int = 0
    strand_from_flag: int = 0
    spliced_without_ts: int = 0
    strand_counts: Dict[str, int] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, int]:
        """Flat dictionary of the counters, for logging or a run manifest."""
        out = {k: v for k, v in self.__dict__.items() if k != "strand_counts"}
        out.update({f"strand_{k}": v for k, v in self.strand_counts.items()})
        return out


def _parse_cigar_text(cigar: str) -> List[Tuple[str, int]]:
    """Parse a SAM CIGAR string into (operation, length) pairs."""
    if cigar == "*":
        return []
    ops: List[Tuple[str, int]] = []
    num = ""
    for char in cigar:
        if char.isdigit():
            num += char
        else:
            if not num or char not in _CIGAR_OPS:
                raise ValueError(f"malformed CIGAR: {cigar!r}")
            ops.append((char, int(num)))
            num = ""
    if num:
        raise ValueError(f"malformed CIGAR: {cigar!r}")
    return ops


def _iter_sam(handle: BinaryIO) -> Iterator[AlignmentRecord]:
    """Yield records from an uncompressed or decompressed SAM stream."""
    for raw in handle:
        line = raw.decode("utf8", "replace")
        if not line or line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            continue
        ts_tag = None
        for extra in fields[11:]:
            if extra.startswith("ts:A:"):
                ts_tag = extra[5:].strip()
                break
        yield AlignmentRecord(
            query_name=fields[0],
            flag=int(fields[1]),
            reference_name=fields[2],
            position=int(fields[3]),
            cigar=_parse_cigar_text(fields[5]),
            ts_tag=ts_tag,
        )


def _skip_bam_tags(block: bytes, offset: int) -> Optional[str]:
    """Walk the tag area of one BAM record and return the ``ts`` value, if any."""
    ts_tag = None
    end = len(block)
    while offset + 3 <= end:
        tag = block[offset : offset + 2].decode("ascii", "replace")
        val_type = chr(block[offset + 2])
        offset += 3
        if val_type == "A":
            value = chr(block[offset])
            offset += 1
            if tag == "ts":
                ts_tag = value
        elif val_type in "cC":
            offset += 1
        elif val_type in "sS":
            offset += 2
        elif val_type in "iIf":
            offset += 4
        elif val_type in "ZH":
            nul = block.index(b"\x00", offset)
            offset = nul + 1
        elif val_type == "B":
            sub = chr(block[offset])
            count = struct.unpack_from("<i", block, offset + 1)[0]
            width = {"c": 1, "C": 1, "s": 2, "S": 2, "i": 4, "I": 4, "f": 4}[sub]
            offset += 5 + count * width
        else:  # unknown type: the rest of the record cannot be walked safely
            break
    return ts_tag


def _iter_bam(handle: BinaryIO) -> Iterator[AlignmentRecord]:
    """Yield records from a BGZF/BAM stream.

    BGZF is a series of standard gzip members, so :mod:`gzip` decompresses it
    transparently and only the BAM binary layout has to be decoded here.
    """
    magic = handle.read(4)
    if magic != b"BAM\x01":
        raise ValueError("not a BAM stream: missing BAM\\1 magic")
    (l_text,) = struct.unpack("<i", handle.read(4))
    handle.read(l_text)
    (n_ref,) = struct.unpack("<i", handle.read(4))
    references: List[str] = []
    for _ in range(n_ref):
        (l_name,) = struct.unpack("<i", handle.read(4))
        name = handle.read(l_name)[:-1].decode("utf8", "replace")
        handle.read(4)  # l_ref, unused
        references.append(name)

    while True:
        header = handle.read(4)
        if len(header) < 4:
            return
        (block_size,) = struct.unpack("<i", header)
        block = handle.read(block_size)
        if len(block) < block_size:
            return
        ref_id, pos, l_read_name, _mapq, _bin, n_cigar_op, flag, l_seq = struct.unpack_from(
            "<iiBBHHHi", block, 0
        )
        offset = 32
        query_name = block[offset : offset + l_read_name - 1].decode("utf8", "replace")
        offset += l_read_name
        cigar: List[Tuple[str, int]] = []
        for _ in range(n_cigar_op):
            (packed,) = struct.unpack_from("<I", block, offset)
            offset += 4
            cigar.append((_CIGAR_OPS[packed & 0xF], packed >> 4))
        offset += (l_seq + 1) // 2 + l_seq  # packed sequence + qualities
        yield AlignmentRecord(
            query_name=query_name,
            flag=flag,
            reference_name=references[ref_id] if 0 <= ref_id < len(references) else "*",
            position=pos + 1,  # BAM stores 0-based; SAM/AlignmentRecord use 1-based
            cigar=cigar,
            ts_tag=_skip_bam_tags(block, offset),
        )


def iter_alignments(path: Path) -> Iterator[AlignmentRecord]:
    """Yield alignments from a SAM, gzipped SAM or BAM file.

    The format is detected from the file's content, not its extension.
    """
    path = Path(path)
    with open(path, "rb") as probe:
        head = probe.read(2)
    if head == b"\x1f\x8b":
        with gzip.open(path, "rb") as handle:
            if handle.read(4) == b"BAM\x01":
                handle.seek(0)
                yield from _iter_bam(handle)
            else:
                handle.seek(0)
                yield from _iter_sam(handle)
    else:
        with open(path, "rb") as handle:
            yield from _iter_sam(handle)


def _blocks_from_cigar(cigar: Sequence[Tuple[str, int]]) -> Tuple[List[int], List[int], int]:
    """Split a CIGAR into BED12 blocks.

    Returns:
        (block_starts, block_sizes, reference_span), all relative to the
        alignment start and in reference coordinates.
    """
    block_starts: List[int] = []
    block_sizes: List[int] = []
    span = 0
    block_open = 0
    for op, length in cigar:
        if op in _REF_CONSUMING:
            span += length
        elif op == "N":
            block_starts.append(block_open)
            block_sizes.append(span - block_open)
            span += length
            block_open = span
    block_starts.append(block_open)
    block_sizes.append(span - block_open)
    return block_starts, block_sizes, span


def _format_bed12(record: AlignmentRecord, strand: str, colour: str) -> str:
    """Render one BED12 line."""
    block_starts, block_sizes, span = _blocks_from_cigar(record.cigar)
    start = record.position - 1
    end = start + span
    return (
        f"{record.reference_name}\t{start}\t{end}\t{record.bed_name}\t1000\t{strand}\t"
        f"{start}\t{end}\t{colour}\t{len(block_starts)}\t"
        f"{','.join(str(s) for s in block_sizes)},\t"
        f"{','.join(str(s) for s in block_starts)},\n"
    )


def splice_to_bed(
    alignment_file: Path,
    bed_file: Path,
    on_missing_ts: str = MissingTsPolicy.UNSTRANDED,
    keep_secondary: bool = False,
) -> ConversionStats:
    """Convert spliced SAM/BAM alignments to BED12 with correct transcript strands.

    Args:
        alignment_file: SAM, gzipped SAM or BAM produced by a spliced aligner.
        bed_file: Destination BED12 path.
        on_missing_ts: Behaviour when Minimap2 reported no usable ``ts`` tag; see
            :class:`.transcript_strand.MissingTsPolicy`. Defaults to emitting ``.``.
        keep_secondary: Retain secondary alignments (FLAG 0x100). False matches
            ``paftools.js splice2bed`` without ``-m``. Supplementary alignments
            (FLAG 0x800) are retained either way, also matching ``splice2bed``.

    Returns:
        A :class:`ConversionStats` describing what was read and written.

    Note:
        Records are grouped by consecutive query name, as ``splice2bed`` does, to
        decide the ``itemRgb`` colour. Aligner output is query-grouped; with a
        coordinate-sorted input only that cosmetic colour differs.
    """
    stats = ConversionStats()
    alignment_file = Path(alignment_file)
    bed_file = Path(bed_file)

    def flush(group: List[Tuple[AlignmentRecord, str]], out) -> None:
        primaries = sum(1 for rec, _ in group if not rec.is_secondary)
        for rec, strand in group:
            if rec.is_secondary:
                colour = _COLOUR_SECONDARY
            else:
                colour = _COLOUR_MULTI_PRIMARY if primaries > 1 else _COLOUR_SINGLE_PRIMARY
            out.write(_format_bed12(rec, strand, colour))
            stats.records_written += 1
            stats.strand_counts[strand] = stats.strand_counts.get(strand, 0) + 1

    with open(bed_file, "w", encoding="utf8") as out:
        group: List[Tuple[AlignmentRecord, str]] = []
        group_name: Optional[str] = None
        for record in iter_alignments(alignment_file):
            stats.alignments_read += 1
            if record.is_unmapped or record.reference_name == "*":
                stats.skipped_unmapped += 1
                continue
            if record.is_secondary and not keep_secondary:
                stats.skipped_secondary += 1
                continue
            if not record.cigar:
                stats.skipped_no_cigar += 1
                continue

            has_ts = record.ts_tag in ("+", "-")
            strand = transcript_strand(
                record.ts_tag,
                record.is_reverse,
                on_missing=on_missing_ts,
                read_name=record.query_name,
            )
            if has_ts:
                stats.strand_from_ts += 1
            else:
                if any(op == "N" for op, _ in record.cigar):
                    stats.spliced_without_ts += 1
                if on_missing_ts == MissingTsPolicy.FLAG:
                    stats.strand_from_flag += 1
                else:
                    stats.strand_unstranded += 1

            if record.query_name != group_name:
                flush(group, out)
                group = []
                group_name = record.query_name
            group.append((record, strand))
        flush(group, out)

    if stats.spliced_without_ts:
        logger.warning(
            "%s: %d spliced alignment(s) had no usable Minimap2 ts tag; "
            "their transcript strand was resolved by the '%s' policy",
            alignment_file.name,
            stats.spliced_without_ts,
            on_missing_ts,
        )
    logger.info("%s -> %s: %s", alignment_file.name, bed_file.name, stats.as_dict())
    return stats
