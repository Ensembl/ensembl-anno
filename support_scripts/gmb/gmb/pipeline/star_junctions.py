"""Quantitative splice-junction evidence from per-library STAR ``SJ.out.tab`` files.

This module only *loads, validates and indexes* junction observations. It makes no
decision about gene models; :mod:`gmb.pipeline.junction_diagnostics` reports what
the evidence says about a structure, and future scoring or reconstruction code can
query the same :class:`JunctionIndex`.

Input contract (see ``docs/star_junctions.md``)
-----------------------------------------------
A tab-separated **manifest** with one row per library:

* required: ``library_id``, ``path`` (relative paths resolve against the manifest's
  directory);
* ``role`` -- which evidence set the library belongs to (e.g. ``construction``,
  ``holdout``). Only libraries whose role is selected are read at all, so held-out
  libraries cannot leak into construction evidence. A manifest without a ``role``
  column is only accepted with ``roles=("all",)``;
* optional identifiers: ``sample_id``, ``study``, ``strain``, ``condition``,
  ``group`` (an evidence-independence class; defaults to study|strain|condition,
  then to the library itself);
* optional alignment facts that make *absence* interpretable: ``max_intron``,
  ``min_intron`` (the aligner's limits, bp), ``noncanonical_filtered``
  (``true``/``false``), ``count_unit`` (``mate`` or ``fragment``).

STAR conventions handled here
-----------------------------
* Columns 2-3 are the first and last intronic base, 1-based inclusive. They are
  stored half-open 0-based, ``(start0, end) = (col2 - 1, col3)``, which is exactly
  ``(previous_exon.End, next_exon.Start)`` in GMB's pyranges coordinates.
* Column 4 strand (0 undefined, 1 +, 2 -) is derived by STAR from the intron
  motif, not from the library; it is kept as reported.
* Column 6 ("annotated") is read and ignored: in two-pass mode it marks junctions
  found in pass 1, so it is not independent evidence.
* Columns 7/8 are read counts in the library's ``count_unit``; nothing here
  converts mates to fragments.
"""

from __future__ import annotations

import bisect
import csv
import gzip
import os
import time
from dataclasses import dataclass, field

import numpy as np

STAR_STRAND = {0: ".", 1: "+", 2: "-"}
STAR_MOTIF = {0: "non-canonical", 1: "GT/AG", 2: "CT/AC", 3: "GC/AG", 4: "CT/GC", 5: "AT/AC", 6: "GT/AT"}
# Library load outcomes. Only LOADED libraries contribute evidence.
LOADED = "LOADED"
FAILED_EMPTY = "FAILED_EMPTY"
MISSING_FILE = "MISSING_FILE"
MALFORMED = "MALFORMED"
NO_KNOWN_SEQIDS = "NO_KNOWN_SEQIDS"
NOT_SELECTED = "NOT_SELECTED"

_TRUE = {"true", "yes", "1", "y", "t"}
_FALSE = {"false", "no", "0", "n", "f"}


class JunctionManifestError(ValueError):
    """The manifest itself cannot be used (as opposed to one bad library)."""


@dataclass
class Library:
    library_id: str
    path: str
    role: str = ""
    sample_id: str = ""
    study: str = ""
    strain: str = ""
    condition: str = ""
    group: str = ""
    count_unit: str = "unknown"
    max_intron: int | None = None
    min_intron: int | None = None
    noncanonical_filtered: bool | None = None
    status: str = ""
    message: str = ""
    n_junctions: int = 0
    n_rows_unknown_seqid: int = 0
    sha256: str = ""
    load_seconds: float = 0.0

    def as_row(self) -> dict:
        return {
            "library_id": self.library_id, "role": self.role, "status": self.status,
            "n_junctions": self.n_junctions, "n_rows_unknown_seqid": self.n_rows_unknown_seqid,
            "study": self.study, "strain": self.strain, "condition": self.condition,
            "group": self.group, "sample_id": self.sample_id, "count_unit": self.count_unit,
            "max_intron": "" if self.max_intron is None else self.max_intron,
            "min_intron": "" if self.min_intron is None else self.min_intron,
            "noncanonical_filtered": "" if self.noncanonical_filtered is None
            else str(self.noncanonical_filtered).lower(),
            "path": self.path, "sha256": self.sha256, "message": self.message,
        }


def _opt_int(value: str, column: str, row: int) -> int | None:
    value = (value or "").strip()
    if not value or value.lower() in ("na", "unknown", "none"):
        return None
    try:
        out = int(value)
    except ValueError:
        raise JunctionManifestError(f"manifest row {row}: {column}={value!r} is not an integer") from None
    if out < 0:
        raise JunctionManifestError(f"manifest row {row}: {column} must be >= 0")
    return out


def _opt_bool(value: str, column: str, row: int) -> bool | None:
    v = (value or "").strip().lower()
    if not v or v in ("na", "unknown", "none"):
        return None
    if v in _TRUE:
        return True
    if v in _FALSE:
        return False
    raise JunctionManifestError(f"manifest row {row}: {column}={value!r} is not true/false")


def read_manifest(path: str, roles=("construction",)) -> list[Library]:
    """Parse a junction manifest. Returns every library; ``status`` is
    ``NOT_SELECTED`` for libraries whose role is not in *roles* (``"all"`` selects
    every library)."""
    if not os.path.exists(path):
        raise JunctionManifestError(f"junction manifest not found: {path}")
    roles = tuple(r.strip() for r in roles if r and r.strip())
    if not roles:
        raise JunctionManifestError("no junction library roles selected")
    base = os.path.dirname(os.path.abspath(path))
    with open(path, newline="") as fh:
        rows = [r for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")]
    if not rows:
        raise JunctionManifestError(f"junction manifest has no libraries: {path}")
    cols = set(rows[0])
    missing = {"library_id", "path"} - cols
    if missing:
        raise JunctionManifestError(f"junction manifest lacks required column(s): {sorted(missing)}")
    if "role" not in cols and "all" not in roles:
        raise JunctionManifestError(
            "junction manifest has no 'role' column; add one (e.g. construction/holdout) or "
            "select every library explicitly with roles=all")
    libs, seen = [], set()
    for i, r in enumerate(rows, start=2):
        lid = (r.get("library_id") or "").strip()
        p = (r.get("path") or "").strip()
        if not lid or not p:
            raise JunctionManifestError(f"manifest row {i}: empty library_id or path")
        if lid in seen:
            raise JunctionManifestError(f"manifest row {i}: duplicate library_id {lid!r}")
        seen.add(lid)
        def get(k, row=r):
            return (row.get(k) or "").strip()

        study, strain, condition = get("study"), get("strain"), get("condition")
        group = get("group") or ("|".join((study, strain, condition))
                                 if (study or strain or condition) else lid)
        lib = Library(
            library_id=lid, path=p if os.path.isabs(p) else os.path.join(base, p),
            role=get("role"), sample_id=get("sample_id"), study=study or group, strain=strain,
            condition=condition, group=group, count_unit=get("count_unit") or "unknown",
            max_intron=_opt_int(r.get("max_intron"), "max_intron", i),
            min_intron=_opt_int(r.get("min_intron"), "min_intron", i),
            noncanonical_filtered=_opt_bool(r.get("noncanonical_filtered"), "noncanonical_filtered", i),
        )
        if "all" not in roles and lib.role not in roles:
            lib.status = NOT_SELECTED
        libs.append(lib)
    return libs


def _sha256(path: str) -> str:
    import hashlib

    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


_SJ_COLUMNS = ["chrom", "s1", "e1", "strand", "motif", "annotated", "unique", "multi", "overhang"]


def parse_sj_file(path: str, seq_lengths: dict | None = None):
    """Parse one STAR ``SJ.out.tab``-format file (plain or gzip, detected by magic).

    Returns ``(status, message, frame, n_unknown_seqid)``. *frame* has columns
    ``chrom, start0, end, strand, motif, unique, multi, overhang`` (half-open
    0-based intron; strand/motif as STAR codes). Any malformed line rejects the
    whole file: a partially parsed library would bias its counts silently.
    """
    import pandas as pd

    empty = pd.DataFrame(columns=["chrom", "start0", "end", "strand", "motif", "unique", "multi", "overhang"])
    if not os.path.exists(path):
        return MISSING_FILE, "file not found", empty, 0
    if os.path.getsize(path) == 0:
        return FAILED_EMPTY, "zero-byte file: failed or missing aligner output (not zero junctions)", empty, 0
    with open(path, "rb") as fh:
        magic = fh.read(2)
    compression = "gzip" if magic == b"\x1f\x8b" else None
    try:
        if compression:
            with gzip.open(path, "rb") as fh:
                last = fh.read()[-1:]
        else:
            with open(path, "rb") as fh:
                fh.seek(-1, os.SEEK_END)
                last = fh.read(1)
    except (OSError, EOFError) as exc:
        return MALFORMED, f"unreadable: {exc}", empty, 0
    if last != b"\n":
        return MALFORMED, "last line not newline-terminated (possible truncation)", empty, 0
    read = dict(sep="\t", header=None, comment="#", compression=compression)
    typed = {0: str, **{i: "int64" for i in range(1, 9)}}
    try:
        df = pd.read_csv(path, dtype=typed, **read)
        if df.shape[1] != 9 or df.isna().any().any():
            raise ValueError("shape")
    except pd.errors.EmptyDataError:
        return MALFORMED, "no junction rows (only comments or blank lines)", empty, 0
    except (ValueError, TypeError, pd.errors.ParserError, UnicodeDecodeError, EOFError, OSError):
        # Typed parsing failed: re-read as text only to say *where* it is malformed.
        try:
            df = pd.read_csv(path, dtype=str, keep_default_na=False, na_filter=False, **read)
        except (pd.errors.ParserError, UnicodeDecodeError, EOFError, OSError) as exc:
            return MALFORMED, f"unparseable: {str(exc).splitlines()[0]}", empty, 0
        if df.shape[1] != 9:
            return MALFORMED, f"{df.shape[1]} columns, expected 9", empty, 0
        for i, c in enumerate(_SJ_COLUMNS[1:], start=1):
            v = pd.to_numeric(df[i], errors="coerce")
            bad = v.isna() | (v != v.round())
            if bad.any():
                return MALFORMED, f"line {int(bad.idxmax()) + 1}: non-integer {c}", empty, 0
        return MALFORMED, "unparseable rows", empty, 0
    df.columns = _SJ_COLUMNS
    num = {c: df[c] for c in _SJ_COLUMNS[1:]}
    bad = ((num["s1"] < 1) | (num["e1"] < num["s1"]) | ~num["strand"].isin(list(STAR_STRAND))
           | ~num["motif"].isin(list(STAR_MOTIF)) | (num["unique"] < 0) | (num["multi"] < 0)
           | (num["overhang"] < 0))
    if bad.any():
        return MALFORMED, f"line {int(bad.idxmax()) + 1}: value out of range", empty, 0
    out = pd.DataFrame({"chrom": df["chrom"], "start0": num["s1"] - 1, "end": num["e1"],
                        "strand": num["strand"].astype("int8"), "motif": num["motif"].astype("int8"),
                        "unique": num["unique"], "multi": num["multi"],
                        "overhang": num["overhang"].astype("int32")})
    dup = out.duplicated(["chrom", "start0", "end"])
    if dup.any():
        r = out.loc[dup.idxmax()]
        return MALFORMED, f"line {int(dup.idxmax()) + 1}: duplicate junction {r.chrom}:{r.start0 + 1}-{r.end}", empty, 0
    unknown = 0
    if seq_lengths is not None:
        known = out["chrom"].isin(set(seq_lengths))
        unknown = int((~known).sum())
        out = out[known]
        if out.empty and unknown:
            return NO_KNOWN_SEQIDS, f"none of {unknown} junction rows is on a genome sequence", empty, unknown
        lens = out["chrom"].map(seq_lengths)
        beyond = out["end"] > lens
        if beyond.any():
            r = out[beyond].iloc[0]
            return MALFORMED, f"intron {r.chrom}:{r.start0 + 1}-{r.end} beyond sequence end", empty, 0
    return LOADED, "", out.reset_index(drop=True), unknown


@dataclass
class JunctionSupport:
    """Aggregated observation of one junction across the loaded libraries."""

    jid: int
    chrom: str
    start0: int
    end: int
    strand: str
    motif: str
    unique_reads: int
    multi_reads: int
    n_libraries: int          # libraries with >= 1 unique read
    n_libraries_any: int      # libraries with >= 1 unique or multimapping read
    n_groups: int
    n_studies: int
    max_overhang: int

    @property
    def start1(self) -> int:
        return self.start0 + 1


@dataclass
class JunctionIndex:
    """All junctions of the loaded libraries, indexed for exact, site and range queries.

    Per-junction aggregates are numpy arrays; per-library counts are kept in CSR
    form (``lib_ptr``/``lib_idx``/``lib_unique``/``lib_multi``) so future weighting
    or independence policies can be computed without re-reading files.
    """

    libraries: list[Library] = field(default_factory=list)
    chroms: list[str] = field(default_factory=list)
    stats: dict = field(default_factory=dict)

    # ---------------------------------------------------------------- build
    @classmethod
    def from_manifest(cls, manifest: str, roles=("construction",), seq_lengths: dict | None = None,
                      checksum: bool = True) -> JunctionIndex:
        t0 = time.perf_counter()
        libs = read_manifest(manifest, roles)
        idx = cls(libraries=libs)
        idx._load([lib for lib in libs if lib.status != NOT_SELECTED], seq_lengths, checksum)
        idx.stats.update(manifest=os.path.abspath(manifest), roles=list(roles),
                         load_seconds=round(time.perf_counter() - t0, 3))
        return idx

    def _load(self, libs, seq_lengths, checksum):
        import pandas as pd

        parts, loaded, chrom_code = [], [], {}
        for lib in libs:
            t = time.perf_counter()
            status, msg, frame, unknown = parse_sj_file(lib.path, seq_lengths)
            lib.status, lib.message, lib.n_rows_unknown_seqid = status, msg, unknown
            if status == LOADED:
                lib.n_junctions = len(frame)
                if checksum:
                    lib.sha256 = _sha256(lib.path)
                local, names = pd.factorize(frame["chrom"])
                glob = np.asarray([chrom_code.setdefault(c, len(chrom_code)) for c in names],
                                  dtype=np.int32)
                n = len(frame)
                parts.append((glob[local] if n else np.zeros(0, np.int32), frame["start0"].to_numpy(np.int64),
                              frame["end"].to_numpy(np.int64), frame["strand"].to_numpy(np.int8),
                              frame["motif"].to_numpy(np.int8), frame["unique"].to_numpy(np.int64),
                              frame["multi"].to_numpy(np.int64), frame["overhang"].to_numpy(np.int32),
                              np.full(n, len(loaded), dtype=np.int32)))
                loaded.append(lib)
            lib.load_seconds = round(time.perf_counter() - t, 4)
        self.loaded = loaded
        # Chromosome codes in sorted-name order, so every array is deterministic.
        self.chroms = sorted(chrom_code)
        remap = np.zeros(max(1, len(chrom_code)), dtype=np.int32)
        for name, code in chrom_code.items():
            remap[code] = self.chroms.index(name)
        cat = [np.concatenate([p[i] for p in parts]) if parts else np.zeros(0, dt)
               for i, dt in enumerate((np.int32, np.int64, np.int64, np.int8, np.int8, np.int64,
                                       np.int64, np.int32, np.int32))]
        ch, s0, en, st, mo, un, mu, ov, li = cat
        ch = remap[ch] if len(ch) else ch
        del parts, cat
        order = np.lexsort((li, en, s0, ch))  # junction key, then library
        ch, s0, en, st, mo, un, mu, ov, li = (a[order] for a in (ch, s0, en, st, mo, un, mu, ov, li))
        del order
        first = np.ones(len(ch), dtype=bool)
        if len(ch) > 1:
            first[1:] = (ch[1:] != ch[:-1]) | (s0[1:] != s0[:-1]) | (en[1:] != en[:-1])
        jid_of_obs = np.cumsum(first) - 1
        J = int(first.sum())
        n_obs = len(ch)
        self.lib_ptr = np.append(np.flatnonzero(first), n_obs).astype(np.int64)
        self.lib_idx, self.lib_unique, self.lib_multi, self.lib_overhang = li, un, mu, ov
        self.j_chrom, self.j_start, self.j_end = ch[first], s0[first], en[first]
        self.j_strand, self.j_motif = st[first], mo[first]
        del s0, en
        incons = np.zeros(J, dtype=bool)
        if J:
            np.logical_or.at(incons, jid_of_obs, (st != self.j_strand[jid_of_obs]) | (mo != self.j_motif[jid_of_obs]))
        inconsistent = int(incons.sum())
        starts = self.lib_ptr[:-1]
        if J:
            self.j_unique = np.add.reduceat(self.lib_unique, starts)
            self.j_multi = np.add.reduceat(self.lib_multi, starts)
            self.j_overhang = np.maximum.reduceat(self.lib_overhang, starts)
            supported = self.lib_unique > 0
            self.j_nlib = np.add.reduceat(supported.astype(np.int32), starts).astype(np.int32)
            self.j_nlib_any = np.diff(self.lib_ptr).astype(np.int32)
        else:
            self.j_unique = self.j_multi = np.zeros(0, np.int64)
            self.j_overhang = self.j_nlib = self.j_nlib_any = np.zeros(0, np.int32)
            supported = np.zeros(0, dtype=bool)
        group_codes = {g: i for i, g in enumerate(sorted({lib.group for lib in loaded}))}
        study_codes = {s: i for i, s in enumerate(sorted({lib.study for lib in loaded}))}
        self.lib_group_code = np.asarray([group_codes[lib.group] for lib in loaded], dtype=np.int64)
        self.lib_study_code = np.asarray([study_codes[lib.study] for lib in loaded], dtype=np.int64)
        self.n_groups_total, self.n_studies_total = len(group_codes), len(study_codes)

        def _distinct(codes, width):
            out = np.zeros(J, dtype=np.int32)
            if J and supported.any():
                pair = np.unique(jid_of_obs[supported] * width + codes[self.lib_idx[supported]])
                np.add.at(out, pair // width, 1)
            return out

        self.j_ngroup = _distinct(self.lib_group_code, max(1, len(group_codes)))
        self.j_nstudy = _distinct(self.lib_study_code, max(1, len(study_codes)))
        keys = list(zip((self.chroms[c] for c in self.j_chrom.tolist()), self.j_start.tolist(),
                        self.j_end.tolist()))
        del st, mo, ch
        strand = self.j_strand
        # Exact lookup and splice-site maps. A site key is (chrom, strand, position);
        # the donor is the intron's 5' end on its strand. Unstranded junctions (".")
        # are indexed under "." and offered to queries on either strand.
        self._key = {k: j for j, k in enumerate(keys)}
        self._donor: dict = {}
        self._acceptor: dict = {}
        for j, (c, s, e) in enumerate(keys):
            st = STAR_STRAND[int(strand[j])]
            d, a = (s, e) if st != "-" else (e, s)
            self._donor.setdefault((c, st, d), []).append(j)
            self._acceptor.setdefault((c, st, a), []).append(j)
        # Range index per chromosome: starts sorted (keys already sorted per chrom).
        self._by_chrom = {}
        for ci, c in enumerate(self.chroms):
            sel = np.flatnonzero(self.j_chrom == ci)
            self._by_chrom[c] = (self.j_start[sel], sel)
        self.max_length = int((self.j_end - self.j_start).max()) if J else 0
        # Alignment limits that bound what could have been observed.
        mx = [lib.max_intron for lib in loaded]
        mn = [lib.min_intron for lib in loaded]
        nc = [lib.noncanonical_filtered for lib in loaded]
        self.max_intron_any = None if (not loaded or any(v is None for v in mx)) else max(mx)
        self.min_intron_all = None if (not loaded or any(v is None for v in mn)) else min(mn)
        self.noncanonical_filtered_all = bool(loaded) and all(v is True for v in nc)
        self.stats.update(
            libraries_in_manifest=len(self.libraries), libraries_selected=len(libs),
            libraries_loaded=len(loaded),
            libraries_by_status=_count(lib.status for lib in self.libraries),
            junctions=J, observations=n_obs, groups=self.n_groups_total,
            studies=self.n_studies_total, junctions_inconsistent_strand_or_motif=inconsistent,
            count_units=sorted({lib.count_unit for lib in loaded}),
            max_intron_any_library=self.max_intron_any, min_intron_all_libraries=self.min_intron_all,
            noncanonical_filtered_all_libraries=self.noncanonical_filtered_all,
            approx_index_bytes=int(sum(a.nbytes for a in (
                self.lib_ptr, self.lib_idx, self.lib_unique, self.lib_multi, self.lib_overhang,
                self.j_chrom, self.j_start, self.j_end, self.j_strand, self.j_motif, self.j_unique,
                self.j_multi, self.j_overhang, self.j_nlib, self.j_nlib_any, self.j_ngroup,
                self.j_nstudy))),
        )

    # -------------------------------------------------------------- queries
    def __len__(self) -> int:
        return len(self.j_start)

    def support(self, j: int) -> JunctionSupport:
        return JunctionSupport(
            jid=int(j), chrom=self.chroms[int(self.j_chrom[j])], start0=int(self.j_start[j]),
            end=int(self.j_end[j]), strand=STAR_STRAND[int(self.j_strand[j])],
            motif=STAR_MOTIF[int(self.j_motif[j])], unique_reads=int(self.j_unique[j]),
            multi_reads=int(self.j_multi[j]), n_libraries=int(self.j_nlib[j]),
            n_libraries_any=int(self.j_nlib_any[j]), n_groups=int(self.j_ngroup[j]),
            n_studies=int(self.j_nstudy[j]), max_overhang=int(self.j_overhang[j]))

    def find(self, chrom: str, start0: int, end: int) -> int | None:
        """Junction id with exactly these half-open intron coordinates (any strand)."""
        return self._key.get((chrom, int(start0), int(end)))

    def lookup(self, chrom: str, start0: int, end: int) -> JunctionSupport | None:
        j = self.find(chrom, start0, end)
        return None if j is None else self.support(j)

    def site_partners(self, chrom: str, strand: str, start0: int, end: int) -> dict:
        """Observed junctions sharing this intron's donor or acceptor on *strand*
        (unstranded junctions included), excluding the intron itself.
        Returns ``{"donor": [jid...], "acceptor": [jid...]}``."""
        d, a = (start0, end) if strand != "-" else (end, start0)
        me = self.find(chrom, start0, end)
        out = {}
        for side, table, pos in (("donor", self._donor, d), ("acceptor", self._acceptor, a)):
            ids = list(table.get((chrom, strand, pos), ()))
            if strand in ("+", "-"):
                # Unstranded junctions are stored in + orientation (donor = start), so
                # a minus-strand query reads them with donor and acceptor swapped.
                flat = table if strand == "+" else (self._acceptor if side == "donor" else self._donor)
                ids += flat.get((chrom, ".", pos), ())
            out[side] = sorted(j for j in set(ids) if j != me)
        return out

    def overlapping(self, chrom: str, lo: int, hi: int):
        """Junction ids whose intron overlaps the half-open interval [lo, hi)."""
        entry = self._by_chrom.get(chrom)
        if entry is None:
            return []
        starts, ids = entry
        i0 = bisect.bisect_left(starts, lo - self.max_length)
        i1 = bisect.bisect_left(starts, hi)
        return [int(ids[i]) for i in range(i0, i1) if self.j_end[ids[i]] > lo]

    def library_counts(self, j: int) -> list[tuple[str, int, int]]:
        """(library_id, unique, multi) for every library observing junction *j*."""
        a, b = int(self.lib_ptr[j]), int(self.lib_ptr[j + 1])
        return [(self.loaded[int(self.lib_idx[x])].library_id, int(self.lib_unique[x]),
                 int(self.lib_multi[x])) for x in range(a, b)]

    def supporting_libraries(self, j: int) -> set[int]:
        a, b = int(self.lib_ptr[j]), int(self.lib_ptr[j + 1])
        return {int(self.lib_idx[x]) for x in range(a, b) if self.lib_unique[x] > 0}

    def library_table(self) -> list[dict]:
        return [lib.as_row() for lib in self.libraries]


def _count(values) -> dict:
    out: dict = {}
    for v in values:
        out[v] = out.get(v, 0) + 1
    return dict(sorted(out.items()))


def preflight_checks(manifest: str, roles, seq_lengths: dict | None) -> list[tuple]:
    """Checks for ``gmb-preflight --star-junctions``: ``[(name, verdict, message, detail)]``.

    FAIL when the manifest is unusable or no selected library loads; WARN for
    individual failed libraries and for facts that limit interpretation."""
    try:
        idx = JunctionIndex.from_manifest(manifest, roles, seq_lengths, checksum=False)
    except JunctionManifestError as exc:
        return [("star_junction_manifest", "FAIL", str(exc), {})]
    st = idx.stats
    out = []
    if st["libraries_selected"] == 0:
        return [("star_junction_manifest", "FAIL",
                 f"no library has role(s) {','.join(roles)}", {"by_status": st["libraries_by_status"]})]
    bad = [lib for lib in idx.libraries if lib.status not in (LOADED, NOT_SELECTED)]
    if not st["libraries_loaded"]:
        out.append(("star_junction_libraries", "FAIL",
                    f"none of {st['libraries_selected']} selected libraries loaded",
                    {"by_status": st["libraries_by_status"]}))
    elif bad:
        out.append(("star_junction_libraries", "WARN",
                    f"{len(bad)} selected librar(ies) not loaded and ignored: "
                    + ", ".join(f"{lib.library_id} {lib.status}" for lib in bad[:5]),
                    {"by_status": st["libraries_by_status"]}))
    else:
        out.append(("star_junction_libraries", "PASS",
                     f"{st['libraries_loaded']} librar(ies), {st['junctions']:,} junctions, "
                     f"{st['groups']} group(s), {st['studies']} stud(ies); "
                     f"{st['libraries_by_status'].get(NOT_SELECTED, 0)} not selected", {}))
    unknown = sum(lib.n_rows_unknown_seqid for lib in idx.loaded)
    if unknown:
        out.append(("star_junction_seqids", "WARN",
                    f"{unknown:,} junction row(s) on sequences absent from the genome were skipped", {}))
    if st["count_units"] != ["fragment"]:
        out.append(("star_junction_count_unit", "WARN",
                    f"count unit(s) {st['count_units']}: read counts are not fragment counts", {}))
    if st["max_intron_any_library"] is None:
        out.append(("star_junction_limits", "WARN",
                    "max_intron not given for every library: absence of long introns cannot be "
                    "interpreted", {}))
    return out
