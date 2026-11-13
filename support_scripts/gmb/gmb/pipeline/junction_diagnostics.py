"""Per-intron and per-model splice-junction diagnostics.

Describes what a :class:`gmb.pipeline.star_junctions.JunctionIndex` says about
the introns of any set of structures -- Helixer models, transcript-derived
candidates, GMB output, a reference annotation or future reconstructed models --
with one representation, so they can be compared directly.

These are *diagnostics*. They do not select or change models, and agreement with
junctions that also fed the transcript assemblies is not independent validation.

Intron status (``status`` / ``evidence_class``)
----------------------------------------------
``SUPPORTED``                    observed with >= 1 unique read (``supported``)
``MULTIMAP_ONLY``                observed, multimapping reads only (``weak``)
``OPPOSITE_STRAND``              these coordinates were observed with the other
                                 strand's motif (``contradicted``)
``ALTERNATIVE_SITE_OBSERVED``    not observed, but another junction sharing its
                                 donor or acceptor is (``contradicted``)
``UNOBSERVABLE_LONG``            longer than every library's aligner maximum
``UNOBSERVABLE_SHORT``           shorter than every library's aligner minimum
``UNOBSERVABLE_NONCANONICAL``    non-canonical motif and every library filtered
                                 non-canonical junctions (all ``unobservable``)
``UNOBSERVED_EXPRESSED_LOCUS``   not observed; other junctions overlapping the
                                 model's CDS span are observed on its strand (``unknown``)
``UNOBSERVED_NO_LOCAL_EVIDENCE`` not observed; no junction overlapping the model's CDS span
                                 (``unknown``: not expressed in these libraries?)

Absence is never reported as contradiction: only an observed competing splice
(alternative site, opposite strand) is. Precedence: opposite strand > aligner length
limits > alternative site > non-canonical > local expression.
"""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass, field

from gmb.pipeline.star_junctions import JunctionIndex

CANONICAL_MOTIFS = {
    "+": {("GT", "AG"), ("GC", "AG"), ("AT", "AC")},
    "-": {("CT", "AC"), ("CT", "GC"), ("GT", "AT")},
}

INTRON_COLUMNS = [
    "model_id",
    "source",
    "gene_id",
    "chrom",
    "strand",
    "intron_number",
    "intron_start",
    "intron_end",
    "length",
    "genome_motif",
    "status",
    "evidence_class",
    "unique_reads",
    "multi_reads",
    "n_libraries",
    "n_libraries_any",
    "n_groups",
    "n_studies",
    "library_fraction",
    "group_fraction",
    "max_overhang",
    "star_strand",
    "star_motif",
    "donor_usage",
    "acceptor_usage",
    "n_alternative_junctions",
    "best_alternative",
    "best_alternative_side",
    "best_alternative_unique_reads",
    "best_alternative_n_libraries",
    "best_alternative_n_groups",
    "site_conflict",
    "reasons",
]
MODEL_COLUMNS = [
    "model_id",
    "source",
    "gene_id",
    "chrom",
    "strand",
    "start",
    "end",
    "n_cds_segments",
    "n_cds_introns",
    "n_supported",
    "n_multimap_only",
    "n_contradicted",
    "n_unknown",
    "n_unobservable",
    "n_minor_site",
    "min_unique_reads",
    "min_n_libraries",
    "min_n_groups",
    "libraries_supporting_all_introns",
    "groups_supporting_all_introns",
    "observed_junctions_in_span",
    "model_class",
    "unsupported_introns",
    "contradicted_introns",
    "extra",
]


@dataclass
class ModelSpec:
    """One structure to annotate. ``segments`` are half-open 0-based (pyranges)."""

    model_id: str
    source: str
    chrom: str
    strand: str
    segments: list
    gene_id: str = ""
    extra: dict = field(default_factory=dict)


def introns_of(segments) -> list[tuple[int, int]]:
    s = sorted(segments)
    return [(s[i][1], s[i + 1][0]) for i in range(len(s) - 1) if s[i + 1][0] > s[i][1]]


def genome_motif(seq: str | None, start0: int, end: int, strand: str) -> str:
    """Donor/acceptor dinucleotides read on *strand*, e.g. ``GT/AG``; ``""`` if unknown."""
    if not seq or end - start0 < 4 or end > len(seq):
        return ""
    a, b = seq[start0 : start0 + 2].upper(), seq[end - 2 : end].upper()
    if strand == "-":
        comp = str.maketrans("ACGTN", "TGCAN")
        a, b = b.translate(comp)[::-1], a.translate(comp)[::-1]
    return f"{a}/{b}"


def _is_canonical(motif: str) -> bool | None:
    if not motif:
        return None
    return tuple(motif.split("/")) in CANONICAL_MOTIFS["+"]


def _fmt(x) -> str:
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.4f}"
    return str(x)


class JunctionAnnotator:
    """Annotates models against one index. Reuse it for many model sets."""

    def __init__(self, index: JunctionIndex, genome: dict | None = None):
        self.index = index
        self.genome = genome or {}
        self.n_lib = max(1, len(index.loaded))
        self.n_grp = max(1, index.n_groups_total)

    # ------------------------------------------------------------ one intron
    def intron(self, chrom: str, strand: str, start0: int, end: int, span=None) -> dict:
        ix = self.index
        length = end - start0
        motif = genome_motif(self.genome.get(chrom), start0, end, strand)
        row = {
            "intron_start": start0 + 1,
            "intron_end": end,
            "length": length,
            "genome_motif": motif,
            "unique_reads": 0,
            "multi_reads": 0,
            "n_libraries": 0,
            "n_libraries_any": 0,
            "n_groups": 0,
            "n_studies": 0,
            "max_overhang": "",
            "star_strand": "",
            "star_motif": "",
            "donor_usage": None,
            "acceptor_usage": None,
            "site_conflict": "",
            "reasons": [],
        }
        j = ix.find(chrom, start0, end)
        sup = ix.support(j) if j is not None else None
        same_strand = sup is not None and (sup.strand == strand or sup.strand == ".")
        partners = ix.site_partners(chrom, strand, start0, end)
        my_reads = sup.unique_reads if same_strand else 0
        best = None
        for side in ("donor", "acceptor"):
            others = [ix.support(p) for p in partners[side]]
            obs = [o for o in others if o.unique_reads > 0]
            total = my_reads + sum(o.unique_reads for o in obs)
            row[f"{side}_usage"] = (my_reads / total) if total else None
            for o in obs:
                if best is None or (o.unique_reads, -o.start0) > (
                    best[0].unique_reads,
                    -best[0].start0,
                ):
                    best = (o, side)
        n_alt = len(set(partners["donor"]) | set(partners["acceptor"]))
        row["n_alternative_junctions"] = n_alt
        if best:
            o, side = best
            row.update(
                best_alternative=f"{o.start1}-{o.end}",
                best_alternative_side=side,
                best_alternative_unique_reads=o.unique_reads,
                best_alternative_n_libraries=o.n_libraries,
                best_alternative_n_groups=o.n_groups,
            )
        if sup is not None:
            row.update(star_strand=sup.strand, star_motif=sup.motif, max_overhang=sup.max_overhang)
        if same_strand:
            row.update(
                unique_reads=sup.unique_reads,
                multi_reads=sup.multi_reads,
                n_libraries=sup.n_libraries,
                n_libraries_any=sup.n_libraries_any,
                n_groups=sup.n_groups,
                n_studies=sup.n_studies,
            )
            if sup.unique_reads > 0:
                status, cls = "SUPPORTED", "supported"
                if best and best[0].unique_reads > sup.unique_reads:
                    row["site_conflict"] = "minor_site"
                    row["reasons"].append(
                        f"a junction sharing the {best[1]} has more unique reads"
                    )
            else:
                status, cls = "MULTIMAP_ONLY", "weak"
                row["reasons"].append("only multimapping reads")
            if sup.strand == ".":
                row["reasons"].append("STAR strand undefined (non-canonical motif)")
        else:
            status, cls = self._absent(
                chrom, strand, start0, end, length, motif, best, sup, span, row
            )
        row["status"], row["evidence_class"] = status, cls
        row["library_fraction"] = row["n_libraries"] / self.n_lib
        row["group_fraction"] = row["n_groups"] / self.n_grp
        if (
            length > 0
            and ix.max_intron_any is not None
            and length > ix.max_intron_any
            and status not in ("UNOBSERVABLE_LONG",)
        ):
            row["reasons"].append("longer than the aligner maximum")
        row["reasons"] = "; ".join(row["reasons"])
        return row

    def _absent(self, chrom, strand, start0, end, length, motif, best, sup, span, row):
        # Order matters: an observed splice (same coordinates on the other strand,
        # or a competitor sharing a site) is positive evidence and outranks the
        # "could not be observed" rules -- except the aligner's length limits, so a
        # junction longer than every library's maximum is never penalised.
        ix = self.index
        if sup is not None:
            row["reasons"].append(
                f"observed only with STAR strand {sup.strand} (motif {sup.motif})"
            )
            return "OPPOSITE_STRAND", "contradicted"
        if ix.max_intron_any is not None and length > ix.max_intron_any:
            row["reasons"].append(
                f"{length} bp > aligner max_intron {ix.max_intron_any}: absence "
                "carries no information"
            )
            return "UNOBSERVABLE_LONG", "unobservable"
        if ix.min_intron_all is not None and length < ix.min_intron_all:
            row["reasons"].append(f"{length} bp < aligner min_intron {ix.min_intron_all}")
            return "UNOBSERVABLE_SHORT", "unobservable"
        noncanonical = ix.noncanonical_filtered_all and _is_canonical(motif) is False
        if best is not None:
            row["reasons"].append(
                f"{best[1]} shared with observed junction {best[0].start1}-{best[0].end}"
            )
            if noncanonical:
                row["reasons"].append(f"own motif {motif} is non-canonical")
            row["site_conflict"] = "alternative_dominant"
            return "ALTERNATIVE_SITE_OBSERVED", "contradicted"
        if noncanonical:
            row["reasons"].append(
                f"non-canonical motif {motif}; aligner filtered non-canonical junctions"
            )
            return "UNOBSERVABLE_NONCANONICAL", "unobservable"
        lo, hi = span if span else (start0, end)
        local = [
            j
            for j in ix.overlapping(chrom, lo, hi)
            if ix.j_unique[j] > 0
            and ix.j_strand[j] in ((1,) if strand == "+" else (2,) if strand == "-" else (1, 2))
        ]
        if local:
            row["reasons"].append(
                f"{len(local)} other junction(s) observed within the model; none shares a site"
            )
            return "UNOBSERVED_EXPRESSED_LOCUS", "unknown"
        row["reasons"].append("no junction observed within the model on this strand")
        return "UNOBSERVED_NO_LOCAL_EVIDENCE", "unknown"

    # ------------------------------------------------------------- one model
    def model(self, m: ModelSpec) -> tuple[list[dict], dict]:
        segs = sorted((int(s), int(e)) for s, e in m.segments)
        ints = introns_of(segs)
        lo, hi = (segs[0][0], segs[-1][1]) if segs else (0, 0)
        if m.strand == "-":
            ints = ints[::-1]  # number introns 5' -> 3'
        rows = []
        for k, (s, e) in enumerate(ints, start=1):
            r = self.intron(m.chrom, m.strand, s, e, span=(lo, hi))
            r.update(
                model_id=m.model_id,
                source=m.source,
                gene_id=m.gene_id,
                chrom=m.chrom,
                strand=m.strand,
                intron_number=k,
            )
            rows.append(r)
        ix = self.index
        cls = [r["evidence_class"] for r in rows]
        sup_rows = [r for r in rows if r["evidence_class"] == "supported"]
        libs = groups = None
        if rows and len(sup_rows) == len(rows):
            sets = [
                ix.supporting_libraries(ix.find(m.chrom, r["intron_start"] - 1, r["intron_end"]))
                for r in rows
            ]
            common = set.intersection(*sets)
            libs = len(common)
            groups = len({int(ix.lib_group_code[x]) for x in common})
        observed_in_span = sum(
            1
            for j in ix.overlapping(m.chrom, lo, hi)
            if ix.j_unique[j] > 0
            and STAR_STRAND_OK(ix.j_strand[j], m.strand)
            and ix.j_start[j] >= lo
            and ix.j_end[j] <= hi
        )
        if not rows:
            label = "no_cds_introns"
        elif "contradicted" in cls:
            label = "contradicted"
        elif len(sup_rows) == len(rows):
            label = "all_supported"
        elif sup_rows or "weak" in cls:
            label = "partially_supported"
        else:
            label = "no_junction_evidence"
        summary = {
            "model_id": m.model_id,
            "source": m.source,
            "gene_id": m.gene_id,
            "chrom": m.chrom,
            "strand": m.strand,
            "start": lo + 1,
            "end": hi,
            "n_cds_segments": len(segs),
            "n_cds_introns": len(rows),
            "n_supported": len(sup_rows),
            "n_multimap_only": cls.count("weak"),
            "n_contradicted": cls.count("contradicted"),
            "n_unknown": cls.count("unknown"),
            "n_unobservable": cls.count("unobservable"),
            "n_minor_site": sum(1 for r in rows if r["site_conflict"] == "minor_site"),
            "min_unique_reads": min((r["unique_reads"] for r in rows), default=""),
            "min_n_libraries": min((r["n_libraries"] for r in rows), default=""),
            "min_n_groups": min((r["n_groups"] for r in rows), default=""),
            "libraries_supporting_all_introns": "" if libs is None else libs,
            "groups_supporting_all_introns": "" if groups is None else groups,
            "observed_junctions_in_span": observed_in_span,
            "model_class": label,
            "unsupported_introns": ",".join(
                f"{r['intron_start']}-{r['intron_end']}:{r['status']}"
                for r in rows
                if r["evidence_class"] != "supported"
            ),
            "contradicted_introns": ",".join(
                f"{r['intron_start']}-{r['intron_end']}"
                for r in rows
                if r["evidence_class"] == "contradicted"
            ),
            "extra": ";".join(f"{k}={v}" for k, v in sorted(m.extra.items())),
        }
        return rows, summary

    def annotate(self, models) -> tuple[list[dict], list[dict]]:
        intron_rows, model_rows = [], []
        for m in sorted(
            models,
            key=lambda x: (
                x.chrom,
                min((s for s, _ in x.segments), default=0),
                x.source,
                x.model_id,
            ),
        ):
            r, s = self.model(m)
            intron_rows.extend(r)
            model_rows.append(s)
        return intron_rows, model_rows


def STAR_STRAND_OK(code, strand: str) -> bool:
    return code == 0 or (code == 1 and strand == "+") or (code == 2 and strand == "-")


def write_tsv(path: str, rows: list[dict], columns: list[str]) -> None:
    with open(path, "w", newline="") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(_fmt(r.get(c, "")) for c in columns) + "\n")


# ------------------------------------------------------------------ readers
_GTF_ATTR = re.compile(r'(\S+) "([^"]*)"')


def read_models(path: str, source: str | None = None, feature: str = "CDS") -> list[ModelSpec]:
    """Models from a GFF3 or GTF file: one per transcript, segments = its *feature*
    rows (``CDS`` by default; ``exon`` for transcript assemblies), converted to
    half-open 0-based. Transcripts without such rows are skipped."""
    gff3 = not path.endswith((".gtf", ".gtf.gz"))
    opener = __import__("gzip").open if path.endswith(".gz") else open
    segs: dict = {}
    meta: dict = {}
    tx_gene: dict = {}
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9:
                continue
            if gff3:
                a = dict(
                    kv.split("=", 1) for kv in p[8].strip().strip(";").split(";") if "=" in kv
                )
                if p[2] in ("mRNA", "transcript"):
                    tx_gene[a.get("ID", "")] = a.get("Parent", "")
                    continue
                if p[2] != feature:
                    continue
                parents = a.get("Parent", "").split(",")
            else:
                a = dict(_GTF_ATTR.findall(p[8]))
                if p[2] != feature:
                    continue
                parents = [a.get("transcript_id", "")]
                tx_gene.setdefault(parents[0], a.get("gene_id", ""))
            for t in parents:
                if not t:
                    continue
                segs.setdefault(t, []).append((int(p[3]) - 1, int(p[4])))
                meta.setdefault(t, (p[0], p[6], source or p[1]))
    return [
        ModelSpec(
            model_id=t,
            source=meta[t][2],
            chrom=meta[t][0],
            strand=meta[t][1],
            segments=sorted(s),
            gene_id=tx_gene.get(t, ""),
        )
        for t, s in segs.items()
    ]


def write_diagnostics(out_dir: str, prefix: str, annotator: JunctionAnnotator, models) -> dict:
    """Annotate *models* and write ``<prefix>_introns.tsv`` and ``<prefix>_models.tsv``."""
    os.makedirs(out_dir, exist_ok=True)
    introns, summaries = annotator.annotate(models)
    write_tsv(os.path.join(out_dir, f"{prefix}_introns.tsv"), introns, INTRON_COLUMNS)
    write_tsv(os.path.join(out_dir, f"{prefix}_models.tsv"), summaries, MODEL_COLUMNS)
    counts: dict = {}
    for s in summaries:
        counts[s["model_class"]] = counts.get(s["model_class"], 0) + 1
    by_status: dict = {}
    for r in introns:
        by_status[r["status"]] = by_status.get(r["status"], 0) + 1
    return {
        "models": len(summaries),
        "introns": len(introns),
        "model_class": dict(sorted(counts.items())),
        "intron_status": dict(sorted(by_status.items())),
    }


def write_library_table(path: str, index: JunctionIndex) -> None:
    rows = index.library_table()
    cols = list(rows[0]) if rows else ["library_id"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


# ------------------------------------------------------------ gmb-build hooks
class BuildJunctionDiagnostics:
    """What ``gmb-build --star-junctions`` does: load once, annotate every
    candidate and every selected transcript, write ``<output-dir>/junction_evidence/``.
    It never returns anything the builder uses for selection."""

    def __init__(self, manifest: str, roles, genome: dict, output_dir: str, index=None):
        import time

        self.t = {}
        t0 = time.perf_counter()
        # An index already loaded for evidence-score selection (same manifest and
        # roles) is reused instead of reading the files twice.
        self.used_for_selection = index is not None
        self.index = index if index is not None else JunctionIndex.from_manifest(
            manifest, roles, {k: len(v) for k, v in genome.items()}
        )
        self.t["load_index_s"] = round(time.perf_counter() - t0, 3)
        self.annotator = JunctionAnnotator(self.index, genome)
        self.out = os.path.join(output_dir, "junction_evidence")
        os.makedirs(self.out, exist_ok=True)
        write_library_table(os.path.join(self.out, "libraries.tsv"), self.index)
        self.reports = {}
        idx = self.index.stats
        print(
            f"  STAR junctions: {idx['libraries_loaded']}/{idx['libraries_selected']} libraries "
            f"(roles {','.join(roles)}), {idx['junctions']:,} junctions, "
            f"{idx['groups']} groups, {self.t['load_index_s']} s; "
            + ("shared with evidence-score selection" if self.used_for_selection else "diagnostics only")
        )
        if idx["libraries_loaded"] < idx["libraries_selected"]:
            print(
                f"  WARNING: junction libraries not loaded: {idx['libraries_by_status']} "
                f"(see {self.out}/libraries.tsv)"
            )

    def candidates(self, candidate_exons, candidate_cds: dict) -> None:
        import time

        t0 = time.perf_counter()
        meta = candidate_exons.drop_duplicates("transcript_id")[
            ["transcript_id", "Chromosome", "Strand", "Source", "gene_id"]
        ]
        models = [
            ModelSpec(
                model_id=str(t),
                source=str(src),
                chrom=str(c),
                strand=str(s),
                segments=list(candidate_cds[t]),
                gene_id=str(g),
            )
            for t, c, s, src, g in meta.itertuples(index=False, name=None)
            if candidate_cds.get(t)
        ]
        self.reports["candidates"] = write_diagnostics(
            self.out, "candidates", self.annotator, models
        )
        self.t["candidates_s"] = round(time.perf_counter() - t0, 3)

    def selected(self, gff3_path: str, output_mrnas: list) -> None:
        import time

        t0 = time.perf_counter()
        models = read_models(gff3_path, source="GMB", feature="CDS")
        info = {m["ID"]: m for m in output_mrnas}
        for m in models:
            row = info.get(m.model_id, {})
            m.extra = {
                "candidate_transcript_id": row.get("candidate_transcript_id", ""),
                "evidence": row.get("Evidence", ""),
            }
        self.reports["selected"] = write_diagnostics(self.out, "selected", self.annotator, models)
        self.t["selected_s"] = round(time.perf_counter() - t0, 3)

    def finish(self) -> None:
        import json
        import resource
        import sys

        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        from gmb.cli.junction_evidence import junction_caveats

        summary = {
            "index": self.index.stats,
            "timing": self.t,
            "reports": self.reports,
            "process_peak_rss_mb_at_finish": round(
                rss / (1 << 20) if sys.platform == "darwin" else rss / 1024, 1
            ),
            "notes": junction_caveats(self.index, self.used_for_selection),
        }
        with open(os.path.join(self.out, "junction_evidence_summary.json"), "w") as fh:
            json.dump(summary, fh, indent=2, sort_keys=True)
        print(f"  Junction diagnostics: {self.out}")
