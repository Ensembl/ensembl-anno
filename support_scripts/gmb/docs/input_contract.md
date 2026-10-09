# GMB input contract

What an upstream pipeline must hand to Gene Model Builder.

Every input is a file path. GMB reads them; it does not produce them (the one exception is
`gmb-longread-consensus`, which collapses long-read alignments another tool produced).

---

## Summary

| input | required | format | role | CLI flag |
|---|---|---|---|---|
| genome FASTA | **yes** | FASTA | coordinate system | `--genome` |
| ab initio backbone | **yes** in practice | GFF3 or GTF | `backbone` | `--helixer` / `--tiberius` / `--backbone` (any predictor) |
| short-read transcript models | recommended | GTF | `short_read_transcriptomic` | `--scallop`, `--stringtie` |
| long-read transcript models | optional | GTF | `long_read_transcriptomic` | `--minimap2` |
| protein-to-genome alignments | recommended | GTF | `protein_alignment` | `--orthodb`, `--uniprot`, `--genblast` |
| DIAMOND protein DB | optional | `.dmnd` | `protein_validation` | config only |
| reference annotation | **never a production input** | GFF3 | evaluation only | `annotation-qc pairwise-compare --reference` (ensembl-genes) |
| seqname map | **never a production input** | TSV | evaluation only | `annotation-qc pairwise-compare --seqname-map` |

The CLI flags are named after the tools that historically produced each file, but **the flag
does not determine behaviour** — the evidence *role* does. Passing an IsoQuant GTF to
`--minimap2` is fine as long as `scoring.longread_label` names the label GMB will apply.

---

## Rules that apply to every input

### Coordinate convention

GFF3/GTF as specified: **1-based, inclusive** start and end. GMB converts internally to
0-based half-open (pyranges) and converts back on output. You never see the internal form.

### Sequence names are authoritative from the genome FASTA

Every evidence file must use the **same seqids as the genome FASTA**. Preflight fails a track
where more than `preflight.max_unknown_seqid_fraction` (default 5%) of its sequence names are
absent from the genome. Remap upstream — GMB does not guess.

Helixer and other tools run on an NCBI download name sequences by GenBank accession
(`CM001196.1`), while Ensembl genomes use names such as `1`. Rename with
`tools/remap_helixer.py --assembly-report <NCBI assembly report>`, which maps the
GenBank-Accn column to the Assigned-Molecule column (including `##sequence-region` headers)
and fails on any sequence the report does not cover. `gmb-build --assembly-report` /
`--seqname-map` can rename at build time, but `gmb-preflight` and `run_gene_model_builder`
cannot, so renamed files are the supported route.

### Strand

`+` or `-`. Features with `.` or missing strand are **excluded from selection**; preflight
warns below 5% and fails above it.

### Transcript identity

Exons are grouped into transcripts by:
- **GTF**: `transcript_id`
- **GFF3**: `Parent` on the exon, which must name an `mRNA`/`transcript`

IDs need only be unique **within a file**; GMB prefixes every ID with its source label. IDs
reused **across sequences** in the same file are namespaced automatically with their seqid.

> **Why that matters.** A predictor that numbers genes per sequence (Tiberius emits `g1`,
> `g2`… restarting on every contig) produces the same `transcript_id` on many sequences.
> Grouping by `transcript_id` alone then fuses unrelated genes into chimeras — in one observed
> case across opposite strands of different chromosomes, with a fabricated 40–96 kb intron and
> internal stop codons starting exactly at the junction. GMB now namespaces these; preflight
> still reports the count, because a large number tells you something real about the upstream
> tool.

### Failure behaviour

| situation | behaviour |
|---|---|
| file missing/unreadable | preflight FAIL; `gmb-build` raises |
| empty file, or no exon/CDS rows with a transcript ID (not GTF/GFF3) | preflight FAIL (`file_readable`) |
| start > end, start < 1, or end beyond the sequence | preflight FAIL (`coordinates_valid`) |
| backbone / transcript track with no `exon` rows | preflight FAIL (`exon_rows_present`); `gmb-build` builds candidates from exon rows only |
| transcripts with CDS rows but no exon rows | preflight WARN; ignored by `gmb-build` |
| backbone with no CDS | preflight WARN; ORFs are inferred from exons instead |
| seqids not in genome | preflight FAIL above 5%, WARN below; unmatched models are unusable |
| unstranded features | preflight FAIL above 5%, WARN below; those features are excluded |
| one ID used for several rows on the same sequence | grouped into one model: exact repeats duplicate its exons, distinct alignments fuse into one span (protein spans > `protein_filter.max_span_bp` are then dropped) |
| IDs reused across seqids | namespaced automatically; preflight WARNs |
| poor splice quality | preflight WARN/FAIL by role (see below) |
| per-read long-read alignments | preflight WARN (`longread_collapsed`) |
| source in no evidence role | preflight **FAIL** — it would get the unknown-source weight |

---

## Per-input detail

### Genome FASTA — required

Plain or gzipped. Defines the coordinate system and supplies the sequence for ORF inference,
splice-site checking and FASTA regeneration. Sequence names are taken up to the first
whitespace.

### Ab initio backbone — required for the standard workflow

`--helixer` (GFF3), `--tiberius` (GTF) or the generic `--backbone` (GTF or GFF3, any
predictor; label from its source column or `--backbone-label`). **Mutually exclusive:
exactly one backbone.** It must have `exon` rows; its `CDS` rows are used as the coding
sequence where present (Helixer: all models, stop codon included in the CDS).

Whichever flag you use sets `scoring.backbone_label`, which is how the backbone role is
resolved. GMB will run without a backbone, but the standard production workflow expects one
and preflight warns.

The backbone's structural quality drives the `backbone_intron_rescue` applicability gate, so
it is measured and reported by preflight.

### Short-read transcript models — 0..N, strongly recommended

GTF with `exon` features and `transcript_id`. Any assembler; list its label under
`scoring.shortread_labels`. Two flags ship (`--scallop`, `--stringtie`) but the role, not the
flag, determines behaviour.

Assembled transcripts are filtered before scoring: any intron longer than
`transcriptomic_filter.max_intron_length`, or a genomic span longer than
`transcriptomic_filter.max_transcript_length` (fungi: 3 kb and 20 kb), removes the
transcript. In compact genomes read-through assemblies join neighbouring genes across short
intergenic gaps with no long intron, so the span rule is the one that catches them.

### Long-read transcript models — optional

Supply **collapsed consensus models**, not raw read alignments. `gmb-longread-consensus` does
the collapse if you have alignments. GMB scores every model as a candidate, so a per-read
track (2.05 M reads for *Z. tritici*) makes a build impractical; preflight warns
(`longread_collapsed`) above `preflight.max_longread_models_per_mb`.

**Strand is taken as given.** GMB uses each transcript's GTF strand; it never re-derives
strand from reads. Upstream, the Ensembl anno Minimap2 path sets strand from Minimap2's `ts`
tag combined with the alignment orientation
(`src/python/ensembl/tools/anno/transcriptomic_annotation/transcript_strand.py`); the SAM
FLAG alone is not a transcript strand for an unstranded library. Preflight's splice check is
computed on the declared strand, so a mis-stranded track shows up as a low canonical fraction.

Long-read evidence is genuinely optional: with none supplied there is no special case, no
penalty, and the long-read guard cannot fire.

> **Splice quality is checked.** A long-read track asserts it *observed* the splice structure,
> so it is judged most strictly (`preflight.longread_splice_fail`, default 0.85). A
> P. falciparum long-read consensus measured **16.4% canonical** splice sites against 99–100%
> for every other track; locus-by-locus it worsened 233 models and improved 228. If your track
> fails, set `scoring.longread_disposition` to `support_only` or `reject`.

### Protein-to-genome alignments — 0..N, recommended

GTF, typically CDS features. **Protein alignments are support and veto evidence only. They
never become candidate gene structures**, in any configuration. They can raise a candidate's
score, satisfy a retention gate, and appear in attribution.

Because they are never candidates, no ORF is ever inferred for them. CDS-only GTFs are
accepted for this role (and only this role).

### DIAMOND protein database — optional

Used with Psauron to score models after they exist. **Never contributes structure**: it can
change which isoform wins, not which structures are possible. Set
`protein_validation.enabled: false` if unavailable.

GMB does **not** install DIAMOND or Psauron; you supply absolute paths.

### Reference annotation — NEVER a production input

> `gmb-preflight`, `gmb-build` and `gmb-finalise` expose **no option** that accepts a reference
> annotation, and no resolved config can reference one. The reference enters only through
> external evaluation (`annotation-qc` in ensembl-genes), after the annotation exists.

Using a reference to choose a production configuration is not validation. Freeze the config
first, then evaluate.

---

## Minimum viable bundle

```
genome.fa            required
backbone.gff3        required in practice
```

## Recommended production bundle

```
genome.fa
backbone.gff3        Helixer or Tiberius
scallop.gtf          short-read assembly
stringtie.gtf        short-read assembly
orthodb.gtf          protein alignments
proteins.dmnd        protein validation
```
