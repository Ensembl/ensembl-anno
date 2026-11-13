# STAR splice-junction evidence (optional, diagnostics only)

GMB can read per-library STAR `SJ.out.tab` files and report, for every candidate and every
selected transcript, how strongly each CDS intron is supported. **In this release the
junctions never change which models GMB builds or selects.** With or without
`--star-junctions`, `consensus.gff3`, the FASTA files, `evidence_attribution.tsv`,
`summary.*` and `resolved_config.yaml` are byte-identical. Only the run manifest
records the extra command-line arguments.

```bash
gmb-build  ... --star-junctions junctions.tsv [--star-junction-roles construction]
gmb-run    ... --star-junctions junctions.tsv          # passed to preflight and build
gmb-preflight ... --star-junctions junctions.tsv       # validates the manifest and files
gmb-junction-evidence --star-junctions junctions.tsv --genome genome.fa \
    --annotation helixer=helixer.gff3 --annotation reference=ref.gff3 --output-dir report/
```

(`gmb-junction-evidence` is also `python -m gmb.cli.junction_evidence`.)

## Manifest

A tab-separated file with one row per library. Lines starting with `#` are comments.
Relative paths resolve against the manifest's directory.

| column | required | meaning |
|---|---|---|
| `library_id` | yes | unique library (run) identifier |
| `path` | yes | STAR `SJ.out.tab` (plain or gzip, detected by content; the name does not matter) |
| `role` | yes, unless `--star-junction-roles all` | evidence set, e.g. `construction`, `holdout`, `auxiliary`, `excluded`. **Only libraries with a selected role are opened.** The default selection is `construction` |
| `sample_id`, `study`, `strain`, `condition` | no | identifiers carried into reports |
| `group` | no | evidence-independence class (libraries from one biological condition). Default: `study\|strain\|condition`, or the library itself if all three are empty |
| `count_unit` | no | `mate` or `fragment`. Use `mate` when paired reads were aligned as independent reads |
| `max_intron`, `min_intron` | no | the aligner's intron-length limits (bp). Absence of an intron outside them for **every** library is reported as unobservable |
| `noncanonical_filtered` | no | `true` when the aligner discarded non-canonical junctions, so a non-canonical intron's absence carries no information |

Example (one Z. tritici construction library):

```
library_id	path	role	sample_id	study	strain	condition	count_unit	max_intron	min_intron	noncanonical_filtered
SRR17893325	SRR17893325.sj.tab	construction	SAMN25686266	PRJNA803733	IPO323	IPO323_FLUX	mate	1000	21	true
```

## How the files are read

* **Coordinates.** STAR columns 2–3 are the first and last intronic base, 1-based inclusive. They are stored half-open
  0-based, `(col2 - 1, col3)`, which equals `(previous_exon.End, next_exon.Start)` in GMB's
  internal (pyranges) coordinates. Reports print the intron 1-based inclusive again.
* **Strand** comes from column 4 (STAR derives it from the intron motif, not from the library). A
  junction with an undefined strand (0) is offered to queries on either strand and flagged.
* **Counts.** Columns 7/8 are unique and multimapping reads, in the library's `count_unit`. Nothing converts
  mates to fragments.
* **Column 6 ("annotated") is ignored.** In two-pass STAR it marks junctions found in pass 1, so it is
  not independent evidence.
* **Validation, per library.** A zero-byte file is `FAILED_EMPTY` (a failed job, *not* zero
  junctions). A missing file is `MISSING_FILE`. Any malformed line rejects the whole library as `MALFORMED`:
  wrong column count, a non-integer, a value out of range, an unterminated last line, a duplicate junction, or an
  end beyond the sequence. Junctions on sequences absent from the genome are skipped and counted; if none is on
  the genome, the status is `NO_KNOWN_SEQIDS`. Rejected libraries contribute nothing and are listed in
  `libraries.tsv`. Preflight warns about them.

## Index

`gmb.pipeline.star_junctions.JunctionIndex` loads the selected libraries once and supports:
* exact intron lookups;
* donor/acceptor partner lookups (shared 5′ or 3′ site on the intron's strand);
* range queries.

Per junction it keeps:
* unique and multimapping reads;
* supporting libraries (≥1 unique read);
* distinct groups and studies;
* maximum overhang, motif and strand.

Per-library counts stay in CSR arrays, so weighting or independence policies can be added later without
re-reading the files. The library count is **not** a confidence score: replicate libraries of one condition share a
`group`.

Measured on 73 Z. tritici libraries (1.67 M observations, 125,546 junctions):
* 2.3 s to load;
* a 46 MB index;
* about 3 µs per exact lookup plus site query.

## Outputs (`<build-dir>/junction_evidence/`)

| file | content |
|---|---|
| `libraries.tsv` | every manifest library with role, load status, junction count, skipped rows and checksum |
| `candidates_introns.tsv` / `candidates_models.tsv` | every candidate with a CDS (Helixer, StringTie, Scallop, long-read…), before selection |
| `selected_introns.tsv` / `selected_models.tsv` | every transcript in `consensus.gff3` (`extra` = candidate id and evidence) |
| `junction_evidence_summary.json` | index statistics, timings, peak memory, caveats |

**Per intron:**
* coordinates (1-based), length, and the motif read from the genome on the model's strand;
* `status` and `evidence_class` (below), unique/multimap reads, supporting libraries, groups and studies, and both as fractions;
* maximum overhang, STAR strand and motif;
* donor and acceptor usage (this intron's unique reads ÷ all junctions sharing that site);
* the strongest alternative junction;
* `site_conflict` (`minor_site` when a junction sharing a site has more reads, `alternative_dominant` when this intron is unobserved but a site partner is observed);
* reasons.

| status | class | meaning |
|---|---|---|
| `SUPPORTED` | supported | ≥1 unique read |
| `MULTIMAP_ONLY` | weak | multimapping reads only |
| `OPPOSITE_STRAND` | contradicted | these coordinates were observed with the other strand's motif |
| `ALTERNATIVE_SITE_OBSERVED` | contradicted | not observed, but a junction sharing its donor or acceptor is |
| `UNOBSERVABLE_LONG` / `UNOBSERVABLE_SHORT` | unobservable | outside every library's aligner intron limits |
| `UNOBSERVABLE_NONCANONICAL` | unobservable | non-canonical motif, and every library filtered such junctions |
| `UNOBSERVED_EXPRESSED_LOCUS` | unknown | not observed, while other junctions (unique reads, same strand) overlap the model's CDS span |
| `UNOBSERVED_NO_LOCAL_EVIDENCE` | unknown | no junction on its strand overlaps the model's CDS span |

Precedence runs from top to bottom: opposite strand, then length limits, then alternative site, then non-canonical, then local expression. **Absence alone is never
contradiction.**

**Per model:**
* CDS intron count, and the count per class (supported, multimap-only, contradicted, unknown, unobservable);
* minor-site introns;
* the minimum reads, libraries and groups across its introns;
* the libraries and groups that observe *every* intron;
* the number of junctions observed inside the model;
* `model_class`: `all_supported`, `partially_supported`, `contradicted`, `no_junction_evidence` or `no_cds_introns`;
* the unsupported and contradicted introns.

## Interpreting the diagnostics

* **They are not independent validation.** If the junctions come from the alignments that also built the
  StringTie/Scallop tracks, agreement between them is expected. Keep held-out libraries in a separate role
  and evaluate them with `--star-junction-roles holdout` in a separate run.
* **Missing evidence is not negative evidence.** An unexpressed gene, an intron beyond the aligner limit
  or a filtered motif leaves no junction.

## Use in selection

With `scoring.primary_selection: evidence_score` the same index (construction libraries only)
also chooses primaries; see [evidence_selection.md](evidence_selection.md). Every other mode
uses junctions for diagnostics only.
