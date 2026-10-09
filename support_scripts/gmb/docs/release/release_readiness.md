# GMB 2.0.0 — release-readiness report

Branch `gmb/release-prep`, 2026-10-09. Nothing has been committed, pushed, merged or tagged;
all changes are in the working tree for review. Detail: [`z_tritici_validation.md`](z_tritici_validation.md)
(fungal evidence) and [`../known_issues.md`](../known_issues.md) (open issues by severity).

## Decisions

| decision | recommendation |
|---|---|
| **Start a production-scale fungal pipeline test** | **GO** |
| **Create a stable tagged GMB release** | **NO-GO now** — becomes GO when the conditions below are met |

**Why GO for the test.** The production path installs cleanly as a wheel (four runtime
dependencies) and runs preflight → build → finalise from any directory with no local files.
It passes 862 tests, and its input checks now catch the problems the real *Z. tritici* bundle
contains. The fungal configuration does what its documentation says, and on 22% of the genome
it keeps exact-CDS agreement while removing most of the structural damage the previous
default did (gene merges −47%, splits −32%). Run it with:

1. the Helixer GFF3 renamed to genome sequence names (`tools/remap_helixer.py`);
2. no raw per-read long-read track (collapse it, or omit long reads);
3. `--preset fungi` on both preflight and build; protein validation as you intend to ship it;
4. timing and RSS recorded;
5. comparisons with `annotation-qc` of the output (all isoforms *and* canonical-only) against
   both the reference and Helixer alone, reported separately for core and accessory chromosomes.

**Why not yet a stable tag.** Conditions, in order:

1. the production-scale test confirms the 2.0.0 fungal behaviour genome-wide (it has been
   validated on nine windows only);
2. a decision on alternate-isoform policy for fungi (alternates are the source of nearly all
   remaining splits and merges and of GMB's whole gain over Helixer);
3. the apicomplexan presets are re-validated, or the tag is explicitly scoped to fungi
   (correctness fixes in 2.0.0 change their output somewhat);
4. CI lint green (failing before this work; no new failures added);
5. one networked clean install (`pip install`) with dependency resolution.

## 1. Whole-genome Helixer vs reference (summary)

R = 10,931 reference protein-coding genes; Q = 14,757 Helixer genes.

| | all | core 1–13 | accessory 14–21 |
|---|---|---|---|
| locus recovered (same strand) | 96.7% | 97.5% | 84.3% |
| CDS coordinate-exact | 36.2% | 37.6% | 13.5% |
| CDS intron chain (multi-CDS genes) | 42.0% | 43.9% | 10.1% |
| missed entirely | 260 (2.4%) | 182 (1.8%) | 78 (11.9%) |
| Helixer genes without reference counterpart | 3,279 (22.2%) | 19.1% | 49.6% |

The remaining recovered loci are incomplete rather than missing (51.9% overlapped with a
different structure). The 500 kb test region matched the core on locus recovery but sits below
it on exact CDS and cannot represent the accessory chromosomes.

## 2. Evidence and potentially novel genes

Missed reference genes are short (median CDS 269 bp) and 75% have RNA-seq support: a backbone
limitation, not an evidence gap. Of the core Helixer genes absent from the reference, 81.8%
have both protein and RNA-seq overlap and 75.3% have protein covering ≥ 50% of the CDS —
plausible unannotated genes, though some may be transposon-derived. Accessory novel genes are
weakly supported (66.4% lack any strict support).

## 3. Clustering mode — recommendation

**`transcript_linked` for fungi, together with the now-working 20 kb span filter and
detached-isoform removal; `exon_overlap` stays the default for `standard` and `apicomplexa`.**
Evidence (core windows): the pair keeps exact CDS (770 vs 771) and cuts merges 465 → 245 and
splits 484 → 327; the span filter *with* `exon_overlap` lowers exact CDS (761) because the
removed chimeras were bridging introns. The choice rests on structural correctness, not on
reference gains. Risks: regional evidence only; > 20 kb short-read transcripts are dropped
(1 reference gene exceeds 20 kb); ~11% more selection work per locus. The alternative keeps
fragment scoring of 5.7% of multi-exon candidates.

## 4. Configuration

| key(s) | status |
|---|---|
| `transcriptomic_filter.max_transcript_length` | **implemented** (fungi 20 kb applied; apicomplexa `null` — its 35 kb would remove 18% of *T. gondii* short-read models unvalidated) |
| `transcriptomic_filter.allow_single_exon` | **implemented** (default true = unchanged) |
| `scoring.locus_clustering` | new; fungi `transcript_linked` |
| 15 obsolete keys (`qc.*`, `export.*`, `reporting.formats`, `orf.stop_codon_char`, `orf.partial_prefix`, `protein_filter.min_exon_count_for_short`, `transcriptomic_filter.strand_consistency_check`, `transcript_splitting.split_on_*`) | **deprecated**: accepted, `FutureWarning` when set |
| 6 reserved keys (`orf.allow_partial_5/3`, `orf.allow_non_atg_start`, `utr.min_protein_coding_score_for_utr`, `utr.max_end_extension_bp`, `interpro_resolver.min_coverage_delta_for_replacement`) | **unsupported**: accepted, `UserWarning` when set; fixed behaviour documented |

Shipped presets and example configs set none of the 21; `resolved_config.yaml` reloads do
not warn. Every other fungal setting was traced to the code that reads it. Still open:
hard-coded −5.0 protein-validation penalty; `max_exon_len_reference: reference` unimplemented.

## 5. Handover interfaces

| check | result |
|---|---|
| production modules reachable from CLI / API | yes; every module reachable from a documented entry point; nothing required was removed |
| external comparison needed to build | no; no runtime or test dependency on ensembl-genes; `gmb-compare` prints where comparison moved |
| preflight catches bad inputs | missing, empty or non-GTF/GFF3 files, impossible coordinates, exon-less backbone/transcript tracks, unknown sequence names (e.g. accession-named Helixer), unstranded tracks: FAIL; CDS without exons, CDS-less backbone, per-read long-read tracks: WARN. Real bundle: 0 false alarms |
| formats, strand, coordinates, provenance | documented in `input_contract.md` / `output_contract.md` |
| stability | API signature, CLI flags and `evidence_attribution.tsv` columns pinned by tests; summary counters additive |
| clean non-editable install | wheel in a fresh venv, run from an unrelated directory: preflight, build, finalise succeed; output byte-identical to the editable install |
| documented fungal command | corrected (Helixer rename step); verified end to end on chr21 |
| local files or paths | none in the package |

Known gap: seqname mapping exists only in `gmb-build`, so accession-named backbones must be
renamed upstream (documented; P2).

## 6. Tests

| suite | baseline (start of work) | now |
|---|---|---|
| GMB `tests/` | 789 passed, 15 skipped | **862 passed, 3 skipped** |
| transcript strand (`tests/test_transcript_strand.py`) | 25 passed | 25 passed |
| clade skill | — | 28 passed |

`tests/test_release_prep_regressions.py` (101 tests) covers: imports and dependencies; the
`gmb-compare` stub; one default preset; UTR end support; clustering (real-data fragmentation,
whole-candidate scoring, nested genes, the reference-exact model it recovers); the span and
single-exon filters; detached-isoform removal; inert-key warnings and the documented set;
preflight malformed-input and long-read checks; the remap tool; API/CLI/attribution
stability; generic-backbone in-process build; fixture builds in both modes with a
gene-overlap invariant. The golden fixture was regenerated for the intended 2.0.0 behaviour
(186 genes; exact CDS 33 → 34, merges 25 → 14, splits 29 → 21 against the reference). The 3
skipped tests need a whole-genome data directory that is not in the repository. Not
automated: genome-scale runtime/memory, DIAMOND/Psauron runs, InterProScan.

## 7. Remaining issues

P1: whole-genome confirmation of 2.0.0 fungal behaviour; alternate-isoform policy;
apicomplexan re-validation. P2: seqname mapping only in build; genome-scale runtime unmeasured;
accessory-chromosome support; near-saturated protein support and possible TE-derived models;
hard-coded penalty; non-independent short-read tracks; input issues GMB does not repair. P3:
lint, default-preset fallback, `builder.py` size. Detail and evidence: `known_issues.md`.

## 8. Changes and proposed commits

Public interface changes: `gmb-visualize`, `gmb.compare`, `gmb.pipeline.fasta_export`,
`gmb.pipeline.reporting` removed; `gmb-compare` exits 2 with a pointer; `gmb-preflight` default
preset `standard` → `fungi` (same as `gmb-build`); new config keys
`scoring.locus_clustering`, `preflight.max_longread_models_per_mb`; `load_config(...,
warn_inert=True)`; new preflight checks and summary counters; `matplotlib` and `biopython`
no longer installed; fungi preset behaviour changed as above; apicomplexa
`max_transcript_length` `35000` → `null` (behaviour unchanged).

Suggested reviewable commits, in order. Only the complete set has been tested; run the suite
at each commit when splitting:

1. **Remove standalone comparison code** — `gmb/compare/`, `cli/visualize.py`, `cli/compare.py`
   stub, dead `fasta_export.py` / `reporting.py`, their tests, `pyproject.toml` /
   `requirements*.txt`.
2. **Builder correctness and performance** — UTR end-support index and same-sequence rule,
   discarded selection pass, stable sort, `--list-presets`, shared default preset
   (`builder.py`, `config.py`, `cli/preflight.py`).
3. **Locus clustering option** — `cluster_candidate_loci`, `scoring.locus_clustering`, counter.
4. **Fungal chimera handling** — span and single-exon filters (`evidence_filter.py`),
   `drop_detached_isoforms` (`gff3_validate.py`), preset changes (`fungi.yaml`,
   `apicomplexa.yaml`, `configs/*.yaml`), regenerated golden fixture, updated
   `test_config.py`.
5. **Inert configuration keys** — `INERT_CONFIG_KEYS` and warnings, `standard.yaml` notes,
   `finalise.py` reload.
6. **Preflight input checks** — long-read collapse and malformed-input checks (`checks.py`,
   `config.py`).
7. **Tools and examples** — `tools/remap_helixer.py`, example scripts, root `.gitignore`.
8. **Tests** — `test_annotate_cds_utrs.py` on the fixture, `test_release_prep_regressions.py`.
   (Alternatively fold each test group into commits 1–7.)
9. **Documentation** — README, HANDOVER, `docs/`, skill docs, release reports.
