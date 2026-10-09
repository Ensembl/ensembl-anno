# GMB known issues

Open problems, ordered by severity, then what was fixed for 2.0.0. Each entry gives the
evidence, the current mitigation and the recommended action. Last reviewed 2026-10-09
(whole-genome Helixer assessment of *Z. tritici*; see `release/z_tritici_validation.md`).

Severity: **P1** affects selected gene models in normal production runs or blocks a stable
release; **P2** affects a subset of runs, inputs or configurations; **P3** developer-facing.

---

## P1 — 2.0.0 fungal behaviour is validated regionally, not genome-wide

The fungi preset's effective behaviour changed in 2.0.0: the 20 kb transcript-span filter now
applies, loci are `transcript_linked`, and isoforms detached from their gene are removed. On
nine windows (8.8 Mb, 22% of the genome, six core and three accessory chromosomes) this
release candidate matched the previous default on exact CDS (core 770 vs 771 of 1,865
reference genes; accessory 29 vs 27 of 236) while cutting gene merges 47% (465 → 245) and
splits 32% (484 → 327), with 12% less build time. **Action:** the production-scale fungal run
is the genome-wide validation; tag only after it.

## P1 — Alternate isoforms drive most gene splits and merges

With one isoform per gene, GMB's core-window output equals Helixer's on every reference
metric (719 CDS-exact, 18 merges, 71 splits vs Helixer 69). All of GMB's measured gain over
Helixer comes from short-read **alternate isoforms**, which also cause the extra
splits/merges: in the release candidate 242 of the 255 merged genes have ≥ 2 isoforms and
226 of them are "span-only" (overlapping alternates chain the gene span into a neighbour).
The reference has one isoform per gene, so it cannot say whether those alternates are real.
**Action (decision, not code):** for the production test, evaluate both the all-isoform and
the canonical-only output; then decide `max_isoforms_per_locus` / alternate-admission rules
(`same_gene_overlap_threshold` 0.15 relative to the shorter model) for fungi.

## P1 — Apicomplexan presets need re-validation before a stable tag

Two 2.0.0 correctness fixes apply to every preset: UTR end support no longer accepts ends on
other sequences, and detached isoforms are removed. Apicomplexan outputs will therefore
change somewhat, unmeasured. `max_transcript_length` is set to `null` in the apicomplexa
preset because applying its documented 35 kb would remove 0.2% (*P. falciparum*), 7.8%
(GCA_000006355.3) and 18.3% (*T. gondii*) of surviving short-read models. **Action:** re-run
the apicomplexan validation (or explicitly scope the 2.0.0 tag to fungi).

## P2 — Seqname mapping exists only in `gmb-build`

`--assembly-report` / `--seqname-map` exist in `gmb-build` but not in `gmb-preflight`,
`gmb-finalise` or `run_gene_model_builder`. Helixer output for an NCBI download uses GenBank
accessions (`CM001196.1`), so the delivered *Z. tritici* backbone fails preflight
(21/21 names absent) until renamed. **Mitigation:** rename upstream with
`tools/remap_helixer.py` (documented in `input_contract.md`; now also renames
`##sequence-region` headers and fails on unmapped sequences). **Action:** either add the same
mapping options to preflight and the API, or deprecate build-time mapping.

## P2 — Genome-scale runtime not measured on 2.0.0

Recorded builds took 29 h (*Z. tritici*, with protein validation) and 47.9 h (*T. gondii*). 2.0.0
removed a discarded second selection pass and a per-model genome-wide UTR scan (1.32 s →
0.01 ms per model); regional builds were 12% faster than the previous default. Clustering
with `transcript_linked` costs 0.8 s genome-wide and raises the selection-cost proxy (Σk²) by
11%. **Action:** record wall time and RSS in the production-scale test.

## P2 — Evidence and predictions on accessory chromosomes (14–21)

Helixer recovers 84.3% of accessory reference loci (core 97.5%) and 13.5% CDS-exact (core
37.6%); half of accessory Helixer genes (757 / 1,526) are absent from the reference, and two
thirds of those have neither ≥ 50% protein coverage nor RNA-seq-supported introns. Expect
lower accuracy and more weakly supported calls there; report them separately.

## P2 — Positional protein support is near-saturated; TE-derived models unknown

Same-strand protein overlap covers 98.6% of reference-overlapping Helixer genes and 76.6% of
novel ones, so `protein_overlap_bonus` discriminates little. Protein homology can also come
from transposable-element proteins, which reference gene sets usually exclude; GMB has no
repeat input. **Action:** evaluate `protein_support_mode: cds_span_compatible`; check novel
models against a repeat annotation.

## P2 — Selection-affecting fixed values and unimplemented modes

- With `protein_validation.policy: penalize` (fungi), a model below `min_score` loses a
  hard-coded 5.0 points (more than the fungal backbone weight, 3.1). Only when protein
  validation is enabled. Make it configurable before enabling it in production.
- `validation.max_exon_len_reference: reference` is accepted but falls back to
  `candidates_supported` with a warning.

## P2 — Short-read tracks are not independent evidence

Both Ensembl anno short-read files for *Z. tritici* are StringTie-format and recover almost
identical reference introns (11,456 vs 11,464). Retention and `multi_source_bonus` count
named sources, so their agreement counts twice. Consider counting roles in the retention gate.

## P2 — Input-side issues GMB does not repair

- **Raw long-read alignments** (2.05 M reads for *Z. tritici*) are not valid `--minimap2`
  input; preflight warns (`longread_collapsed`). Collapse with `gmb-longread-consensus`.
- **Protein IDs reused on one sequence.** The OrthoDB genBlastG GTF repeats 33,056 records
  exactly and reuses 9,480 IDs for distinct alignments; GMB groups by ID, so distant reuse
  fuses into one span (up to 5.69 Mb) that `max_span_bp` then drops. Small effect on support;
  de-duplicate upstream or namespace same-sequence reuse in `load_evidence`.
- **Score thresholds** `protein_filter.min_alignment_coverage` / `min_percent_identity` /
  `min_bitscore` need per-exon attributes the genBlastG GTFs lack (0 of 350,980 removed).

## P3 — Interface and developer notes

- Default preset: `gmb-build`, `gmb-preflight` and `load_config()` fall back to `fungi`
  (announced in the log); `run_gene_model_builder()` defaults to `standard`. Name the preset.
- `summary.json` carries the counters under both `summary` and `filtering`.
- 21 config keys are accepted but have no effect; setting one warns (`configuration.md`).
- CI lint was already failing before 2.0.0 work (ruff 44 / black 20 files / isort 6; now
  43 / 20 / 5, no new failures).
- Three slow tests need a whole-genome `z_tritici/` data directory not shipped with the repo.
- `builder.py` `main()` is ~1,000 lines; refactor after the production test.

---

## Fixed in 2.0.0

| issue | fix | evidence |
|---|---|---|
| `transcriptomic_filter.max_transcript_length` (fungi 20 kb) and `allow_single_exon` accepted but ignored | implemented in `filter_chimeras`; apicomplexa set to `null` pending re-validation | 4,962 of 31,529 *Z. tritici* short-read models newly removed; 98% of the > 20 kb survivors spanned ≥ 2 reference genes |
| multi-exon candidates split across loci (5.7% genome-wide) | `scoring.locus_clustering`; `transcript_linked` in fungi | 0 split; with the span filter, exon_overlap lost 10 CDS-exact in core windows, transcript_linked gained 3 |
| genes whose isoforms do not overlap (read-through alternate trimmed by validation) | `drop_detached_isoforms` before gene bounds | 244 of 481 merged genes in windows; 0 after |
| UTR end support counted ends on other sequences; per-model genome-wide scan | same sequence only; prebuilt index | 5.8% / 6.7% of 5′/3′ ends; 1.32 s → 0.01 ms per model |
| second, discarded `select_isoforms` pass | removed | — |
| `gmb-preflight` and `gmb-build` default presets differed | one `DEFAULT_PRESET` | — |
| `gmb-build --list-presets` required `--output-dir` | fixed | — |
| preflight passed empty/unparseable files, impossible coordinates, exon-less backbones, raw long-read tracks | `file_readable`, `coordinates_valid`, `exon_rows_present`, `backbone_cds_present`, `cds_without_exons`, `longread_collapsed` | real bundle: 36 pass, 0 false alarms |
| 23 config keys silently inert | 2 implemented; 21 classified deprecated/unsupported, warn when set; removed from shipped configs | — |
| GFF3 parent-before-child order not guaranteed on ties | stable sort | — |
| `matplotlib`, `biopython` mandatory; comparison code duplicated from ensembl-genes | removed; `gmb-compare` prints where it moved | — |
| example scripts excluded by `.gitignore` (`*_build*`) | re-included | — |
| 10 CDS/UTR tests skipped (needed an absent *Candida* genome) | run on the bundled fixture | — |
