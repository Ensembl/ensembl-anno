# GMB known issues

Open problems, ordered by severity, then what was fixed for 2.0.0. Each entry gives the
evidence, the current mitigation and the recommended action. Last reviewed 2026-10-09
(whole-genome Helixer assessment of *Z. tritici*, `release/z_tritici_validation.md`; evidence
integration, `release/evidence_integration.md`).

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

## P1 — Evidence cannot change a backbone model under the fungi preset

Where Helixer predicts a gene, the selected primary is Helixer's model: a Helixer model with
positional protein overlap scores 5.1, the best structure supported only by short-read
transcripts 4.5. On 17 *Z. tritici* windows the primary transcripts of the fungi build equal
Helixer's on every reference metric; protein evidence alone, one short-read track alone, or
long reads alone change nothing at all. Canonical selection then swaps ~6% of canonical
transcripts to Scallop+StringTie isoforms because the two assemblers count as two sources,
which is net harmful (held-out core CDS-exact: Helixer 710, fungi canonical 695).
**Mitigation:** `configs/fungi_evidence_integration.experimental.yaml`
(`primary_selection: junction_supported`, `isoform_cds_overlap: drop`, `min_cds_bp: 150`,
`canonical_selection.prefer_build_primary`): held-out core canonical CDS-exact 733 (754 with
collapsed long reads), splits/merges at Helixer's level. **Action:** run both configurations
in the production-scale test; adopt the overlay if the genome-wide result matches.

## P1 — Alternate isoforms drive most gene splits and merges

With one isoform per gene, GMB's core-window output equals Helixer's on every reference
metric. In the fungi preset 41% of alternate isoforms share no coding base with their gene's
primary: UTR and read-through overlaps admit a neighbouring ORF as an "isoform", chaining
genes together (all-isoform splits/merges 342/253 on held-out core windows vs Helixer 64/23).
**Mitigation:** `scoring.isoform_cds_overlap: drop` (in the experimental overlay) refuses
them while keeping coding-sharing alternates: splits/merges 65/24, CDS-exact unchanged or
better. **Action (decision):** adopt `drop` for fungi; then decide how many coding-sharing
alternates to ship (`max_isoforms_per_locus`; half of the remaining alternates differ from the
primary only in UTR).

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

- The protein-validation penalty is now `protein_validation.penalty` (default 5.0, unchanged
  behaviour). With the fungal weights (`diamond_weight` 0.7, `min_score` 0.7) every model
  without a DIAMOND hit is penalised whatever its Psauron score (186 of 399 transcripts in one
  window against UniProt eukaryota), and 5.0 exceeds every backbone weight: enabling it
  removes most alternates and changes no primary (dev core: transcripts 3,960 → 2,650,
  CDS-exact 719 → 719, missed 35 → 41). Revisit `min_score`/weights before enabling it.
- `validation.max_exon_len_reference: reference` is accepted but falls back to
  `candidates_supported` with a warning.

## P2 — Transcript-only antisense and lncRNA ORFs become protein-coding genes

About 5% of fungi-preset genes on *Z. tritici* windows are transcript-only models antisense
to a reference gene (held-out core: 123 of 2,496; 188, 7.5%, counting Helixer's own antisense
calls) and a further ~2% are transcript-only intergenic novel models. Their introns are
canonical on the declared strand (genuine antisense transcription), but the ORF is
incidental: median Psauron 0.10, no DIAMOND hit, and almost none agree with a spliced protein
alignment. Positional protein support does not separate them (60% "strong"). Neither the
opposite-strand-overlap test nor protein coverage can tell them from the ~15 real genes per
2,100 reference genes (dev windows) where Helixer called the wrong strand.
**Action:** decide on a coding-potential gate for transcript-only new loci (Psauron < 0.5
would remove ~80% of antisense models and ~65% of the transcript-only models that overlap a
reference gene) or report them in a separate low-confidence class.

## P2 — Long-read evidence must be collapsed upstream, and is the strongest evidence

The *Z. tritici* long-read track is per-read (unusable as `--minimap2`), but collapsed with
`gmb-longread-consensus` (fungal span 20 kb / intron 3 kb) it is the most precise transcript
evidence available (58.6% of its introns are reference introns vs 43.9–56.5% for the
short-read tracks) and the only independent one. With it the experimental overlay reaches
held-out core canonical CDS-exact 754 (Helixer 710). **Action:** collapse the long-read
track in the production pipeline before GMB.

## P2 — Short-read tracks are not independent evidence

Canonical selection's named-source tier counted them twice and swapped canonical transcripts
to them (`canonical_selection.prefer_build_primary` mitigates).
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

## P2 — Reference genes Helixer misses are mostly not recoverable from this evidence

Of 44 core reference genes Helixer misses in six windows, 18 have a short-read transcript
containing the whole CDS, but always inside a long read-through assembly (median span 13–160
kb) whose longest ORF is elsewhere; 13 of those transcripts contain an ATG ORF equal to the
reference CDS. GMB infers one ORF per transcript. Five are reproduced exactly by a protein
alignment, which GMB cannot use as a candidate; gap-filling with protein-defined CDSs would
add 73 models per nine windows of which 4 are reference-exact. **Action:** longer term —
protein-to-genome model building with a repeat annotation; multi-ORF inference only with
coding-potential support.

## P3 — Interface and developer notes

- Canonical selection takes ORF completeness from `protein_validation.tsv`, so without
  protein validation the complete-ORF tier ties every isoform; its `protein_alignment`
  evidence class is never counted (protein tracks are reported in
  `protein_alignment_sources`, not `evidence_sources`), so every gene is labelled
  `LOW_CONFIDENCE_NO_PROTEIN_SUPPORT` when protein validation is off.
- Scores take 22 distinct values on *Z. tritici* (64% of candidates at exactly 3.0 or 5.1);
  24% of multi-isoform genes have an alternate tied with the primary, broken by input order.

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
| `scoring.min_cds_bp` (150) never applied (read a field nothing set) | applied from the candidate CDS; shipped value 0 keeps validated output; 150 in the experimental fungal overlay | dev windows: 150 removes 46 core / 36 accessory genes, nearly all novel; no reference CDS < 150 bp |
| hard-coded −5.0 protein-validation penalty | `protein_validation.penalty` | — |
| second, discarded `select_isoforms` pass | removed | — |
| `gmb-preflight` and `gmb-build` default presets differed | one `DEFAULT_PRESET` | — |
| `gmb-build --list-presets` required `--output-dir` | fixed | — |
| preflight passed empty/unparseable files, impossible coordinates, exon-less backbones, raw long-read tracks | `file_readable`, `coordinates_valid`, `exon_rows_present`, `backbone_cds_present`, `cds_without_exons`, `longread_collapsed` | real bundle: 36 pass, 0 false alarms |
| 23 config keys silently inert | 2 implemented; 21 classified deprecated/unsupported, warn when set; removed from shipped configs | — |
| GFF3 parent-before-child order not guaranteed on ties | stable sort | — |
| `matplotlib`, `biopython` mandatory; comparison code duplicated from ensembl-genes | removed; `gmb-compare` prints where it moved | — |
| example scripts excluded by `.gitignore` (`*_build*`) | re-included | — |
| 10 CDS/UTR tests skipped (needed an absent *Candida* genome) | run on the bundled fixture | — |
