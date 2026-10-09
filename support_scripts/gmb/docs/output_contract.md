# GMB output contract

What GMB produces, which files are the handover, and which are intermediate.

---

## The rule

```
build/      INTERMEDIATE / DEBUG   — do not hand over
finalise/   PRODUCTION HANDOVER    — this is what you ship
```

`gmb-build` writes a complete annotation, but `gmb-finalise` is where cDNA/CDS/protein are
**regenerated from the final GFF3**, where canonical transcripts are chosen, and where the
handover manifest is written.

> **Why this distinction is load-bearing.** Sequences were once captured mid-pipeline while
> the GFF3 continued to be trimmed afterwards, so **5.4% of cDNAs silently disagreed with the
> annotation** — and FASTA QC passed anyway, because it only checked ID coverage. Regenerating
> from the final GFF3 makes that class of drift impossible. Handing over `build/` outputs
> reintroduces the risk.

---

## Handover directory (`finalise/`)

| file | what it is | ship |
|---|---|---|
| `canonical/consensus.canonical_annotated.gff3` | **primary annotation**; canonical transcript tagged `Ensembl_canonical` | **yes** |
| `consensus.gff3` | full annotation including alternative isoforms | **yes** |
| `cdna.fa` | cDNA, regenerated from the final GFF3 | **yes** |
| `cds.fa` | CDS, regenerated from the final GFF3 | **yes** |
| `prot.fa` | protein, regenerated from the final GFF3 | **yes** |
| `canonical/canonical_transcripts.tsv` | which transcript was chosen per gene | **yes** |
| `canonical/transcript_ranking.tsv` | full per-transcript ranking behind that choice | **yes** |
| `canonical/canonical_selection_summary.json` | canonical selection summary | yes |
| `fasta_qc_report.json` | sequence-vs-annotation QC evidence | **yes** |
| `utr_qc_report.json` | UTR invariant QC evidence | **yes** |
| `handover_manifest.json` / `.tsv` | output inventory with sizes and SHA-256 | **yes** |

## Supporting files from `build/`

Ship these alongside the handover — they explain *why* each model was chosen.

| file | what it is | ship |
|---|---|---|
| `evidence_attribution.tsv` | per-transcript evidence sources, protein support, selection reason, rescue flag | **yes** |
| `protein_validation.tsv` | DIAMOND/Psauron per-model scores | yes |
| `collapsed_duplicate_transcripts.tsv` | which duplicate structures were collapsed | yes |
| `resolved_config.yaml` + `resolved_config_sha256` | the exact configuration used | **yes** |
| `run_manifest.json` / `.tsv` | full provenance (see `reproducibility.md`) | **yes** |
| `summary.json` / `summary.tsv` | build-stage counts | yes |
| `gmb.log` | full build log | optional, large |
| `consensus.gff3`, `cdna.fa`, `cds.fa`, `prot.fa` | **pre-finalisation copies** | **no — use `finalise/`** |

---

## `evidence_attribution.tsv`

One row per emitted transcript. The audit trail for every selection decision.

| column | meaning |
|---|---|
| `gene_id`, `transcript_id` | emitted IDs |
| `evidence_sources` | structural sources that produced this structure |
| `protein_alignment_sources` | protein tracks supporting it (attribution only) |
| `protein_support_strength` | `strong` / `weak` / `none` |
| `structural_support_sources`, `n_structural_support_sources` | independent structural agreement |
| `backbone_shortread_agreement` | backbone and an assembly produced the identical intron chain |
| `longread_structural_role` | `primary` / `support_only` / `not_longread` |
| `backbone_intron_rescue` | whether the rescue rule applied |
| `selection_reason` | the rule that selected it |
| `exon_count`, `cds_bp`, `utr_5p_bp`, `utr_3p_bp`, spans | structure |
| `gmb_score` | numeric score (tie-break within a ranking tier) |
| `utr_*_supported` / `_action` / `_reason` | UTR retention decisions |

`selection_reason` values, highest ranking tier first:

| value | meaning |
|---|---|
| `backbone_shortread_agreement` | tier 3 — backbone and an assembly agree exactly |
| `multi_source_structural_agreement` | tier 2 — several independent sources agree |
| `backbone_intron_rescue` | tier 1 — a spliced assembly replaced a collapsed backbone CDS |
| `protein_cds_span_support` | tier 0, described by its protein support |
| `longread_demoted_support_only` | a long-read model demoted by the guard or disposition |
| `longread_only_locus` / `single_exon_longread` | only long-read evidence at this locus |
| `best_single_source` | tier 0, no structural corroboration |

> A transcript may carry `backbone_intron_rescue=True` and still report a *different*
> `selection_reason`: tiers are ordered, so a higher tier names the winner. Both fields are
> emitted because they answer different questions.

---

## File formats and naming

These are stable interfaces: a change to any of them is a breaking change and must be listed
in the release notes.

### `consensus.gff3`

- GFF3 (`##gff-version 3`), **1-based, closed** coordinates. Column 2 is always `GMB`.
- Features: `gene` → `mRNA` → `exon`, `CDS`, `five_prime_UTR`, `three_prime_UTR`.
- Rows are ordered by sequence name, then start; a parent row always precedes its children.
- IDs: gene `{gene_prefix}_{n:05d}` (e.g. `ZTGMB_00042`); transcript `{gene_id}.t{i}`;
  children `{transcript_id}.exon{j}`, `.cds{j}`, `.5utr{j}`, `.3utr{j}`. IDs are unique within a
  run but **not stable across runs** — do not use them as stable identifiers.
- `CDS` column 8 is the GFF3 phase (0/1/2), computed in transcript orientation.
- mRNA attributes: `Evidence` (comma-separated structural sources that produced the
  structure), `ProteinEvidence` (supporting protein-alignment tracks; attribution only) and,
  where a duplicate was collapsed into it, `CollapsedFrom`.
- Transcripts whose strand could not be resolved are never written; they are listed in
  `build/unstranded_exclusions.tsv`.

### `cdna.fa`, `cds.fa`, `prot.fa`

- One record per mRNA (`cdna.fa`) or per mRNA with a CDS (`cds.fa`, `prot.fa`); header
  `>{transcript_id}`.
- Spliced sequence in 5′→3′ transcript orientation (reverse-complemented on `-`).
- `cds.fa` **includes** the terminal stop codon when the model has one; `prot.fa` omits the
  terminal `*` (Ensembl convention).

### `build/` sidecar tables

| file | written when | columns |
|---|---|---|
| `evidence_attribution.tsv` | always | see above |
| `protein_validation.tsv` | `protein_validation.enabled` | `gene_id`, `transcript_id`, `candidate_transcript_id`, `protein_sha256`, `diamond_hit`, `diamond_pident`, `diamond_qcov`, `diamond_scov` (0–100), `diamond_bitscore`, `diamond_evalue`, `psauron_score` (0–1), `protein_length`, `orf_label`, `is_partial_5`, `is_partial_3`, `internal_stop_count`, `protein_coding_score` (0–1), `gmb_score`, `evidence_sources`, `validation_status`, `validation_reason`, `protein_validation_source`, `protein_validation_reused_from`, `diamond_version`, `psauron_version`, `diamond_db` |
| `collapsed_duplicate_transcripts.tsv` | any exact duplicate collapsed | `gene_id`, `retained_transcript_id`, `removed_transcript_id`, `structure_signature`, `duplicate_classification`, `retained_id_reason`, `evidence_merged`, `gmb_score_before`, `gmb_score_recalculated` |
| `unstranded_exclusions.tsv` | any locus with unresolved strand | `gene_id`, `candidate_transcript_id`, `chrom`, `strand`, `evidence_sources`, `rejection_reason` |
| `subset_regions.tsv` | `--seqname` / `--region` / `--regions-file` / `--sample-loci` | `# seed` line, then `seqname`, `start`, `end` |
| `summary.json` / `summary.tsv` | always | run counters (see below) |

`summary.json` has three keys: `summary` and `filtering` (both the full counter dictionary —
kept duplicated for backward compatibility) and `utr` (`transcripts_dropped`). `summary.tsv` is
the same counters as `Metric<TAB>Value` rows. Counters are additive: new ones may appear in a
minor release, existing ones are not renamed. Notable counters: `total_loci` (genes written),
`protein_supported_candidates`, `backbone_intron_rescue` (the gate decision),
`candidates_split_across_loci` (0 under `transcript_linked`), `chimeras_large_intron`,
`chimeras_long_span`, `single_exon_removed` (transcriptomic filter),
`genes_with_detached_isoforms` / `detached_isoforms_removed` (isoforms that no longer overlap
their gene after validation), `validation_*`, `dedup_*`, `duplicate_collapse_*`.

---

## Guarantees

Every handover annotation satisfies, by construction:

- gene coordinates equal the union of that gene's surviving transcripts;
- every CDS/UTR segment lies inside an exon of its own transcript;
- cDNA/CDS/protein sequences are derived from the final GFF3 and match it;
- no transcript contains an unintended internal stop codon;
- exactly one canonical transcript per gene;
- no transcript's parent gene is on a different sequence or strand;
- no two transcripts of a gene share an identical structure;
- every transcript of a gene overlaps the gene's primary transcript, directly or through a
  chain of overlapping transcripts.

These are verified — see `qc.md`. A build that violates one is not a handover candidate.
