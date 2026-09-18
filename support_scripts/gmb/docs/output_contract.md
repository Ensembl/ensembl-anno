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

## Guarantees

Every handover annotation satisfies, by construction:

- gene coordinates equal the union of that gene's surviving transcripts;
- every CDS/UTR segment lies inside an exon of its own transcript;
- cDNA/CDS/protein sequences are derived from the final GFF3 and match it;
- no transcript contains an unintended internal stop codon;
- exactly one canonical transcript per gene;
- no transcript's parent gene is on a different sequence or strand;
- no two transcripts of a gene share an identical structure.

These are verified — see `qc.md`. A build that violates one is not a handover candidate.
