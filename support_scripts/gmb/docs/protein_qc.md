# Protein validation and post-build protein QC

GMB can run DIAMOND and PSAURON in two separate places:

| | in-build protein validation | post-build protein QC |
|---|---|---|
| command | `gmb-build` with `protein_validation.enabled: true` | `gmb-protein-qc` |
| changes the annotation | **yes**: penalises or drops models during selection | **no** |
| default | off in every shipped preset | run on request |
| proteins scored | backbone/transcript candidates (not generated secondary ORFs) | every final protein, secondary ORFs included |
| output | `build/protein_validation.tsv`, read automatically by `gmb-finalise` | `protein_qc.tsv` and `protein_qc_summary.json`, read by nothing in GMB |
| recommended for fungi | **no** (see below) | yes, as a QC report |

Both call the same code (`gmb.pipeline.protein_validation.batch_score_proteins`): one DIAMOND
`blastp` run and one PSAURON run on the unique sequences, with SHA-256 reuse and recorded tool
versions and database.

## Why in-build protein validation is not recommended for fungi

On the full *U. maydis* 521 genome (GCA_000328475.2, Ensembl Fungi reference used only for
evaluation), adding protein validation with the fungi preset settings to P4d with secondary ORFs
**lowered** exact reference CDS agreement of build primaries from 5,148 to 5,080 (+33 / −101
genes) and removed 37% of transcripts (12,631 → 7,892). An *A. fumigatus* Af293 pilot region
showed the same direction (733 → 729).

The mechanism:

- With the fungi weights (`diamond_weight` 0.7 on `min(1, bitscore/1000)`, `psauron_weight`
  0.3, `min_score` 0.7) almost no protein reaches `min_score`: 86% of final *U. maydis*
  transcripts are penalised, including 91% of reference-exact primaries.
- The fixed 5-point penalty pushes almost every alternate below `scoring.min_alternate_score`,
  so the junction- and evidence-supported promotions lose the candidates they would promote.
- Primaries are penalised, not removed, so low-coding-potential novel primaries are kept anyway
  (650 of 653 with PSAURON < 0.5).

PSAURON on its own separates likely non-coding models well (novel primaries without reference
overlap: median 0.21, 70% below 0.5; reference-exact primaries: 45 of 5,148 below 0.5). That
signal is what `gmb-protein-qc` reports. Using it to filter genes would be a new policy and
needs its own evaluation.

Details: `gmb_evidence_workspace/experiments/20261117_fungal_multigenome_full_validation/`
(`FUNGAL_BENCHMARK_REPORT.md`, `PROTEIN_VALIDATION_AUDIT.md`).

## Running post-build protein QC

```bash
gmb-protein-qc \
  --proteins    out/finalise/prot.fa \
  --gff3        out/finalise/canonical/consensus.canonical_annotated.gff3 \
  --build-dir   out/build \
  --diamond-db  /data/uniprot_fungi_excl_target.dmnd \
  --diamond-path /path/to/diamond --psauron-path /path/to/psauron --psauron-cpu \
  --output-dir  out/protein_qc
```

- An editable install made before this command existed has no `gmb-protein-qc` script: re-run
  `pip install -e .`, or call `python -m gmb.cli.protein_qc` with the same arguments.
- `--proteins` and `--gff3` must come from the same annotation. Every FASTA ID must be an mRNA of
  the GFF3, or the tool stops. With the canonical-annotated GFF3, the `canonical` column and the
  `canonical` summary group are filled in.
- `--build-dir` reads tool paths, DIAMOND coverage gates and the score weights from the build's
  `resolved_config.yaml`. Without it, `--preset` / `--config` are used. The flags override either.
  The build's own configuration is never modified: validation is switched on in a private copy.
- The tool refuses to overwrite an existing `protein_qc.tsv` unless `--overwrite` is given.
- Use a DIAMOND database **without the target species**. Otherwise hits are self-hits, not
  independent evidence. The database path, size, modification time and `diamond dbinfo` counts
  are recorded in the summary.
- Runtime: *U. maydis* chromosome 1 (1,618 transcripts, 1,212 unique proteins) took 43 s wall time
  (103 s CPU) against an 896,208-sequence UniProt eukaryote database, with PSAURON on CPU.

### Outputs

`protein_qc.tsv`, one row per mRNA of the GFF3:

| column | meaning |
|---|---|
| `gene_id`, `transcript_id`, `evidence` | from the GFF3 |
| `candidate_origin` | `existing` or `secondary_orf` (GFF3 `CandidateOrigin`) |
| `canonical` | `True`/`False` from `tag=Ensembl_canonical`; empty if the GFF3 has no canonical tags |
| `protein_sha256`, `protein_length` | as in `protein_validation.tsv` |
| `diamond_hit`, `diamond_pident`, `diamond_qcov`, `diamond_scov` (0–100), `diamond_bitscore`, `diamond_evalue` | best DIAMOND hit passing the configured coverage gates |
| `psauron_score` (0–1) | PSAURON in-frame score |
| `psauron_is_protein` | PSAURON's own call at its built-in threshold; GMB adds no threshold |
| `protein_coding_score` | the weighted DIAMOND/PSAURON combination under the configured weights, for comparison with in-build validation; not a pass/fail call |
| `validation_status`, `validation_reason`, `protein_validation_source`, `protein_validation_reused_from`, `diamond_version`, `psauron_version`, `diamond_db` | provenance, as in `protein_validation.tsv` |

`protein_qc_summary.json` holds input paths and SHA-256s, the GMB version and Git state, tool
versions, database details, status counts, and per group (`all`, `origin:existing`,
`origin:secondary_orf`, `canonical`) the DIAMOND hit fraction, PSAURON quartiles and
`psauron_is_protein` counts.

`protein_qc.tsv` carries the columns `gmb-interpro-review` reads from `protein_validation.tsv`
(`transcript_id`, `protein_sha256`, DIAMOND and PSAURON fields), so it can be passed explicitly as
`--protein-validation protein_qc.tsv` to review a build made without in-build validation. That
path has not been tested end to end, and a subsequent `gmb-interpro-resolve` can change canonical
choices, so treat it as a deliberate, separate step.

**Not consumed.** No GMB module reads `protein_qc*` files (enforced by
`tests/test_protein_qc.py`). In particular `gmb-finalise` only reads
`build/protein_validation.tsv`, so running the QC never changes canonical selection.

## Known limitations of in-build validation

These apply only with `protein_validation.enabled: true`; none affect the shipped presets.

- **Fixed in this release (D1).** Merging identical exon structures used to copy each later
  member's `protein_coding_score` onto the representative, so a structure could be judged on
  another model's protein (*U. maydis* chromosome 1: 141 of 2,414 structures, 18 pass/fail
  flips). Each model now keeps its own score.
- **D2: generated secondary ORFs are not validated in-build.** They carry
  `protein_coding_score: None` and bypass the penalty. `gmb-protein-qc` scores them.
- **D3: `policy: drop` judges a CDS structure by every carrier.** In evidence-score promotion
  (`_promote_evidence_supported`), a structure is rejected if any transcript carrying it is below
  `min_score`, and transcripts without a score count as 0.0. Not reachable with the shipped
  `penalize` policy.
- **D4: score scale.** With the fungi weights 78–87% of fungal transcripts fall below `min_score`
  (bitscore/1000 rarely approaches 1). The weights are historical opt-in settings, unchanged.
- **Canonical selection depends on `protein_validation.tsv` (F1).** The complete-ORF / internal-stop
  tier reads ORF fields only from that file, so without in-build validation every isoform ties at
  that tier (`evidence_attribution.tsv` has the same fields but is not used there). The protein
  tier ranks above evidence breadth and `prefer_build_primary`, and without the file every gene is
  labelled `LOW_CONFIDENCE_NO_PROTEIN_SUPPORT` (the deciding tier is in `deciding_tier`).
  Changing this would change default canonical outputs, so it is recorded, not changed.
- `utr.end_support_mode: protein_validated` (or `either`) needs validation enabled; with it off,
  every UTR end is treated as unsupported without a warning.
