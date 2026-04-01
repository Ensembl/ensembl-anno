# Known Issues

This file documents known bugs and limitations discovered during the production-readiness
audit.  **None of these are fixed in this PR** — they are catalogued here so they can be
addressed in a subsequent, dedicated PR with targeted test coverage.

---

## KI-001 · Evidence attribution TSV may not be written on all code paths

**Severity:** Low
**Affected file:** `gene_model_builder.py` / `reporting.py`
**Observed:** The `evidence_attribution.tsv` file is listed as an expected output in
documentation and test fixtures, but the file is not guaranteed to be written in all
configurations. When protein validation is disabled, the attribution TSV is sometimes
absent.
**Impact:** Tests that assert the file exists may skip or fail depending on config.
**Workaround:** Run with `--preset fungi_default` and protein validation disabled; check
`summary.json` for attribution data instead.

---

## KI-002 · Seqname mapping is applied globally before any subsetting

**Severity:** Low / informational
**Affected file:** `subset_utils.py`
**Observed:** The `--assembly-report` / `--seqname-map` remapping is applied to **all**
input DataFrames before `--seqname` / `--region` filtering. This is correct behaviour,
but it means the user must always specify `--seqname` using the *mapped* (post-remap)
chromosome name, not the raw accession from the source GTF/GFF3.
**Documentation gap:** The README fast-test section documents this correctly, but the
`--help` text for `--seqname` does not mention that the name should be post-mapping.
**Workaround:** Pass the mapped name (e.g. `--seqname 1`) not the raw accession
(e.g. `CM001642.1`).

---

## KI-003 · `compare_annotations.py` may produce a divide-by-zero when reference has 0 genes on a requested seqname

**Severity:** Medium
**Affected file:** `compare_annotations.py`
**Observed:** When `--seqname` is used and the reference GFF3 contains no genes on that
chromosome, sensitivity metrics (computed as `TP / (TP + FN)`) produce a division by
zero because `TP + FN == 0`. The script currently catches this silently with `nan`, but
the downstream TSV/JSON contains `NaN` entries which confuse downstream consumers.
**Impact:** `NaN` rows in comparison reports when subsetting to a chromosome that the
reference annotation does not cover.
**Workaround:** Only subset to chromosomes that appear in both the consensus and reference
annotations.

---

## KI-004 · Protein validation score not available without DIAMOND + Psauron

**Severity:** Low / by design
**Affected file:** `protein_validation.py`, `gene_model_builder.py`
**Observed:** The `protein_coding_score` field (used by the scoring module when
`protein_validation.enabled: true`) is only populated if DIAMOND and Psauron are
installed and the `--check-deps` gate passes. When these tools are absent, the field
defaults to `0.0`, which can silently demote otherwise high-quality models.
**Impact:** Pipeline outputs differ between environments with/without DIAMOND+Psauron,
even for the same input data.
**Workaround:** Set `protein_validation.enabled: false` (the default) when running
without these tools, or ensure both are on `$PATH`.

---

## KI-005 · `subset_regions.tsv` written before pipeline stages complete

**Severity:** Cosmetic
**Affected file:** `subset_utils.py` / `gene_model_builder.py`
**Observed:** `subset_regions.tsv` is written at the start of the run (immediately after
argument parsing), so it exists even if the pipeline subsequently fails. This means
checking for its existence does not confirm a successful pipeline run.
**Impact:** Tests that use `subset_regions.tsv` as a success indicator may pass even on a
crashed run.
**Workaround:** Use `summary.json` or `consensus.gff3` as success indicators.

---

## KI-006 · Single-exon genes with no protein support may be filtered inconsistently

**Severity:** Low
**Affected file:** `evidence_filter.py` (`filter_helixer_models`)
**Observed:** The Helixer single-exon filter (`require_protein_support_for_single_exon`)
applies to Helixer-sourced models but not to single-exon transcriptomic models
(Scallop/StringTie). This asymmetry is intentional per the current config design, but
is not documented, leading to confusion about why single-exon Scallop models appear in
the output when single-exon Helixer models are filtered.
**Documentation gap:** The filter logic comment says "Helixer single-exon filter" but
does not clarify the deliberate asymmetry.
**Workaround:** Increase `scoring.min_combined_evidence` to indirectly penalise
single-exon transcriptomic models without protein support.
