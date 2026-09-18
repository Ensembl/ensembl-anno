# Creating a configuration for a new clade

Short and practical. The goal is a frozen, defensible configuration — not the best possible
score on one genome.

**The rule that governs everything below:**

> A reference annotation may be used to **evaluate a completed build**. It must never be used
> to **choose the production configuration** for that genome. The moment you adjust a setting
> because of what the reference said, you have stopped validating and started fitting, and the
> resulting number no longer predicts anything about the next genome.

---

## Step 1 — Run the upstream evidence modules

Produce the bundle in `input_contract.md`: genome FASTA, ab initio backbone, short-read
assemblies, optionally long-read transcript models, protein alignments. GMB does not generate
these.

## Step 2 — Run GMB preflight

```bash
gmb-preflight --preset standard \
  --genome "$GENOME" --helixer "$BACKBONE" \
  --scallop "$SCALLOP" --stringtie "$STRINGTIE" --orthodb "$PROTEINS" \
  --output-dir "$OUT/preflight"
```

Reference-free. Fix anything it FAILs before going further.

## Step 3 — Read the numbers that decide the configuration

From `preflight_report.txt`:

| what to look at | what it tells you |
|---|---|
| **model counts per track** | whether a track has enough evidence to matter |
| **single vs multi-exon fraction per track** | the single most important comparison — see below |
| **canonical splice fraction per track** | whether a track's splice structures are trustworthy at all |
| **resolved role and weight per track** | that nothing fell through to the unknown-source weight |
| **source overlap** | how much independent corroboration is available |
| **backbone vs assembled-transcript complexity** | which source is the better structural authority |

### The comparison that matters most

Compare the **backbone's multi-exon fraction** with the **assembled transcripts'**.

| observed | backbone | assembled | ratio | what it means |
|---|---|---|---|---|
| P. falciparum (Tiberius) | 17.3% | 98.7% | 5.7× | backbone collapses introns; assemblies are the better structure |
| T. gondii (Tiberius) | 38.7% | 95.4% | 2.5× | same state, less extreme |
| Z. tritici (Helixer) | 71.9% | 92.1% | 1.3× | backbone is sound; assemblies are *not* better |

This single comparison predicted the outcome in all three genomes, **before any reference was
consulted**. It is what the `backbone_intron_rescue: auto` gate automates.

## Step 4 — Start from the neutral preset

```bash
cp configs/new_clade_template.yaml my_clade.yaml
```

Fill in the evidence role labels. That is often the only edit you need.

**Do not start from `apicomplexa.yaml` or `fungi.yaml`.** They encode evidence states that may
not be yours — and those two states are opposites of each other.

## Step 5 — Change a biological policy only for an evidence-derived reason

For each policy, the reason must be visible in preflight output:

| policy | enable when |
|---|---|
| `backbone_intron_rescue: auto` | backbone multi-exon fraction is low and assembled transcripts are both substantially more spliced and canonically spliced. `auto` checks this itself — prefer it to `on` |
| `structural_corroboration: true` | several credible independent structural tracks exist |
| `protein_support_mode: cds_span_compatible` | protein support predicts structural correctness *monotonically* for your evidence |
| `longread_structural_guard: true` | you supplied a long-read track whose splice quality preflight did not fully vouch for |
| `longread_disposition: support_only` / `reject` | preflight FAILed the long-read track's splice-quality check |
| `weights.*` | the multi-exon and splice comparisons say one role is the better structural authority |

Leave everything else alone. An unchanged neutral preset is a perfectly respectable
production configuration.

## Step 6 — Validate against a reference, if one exists

```bash
gmb-build    --preset standard --config my_clade.yaml ... --output-dir "$OUT/build"
gmb-finalise --preset standard --config my_clade.yaml ... --output-dir "$OUT/finalise"
gmb-compare  --query "$OUT/finalise/consensus.gff3" --reference "$REFERENCE" \
             --reference-fasta "$GENOME" --evaluation-mode protein_coding \
             --output-dir "$OUT/comparison"
```

Report **CDS exact** and **CDS exact — multi-exon** as the headline. Do not report
Exact Match as coding accuracy: it includes UTR extent, and GMB emits fewer UTRs than a curated
reference (in one build, 1,206 reference genes had a byte-identical CDS yet failed Exact Match
on UTR extent alone).

Run **both** arms — your candidate config *and* the plain neutral preset. If the candidate is
not better, keep the neutral one. That comparison is the entire point, and it is how the
fungal default came to be the baseline.

## Step 7 — Freeze the config

Stop editing. Record the SHA-256 (the run manifest does this automatically). A config that
keeps changing has not been validated — only its latest version has, on one genome.

## Step 8 — Test the frozen config on an independent genome before making it a clade default

This is the step that distinguishes a clade preset from a one-genome tuning.

Run the **byte-identical** config on a second genome of the same clade and compare against
that genome's reference. Only then propose it as a preset.

Two outcomes from doing exactly this:

- The Apicomplexa policy, frozen on P. falciparum, transferred to T. gondii — so it became a
  preset.
- The same policy, applied unchanged to Z. tritici, made the annotation **worse**
  (improvements : regressions = 249 : 555) — so fungi kept the baseline instead.

Without step 8, the second case would have shipped as a "validated" fungal default.

---

## Checklist before proposing a preset

- [ ] preflight PASSes (or every WARN is understood and documented)
- [ ] every policy change traces to a number in the preflight report
- [ ] the config was frozen before the reference was consulted
- [ ] candidate beats the plain neutral preset on CDS exact and multi-exon CDS exact
- [ ] the byte-identical config was tested on a second, independent genome
- [ ] the preset file records the evidence for each non-obvious value
- [ ] `tests/test_production_contract.py` has a case pinning its resolved values
