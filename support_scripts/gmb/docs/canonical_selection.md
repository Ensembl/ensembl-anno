# Canonical transcript selection

`gmb-canonical-selection` is a standalone, non-destructive post-build step
that ranks isoforms and selects one representative canonical transcript per
gene. It reads a completed build's output files and never modifies them.

---

## Usage

```bash
gmb-canonical-selection \
    --consensus-gff3        output/consensus.gff3 \
    --evidence-attribution  output/evidence_attribution.tsv \
    --protein-validation    output/protein_validation.tsv \    # optional
    --output-dir            output/canonical_selection/
```

`--protein-validation` is optional. Without it, the protein-validation
evidence class is absent from scoring (treated as having no protein support).

---

## Outputs

| File | Description |
| :--- | :---------- |
| `canonical_transcripts.tsv` | One row per gene: the selected canonical transcript ID and key scores |
| `transcript_ranking.tsv` | All isoforms per gene with their rank key components |
| `canonical_selection_summary.json` | Run statistics (genes processed, single-isoform fraction, etc.) |

---

## Selection policy

Canonical selection uses a deterministic lexicographic priority order:

1. **Complete ORF** — has both start and stop codon (preferred). ORF status comes from
   `protein_validation.tsv`, so without protein validation every isoform ties here.
2. **No internal stop codon**
3. **Protein validation support** — DIAMOND balanced coverage, identity and Psauron
4. **Biological evidence-class breadth** — number of distinct evidence
   *classes* represented (not raw named-source count):
   - `backbone` (Helixer / Tiberius)
   - `short_read_transcriptomic` (Scallop, StringTie — counted as one class
     even though both are present)
   - `long_read_transcriptomic` (Minimap2Consensus)
   - `protein_alignment` (only when named in `evidence_sources`; the build reports
     protein alignments separately, so in practice this class is not counted)
   - `protein_validation` (DIAMOND / Psauron result)
5. **Build primary** — only with `canonical_selection.prefer_build_primary: true`: the
   build's primary transcript (`<gene>.t1`) wins here
6. **Named-source count** — Scallop and StringTie count as two
7. **GMB continuous score** — the weighted sum from `gmb-build`
8. **CDS length** (longer preferred)
9. **Transcript ID** (lexicographic tiebreaker for reproducibility)

Without `prefer_build_primary`, tier 6 lets a Scallop+StringTie isoform replace a backbone
primary whenever the earlier tiers tie, which they do when protein validation is off. On
*Z. tritici* that changed 199 of ~3,100 canonical transcripts in nine windows and lost more
reference-exact models than it gained (`release/evidence_integration.md`).

**What can still override the build primary when `prefer_build_primary` is on.** Only the
tiers above it: an isoform with a complete ORF where the primary's is partial, no internal stop
where the primary has one (both only when protein validation ran, see tier 1), a better
protein-validation result, or more evidence *classes*. In practice the last means a long-read +
short-read isoform against a backbone-only or short-read-only primary. On the 17
*Z. tritici* evaluation windows the experimental configuration overrode no primary without long
reads, and 243 of 6,095 with collapsed long reads, all by evidence-class breadth (21
reference-exact canonicals gained, 14 lost). With protein validation on, 14 of 5,988 were
overridden, all by the protein-validation tier.

**Reading the result.** `selection_reason` reports `LOW_CONFIDENCE_NO_PROTEIN_SUPPORT` for any
canonical without a protein-validation result, which hides the tier that decided. Use
`deciding_tier` and `canonical_is_build_primary` in `canonical_transcripts.tsv`, and
`canonical_not_build_primary_by_tier` in `canonical_selection_summary.json`.

**Configuration.** `gmb-finalise` uses the build's `resolved_config.yaml`, so a build made with
the experimental overlay is finalised with `prefer_build_primary: true` automatically. The
standalone `gmb-canonical-selection` command takes its own `--preset` / `--config`; pass the
same overlay there.

The `canonical_total_score` column in `canonical_transcripts.tsv` is the
continuous weighted score. It reflects *how much better* the canonical is
relative to alternatives — it is **not** what picked it (the priority order
above is).

---

## Evidence-class breadth terminology

The scoring reports "biological evidence-class breadth" rather than
"independent evidence-source count" because Scallop and StringTie are run
over the same RNA-seq libraries and do not constitute two independent pieces
of evidence. Grouping named sources into evidence *classes* avoids inflating
the apparent support for transcripts backed only by short-read assembly.

Unknown or unrecognised source names are mapped to an `other` class and
reported in `canonical_selection_summary.json` as a warning (once per run,
not once per transcript).

---

## InterPro resolver (second stage)

For genes where canonical selection could not distinguish between isoforms
confidently, `gmb-interpro-review` + `gmb-interpro-resolve` provides a
second-stage resolver backed by InterProScan domain evidence. This is
disabled by default and is not required for a normal build.

See **[docs/interpro_resolver.md](interpro_resolver.md)**.
