# InterProScan resolver for ambiguous canonical transcripts

InterProScan is used here **only** as a second-stage resolver for canonical
transcript choices GMB could not make confidently. It is never part of a normal
GMB run, is never required for GMB to complete, and never triggers a second
DIAMOND/psauron pass.

This is explicitly **not** a geneset QC tool: it does not compute whole-proteome
domain coverage, reference-vs-GMB comparisons, or Pfam coverage statistics.

## Master switch: disabled by default

All of this lives under `canonical_selection.interpro_resolver:` in config
(moved from a top-level `interpro_review:` section):

```yaml
canonical_selection:
  interpro_resolver:
    enabled: false            # master switch -- InterPro involvement is opt-in
    run_interproscan: false   # Mode A (true) vs Mode B (false, default)
    apply_replacements: true  # whether a passed-safeguard replacement actually happens
```

`enabled: false` (the default) means none of this runs, and a first-pass GMB
build completes exactly as if this module did not exist -- confirmed
architecturally, not just by config default: `gmb.pipeline.builder`/`gmb.cli.build`
never import `interpro_review` or `interpro_resolver` at all; both are
separate, explicitly-invoked CLIs (`gmb-interpro-review`, `gmb-interpro-resolve`).

`apply_replacements` and `run_interproscan` are independent of `enabled` and
of each other:

| `enabled` | `run_interproscan` | `apply_replacements` | Behaviour |
| :-: | :-: | :-: | :--- |
| false | – | – | Nothing runs. |
| true | false (Mode B) | true (default) | Prepare + consume an already-completed InterProScan run; safeguarded replacements are applied. |
| true | false | false | Prepare + consume, but report-only -- verdicts/reason codes recorded, canonical never replaced. Useful for evaluating the evidence model on a new dataset first. |
| true | true (Mode A) | true/false | GMB also launches InterProScan 6 itself via Nextflow (see below), then behaves as above. |

## Where it sits

```
candidate models
    ↓
translated proteins
    ↓
DIAMOND + Psauron  ← runs ONCE per build, on deduplicated sequences
    ↓
protein_validation.tsv   (normalised, checksummed, provenance-stamped;
                          records calculated-vs-reused per transcript)
    ├── model scoring / filtering policy   (gmb.pipeline.scoring)
    └── initial canonical selection        (gmb-canonical-selection)
             ↓
       confidence / ambiguity assessment
             ↓
       if interpro_resolver.enabled=false: final = initial, unconditionally
             ↓
       if enabled=true:
           ambiguous genes only            (gmb-interpro-review)
                ↓
           InterProScan  ← Mode A: GMB launches it (Nextflow); or
                            Mode B: already-completed run, any execution profile
                ↓
           conservative replacement policy (gmb-interpro-resolve)
                ↓
           final canonical + full attribution (canonical_decisions.tsv,
           consensus.final_canonical_annotated.gff3)
```

Canonical selection reads `protein_validation.tsv` and never re-invokes DIAMOND
or Psauron -- this is enforced by a test
(`TestCanonicalSelectionDoesNotInvokeExternalTools`).

## Step 1 -- prepare review candidates

```bash
gmb-interpro-review \
    --canonical-transcripts out/canonical_selection/canonical_transcripts.tsv \
    --transcript-ranking    out/canonical_selection/transcript_ranking.tsv \
    --protein-validation    out/protein_validation.tsv \
    --prot-fa               out/prot.fa \
    --config                my_config.yaml \
    --output-dir            out/interpro_review
```

Writes:

| File | Contents |
| :--- | :--- |
| `interpro_review_candidates.faa` | One batched FASTA. Headers are protein SHA-256 checksums, so identical proteins are submitted once. |
| `interpro_review_manifest.tsv` | Gene/transcript/checksum mapping plus the DIAMOND, Psauron, evidence-class and confidence context behind the review. |
| `interpro_review_summary.json` | Counts, trigger breakdown, and the policy settings used. |

Ambiguity triggers are configured under `canonical_selection.interpro_resolver:`
(see `gmb/pipeline/config.py`'s `InterProResolverConfig`). Single-transcript
genes, genes whose top two candidates encode identical proteins, and confident
well-separated choices are never submitted.

## Step 2 -- get InterProScan output (Mode A or Mode B)

### Mode B (default): consume an already-completed run

InterProScan 6 is already a Nextflow workflow; by default GMB does not wrap
it. The resolver only needs the output files, so **any** execution profile
works.

**Laptop development (Docker)** -- matches the shape of the verified local run
(InterProScan 6.0.1, Nextflow 26.04.6, InterPro 109.0):

```bash
nextflow run ebi-pf-team/interproscan6 -r 6.0.1 \
    -profile docker \
    --input   out/interpro_review/interpro_review_candidates.faa \
    --datadir /path/to/interproscan6/data \
    --interpro latest \
    --outdir  out/interpro_run/results \
    -work-dir out/interpro_run/work
```

**Cluster (Slurm + Singularity/Apptainer)** -- the workflow already ships
`slurm`, `singularity` and `apptainer` profiles, so no new workflow is needed,
only an external config supplying site values. Do **not** assume Docker is
available on the cluster:

```bash
nextflow run ebi-pf-team/interproscan6 -r 6.0.1 \
    -profile slurm,singularity \
    -c        docs/interproscan6_cluster.example.config \
    --input   out/interpro_review/interpro_review_candidates.faa \
    --datadir "$SHARED_INTERPRO_DATA" \
    --interpro latest \
    --outdir  "$RESULTS_DIR" \
    -work-dir "$WORK_DIR"
```

An example config is provided at
[`docs/interproscan6_cluster.example.config`](interproscan6_cluster.example.config).
Every site-specific value in it is a clearly marked placeholder.

### Mode A: GMB launches InterProScan 6 itself

Set `run_interproscan: true` and fill in
`canonical_selection.interpro_resolver.nextflow:`:

```yaml
canonical_selection:
  interpro_resolver:
    enabled: true
    run_interproscan: true
    nextflow:
      nextflow_executable: nextflow
      workflow: ebi-pf-team/interproscan6
      revision: "6.0.1"
      profile: docker                 # or "slurm,singularity" / "slurm,apptainer"
      config_file: null                # -c <file>, e.g. the cluster example above
      data_dir: /path/to/interproscan6/data   # REQUIRED for Mode A -- no default
      interpro_release: latest
      applications: null                # e.g. "Pfam,PANTHER"; null = workflow default set
      work_dir: null                    # defaults to <output-dir>/interpro_run/work
      output_dir: null                  # defaults to <output-dir>/interpro_run/results
      cpus: null
      max_workers: null
      extra_args: []
```

`gmb.pipeline.interpro_review.build_interproscan_nextflow_command()` builds
the exact command line from these fields alone (independently testable
without Nextflow installed); `run_interproscan_workflow()` then runs it via
`subprocess` and raises on a non-zero exit rather than silently producing an
empty result. **No Docker image name, Slurm partition, Singularity/Apptainer
cache path, or any other site-specific value is ever hardcoded in GMB Python
code** -- every one of them is a field above, defaulting to `null`/a generic
placeholder rather than a guessed real value.

`gmb-interpro-review` calls this automatically when `run_interproscan: true`
and there is at least one review candidate.

## Step 3 -- resolve

```bash
gmb-interpro-resolve \
    --manifest             out/interpro_review/interpro_review_manifest.tsv \
    --interpro-jsonl        out/interpro_run/results/*.jsonl \
    --config                my_config.yaml \
    --canonical-transcripts out/canonical_selection/canonical_transcripts.tsv \
    --consensus-gff3        out/consensus.gff3 \
    --output-dir            out/interpro_review
```

`--interpro-gff3` is accepted as an alternative to `--interpro-jsonl`.
**JSONL is preferred**: it is the only emitted format carrying the
`representative` flag, the integrated InterPro entry (with its type), and the
member-database release version together. TSV lacks the representative flag
entirely.

`--canonical-transcripts`/`--consensus-gff3` are optional: without them you
get the resolver's own verdict/reason-code report only (no attribution files);
with them you additionally get `canonical_decisions.tsv` and
`consensus.final_canonical_annotated.gff3` (see Attribution below).

GMB takes explicit output paths and is indifferent to how they were produced --
no container runtime, image name, work directory, or site path appears anywhere
in the resolver code (Mode B), and Mode A's own launch is fully parameterised
(see above) rather than assuming any of that either.

## Evidence model

Raw match count is deliberately **not** scored. In the reference run, 3,864
locations reduced to 424 representative ones -- counting matches would reward
redundant overlapping signatures rather than architecture. Scoring uses:

- representative locations only (merged, so overlapping signatures count once);
- integrated InterPro accessions, and whether multiple member databases
  corroborate the same entry;
- the fraction of the protein covered by representative locations.

Feature types are not equally weighted: `Family`/`Domain`/`Homologous_superfamily`
are architecture-defining; `Conserved_site`/`Repeat` are supporting;
`Coiled_coil`/`Region` and the disorder/low-complexity predictors (MobiDB-lite,
COILS) are not functional evidence. **AntiFam** matches are treated as evidence
*against* a model being a real protein.

**Absence of InterPro matches never penalises a candidate.** A lineage-specific
Apicomplexan protein with no signature is an expected outcome, not evidence of a
bad model. (Caveat, not a rule change: this is indistinguishable from "this
protein was never scanned" if the candidate's checksum simply is not in the
provided JSONL/GFF3 at all -- see Limitations below.)

## Resolver behaviour: verdicts, reason codes, and replacement

The resolver first computes one of four generic verdicts (`compare_candidates`):
`supports_current`, `supports_runner_up`, `uninformative`, `conflicting`. It
then classifies *why*, as one of ten specific reason codes
(`classify_reason_code`):

| Reason code | Meaning |
| :--- | :--- |
| `INTERPRO_SUPPORTS_INITIAL` | Verdict was `supports_current` -- initial choice stands. |
| `INTERPRO_EQUIVALENT` | Both candidates have domains, architecture is a tie. |
| `INTERPRO_UNINFORMATIVE` | Neither candidate has any representative domain evidence. |
| `INTERPRO_CONFLICTING` | Each candidate has integrated entries the other lacks -- a human should look. |
| `INTERPRO_ANTIFAM_AVOIDED` | Initial choice matches an AntiFam (spurious-ORF) signature the runner-up does not. |
| `INTERPRO_COMPLETE_DOMAIN_RESTORED` | Initial choice has no domain evidence at all; runner-up has a real one. |
| `INTERPRO_N_TERMINAL_TRUNCATION_RESOLVED` | Runner-up's architecture extends further at the N-terminus only. |
| `INTERPRO_C_TERMINAL_TRUNCATION_RESOLVED` | Runner-up's architecture extends further at the C-terminus only. |
| `INTERPRO_ARCHITECTURE_RESTORED` | A genuine architecture improvement that isn't cleanly localised to one terminus. |
| `INTERPRO_FUSION_AVOIDED` | Initial choice is much longer and carries a strong-domain region entirely disjoint from the runner-up's own domain -- a suspected fusion/read-through model (`detect_fusion`, checked *before* the generic verdict, since a fusion can have MORE integrated entries and would otherwise look "stronger"). |

Only the last five are ever replacement *candidates*
(`INTERPRO_SUPPORTS_INITIAL`/`EQUIVALENT`/`UNINFORMATIVE`/`CONFLICTING` never
replace). Every replacement candidate is still gated by `check_safeguards()`
before it takes effect:

- **no worse structural class** -- never replace with an incomplete-ORF
  candidate when the initial choice was complete;
- **no new internal stops** -- never replace with a candidate that introduces
  one;
- **materiality / no silent override of much-stronger protein-validation
  evidence** -- a DIAMOND-balanced-coverage or psauron-score gap above
  `max_protein_validation_override_gap` (default 0.15) in the initial choice's
  favour blocks the replacement, **except** for `INTERPRO_ANTIFAM_AVOIDED` and
  `INTERPRO_FUSION_AVOIDED`, which are safety-critical overrides: a spurious
  ORF or a fusion/read-through model can still score deceptively well on
  DIAMOND/psauron precisely because it contains the real protein's sequence,
  so this specific safeguard must not be allowed to protect them.

`apply_replacements: false` keeps the resolver report-only regardless of
reason code/safeguards -- useful for evaluating the policy on a new dataset
before trusting it to touch anything.

`final = initial` in every case where the reason code is not a replacement
candidate, a safeguard fails, or `apply_replacements` is false.

## Attribution outputs (initial vs recommendation vs final)

Every reviewed gene's decision is distinguished three ways, never
conflated: `initial_canonical_transcript_id` (what canonical_selection alone
picked), `interpro_recommended_transcript_id` (what the resolver's evidence
suggests, whether or not it was applied), and `final_canonical_transcript_id`
(what actually stands after safeguards/`apply_replacements`).

- `interpro_resolver_report.tsv` -- one row per *reviewed* gene, full
  verdict/reason-code/safeguard detail.
- `canonical_decisions.tsv` (when `--canonical-transcripts`/`--consensus-gff3`
  are given) -- one row per **gene in the whole build**, not just reviewed
  ones. Non-reviewed genes get a trivial row
  (`canonical_selection_stage="initial"`, `final == initial`). This is what
  makes "exactly one final canonical per gene" checkable across the entire
  geneset, not just the ambiguous subset.
- `consensus.final_canonical_annotated.gff3` -- a COPY of `consensus.gff3`
  (never modified in place) with `canonical=1/0` (based on the FINAL
  canonical), `canonical_selection_stage=initial|interpro`,
  `initial_canonical=<ID>`, `canonical_reason=<code>` added to every mRNA.
  Kept separate from `canonical_selection`'s own initial-only annotated GFF3
  (`consensus.canonical_annotated.gff3`, from `--annotate-gff3`), so a build
  that never touches this resolver keeps producing that file unchanged.

`interpro_resolver_summary.json`'s `canonical_outputs_modified` is `true` only
when `canonical_decisions.tsv`/`consensus.final_canonical_annotated.gff3` were
written as **new** files this run -- `consensus.gff3` and
`canonical_transcripts.tsv` themselves are never modified in place, by this or
any other GMB module.

## Caching / reuse

The natural cache key is:

```
protein SHA-256 + InterProScan workflow version + InterPro data release
                + selected applications + representative-location setting
```

All five components matter: the same protein scanned against a different
InterPro release, or with a different application set, is not the same result.
Reuse is only safe when every component matches. A local results directory keyed
this way is sufficient; no cross-project cache service is implemented or needed
here.

## Limitations

- **A candidate absent from the InterProScan output is indistinguishable from
  a candidate with no domains at all.** `summarise_architecture()` looks matches
  up by protein checksum and treats "no entry for this checksum" identically to
  "scanned, zero hits" (`has_domains=False`). If the *initial* canonical
  transcript changes (e.g. after a scoring-policy change) to a candidate that
  was never submitted for scanning, the resolver can report
  `INTERPRO_COMPLETE_DOMAIN_RESTORED` purely because of the missing data, not
  because of a genuine biological difference. This was observed concretely on
  the chr1 validation run (see the project report) -- always check
  `current_n_integrated_entries`/`current_representative_coverage` for a
  suspicious `0`/`0.0` before trusting a verdict that recommends moving away
  from the initial choice.
- **`INTERPRO_FUSION_AVOIDED`'s detection heuristic (`detect_fusion`) is
  conservative but not exhaustive**: it looks for a strong-domain interval
  entirely disjoint from the runner-up's own span on a much-longer candidate.
  A fusion whose extra region partially overlaps the runner-up's domain, or
  whose length ratio falls under `length_ratio_flag`, will not be flagged this
  way and will instead be scored by the generic architecture comparison
  (which can, in that case, still favour the fusion candidate on raw entry
  count -- exactly the failure mode this heuristic exists to catch in the
  clear-cut case, not every case).

## Still to confirm on the cluster

None of the following could be verified from this laptop:

- Nextflow version available (local dev used 26.04.6);
- Java version (Nextflow requires a supported JDK);
- Slurm availability, and the correct partition/queue name;
- Singularity or Apptainer version, and which of the two is installed;
- whether the InterProScan container image can be pulled on compute nodes, or
  must be pre-staged to a shared location;
- shared InterPro data directory path and its release version;
- CPU architecture (the local run was arm64; cluster is likely x86-64 -- the
  image must match);
- storage quota for the Nextflow work directory and the InterPro data (the
  local `data/` directory spans 16 member databases);
- required bind mounts for the data and work directories;
- per-job CPU/memory/walltime limits.
