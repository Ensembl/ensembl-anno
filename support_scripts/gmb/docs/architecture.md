# GMB architecture

## Production boundary — what GMB owns

GMB is the **evidence-integration and final gene-model selection stage**. It consumes
evidence another module produced and emits a finished annotation.

```
   upstream evidence-producing modules
                 |
                 v
     well-defined GMB input contract        (docs/input_contract.md)
                 |
                 v
   +-------------------------------+
   |       GENE MODEL BUILDER      |
   |  preflight -> build -> finalise |
   +-------------------------------+
                 |
                 v
   annotation + FASTAs + attribution + QC   (docs/output_contract.md)
```

### GMB owns

| responsibility | where |
|---|---|
| input validation | `gmb/preflight/` |
| evidence normalisation | `gmb/pipeline/builder.py::load_evidence`, `evidence_filter.py` |
| candidate integration | `gmb/pipeline/builder.py` |
| ORF / CDS handling | `gmb/pipeline/annotate_cds_utrs.py` |
| candidate retention | `gmb/pipeline/scoring.py` (retention gate) |
| candidate ranking | `gmb/pipeline/scoring.py` (`rank_tier`) |
| gene construction | `gmb/pipeline/dedup_genes.py` |
| duplicate collapse | `gmb/pipeline/duplicate_transcript_collapse.py` |
| gene-boundary recomputation | `gmb/pipeline/gff3_validate.py::recompute_gene_bounds` |
| sequence regeneration | `gmb/pipeline/builder.py::regenerate_final_fasta` |
| QC | `gmb/pipeline/fasta_qc.py`, `utr_validator.py`, `gff3_validate.py` |
| protein validation | `gmb/pipeline/protein_validation.py` (DIAMOND + Psauron) |
| canonical transcript selection | `gmb/pipeline/canonical_selection.py` |
| handover outputs | `gmb/pipeline/finalise.py` |
| provenance / manifests | `gmb/provenance/` |

### GMB does NOT own

Running the evidence generators. GMB never invokes Helixer, Tiberius, Scallop, StringTie,
Minimap2, GenBlast or an aligner, and does not install them. The one exception is
`gmb-longread-consensus`, which collapses long-read alignments that another tool produced.

DIAMOND and Psauron are *invoked* by GMB for protein validation but are **not installed** by
it; you supply the paths.

### Reference annotations

A reference annotation is **evaluation-only**. `gmb-build`, `gmb-finalise` and `gmb-preflight`
expose no option that accepts one, and no resolved config can reference one. The reference
enters only through `gmb-compare`, which runs after the annotation exists and cannot
influence it.

---

## Module map

The package is `gmb`. Responsibilities are expressed by module, and the existing layout
already separates them; it has been documented rather than reorganised.

```
gmb/
  cli/                 command-line entry points — one module per command
    preflight.py         gmb-preflight    input validation
    build.py             gmb-build        evidence integration + selection
    finalise.py          gmb-finalise     FASTA regeneration, QC, canonical, handover
    compare.py           gmb-compare      EVALUATION ONLY
    longread_consensus.py, canonical_selection.py, interpro_*.py, visualize.py

  configs/             shipped configuration
    standard.yaml        neutral base, always loaded
    apicomplexa.yaml     validated clade preset
    fungi.yaml           validated clade preset
    longread_consensus/  presets for the long-read collapse step

  preflight/           PRE-BUILD INPUT VALIDATION
    checks.py            per-track summaries, splice quality, verdicts

  pipeline/            CONSTRUCTION AND SELECTION
    config.py            schema, layering, deprecated aliases, strict YAML loading
    canonical_evidence.py  EvidenceRoles — the ONLY source-label -> role mapping
    applicability.py     reference-free gates (backbone_intron_rescue auto)
    builder.py           orchestration: load -> filter -> ORF -> select -> emit
    evidence_filter.py   normalisation and filtering of loaded evidence
    scoring.py           score_model, rank_tier, retention gate, select_isoforms
    annotate_cds_utrs.py ORF inference, CDS/UTR annotation, splice-site checks
    dedup_genes.py       gene construction from selected transcripts
    duplicate_transcript_collapse.py
    gff3_validate.py     structural invariants + recompute_gene_bounds
    utr_validator.py     UTR invariant checks
    fasta_export.py      cDNA/CDS/protein emission
    fasta_qc.py          sequence-vs-annotation acceptance checks
    canonical_selection.py  canonical transcript ranking
    finalise.py          handover assembly
    protein_validation.py   DIAMOND + Psauron
    longread/            long-read consensus collapse

  provenance/          RUN MANIFEST
    manifest.py          versions, input hashes, resolved policy, runtime

  compare/             EVALUATION ONLY — not part of the production path
  utils/               fasta, gff, intervals, io, logging
```

### Why this layout rather than the notional `evidence/ validation/ qc/ canonical/ io/`

Those responsibilities exist and are separated — they live as named modules inside
`pipeline/` rather than as sibling packages. Splitting them into new top-level packages would
move ~20 files and rewrite every import for no behavioural gain, against the explicit
instruction not to reorganise for its own sake. Two genuinely new capabilities did get their
own packages, because nothing existed to put them in: `preflight/` and `provenance/`.

---

## The three stages

### 1. `gmb-preflight` — validate the evidence bundle

Reference-free. Reports seqid compatibility, ID collisions, strand completeness, structural
distributions, **per-track canonical splice fraction**, resolved roles and weights, and which
optional policies will actually fire. Exits non-zero on FAIL.

### 2. `gmb-build` — integrate evidence and select models

```
load evidence  ->  namespace colliding IDs per seqid
               ->  filter / split
               ->  infer ORFs, annotate CDS and UTRs
               ->  resolve applicability gates (once per run)
               ->  cluster loci
               ->  score, retain, rank, select
               ->  construct genes, dedup, collapse duplicates
               ->  recompute gene bounds from surviving transcripts
               ->  regenerate FASTA from the final GFF3
               ->  FASTA QC
```

### 3. `gmb-finalise` — produce the handover

Regenerates cDNA/CDS/protein **from the final GFF3** so sequence and annotation cannot drift,
re-runs QC, selects one canonical transcript per gene, writes the handover manifest.

**`finalise/` is the production output. `build/` is intermediate.**

---

## Two design rules worth knowing before changing anything

### Evidence roles, not tool names

Selection logic never tests a literal tool name. `EvidenceRoles` is the single mapping from
source label to role (`backbone`, `short_read_transcriptomic`, `long_read_transcriptomic`,
`protein_alignment`, `protein_validation`), and everything downstream — ranking, retention,
corroboration counting, **and numeric weights** — acts on the role.

This is why an assembler called `IsoQuant` or `AssemblerX` behaves exactly like `StringTie`
once listed under `shortread_labels`. `tests/test_production_contract.py` enforces it using
invented tool names.

### Correctness is not policy

Two categories, deliberately never mixed in the same flag:

**Always on, not configurable** — strand/end-aware UTR correction, exon reconstruction,
FASTA regeneration from the final GFF3, strict sequence QC, exact duplicate collapse,
same-strand protein support, protein evidence attribution, cross-seqid ID namespacing,
gene-boundary recomputation, gene/transcript containment validation, one canonical transcript
per gene, run provenance.

**Explicit and configurable** — structural corroboration, CDS-compatible protein support,
long-read structural guard, long-read disposition, backbone intron rescue, numeric weights.

A correctness property is never expressed as a tuning switch, and a tuning switch never
silently changes a correctness property.
