# Output contracts

All output files written by `gmb-build` are listed here with their format,
coverage guarantees, and consistency rules.

---

## Primary outputs

### `consensus.gff3`

GFF3 file with features: `gene`, `mRNA`, `exon`, `CDS`, `five_prime_UTR`,
`three_prime_UTR`.

- One `gene` feature per locus.
- One or more `mRNA` children per gene (up to `scoring.max_isoforms_per_locus`).
- `exon` features cover the spliced transcript region (inclusive of UTR sequence).
- `CDS` features are present only when a coding sequence was inferred. A transcript
  with no `CDS` child is a non-coding or UTR-only model.
- IDs are stable within a run: `{gene_prefix}_G{n}` for genes,
  `{gene_prefix}_T{n}.{isoform}` for transcripts.

### `cdna.fa`

- One record per `mRNA` in `consensus.gff3`.
- Sequence is the spliced cDNA: exonic genomic sequence concatenated in 5′→3′
  transcript order, reverse-complemented for minus-strand transcripts.
- FASTA header: `>{transcript_id}` matching the `ID=` attribute of the mRNA row.

### `cds.fa`

- One record per `mRNA` that has at least one `CDS` child in `consensus.gff3`.
- Sequence is the concatenated coding sequence. Terminal stop codons are stripped.
- Absent when no CDS features were predicted in the run.

### `prot.fa`

- One record per `mRNA` that has at least one `CDS` child in `consensus.gff3`.
- Sequence is the translated protein from the CDS interval(s).
- FASTA header: `>{transcript_id}` — same stable ID as the mRNA row.
- UTR-only transcripts (no CDS rows) are excluded.

**ID guarantee:** Every `>id` in `prot.fa` and `cdna.fa` maps to exactly one
`mRNA` row in `consensus.gff3`. No transcript appears in the FASTA files unless
it survived all GFF post-processing steps (structural validation, deduplication).

---

## Tabular outputs

### `summary.json` / `summary.tsv`

Pipeline metrics: gene counts, filtering statistics, per-evidence-track record
counts. TSV is the transposed form of the JSON.

### `evidence_attribution.tsv`

Per-transcript evidence provenance written immediately after the build.

| Column | Description |
| :----- | :---------- |
| `transcript_id` | Transcript identifier |
| `evidence_sources` | Comma-separated list of named source tracks |
| `gmb_score` | Weighted scoring sum |
| `protein_coding_score` | DIAMOND/Psauron combined score (null when disabled) |

### `protein_validation.tsv`

Written when `protein_validation.enabled: true`.

| Column | Description |
| :----- | :---------- |
| `transcript_id` | Transcript identifier |
| `diamond_hit_id` | Best DIAMOND hit accession (null if no hit) |
| `diamond_bitscore` | DIAMOND bitscore |
| `diamond_qcov` | Query coverage (0–100) |
| `diamond_scov` | Target coverage (0–100) |
| `psauron_score` | Psauron protein-coding score (0–1) |
| `protein_coding_score` | Combined weighted score |
| `calculated` | `true` if computed this run; `false` if reused from cache |

### `collapsed_duplicate_transcripts.tsv`

Written when any exact-duplicate transcripts are collapsed
(`duplicate_transcript_collapse.collapse_exact_duplicates: true`).

| Column | Description |
| :----- | :---------- |
| `gene_id` | Gene containing the duplicates |
| `retained_transcript_id` | Transcript ID kept |
| `collapsed_transcript_ids` | Comma-separated IDs that were collapsed into the retained one |

### `subset_regions.tsv`

Written when `--seqname`, `--region`, or `--regions-file` is used.

| Column | Description |
| :----- | :---------- |
| `seqname` | Selected sequence name |
| `start` | Region start (1-based; full contig if no region was specified) |
| `end` | Region end |

---

## FASTA QC report

`fasta_qc_report.json` is written when `--validate-fasta` is used or
`python -m gmb.pipeline.fasta_qc output/` is run.

It reports:
- Transcript / protein / cDNA record counts.
- Missing or extra FASTA records relative to GFF3 (should both be zero).
- Duplicate FASTA headers.
- Protein and cDNA sequence reconstruction mismatches (when `--genome` is
  supplied).

---

## Canonical selection outputs

Written by `gmb-canonical-selection` to its own `--output-dir`.
See **[docs/canonical_selection.md](docs/canonical_selection.md)**.

## InterPro resolver outputs

Written by `gmb-interpro-resolve` to its own `--output-dir`.
See **[docs/interpro_resolver.md](docs/interpro_resolver.md)**.
