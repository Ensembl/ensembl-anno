# Gene Model Builder — pipeline handover

This is the one document to read before running GMB in the Ensembl annotation pipeline. It
assumes no knowledge of the development work. Deeper references are linked where needed.

## 1. What GMB does

GMB is the **evidence-integration and final gene-model-selection stage**. Given a genome, an ab
initio backbone (Helixer), assembled transcripts and protein alignments, it chooses the
best-supported transcript structures at each locus and writes a finished annotation with
regenerated sequences, per-model evidence attribution, QC and provenance. It does not generate
evidence and does not install the tools that do. A reference annotation is **evaluation-only**
and can never influence model selection.

```
gmb-longread-consensus      gmb-preflight     →     gmb-build     →     gmb-finalise
(only if long reads;        validate inputs        integrate          regenerate FASTA from the
 collapse per-read          (exit 1 on FAIL)       evidence, select   final GFF3, QC, canonical
 alignments)                                       models             transcript, handover
```

## 2. Install

```bash
mamba create -n gmb -c conda-forge python=3.12 'pandas>=2,<3' 'pyranges<=0.1.4' pyyaml numpy
mamba activate gmb
pip install /path/to/ensembl-anno/support_scripts/gmb      # or: pip install -e for development
gmb-build --help && gmb-build --list-presets               # apicomplexa, fungi
```

`pyranges` must be `<=0.1.4`. Record the commit you install (`git rev-parse HEAD`). The
experimental overlay (§3) is a file in the repository, `configs/`, **not** part of the
installed package: reference it from the checkout or copy it into the pipeline's configuration.

## 3. Two supported fungal configurations

| | **A. Current fungal baseline** | **B. Experimental evidence-integrated** |
|---|---|---|
| options | `--preset fungi` | `--preset fungi --config configs/fungi_evidence_integration.experimental.yaml` |
| status | validated default, unchanged by the evidence-tuning work | opt-in; validated on 17 *Z. tritici* windows; to be decided by the production-scale test |

Configuration B changes four settings, and nothing else:

| setting (B) | what changes, scientifically |
|---|---|
| `scoring.primary_selection: junction_supported` | Under A, where Helixer predicts a gene its model is always the primary transcript: no structure supported only by transcripts can outscore it. B promotes an RNA-seq-derived alternate when Helixer's model has an intron that **no** transcript observes and the alternate is a complete ORF covering ≥ 90% of Helixer's CDS at the same coding locus, with spliced-protein agreement no worse. Helixer's model is kept as an alternate. |
| `scoring.isoform_cds_overlap: drop` | Under A, a transcript joins a gene by sharing an intron or > 15% of span, so a neighbouring ORF reached through a UTR or read-through becomes an "isoform"; 41% of A's alternates share no coding base with their primary, causing most gene merges. B discards such structures (listed in `build/excluded_candidates.tsv`) and never discards a Helixer model. |
| `canonical_selection.prefer_build_primary: true` | Under A, canonical selection counts Scallop and StringTie (two assemblies of the same reads) as two sources and overrides the build's choice; on held-out windows this made the canonical set *less* accurate than Helixer alone. B keeps the build primary unless ORF integrity, protein validation or evidence-class breadth says otherwise. |
| `scoring.min_cds_bp: 150` | Drops models whose CDS is < 150 bp. No *Z. tritici* reference CDS is shorter; the models removed are almost all unsupported novel loci. (Under A the setting is 0, off.) |

Everything else — filters, weights, clustering, UTR rules, QC — is identical. Results on held-out
*Z. tritici* windows, canonical transcripts, core chromosomes (1,874 reference genes):

| run | exact CDS | exact CDS intron chain | splits / merges | F1 |
|---|---|---|---|---|
| Helixer alone | 710 | 599 | 64 / 23 | 0.336 |
| A | 695 | 589 | 88 / 22 | 0.316 |
| B | 733 | 636 | 63 / 23 | 0.335 |
| B + collapsed long reads | 754 | 668 | 63 / 20 | 0.342 |

Full evaluation: [`docs/release/evidence_integration.md`](docs/release/evidence_integration.md).

## 4. Inputs and their requirements

Every file must use the **same sequence names as the genome FASTA**. Full contract:
[`docs/input_contract.md`](docs/input_contract.md).

| input | flag | requirement |
|---|---|---|
| genome FASTA | `--genome` | required; sequence names as in every evidence file |
| Helixer GFF3 | `--helixer` | required; must carry CDS. Helixer output for an NCBI assembly names sequences by GenBank accession (`CM001196.1`): rename first with `tools/remap_helixer.py` (§5); preflight FAILs the track otherwise |
| protein alignments (genBlastG GTF) | `--orthodb`, `--uniprot` | recommended. The Ensembl anno `*_genblast_annotation.gtf` is UniProt: pass it as `--uniprot`. Proteins support and veto candidates; they never create gene models |
| short-read assemblies | `--scallop`, `--stringtie` | recommended; **stranded**; GMB never re-derives strand |
| long-read transcript models | `--minimap2` | optional, but the strongest evidence in the evaluation. **Collapsed consensus models only** (§5). A per-read track triggers `longread_collapsed` in preflight and a `WARNING` in the build; neither refuses it |
| DIAMOND / Psauron protein validation | config | optional and **off** in both configurations; see §8 |

## 5. Commands (verified on the bundled fixture and *Z. tritici* windows)

```bash
GMB=/path/to/ensembl-anno/support_scripts/gmb
EV=/data/GCA_000219625.1                       # Ensembl anno evidence outputs (prefix)
G=/data/genome.fa                               # sequence names 1..21
OUT=/work/gmb

# 1. Helixer: rename sequences to the genome's names
python $GMB/tools/remap_helixer.py --input ${EV}_helixer.gff3 \
  --assembly-report ${EV}_MYCGR_v2.0_assembly_report.txt --output $OUT/helixer.gff3

# 2. Long reads (optional): stranded Minimap2 alignments -> consensus transcript models.
#    Strand must already be correct upstream (Ensembl anno sets it from the Minimap2 ts tag).
gmb-longread-consensus --input ${EV}_minimap2_annotation.gtf --preset fungi \
  --output-dir $OUT/longread
LR=(--minimap2 $OUT/longread/minimap2_consensus.gtf)    # or LR=() without long reads

EVID=(--genome $G --helixer $OUT/helixer.gff3
      --scallop ${EV}_scallop_annotation.gtf --stringtie ${EV}_stringtie_annotation.gtf
      --orthodb ${EV}_orthodb_annotation.gtf --uniprot ${EV}_genblast_annotation.gtf "${LR[@]}")

# Configuration A
gmb-preflight --preset fungi "${EVID[@]}" --output-dir $OUT/A/preflight
gmb-build     --preset fungi "${EVID[@]}" --gene-prefix ZTGMB --output-dir $OUT/A/build --validate-fasta
gmb-finalise  --build-dir $OUT/A/build --genome $G --output-dir $OUT/A/finalise

# Configuration B: identical, plus the overlay
X=$GMB/configs/fungi_evidence_integration.experimental.yaml
gmb-preflight --preset fungi --config $X "${EVID[@]}" --output-dir $OUT/B/preflight
gmb-build     --preset fungi --config $X "${EVID[@]}" --gene-prefix ZTGMB --output-dir $OUT/B/build --validate-fasta
gmb-finalise  --build-dir $OUT/B/build --genome $G --output-dir $OUT/B/finalise
```

- `gmb-preflight` exits 1 on FAIL: fix the input. Expected WARNs on the *Z. tritici* bundle:
  unstranded rows in the short-read tracks (0.4% / 1.7%), reused OrthoDB IDs, very long OrthoDB
  spans and introns. Check that the long-read track's `splice_quality` passes (≥ 85% canonical;
  the collapsed *Z. tritici* track passed: 95.5% GT-AG over the 9 development windows, 98.6% by
  preflight's count on one window).
- Always pass the same `--preset`/`--config` to preflight and build. `gmb-finalise` reads the
  build's `resolved_config.yaml`, so it needs no configuration options.
- `--validate-fasta` makes a QC failure exit non-zero: that is a QC verdict, the output is
  written.
- A smoke test of both configurations on the bundled fixture, with independent output checks:
  `validation/smoke_fixture.sh <outdir> [<baseline checkout>]`.

## 6. Outputs

**Hand over `finalise/`. `build/` is intermediate** — its GFF3/FASTA are pre-finalisation
copies; sequences in `finalise/` are regenerated from the final GFF3.

| file | what it is |
|---|---|
| `finalise/canonical/consensus.canonical_annotated.gff3` | **the annotation**: all genes, all isoforms, one transcript per gene tagged `Ensembl_canonical` (the canonical-only gene set = tagged transcripts) |
| `finalise/consensus.gff3` | the full gene set, all isoforms, untagged |
| `finalise/cdna.fa`, `cds.fa`, `prot.fa` | sequences regenerated from the final GFF3 |
| `finalise/canonical/canonical_transcripts.tsv` | canonical choice per gene; `deciding_tier` says which rule decided |
| `finalise/fasta_qc_report.json`, `utr_qc_report.json` | QC evidence (must pass) |
| `finalise/handover_manifest.json` | output inventory with checksums and the configuration used |
| `build/evidence_attribution.tsv` | per transcript: evidence sources, protein support, selection reason, junction support, ORF integrity, locus context |
| `build/excluded_candidates.tsv` | structures refused by `isoform_cds_overlap: drop` (B) |
| `build/resolved_config.yaml` (+ `_sha256`) | the exact configuration |
| `build/run_manifest.json` | GMB version, input hashes, resolved policy, runtime |
| `build/summary.json` | run counters |

Transcript-only loci (`gene_has_backbone_model = False` in the attribution) and antisense coding
predictions among them (`opposite_strand_cds_overlap = True`) can be assessed separately (§8).
Full contract: [`docs/output_contract.md`](docs/output_contract.md).

## 7. First production-scale test

Run Helixer alone, A, B, B + long reads and A + long reads on the whole *Z. tritici* genome, and
compare them with ensembl-genes `annotation-qc` against fixed promotion criteria:
[`docs/release/production_scale_test_plan.md`](docs/release/production_scale_test_plan.md).

## 8. Known limitations to watch

Details and evidence: [`docs/known_issues.md`](docs/known_issues.md) and
[`docs/release/handover_review.md`](docs/release/handover_review.md).

1. **Validated regionally only.** Both configurations: 17 windows (38% of the genome). A
   genome-wide run has not been made with the current code; runtime and memory at genome scale
   are unmeasured (earlier builds: 29 h with protein validation).
2. **Transcript-only antisense coding predictions.** About 5% of genes are RNA-seq-only models
   antisense to another gene, with incidental ORFs (low coding potential). Neither configuration
   filters them; they are flagged in the attribution for separate assessment.
3. **Proteins never construct models**, and an assembled transcript yields one ORF. Most
   reference genes Helixer misses lie inside read-through assemblies and cannot be recovered.
4. **Protein validation** with the fungi weights requires a DIAMOND hit for every model and
   costs accuracy in configuration B (prunes the alternates the junction rule promotes from).
   Leave it off for the first run.
5. **Isoforms:** about half of B's alternates differ from the primary only in UTR.
6. **Accessory chromosomes** (14–21) have far weaker evidence and accuracy; report them
   separately.

## 9. If something goes wrong

QC failed or unsure what a pass means → [`docs/qc.md`](docs/qc.md). Anything else →
[`docs/troubleshooting.md`](docs/troubleshooting.md). A QC pass proves internal consistency,
not biological accuracy.

Development evidence: `validation/z_tritici_evidence_integration/` reproduces the 17-window
evaluation from the original inputs; `docs/release/` holds the reports.
