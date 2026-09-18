# GMB quickstart

From a prepared evidence bundle to a handover, in three commands.

---

## 1. Install

```bash
mamba create -n gmb -c conda-forge python=3.12 pandas 'pyranges<=0.1.4' biopython pyyaml matplotlib numpy
mamba activate gmb
pip install -e /path/to/ensembl-anno/support_scripts/gmb
gmb-build --help
```

`pyranges` must be `<=0.1.4` — the API changed after that and it is the most common install
failure. GMB does **not** install Helixer, Tiberius, Scallop, StringTie, Minimap2, DIAMOND or
Psauron.

## 2. What you need

```bash
GENOME=/data/genome.fa          # required
BACKBONE=/data/backbone.gff3    # required in practice (Helixer or Tiberius)
SCALLOP=/data/scallop.gtf       # recommended
STRINGTIE=/data/stringtie.gtf   # recommended
ORTHODB=/data/orthodb.gtf       # recommended
OUT=/work/gmb
```

Every evidence file must use the **same sequence names as the genome FASTA**.
Long-read evidence is optional. A reference annotation is **never** an input.

## 3. Pick a preset

| your situation | preset |
|---|---|
| new or unknown clade | `standard` |
| Apicomplexa with a general-purpose ab initio backbone | `apicomplexa` |
| fungi with a strong Helixer-like backbone | `fungi` |

Details and evidence: `presets.md`.

## 4. Run

```bash
# 1. validate the evidence bundle — exits 1 on FAIL
gmb-preflight --preset fungi \
  --genome "$GENOME" --helixer "$BACKBONE" \
  --scallop "$SCALLOP" --stringtie "$STRINGTIE" --orthodb "$ORTHODB" \
  --output-dir "$OUT/preflight"

# 2. integrate evidence and select gene models
gmb-build --preset fungi \
  --genome "$GENOME" --helixer "$BACKBONE" \
  --scallop "$SCALLOP" --stringtie "$STRINGTIE" --orthodb "$ORTHODB" \
  --gene-prefix XXGMB --output-dir "$OUT/build" --validate-fasta

# 3. regenerate FASTA from the final GFF3, QC, canonical, handover
gmb-finalise --build-dir "$OUT/build" --genome "$GENOME" \
  --output-dir "$OUT/finalise"
```

Use `--tiberius` instead of `--helixer` for a Tiberius backbone. Omit any evidence flag you do
not have. Ready-made wrappers: `examples/run_gene_model_builder.sh`, `.py`, `submit_gmb.slurm`.

## 5. What you get

```
$OUT/finalise/canonical/consensus.canonical_annotated.gff3   <- the annotation
$OUT/finalise/consensus.gff3                                 all isoforms
$OUT/finalise/{cdna,cds,prot}.fa                             sequences
$OUT/finalise/{fasta_qc_report,utr_qc_report}.json           QC evidence
$OUT/build/evidence_attribution.tsv                          why each model won
$OUT/build/run_manifest.json                                 full provenance
```

**Hand over `finalise/`, not `build/`** — `finalise/` is where sequences are regenerated from
the final GFF3. Full list: `output_contract.md`.

## 6. Check it

```bash
python - <<'PY'
import json
qc = json.load(open("$OUT/finalise/fasta_qc_report.json"))
print("QC pass:", qc["pass"], "| failed:", qc.get("failed_checks"))
PY
```

All of these must be zero: cDNA / CDS / protein mismatches, internal stops, gene-boundary
violations, UTR violations. See `qc.md`.

> **A non-zero `gmb-build` exit under `--validate-fasta` is a QC verdict, not a crash.** The
> annotation is fully written; read `build/fasta_qc_report.json`.

> **QC PASS does not mean high biological accuracy.** It means the annotation is internally
> consistent. See `qc.md`.

## 7. Optional — evaluate against a reference

```bash
gmb-compare --query "$OUT/finalise/consensus.gff3" \
  --reference "$REFERENCE" --reference-fasta "$GENOME" \
  --evaluation-mode protein_coding --output-dir "$OUT/comparison"
```

Evaluation only — never part of the production path, and never used to choose a configuration.

---

## Next

| I want to… | read |
|---|---|
| know exactly what to supply | `input_contract.md` |
| know exactly what I get | `output_contract.md` |
| understand presets | `presets.md` |
| configure a new clade | `creating_a_preset.md` |
| wire GMB into a pipeline | `pipeline_integration.md` |
| understand QC | `qc.md` |
| fix something | `troubleshooting.md` |
| reproduce a run | `reproducibility.md` |
