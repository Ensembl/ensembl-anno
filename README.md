# Ensembl Anno

Ensembl Anno is a modular genome annotation framework developed by Ensembl.

It provides a single wrapper (`ensembl_anno.py`) that orchestrates a complete
genome annotation workflow, including:

- Repeat annotation
- Simple feature annotation
- Small non-coding RNA annotation
- Transcriptomic evidence generation
- Protein homology annotation
- Gene set finalisation
- Optional loading into an Ensembl Core database

---

## Documentation

Full documentation is available at

https://ensembl.github.io/ensembl-anno/

The documentation contains

- Installation
- Configuration
- Pipeline overview
- Running the wrapper
- Examples
- Individual analysis modules
- API reference

---

## Quick start

```bash
python ensembl_anno.py \
    --genome_file genome.fa \
    --output_dir output \
    --run_full_annotation