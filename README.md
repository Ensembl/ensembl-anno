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
```

## Proposed File Structure

```
.
├── main.nf                       # workflow entry point
├── workflows/                    # workflows
├── subworkflows/                 # subworkflows (eg. long reads, repeats, etc)
├── modules/                      # tools (eg. Dust, minimap2 etc)                 
├── bin/                          # Python code required to run anno (happy to move in future)
├── conf/                         # NF configs
├── tests/                        # Testing directory (useful for CI)
    └──nf/                          # Nextflow tests
        ├── data/                       # nf testing data 
        └── tests/                      # Nextflow tests (system testing)
    └── python                      # Python tests
        ├── data/                       # Any python test inputs
        └── tests/                      # Python tests (unit testing)    
└── docs/                         # All documentation
└── .github/                      # Github actions (CI)
```