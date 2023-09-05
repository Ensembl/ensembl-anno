# anno

Toolkit for annotation

## In this repo

### ensembl_anno.py

### support_files

- alignscore.txt (used in run_genblast_align)

### support_classes (should be moved to https://github.com/Ensembl/ensembl-genes-api)

- exon.py
- gene.py
- gtf_adapter.py
- intron.py
- sequence.py
- transcript.py

### support_scripts

- subsample_fastq.py

### support_scripts_perl (should be re-written in python, when API has enough functionality, and added to support_scripts)

- clean_geneset.pl
- clean_utrs_and_lncRNAs.pl
- finalise_geneset.pl
- gtf_to_seq.pl
- load_gtf_ensembl.pl
- select_best_transcripts.pl

# Dependencie: Python EnsEMBL repositories 

| Repository name | branch | URL|
|-----------------|--------|----|
| ensembl-genes | default | https://github.com/Ensembl/ensembl-genes.git |


# Dependencies: Perl EnsEMBL repositories

| Repository name | branch | URL|
|-----------------|--------|----|
| ensembl | default | https://github.com/Ensembl/ensembl.git |
| ensembl-analysis | experimental/gbiab | https://github.com/Ensembl/ensembl-analysis.git | (need to make sure depencies are on main and update this to main/default for branch)
| ensembl-variation | default | https://github.com/Ensembl/ensembl-variation.git |


### Python virtual environment

# Dependencies

## Included Softwares

### DustMasker

DustMasker is a program that identifies and masks out low complexity parts of a genome using an improved DUST algorithm.
The main advantages of the new algorithm are symmetry with respect to taking reverse complements, context insensitivity, and much better performance.

**Citation:** [Morgulis et al., "A Fast and Symmetric DUST Implementation to Mask Low-Complexity DNA Sequences"](https://pubmed.ncbi.nlm.nih.gov/16796549/)

### RepeatMasker

RepeatMasker screens DNA sequences for interspersed repeats and low complexity DNA sequences.

**Citation:** [Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-4.0. 2013-2015](http://www.repeatmasker.org)

### Tandem Repeats Finder

Tandem Repeats Finder locates and displays tandem repeats in DNA sequences.

**Citation:** [Benson, "Tandem repeats finder: a program to analyze DNA sequences"](https://pubmed.ncbi.nlm.nih.gov/9862982/)

### Red

Red is a repeat-detection tool capable of labeling its training data and training itself automatically on an entire genome.

**Citation:** [Girgis, "Red: an intelligent, rapid, accurate tool for detecting repeats de-novo on the genomic scale"](https://doi.org/10.1186/s12859-015-0654-5)

### Eponine

Eponine is a probabilistic method for detecting transcription start sites (TSS) in mammalian genomic sequences.

**Citation:** [Down & Hubbard, "Computational detection and location of transcription start sites in mammalian genomic DNA"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC155284/)

### CpG

CpG is a set of discriminant functions that can recognize structural and compositional features
such as CpG islands, promoter regions and first splice-donor sites.

**Citation:** [Davuluri et al., "Computational identification of promoters and first exons in the human genome"](https://pubmed.ncbi.nlm.nih.gov/11726928/)

### tRNAscan-SE

tRNAscan-SE identifies transfer RNA genes in DNA sequences.

**Citation:** [Lowe & Eddy, "tRNAscan-SE: a program for improved detection of transfer RNA genes in genomic sequence"](https://pubmed.ncbi.nlm.nih.gov/9023104/)

### cmsearch

Infernal and its "cmsearch" tool are used for detecting sncRNAs in sequence databases.

**Citation:** [Nawrocki et al., "Infernal 1.0: inference of RNA alignments"](https://academic.oup.com/bioinformatics/article/25/10/1335/270663)

### GenBlast

GenBlast identifies homologous gene sequences in genomic databases. One of the key features of GenBlast is its flexibility to handle comparative genomics tasks and accurately identify homologs even when the sequences have undergone significant evolutionary changes. This capability makes it a valuable resource for researchers studying gene evolution, gene families, and gene function across diverse species.

**Citation:** [She et al., "GenBlastA: enabling BLAST to identify homologous gene sequences"](https://pubmed.ncbi.nlm.nih.gov/18838612/)

### STAR

STAR-Spliced Transcripts Alignment to a Reference is an ultrafast universal RNA-seq aligner for aligning RNA-seq data to a reference genome.

**Citation:** [Dobin et al., "STAR: ultrafast universal RNA-seq aligner"](https://pubmed.ncbi.nlm.nih.gov/23104886/)

### StringTie

StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. It uses a novel network flow algorithm as well as an optional de novo assembly step to assemble and quantitate full-length transcripts representing multiple splice variants for each gene locus.

**Citation:** [Pertea M et al., "StringTie enables improved reconstruction of a transcriptome from RNA-seq reads"](https://www.nature.com/articles/nbt.3122)

### Scallop

Scallop is a high-performance tool designed for the accurate and efficient quantification of transcriptome assembly. It's capable of handling large-scale transcriptomic data while providing precise estimates of transcript abundances. Scallop's algorithmic approach allows it to efficiently reconstruct transcript structures and quantify their expression levels, making it a valuable resource for studying gene expression and transcriptome analysis.

**Citation:** [Shao et al., Accurate assembly of transcripts through phase-preserving graph decomposition.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5722698/)


## Usage

Please refer to the individual software modules for specific instructions on how to use each tool. [Link to Documentation](docs/build/index.html)






## Running Anno

```
python3 ensembl_anno.py -h

        [-h] [--output_dir OUTPUT_DIR] --genome_file
        GENOME_FILE [--num_threads NUM_THREADS] [--run_masking]
        [--red_path RED_PATH] [--genblast_path GENBLAST_PATH]
        [--convert2blastmask_path CONVERT2BLASTMASK_PATH]
        [--makeblastdb_path MAKEBLASTDB_PATH] [--run_genblast]
        [--run_busco] [--protein_file PROTEIN_FILE]
        [--busco_protein_file BUSCO_PROTEIN_FILE]
        [--rfam_accessions_file RFAM_ACCESSIONS_FILE]
        [--run_star] [--star_path STAR_PATH]
        [--max_reads_per_sample [MAX_READS_PER_SAMPLE]]
        [--max_total_reads [MAX_TOTAL_READS]]
        [--short_read_fastq_dir SHORT_READ_FASTQ_DIR]
        [--max_intron_length [MAX_INTRON_LENGTH]]
        [--run_minimap2] [--minimap2_path]
        [--paftools_path PAFTOOLS_PATH]
        [--long_read_fastq_dir LONG_READ_FASTQ_DIR]
        [--run_augustus] [--augustus_path AUGUSTUS_PATH]
        [--run_stringtie] [--run_scallop]
        [--stringtie_path STRINGTIE_PATH]
        [--scallop_path SCALLOP_PATH]
        [--subsample_script_path SUBSAMPLE_SCRIPT_PATH]
        [--samtools_path SAMTOOLS_PATH] [--finalise_geneset]
        [--db_details DB_DETAILS] [--run_cmsearch] [--run_trf]
        [--trf_path TRF_PATH] [--run_dust]
        [--dust_path DUST_PATH] [--run_repeatmasker]
        [--repeatmasker_path] [--run_trnascan]
        [--trnascan_path TRNASCAN_PATH]
        [--trnascan_filter_path TRNASCAN_FILTER_PATH]
        [--run_cpg] [--cpg_path CPG_PATH] [--run_eponine]
        [--eponine_path EPONINE_PATH] [--java_path JAVA_PATH]
        [--run_full_annotation] [--run_repeats]
        [--run_simple_features] [--run_sncrnas]
        [--run_transcriptomic] [--run_proteins]
        [--diamond_validation_db DIAMOND_VALIDATION_DB]
        [--validation_type VALIDATION_TYPE]
        [--load_to_ensembl_db] [--trim_fastq]
        [--delete_pre_trim_fastq]
        [--repeatmasker_library REPEATMASKER_LIBRARY]
        [--repeatmasker_species REPEATMASKER_SPECIES]
```

## Anno Run Options

### Run full annotation
```
--run_full_annotation
```

Includes:

```
--run_repeats
```
| Analysis | Software | Input | Output |
|----------|----------|-------|--------|
| --run_masking | red |  genome_file, red_path, work_dir | <work_dir>/redt_output/mask_output/<genome_file_name>.msk |
| --run_dust | dust | genome_file, dust_path, work_dir, num_threads | <work_dir>/dust_output/<seq_region_name>.dust.gtf (per region) |
| --run_trf | trf | genome_file, trf_path, work_dir, num_threads | <work_dir>/trf_output/<seq_region_name>.trf.gtf (per region) |

```
--run_simple_features
```
| Analysis | Software | Input | Output |
|----------|----------|-------|--------|
| --run_cpg | cpg | genome_file, cpg_path, work_dir, num_threads | <work_dir>/cpg_output/<seq_region_name>.cpg.gtf (per region) |
| --run_eponine | eponine | genome_file, java_path, eponine_path, work_dir, num_threads | <work_dir>/cpg_output/<seq_region_name>.epo.gtf (per region) |

```
--run_sncrnas
```
| Analysis | Software | Input | Output |
|----------|----------|-------|--------|
| --run_trnascan | trnascan | genome_file, trnascan_path, trnascan_filter_path, work_dir, num_threads | <work_dir>/trnascan_output/<seq_region_name>.trna.gtf (per region) |
| --run_cmsearch* | cmsearch | genome_file, cmsearch_path, rfam_cm_db_path, rfam_seeds_file_path, rfam_accession_file, work_dir, num_threads | <work_dir>/rfam_output/<seq_region_name>.rfam.gtf (per region) |

*requires --rfam_accessions_file, will skip if missing

```
--run_transcriptomic
```
**short-read data**

| Analysis | Software | Input | Output |
|----------|----------|-------|--------|

**long-read data**

| Analysis | Software | Input | Output |
|----------|----------|-------|--------|

```
--run_proteins
```
| Analysis | Software | Input | Output |
|----------|----------|-------|--------|

```
--finalise_geneset
```
| Analysis | Software | Input | Output |
|----------|----------|-------|--------|

### Other analyses

```
--run_repeatmasker
```

