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

## Dependencies

### Software

| Software | Version | URL| Required |
|----------|---------|----|----------|


### Python EnsEMBL repositories you need to have

| Repository name | branch | URL|
|-----------------|--------|----|
| ensembl-genes | default | https://github.com/Ensembl/ensembl-genes.git |


### Perl EnsEMBL repositories you need to have

| Repository name | branch | URL|
|-----------------|--------|----|
| ensembl | default | https://github.com/Ensembl/ensembl.git |
| ensembl-analysis | experimental/gbiab | https://github.com/Ensembl/ensembl-analysis.git | (need to make sure dependencies are on main and update this to main/default for branch)
| ensembl-io | default | https://github.com/Ensembl/ensembl-io.git |
| ensembl-taxonomy | default | https://github.com/Ensembl/ensembl-taxonomy.git |
| ensembl-variation | default | https://github.com/Ensembl/ensembl-variation.git |


### Python virtual environment


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

