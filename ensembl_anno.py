# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# standard library
import argparse
import gc
import glob
import io
import math
import multiprocessing
import os
import pathlib
import random
import re
import shutil
import signal
import subprocess
import tempfile

from typing import Union

# project imports

from repeatmasking_utils import (
    run_repeatmasker_regions,
    multiprocess_repeatmasker,
    create_repeatmasker_gtf,
    run_dust_regions,
    multiprocess_dust,
    create_dust_gtf,
    run_trf_repeats,
    multiprocess_trf,
    create_trf_gtf,
    run_red,
    create_red_gtf,
    )

from simple_features_utils import (
    run_eponine_regions,
    multiprocess_eponine,
    create_eponine_gtf,
    run_cpg_regions,
    multiprocess_cpg,
    create_cpg_gtf,
    )

from sncRNA_utils import (
    run_trnascan_regions,
    multiprocess_trnascan,
    create_trnascan_gtf,
    run_cmsearch_regions,
    multiprocess_cmsearch,
    get_rfam_seed_descriptions,
    extract_rfam_metrics,
    parse_rfam_tblout,
    remove_rfam_overlap,
    filter_rfam_results,
    create_rfam_gtf,
    check_rnafold_structure,
)


from  protein_utils import (
    run_genblast_align,
    multiprocess_genblast,
    generate_genblast_gtf,
    split_protein_file,
    run_convert2blastmask,
    run_makeblastdb,
    )


from transcriptomic_utils import (
    run_trimming,
    multiprocess_trim_galore,
    create_paired_paths,
    run_star_align,
    run_subsample_script,
    check_for_fastq_subsamples,
    run_minimap2_align,
    bed_to_gtf,
    bed_to_gff,
    bed_to_exons,
    check_transcriptomic_output,
    augustus_output_to_gtf,
    run_augustus_predict,
    generate_hints,
    multiprocess_augustus_hints,
    multiprocess_augustus_id,
    create_slice_hints_file,
    run_stringtie_assemble,
    run_scallop_assemble,
    model_builder,
)
from utils import (
    add_log_file_handler,
    check_exe,
    check_file,
    check_gtf_content,
    create_dir,
    get_seq_region_lengths,
    logger,
    prlimit_command,
    load_results_to_ensembl_db,
    generic_load_records_to_ensembl_db,
    multiprocess_load_records_to_ensembl_db,
    batch_gtf_records,
    run_find_orfs,
    find_orf_phased_region,
    slice_output_to_gtf,
    convert_gff_to_gtf,
    set_attributes,
    create_slice_ids,
    update_gtf_genes,
    read_gtf_genes,
    fasta_to_dict,
    splice_junction_to_gff,
    split_genome,
    multiprocess_generic,
    reverse_complement,
    get_seq_region_names,
    slice_genome,
    subprocess_run_and_log,
    get_sequence,
    seq_region_names,
#    coallate_results,
)

from finalisation_utils import (
    run_finalise_geneset,
    validate_coding_transcripts,
    diamond_validation,
    multiprocess_diamond,
    read_rnasamba_results,
    read_cpc2_results,
    read_diamond_results,
    combine_results,
    merge_finalise_output_files,
    multiprocess_finalise_geneset,
)


def main():
    """
    main function
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output_dir",
        type=str,
        help="Path where the output and temp files will write to. Uses current dir by default",
    )
    parser.add_argument(
        "--genome_file", type=str, required=True, help="Path to the fasta genome file"
    )
    parser.add_argument(
        "--num_threads", type=int, default=1, help="Number of threads to use"
    )
    parser.add_argument(
        "--run_masking",
        action="store_true",
        help="Run Red to find repeats and softmask the genome. Otherwise provide a softmasked genome",
    )
    parser.add_argument(
        "--red_path",
        type=str,
        help="Path to Red executable. See http://toolsmith.ens.utulsa.edu",
    )
    parser.add_argument(
        "--genblast_path",
        type=str,
        help="Path to GenBlast executable. See http://genome.sfu.ca/genblast/download.html",
    )
    parser.add_argument(
        "--convert2blastmask_path",
        type=str,
        help="Path to convert2blastmask executable",
    )
    parser.add_argument(
        "--makeblastdb_path", type=str, help="Path to makeblastdb executable"
    )
    parser.add_argument(
        "--run_genblast",
        action="store_true",
        help="Run GenBlast to align protein sequences",
    )
    parser.add_argument(
        "--genblast_timeout",
        type=int,
        default=10800,
        help="GenBlast timeout in seconds",
    )
    parser.add_argument(
        "--run_busco",
        action="store_true",
        help="Run GenBlast to align BUSCO protein sequences",
    )
    parser.add_argument(
        "--protein_file",
        type=str,
        help="Path to a fasta file with protein sequences",
    )
    parser.add_argument(
        "--busco_protein_file",
        type=str,
        help="Path to a fasta file with BUSCO protein sequences",
    )
    parser.add_argument(
        "--rfam_accessions_file",
        type=str,
        help="Path to a file with Rfam CM accessions, one accession per line, to use with cmsearch",
    )
    parser.add_argument(
        "--run_star", action="store_true", help="Run STAR for short read alignment"
    )
    parser.add_argument(
        "--star_path", type=str, help="Path to STAR for short read alignment"
    )
    parser.add_argument(
        "--max_reads_per_sample",
        nargs="?",
        type=int,
        default=0,
        help="The maximum number of reads to use per sample. Default=0 (unlimited)",
    )
    parser.add_argument(
        "--max_total_reads",
        nargs="?",
        type=int,
        default=0,
        help="The maximum total number of reads. Default=0 (unlimited)",
    )
    parser.add_argument(
        "--short_read_fastq_dir",
        help="Path to short read fastq dir for running with STAR",
    )
    parser.add_argument(
        "--max_intron_length",
        nargs="?",
        type=int,
        default=100_000,
        help="The maximum intron size for alignments. Default=100,000",
    )
    parser.add_argument(
        "--run_minimap2",
        action="store_true",
        help="Run minimap2 for long read alignment",
    )
    parser.add_argument(
        "--minimap2_path",
        type=str,
        help="Path to minimap2 for long read alignment",
    )
    parser.add_argument(
        "--paftools_path",
        type=str,
        help="Path to paftools for SAM to BED conversion",
    )
    parser.add_argument(
        "--long_read_fastq_dir",
        type=str,
        help="Path to long read fastq dir for running with minimap2",
    )
    parser.add_argument(
        "--run_augustus",
        action="store_true",
        help="Run Augustus with hints for gene/transcript prediction",
    )
    parser.add_argument("--augustus_path", type=str, help="Path to Augustus")
    parser.add_argument(
        "--run_stringtie",
        action="store_true",
        help="Run Stringtie on the results from the STAR alignments",
    )
    parser.add_argument(
        "--run_scallop",
        action="store_true",
        help="Run Scallop on the results from the STAR alignments",
    )
    parser.add_argument("--stringtie_path", type=str, help="Path to Stringtie")
    parser.add_argument("--scallop_path", type=str, help="Path to Scallop")
    parser.add_argument(
        "--subsample_script_path",
        type=str,
        help="Path to ensembl-anno subsampling script",
    )
    parser.add_argument("--samtools_path", type=str, help="Path to subsampling script")
    parser.add_argument(
        "--finalise_geneset",
        action="store_true",
        help="Used to finalise the gene set from the various GTF files generated",
    )
    parser.add_argument(
        "--db_details",
        type=str,
        help="A comma separated string of dbname,host,port,user,pass",
    )
    parser.add_argument(
        "--run_cmsearch",
        action="store_true",
        help="Search for sncRNA structures using Rfam and cmsearch",
    )
    parser.add_argument(
        "--run_trf", action="store_true", help="Run TRF to find tandem repeats"
    )
    parser.add_argument("--trf_path", type=str, help="Path to TRF")
    parser.add_argument(
        "--run_dust",
        action="store_true",
        help="Run Dust to find low complexity regions",
    )
    parser.add_argument("--dust_path", type=str, help="Path to Dust")
    parser.add_argument(
        "--run_repeatmasker",
        action="store_true",
        help="Run RepeatMasker to find repeat regions",
    )
    parser.add_argument(
        "--repeatmasker_path", action="store_true", help="Path to RepeatMasker"
    )
    parser.add_argument(
        "--run_trnascan", action="store_true", help="Run tRNAscan-SE to find tRNAs"
    )
    parser.add_argument("--trnascan_path", type=str, help="Path to tRNAscan-SE")
    parser.add_argument(
        "--trnascan_filter_path",
        type=str,
        help="Path to tRNAscan-SE high confidence filter",
    )
    parser.add_argument(
        "--run_cpg", action="store_true", help="Run cpg_lh to find CpG islands"
    )
    parser.add_argument("--cpg_path", type=str, help="Path to cpg_lh")
    parser.add_argument(
        "--run_eponine",
        action="store_true",
        help="Run Eponine to find transcription start sites",
    )
    parser.add_argument("--eponine_path", type=str, help="Path to Eponine jar file")
    parser.add_argument(
        "--java_path", type=str, help="Path to Java for use with Eponine"
    )
    parser.add_argument(
        "--run_full_annotation",
        action="store_true",
        help="Run a full annotation, will automatically check for input data and run tools based on that",
    )
    parser.add_argument("--run_repeats", action="store_true", help="Run Red, Dust, TRF")
    parser.add_argument(
        "--run_simple_features", action="store_true", help="Run CpG, Eponine"
    )
    parser.add_argument(
        "--run_sncrnas", action="store_true", help="Run Rfam, tRNAscan-SE"
    )
    parser.add_argument(
        "--run_transcriptomic",
        action="store_true",
        help="Run STAR, Stringtie2, Scallop, minimap2 (if short_read_fastq_dir and/or long_read_fastq_dir are provided)",
    )
    parser.add_argument(
        "--run_proteins",
        action="store_true",
        help="Run GenBlast if protein_file and/or busco_protein_file have also been set",
    )
    parser.add_argument(
        "--diamond_validation_db",
        type=str,
        help="Use a Diamond db with blastp mode to help validate cds sequences",
    )
    parser.add_argument(
        "--validation_type",
        type=str,
        help='The strength of evidence needed to validate an ORF as protein coding, can be "relaxed" or "moderate"',
    )
    parser.add_argument(
        "--load_to_ensembl_db",
        # TODO
        # resolve boolean vs string argument type
        action="store_true",
        help="Load results to an Ensembl db, must also provide the db_details flag",
    )
    parser.add_argument(
        "--trim_fastq",
        action="store_true",
        help="Trim the short read files using Trim Galore",
    )
    parser.add_argument(
        "--delete_pre_trim_fastq",
        action="store_true",
        help="Delete the original fastq files after trimming",
    )
    parser.add_argument(
        "--repeatmasker_library", type=str, help="Specify library for repeatmasker"
    )
    parser.add_argument(
        "--repeatmasker_species",
        type=str,
        # TODO
        # does a default value have to been defined here?
        # run_repeatmasker_regions function doesn't assume an existing value for species
        # default="homo",
        help="Specify species for repeatmasker (default homo)",
    )
    args = parser.parse_args()

    work_dir = args.output_dir
    genome_file = args.genome_file
    num_threads = args.num_threads
    masked_genome_file = genome_file  # This will be updated later if Red is run
    run_masking = args.run_masking
    red_path = args.red_path
    genblast_path = args.genblast_path
    convert2blastmask_path = args.convert2blastmask_path
    makeblastdb_path = args.makeblastdb_path
    run_genblast = args.run_genblast
    genblast_timeout = args.genblast_timeout
    run_busco = args.run_busco
    protein_file = args.protein_file
    busco_protein_file = args.busco_protein_file
    rfam_accessions_file = args.rfam_accessions_file
    run_star = args.run_star
    star_path = args.star_path
    short_read_fastq_dir = args.short_read_fastq_dir
    max_intron_length = args.max_intron_length
    max_reads_per_sample = args.max_reads_per_sample
    max_total_reads = args.max_total_reads
    run_minimap2 = args.run_minimap2
    minimap2_path = args.minimap2_path
    paftools_path = args.paftools_path
    long_read_fastq_dir = args.long_read_fastq_dir
    run_augustus = args.run_augustus
    augustus_path = args.augustus_path
    run_stringtie = args.run_stringtie
    run_scallop = args.run_scallop
    stringtie_path = args.stringtie_path
    scallop_path = args.scallop_path
    subsample_script_path = args.subsample_script_path
    samtools_path = args.samtools_path
    finalise_geneset = args.finalise_geneset
    db_details = args.db_details
    run_cmsearch = args.run_cmsearch
    run_trf = args.run_trf
    trf_path = args.trf_path
    run_dust = args.run_dust
    dust_path = args.dust_path
    run_trnascan = args.run_trnascan
    trnascan_path = args.trnascan_path
    trnascan_filter_path = args.trnascan_filter_path
    run_cpg = args.run_cpg
    cpg_path = args.cpg_path
    run_eponine = args.run_eponine
    eponine_path = args.eponine_path
    java_path = args.java_path
    run_repeatmasker = args.run_repeatmasker
    repeatmasker_path = args.repeatmasker_path
    run_full_annotation = args.run_full_annotation
    run_repeats = args.run_repeats
    run_simple_features = args.run_simple_features
    run_sncrnas = args.run_sncrnas
    run_transcriptomic = args.run_transcriptomic
    run_proteins = args.run_proteins
    diamond_validation_db = args.diamond_validation_db
    validation_type = args.validation_type
    load_to_ensembl_db = args.load_to_ensembl_db
    trim_fastq = args.trim_fastq
    delete_pre_trim_fastq = args.delete_pre_trim_fastq
    library = args.repeatmasker_library
    species = args.repeatmasker_species

    main_script_dir = pathlib.Path(os.path.realpath(__file__)).parent

    genome_file = pathlib.Path(genome_file)
    if not genome_file.exists():
        raise FileNotFoundError("File does not exist: %s" % genome_file)

    if not work_dir:
        work_dir = pathlib.Path.cwd()
    work_dir = pathlib.Path(work_dir)

    # save log to a file
    log_file_path = work_dir / "ensembl_anno.log"
    add_log_file_handler(logger, log_file_path)

    logger.info("work directory: %s" % work_dir)

    if not work_dir.exists():
        logger.info("Work dir does not exist, creating:")
        create_dir(work_dir)

    if num_threads == 1:
        logger.warning(
            "Thread count is set to the default value 1; this might be slow."
        )

    # If the run_full_annotation flag is set then we want to set a standardised set of analyses
    if run_full_annotation:
        run_repeats = True
        run_simple_features = True
        run_sncrnas = True
        run_transcriptomic = True
        run_proteins = True
        finalise_geneset = True

    # These are subsets of the analyses that can be run, group by type
    if run_repeats:
        run_masking = True
        run_dust = True
        run_trf = True

    if run_simple_features:
        run_cpg = True
        run_eponine = True

    if run_sncrnas:
        if rfam_accessions_file:
            run_cmsearch = True
        run_trnascan = True

    if run_transcriptomic:
        if short_read_fastq_dir:
            run_star = True
            run_scallop = True
            run_stringtie = True
        if long_read_fastq_dir:
            run_minimap2 = True

    if run_proteins:
        if protein_file:
            run_genblast = True
        if busco_protein_file:
            run_busco = True

    # Collect a list of seq region names, most useful for multiprocessing regions
    seq_region_names = get_seq_region_names(genome_file)
    logger.debug("seq region names:\n%s" % seq_region_names)

    # Repeat analyses
    ############################################################################
    if run_masking:
        logger.info("Running masking via Red")
        masked_genome_file = run_red(red_path, genome_file, work_dir)
        logger.info("Masked genome file: %s" % masked_genome_file)
    else:
        logger.info("Not running masking, presuming the genome file is softmasked")

    if run_dust:
        logger.info("Annotating low complexity regions")
        run_dust_regions(genome_file, dust_path, work_dir, num_threads)

    if run_trf:
        logger.info("Annotating tandem repeats")
        run_trf_repeats(genome_file, trf_path, work_dir, num_threads)

    if run_repeatmasker:
        logger.info("Annotating repeats with RepeatMasker")
        run_repeatmasker_regions(
            genome_file, repeatmasker_path, library, species, work_dir, num_threads
        )
    ############################################################################

    # Simple feature analyses
    ############################################################################
    if run_cpg:
        logger.info("Annotating CpG islands")
        run_cpg_regions(genome_file, cpg_path, work_dir, num_threads)

    if run_eponine:
        logger.info("Running Eponine to find transcription start sites")
        run_eponine_regions(genome_file, java_path, eponine_path, work_dir, num_threads)
    ############################################################################

    # sncRNA analyses
    ############################################################################
    # Search Rfam with cmsearch
    if run_cmsearch:
        logger.info("Annotating sncRNAs")
        run_cmsearch_regions(
             genome_file=genome_file,
             cmsearch_path=None,
             rfam_cm_db_path=None,
             rfam_seeds_file_path=None,
             rfam_accession_file=rfam_accessions_file,
             main_output_dir=work_dir,
             num_threads=num_threads,)

    if run_trnascan:
        logger.info("Annotating tRNAs")
        run_trnascan_regions(
            genome_file, trnascan_path, trnascan_filter_path, work_dir, num_threads
        )
    ############################################################################

    # Transcriptomic analyses
    ############################################################################
    if trim_fastq:
        run_trimming(work_dir, short_read_fastq_dir, delete_pre_trim_fastq, num_threads)

    # Run STAR
    if run_star:
        logger.info("Running STAR")
        run_star_align(
            star_path=star_path,
            trim_fastq=trim_fastq,
            subsample_script_path=subsample_script_path,
            main_output_dir=work_dir,
            short_read_fastq_dir=short_read_fastq_dir,
            genome_file=genome_file,
            max_reads_per_sample=max_reads_per_sample,
            max_total_reads=max_total_reads,
            max_intron_length=max_intron_length,
            num_threads=num_threads,
            main_script_dir=main_script_dir,
        )

    # Run Scallop
    if run_scallop:
        logger.info("Running Scallop")
        run_scallop_assemble(scallop_path, stringtie_path, work_dir)

    # Run Stringtie
    if run_stringtie:
        logger.info("Running StringTie")
        run_stringtie_assemble(stringtie_path, samtools_path, work_dir, num_threads)

    # Run minimap2
    if run_minimap2:
        logger.info("Running minimap2")
        run_minimap2_align(
            minimap2_path,
            paftools_path,
            work_dir,
            long_read_fastq_dir,
            genome_file,
            max_intron_length,
            num_threads,
        )

    if run_transcriptomic:
        check_transcriptomic_output(work_dir)
    ############################################################################

    # Protein analyses
    ############################################################################
    # Run GenBlast
    if run_genblast:
        logger.info("Running GenBlast")
        run_genblast_align(
            genblast_path,
            convert2blastmask_path,
            makeblastdb_path,
            work_dir / "genblast_output",
            protein_file,
            masked_genome_file,
            max_intron_length,
            num_threads,
            genblast_timeout,
            main_script_dir,
        )

    # Run GenBlast on BUSCO set, gives higher priority when creating the final genes in cases where transcriptomic data are missing or fragmented
    if run_busco:
        logger.info("Running GenBlast of BUSCO proteins")
        run_genblast_align(
            genblast_path,
            convert2blastmask_path,
            makeblastdb_path,
            work_dir / "busco_output",
            busco_protein_file,
            masked_genome_file,
            max_intron_length,
            num_threads,
            genblast_timeout,
            main_script_dir,
        )
    ############################################################################

    # Finalisation analyses
    ############################################################################
    # Do some magic
    if finalise_geneset:
        logger.info("Finalise geneset")
        run_finalise_geneset(
            main_script_dir,
            work_dir,
            genome_file,
            seq_region_names,
            validation_type,
            diamond_validation_db,
            num_threads,
        )
    ############################################################################

    # Other analyses
    ############################################################################
    # Run Augustus
    if run_augustus:
        logger.info("Running Augustus")
        run_augustus_predict(augustus_path, work_dir, masked_genome_file, num_threads)
    ############################################################################

    if load_to_ensembl_db:
        load_results_to_ensembl_db(
            main_script_dir,
            load_to_ensembl_db,
            genome_file,
            work_dir,
            db_details,
            num_threads,
        )


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.info("Interrupted with CTRL-C, exiting...")
