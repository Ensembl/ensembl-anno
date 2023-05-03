# pylint: disable=logging-not-lazy, invalid-name, missing-function-docstring, subprocess-run-check, unused-variable, redefined-outer-name, too-many-arguments, too-many-locals, too-many-branches, too-many-statements, unused-argument, no-else-return, undefined-variable, no-else-continue, no-else-raise, missing-docstring, consider-swap-variables, consider-using-in, too-many-lines, unused-import
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
import sys
import argparse
import errno
import gc
import glob
import io
import json
import logging
import logging.config
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

from repeats.repeatmasker import run_repeatmasker_regions
import simple_feature_utils
import utils

with open(os.environ["ENSCODE"] + "/ensembl-anno/config.json", "r") as f:
    config = json.load(f)






def prlimit_command(command_list, virtual_memory_limit):

    """
    Prepend memory limiting arguments to a command list to be run with subprocess.

    This method uses the `prlimit` program to set the memory limit.

    The `virtual_memory_limit` size is in bytes.

    prlimit arguments:
    -v, --as[=limits]
           Address space limit.
    """
    return ["prlimit", f"-v{virtual_memory_limit}"] + command_list






if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output_dir",
        type=str,
        help="Path where the output and temp files will write to. \
        Uses current dir by default",
    )
    parser.add_argument(
        "--genome_file",
        type=str,
        help="Path to the fasta genome file",
        required=True,
    )
    parser.add_argument(
        "--num_threads", type=int, default=1, help="Number of threads to use"
    )
    parser.add_argument(
        "--run_masking",
        action="store_true",
        help="Run Red to find repeats and softmask the genome. \
        Otherwise provide a softmasked genome",
    )
    parser.add_argument(
        "--red_path",
        type=str,
        help="Path to Red executable. See http://toolsmith.ens.utulsa.edu",
    )
    parser.add_argument(
        "--genblast_path",
        type=str,
        help="Path to GenBlast executable. See http://genome.sfu.ca/\
        genblast/download.html",
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
        help="GenBlast timeout in seconds",
        default=10800,
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
        help="Path to a file with Rfam CM accessions, one accession\
        per line, to use with cmsearch",
    )
    parser.add_argument(
        "--run_star",
        action="store_true",
        help="Run Star for short read alignment",
    )
    parser.add_argument(
        "--star_path", type=str, help="Path to Star for short read alignment"
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
        help="Path to short read fastq dir for running with Star",
    )
    parser.add_argument(
        "--max_intron_length",
        nargs="?",
        type=int,
        default=100000,
        help="The maximum intron size for alignments. Default=100000",
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
        "--run_trnascan",
        action="store_true",
        help="Run tRNAscan-SE to find tRNAs",
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
    parser.add_argument("--java_path", type=str, help="Path to Java for use with Eponine")
    parser.add_argument(
        "--run_full_annotation",
        action="store_true",
        help="Run a full annotation, will automatically check \
        for input data and run tools based on that",
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
        help="Run STAR, Stringtie2, Scallop, minimap2 (if short_read_fastq_dir \
        and/or long_read_fastq_dir are provided)",
    )
    parser.add_argument(
        "--run_proteins",
        action="store_true",
        help="Run GenBlast if protein_file and/or busco_protein_file",
    )
    parser.add_argument(
        "--diamond_validation_db",
        type=str,
        help="Use a Diamond db with blastp mode to help validate cds sequences",
    )
    parser.add_argument(
        "--validation_type",
        type=str,
        help='The strength of evidence needed to validate and ORF as \
        protein coding, can be "relaxed" or "moderate"',
    )
    parser.add_argument(
        "--load_to_ensembl_db",
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
        "--repeatmasker_library",
        type=str,
        help="Specify library for repeatmasker",
    )
    parser.add_argument(
        "--repeatmasker_species",
        type=str,
        help="Specify species for repeatmasker (default homo)",
    )
    args = parser.parse_args()

    work_dir = args.output_dir
    genome_file = args.genome_file
    num_threads = args.num_threads
    # masked_genome_file = genome_file  # This will be updated later if Red is run
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

    main_script_dir = os.path.dirname(os.path.realpath(__file__))
    # work_dir=glob.glob(work_dir)
    if not os.path.exists(genome_file):
        raise IOError("File does not exist: %s" % genome_file)

    if not work_dir:
        work_dir = os.getcwd()
        # work_dir=glob.glob(work_dir)

    # set up logger
    log_file_path = pathlib.Path(work_dir) / "ensembl_anno.log"
    loginipath = pathlib.Path(os.environ["ENSCODE"] + "/ensembl-anno/logging.conf")
    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": log_file_path},
        disable_existing_loggers=False,
    )
    logger = logging.getLogger()
    logger.propagate = False

    logger.info("work directory: %s" % work_dir)
    if not os.path.exists(work_dir):
        logger.info("Work dir does not exist, will create")
        utils.create_dir(work_dir, None)

    if num_threads == 1:
        logger.info("Thread count is set to the default value 1; this might be slow.")

    if os.path.exists(
        os.path.join(work_dir, "red_output", "mask_output")
    ) or os.path.join(work_dir, "red_output", "mask_output").endswith(".msk"):
        red_genome_file = [
            f
            for f in os.listdir(os.path.join(work_dir, "red_output", "mask_output"))
            if f.endswith(".msk")
        ]
        logger.info("red_genome_file %s", red_genome_file)
        masked_genome_file = os.path.join(
            work_dir, "red_output", "mask_output", red_genome_file[0]
        )
    else:
        masked_genome_file = genome_file
    logger.info("Masked genome file %s", masked_genome_file)
    # If the run_full_annotation flag is set then we want to set
    # a standardised set of analyses
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
    seq_region_names = seq_region_names(genome_file)
    for i in seq_region_names:
        logger.info(i)

    #################################
    # Repeat analyses
    #################################
    if run_masking:
        logger.info("Running masking via Red")
        masked_genome_file = repeatmasking_utils.run_red(red_path, work_dir, genome_file)
        logger.info("Masked genome file: " + masked_genome_file)

    else:
        logger.info("Not running masking, presuming the genome file is softmasked")

    if run_dust:
        logger.info("Annotating low complexity regions")
        logger.info("run_dust genome file %s", genome_file)
        repeatmasking_utils.run_dust_regions(
            genome_file, dust_path, work_dir, num_threads
        )

    if run_trf:
        logger.info("Annotating tandem repeats")
        logger.info("run_trf genome file %s", genome_file)
        repeatmasking_utils.run_trf_repeats(genome_file, trf_path, work_dir, num_threads)

    if run_repeatmasker:
        logger.info("Annotating repeats with RepeatMasker")
        logger.info("run_repeatmasker genome file %s", genome_file)
        repeatmasker.run_repeatmasker_regions(
            genome_file,
            work_dir,
            config["repeatmasker"]["software"],
            library,
            config["repeatmasker"]["engine"],
            species,
            num_threads,
        )

    #################################
    # Simple feature analyses
    #################################
    if run_cpg:
        logger.info("Annotating CpG islands")
        logger.info("run_cpg genome file %s", genome_file)
        simple_feature_utils.run_cpg_regions(genome_file, cpg_path, work_dir, num_threads)

    if run_eponine:
        logger.info("Running Eponine to find transcription start sites")
        logger.info("run_eponine genome file %s", genome_file)
        simple_feature_utils.run_eponine_regions(
            genome_file, java_path, eponine_path, work_dir, num_threads
        )

    #################################
    # sncRNA analyses
    #################################
    # Search Rfam with cmsearch
    if run_cmsearch:
        logger.info("Annotating sncRNAs")
        logger.info("run_cmsearch genome file %s", genome_file)
        run_cmsearch_regions(
            genome_file,
            None,
            None,
            None,
            rfam_accessions_file,
            work_dir,
            num_threads,
        )

    if run_trnascan:
        logger.info("Annotating tRNAs")
        logger.info("run_trnascan genome file %s", genome_file)
        run_trnascan_regions(
            genome_file,
            trnascan_path,
            trnascan_filter_path,
            work_dir,
            num_threads,
        )

    #################################
    # Transcriptomic analyses
    #################################
    if trim_fastq:
        run_trimming(work_dir, short_read_fastq_dir, delete_pre_trim_fastq, num_threads)

    # Run STAR
    if run_star:
        logger.info("Running Star")
        logger.info("run_star genome file %s", genome_file)
        run_star_align(
            star_path,
            trim_fastq,
            subsample_script_path,
            work_dir,
            short_read_fastq_dir,
            genome_file,
            max_reads_per_sample,
            max_total_reads,
            max_intron_length,
            num_threads,
        )

    # Run Scallop
    if run_scallop:
        logger.info("Running Scallop")
        run_scallop_assemble(scallop_path, stringtie_path, work_dir)

    # Run Stringtie
    if run_stringtie:
        logger.info("Running Stringtie")
        logger.info("run_stringtie genome file %s", genome_file)
        run_stringtie_assemble(
            stringtie_path, samtools_path, work_dir, genome_file, num_threads
        )

    # Run minimap2
    if run_minimap2:
        logger.info("Running minimap2")
        logger.info("run_minimap2 genome file %s", genome_file)
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

    #################################
    # Protein analyses
    #################################
    # Run GenBlast
    if run_genblast:
        logger.info("Running GenBlast")
        logger.info("run_genblast genome file %s", masked_genome_file)
        run_genblast_align(
            genblast_path,
            convert2blastmask_path,
            makeblastdb_path,
            os.path.join(work_dir, "genblast_output"),
            protein_file,
            masked_genome_file,
            max_intron_length,
            num_threads,
            genblast_timeout,
        )

    # Run GenBlast on BUSCO set, gives higher priority when creating thei
    # final genes in cases where transcriptomic data are missing or fragmented
    if run_busco:
        logger.info("Running GenBlast of BUSCO proteins")
        logger.info("run_busco genome file %s", masked_genome_file)
        run_genblast_align(
            genblast_path,
            convert2blastmask_path,
            makeblastdb_path,
            os.path.join(work_dir, "busco_output"),
            busco_protein_file,
            masked_genome_file,
            max_intron_length,
            num_threads,
            genblast_timeout,
        )

    #################################
    # Finalisation analyses
    #################################
    # Do some magic
    if finalise_geneset:
        logger.info("Finalise geneset")
        logger.info("finalise geneset genome file %s", genome_file)
        run_finalise_geneset(
            main_script_dir,
            work_dir,
            genome_file,
            seq_region_names,
            validation_type,
            diamond_validation_db,
            num_threads,
        )

    #################################
    # Other analyses
    #################################
    # Run Augustus
    if run_augustus:
        logger.info("Running Augustus")
        run_augustus_predict(augustus_path, work_dir, masked_genome_file, num_threads)

    if load_to_ensembl_db:
        load_results_to_ensembl_db(
            main_script_dir,
            load_to_ensembl_db,
            genome_file,
            work_dir,
            db_details,
            num_threads,
        )

#  coallate_results(work_dir)

#  run_find_orfs(genome_file,work_dir)
