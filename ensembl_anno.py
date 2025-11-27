import os
import pathlib
import logging.config
import argparse
# import genebuild modules
from src.python.ensembl.tools.anno.utils import _utils
from src.python.ensembl.tools.anno.repeat_annotation import red
from src.python.ensembl.tools.anno.repeat_annotation import dust
from src.python.ensembl.tools.anno.repeat_annotation import trf
from src.python.ensembl.tools.anno.repeat_annotation import repeatmasker
from src.python.ensembl.tools.anno.simple_feature_annotation import cpg
from src.python.ensembl.tools.anno.simple_feature_annotation import eponine
from src.python.ensembl.tools.anno.snc_rna_annotation import cmsearch
from src.python.ensembl.tools.anno.snc_rna_annotation import trnascan
from src.python.ensembl.tools.anno.transcriptomic_annotation import augustus
from src.python.ensembl.tools.anno.transcriptomic_annotation import minimap
from src.python.ensembl.tools.anno.transcriptomic_annotation import scallop
from src.python.ensembl.tools.anno.transcriptomic_annotation import star
from src.python.ensembl.tools.anno.transcriptomic_annotation import stringtie
from src.python.ensembl.tools.anno.protein_annotation import genblast
from src.python.ensembl.tools.anno.finalise_genset import finalise_geneset_utils

logger = logging.getLogger(__name__)

def set_logging(work_dir:str) -> None:
    log_file_path = pathlib.Path(work_dir) / "ensembl_anno.log"
    loginipath = pathlib.Path(os.environ["ENSCODE"]) / "ensembl-anno/logging.conf"

    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": log_file_path},
        disable_existing_loggers=False,
    )

def configure_analysis_flags(
    run_full_annotation=False,
    run_repeats=False,
    run_simple_features=False,
    run_sncrnas=False,
    run_transcriptomic=False,
    run_proteins=False,
    run_masking = False,
    run_dust = False,
    run_trf = False,
    run_cpg=False,
    run_eponine=False,
    run_cmsearch=False,
    run_trnascan=False,
    run_star=False,
    run_scallop=False,
    run_stringtie = False,
    run_minimap2=False,
    run_genblast=False,
    run_busco=False,
    rfam_accessions_file=None,
    short_read_fastq_dir=None,
    long_read_fastq_dir=None,
    protein_file=None,
    busco_protein_file=None
):
    """
    Sets up analysis flags based on the high-level options and input files.
    Returns a dictionary with all the analysis flags.
    """
    flags = {}

    # If full annotation is requested, enable all analyses
    if run_full_annotation:
        run_repeats = run_simple_features = run_sncrnas = run_transcriptomic = run_proteins = True
        finalise_geneset = True
    else:
        finalise_geneset = False

    # Repeats
    flags['run_masking'] = run_repeats
    flags['run_masking'] = run_masking
    flags['run_dust'] = run_repeats
    flags['run_dust'] = run_dust
    flags['run_trf'] = run_repeats
    flags['run_trf'] = run_trf

    # Simple features
    flags['run_cpg'] = run_simple_features
    flags['run_cpg'] = run_cpg
    flags['run_eponine'] = run_simple_features
    flags['run_eponine'] = run_eponine

    # sncRNAs
    flags['run_cmsearch'] = run_sncrnas and rfam_accessions_file is not None
    flags['run_cmsearch'] = run_cmsearch and rfam_accessions_file is not None
    flags['run_trnascan'] = run_sncrnas
    flags['run_trnascan'] = run_trnascan

    # Transcriptomic
    flags['run_star'] = run_transcriptomic and short_read_fastq_dir is not None
    flags['run_star'] = run_star and short_read_fastq_dir is not None
    flags['run_scallop'] = run_transcriptomic and short_read_fastq_dir is not None
    flags['run_scallop'] = run_scallop and short_read_fastq_dir is not None
    flags['run_stringtie'] = run_transcriptomic and short_read_fastq_dir is not None
    flags['run_stringtie'] = run_stringtie and short_read_fastq_dir is not None
    flags['run_minimap2'] = run_transcriptomic and long_read_fastq_dir is not None
    flags['run_minimap2'] = run_minimap2 and long_read_fastq_dir is not None

    # Proteins
    flags['run_genblast'] = run_proteins and protein_file is not None
    flags['run_genblast'] = run_genblast and protein_file is not None
    flags['run_busco'] = run_proteins and busco_protein_file is not None
    flags['run_busco'] = run_busco and busco_protein_file is not None

    # Finalisation
    flags['finalise_geneset'] = finalise_geneset

    return flags

def parse_args():
    parser = argparse.ArgumentParser(description="Run Ensembl annotation script")
    parser.add_argument("--output_dir", type=str, required=True, help="Path where output and temp files will be written")
    parser.add_argument("--genome_file", type=str, required=True, help="Path to the fasta genome file")
    parser.add_argument("--num_threads", type=int, default=1, help="Number of threads to use")
    parser.add_argument("--run_masking", action="store_true",help="Run Red to find repeats and softmask the genome. Otherwise provide a softmasked genome")
    parser.add_argument("--red_path", type=str, help="Path to Red executable. See http://toolsmith.ens.utulsa.edu")
    parser.add_argument("--genblast_path", type=str, help="Path to GenBlast executable. See http://genome.sfu.ca/genblast/download.html")
    parser.add_argument("--convert2blastmask_path", type=str, help="Path to convert2blastmask executable")
    parser.add_argument("--makeblastdb_path", type=str, help="Path to makeblastdb executable")
    parser.add_argument("--run_genblast", action="store_true", help="Run GenBlast to align protein sequences")
    parser.add_argument("--genblast_timeout", type=int, help="GenBlast timeout in seconds", default=10800)
    parser.add_argument("--run_busco", action="store_true", help="Run GenBlast to align BUSCO (OrthoDB) protein sequences")
    parser.add_argument("--protein_file", type=str, help="Path to a fasta file with protein sequences")
    parser.add_argument("--busco_protein_file", type=str, help="Path to a fasta file with BUSCO (OrthoDB) protein sequences")
    parser.add_argument("--rfam_accessions_file", type=str, help="Path to a file with Rfam CM accessions, one accession per line, to use with cmsearch")
    parser.add_argument("--run_star", action="store_true", help="Run Star for short read alignment")
    parser.add_argument("--star_path", type=str, help="Path to Star for short read alignment")
    parser.add_argument("--max_reads_per_sample", nargs="?", type=int, default=0, help="The maximum number of reads to use per sample. Default=0 (unlimited)")
    parser.add_argument("--short_read_fastq_dir", help="Path to short read fastq dir for running with Star")
    parser.add_argument("--max_intron_length", nargs="?", type=int, default=100000, help="The maximum intron size for alignments. Default=100000")
    parser.add_argument("--run_minimap2", action="store_true", help="Run minimap2 for long read alignment")
    parser.add_argument("--minimap2_path", type=str, help="Path to minimap2 for long read alignment")
    parser.add_argument("--paftools_path", type=str,  help="Path to paftools for SAM to BED conversion")
    parser.add_argument("--long_read_fastq_dir", type=str, help="Path to long read fastq dir for running with minimap2")
    parser.add_argument("--run_augustus", action="store_true", help="Run Augustus with hints for gene/transcript prediction")
    parser.add_argument("--augustus_path", type=str, help="Path to Augustus")
    parser.add_argument("--run_stringtie", action="store_true",  help="Run Stringtie on the results from the STAR alignments")
    parser.add_argument("--run_scallop", action="store_true", help="Run Scallop on the results from the STAR alignments")
    parser.add_argument("--stringtie_path", type=str, help="Path to Stringtie")
    parser.add_argument("--scallop_path", type=str, help="Path to Scallop")
    parser.add_argument("--subsample_script_path", type=str, help="Path to ensembl-anno subsampling script")
    parser.add_argument("--samtools_path", type=str, help="Path to subsampling script")
    parser.add_argument("--finalise_geneset", action="store_true", help="Used to finalise the gene set from the various GTF files generated")
    parser.add_argument("--db_details", type=str, help="A comma separated string of dbname,host,port,user,pass")
    parser.add_argument("--run_cmsearch", action="store_true", help="Search for sncRNA structures using Rfam and cmsearch")
    parser.add_argument("--run_trf", action="store_true", help="Run TRF to find tandem repeats")
    parser.add_argument("--trf_path", type=str, help="Path to TRF")
    parser.add_argument("--run_dust", action="store_true", help="Run Dust to find low complexity regions")
    parser.add_argument("--dust_path", type=str, help="Path to Dust")
    parser.add_argument("--run_repeatmasker", action="store_true", help="Run RepeatMasker to find repeat regions")
    parser.add_argument("--repeatmasker_path", action="store_true", help="Path to RepeatMasker")
    parser.add_argument("--run_trnascan", action="store_true", help="Run tRNAscan-SE to find tRNAs")
    parser.add_argument("--trnascan_path", type=str, help="Path to tRNAscan-SE")
    parser.add_argument("--trnascan_filter_path", type=str, help="Path to tRNAscan-SE high confidence filter")
    parser.add_argument("--run_cpg", action="store_true", help="Run cpg_lh to find CpG islands")
    parser.add_argument("--cpg_path", type=str, help="Path to cpg_lh")
    parser.add_argument("--run_eponine", action="store_true", help="Run Eponine to find transcription start sites")
    parser.add_argument("--eponine_path", type=str, help="Path to Eponine jar file")
    parser.add_argument("--java_path", type=str, help="Path to Java for use with Eponine")
    parser.add_argument("--run_full_annotation",action="store_true",help="Run a full annotation, will automatically check for input data and run tools based on that")
    parser.add_argument("--run_repeats", action="store_true", help="Run Red, Dust, TRF")
    parser.add_argument("--run_simple_features", action="store_true", help="Run CpG, Eponine")
    parser.add_argument("--run_sncrnas", action="store_true", help="Run Rfam, tRNAscan-SE")
    parser.add_argument("--run_transcriptomic", action="store_true", help="Run STAR, Stringtie2, Scallop, minimap2 (if short_read_fastq_dir and/or long_read_fastq_dir are provided)")
    parser.add_argument("--run_proteins", action="store_true", help="Run GenBlast if protein_file and/or busco_protein_file")
    parser.add_argument("--diamond_validation_db", type=str, help="Use a Diamond db with blastp mode to help validate cds sequences")
    parser.add_argument("--validation_type", type=str, help='The strength of evidence needed to validate and ORF as protein coding, can be "relaxed" or "moderate"')
    parser.add_argument("--load_to_ensembl_db", action="store_true", help="Load results to an Ensembl db, must also provide the db_details flag")
    parser.add_argument("--trim_fastq", action="store_true", help="Trim the short read files using Trim Galore")
    parser.add_argument("--delete_pre_trim_fastq", action="store_true", help="Delete the original fastq files after trimming")
    parser.add_argument("--repeatmasker_library", type=str, help="Specify library for repeatmasker")
    parser.add_argument("--repeatmasker_species", type=str, help="Specify species for repeatmasker (default homo sapiens)")
    parser.add_argument("--repeatmasker_analysis", type=str, default="repeatmask_repbase_human",  help="Specify logic name for repeatmasker analysis (default repeatmask_repbase_human)")
    parser.add_argument("--trim_galore_path", type=str, help="Path to trim_galore")

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    work_dir = args.output_dir or os.getcwd()
    genome_file = args.genome_file
    num_threads = args.num_threads
    # masked_genome_file = genome_file  # This will be updated later if Red is run
    red_path = args.red_path
    genblast_path = args.genblast_path
    convert2blastmask_path = args.convert2blastmask_path
    makeblastdb_path = args.makeblastdb_path
    genblast_timeout = args.genblast_timeout
    protein_file = args.protein_file
    busco_protein_file = args.busco_protein_file
    rfam_accessions_file = args.rfam_accessions_file
    star_path = args.star_path
    short_read_fastq_dir = args.short_read_fastq_dir
    max_intron_length = args.max_intron_length
    max_reads_per_sample = args.max_reads_per_sample
    minimap2_path = args.minimap2_path
    paftools_path = args.paftools_path
    long_read_fastq_dir = args.long_read_fastq_dir
    run_augustus = args.run_augustus
    augustus_path = args.augustus_path
    stringtie_path = args.stringtie_path
    scallop_path = args.scallop_path
    subsample_script_path = args.subsample_script_path
    samtools_path = args.samtools_path
    finalise_geneset = args.finalise_geneset
    db_details = args.db_details
    trf_path = args.trf_path
    dust_path = args.dust_path
    trnascan_path = args.trnascan_path
    trnascan_filter_path = args.trnascan_filter_path
    cpg_path = args.cpg_path
    eponine_path = args.eponine_path
    java_path = args.java_path
    repeatmasker_path = args.repeatmasker_path
    run_transcriptomic = args.run_transcriptomic
    diamond_validation_db = args.diamond_validation_db
    validation_type = args.validation_type
    load_to_ensembl_db = args.load_to_ensembl_db
    trim_fastq = args.trim_fastq
    delete_pre_trim_fastq = args.delete_pre_trim_fastq
    repeatmasker_library = args.repeatmasker_library
    species = args.repeatmasker_species
    repeatmasker_analysis = args.repeatmasker_analysis
    trim_galore_path = args.trim_galore_path


    # Initialize logger
    set_logging(work_dir)
    logger.info("Logger initialized")

    logger.info("Working directory is set as: %s" % work_dir)

    # Validate paths and files
    if not os.path.exists(work_dir):
        logger.info("Work dir does not exist, will create")
        _utils.create_dir(work_dir, None)

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
        logger.info('Masking genome file not found. Defaulting to genome file')
    logger.info("Masked genome file %s", masked_genome_file)

    if not os.path.exists(genome_file):
        raise FileNotFoundError(f"Genome file does not exist: {genome_file}")

    if num_threads == 1:
        logger.info("Thread count is set to the default value 1; this might be slow.")


    # Get run settings based on arguments
    # Configure which analyses to run
    analysis_flags = configure_analysis_flags(
        run_full_annotation=args.run_full_annotation,
        run_repeats=args.run_repeats,
        run_simple_features=args.run_simple_features,
        run_sncrnas=args.run_sncrnas,
        run_transcriptomic=args.run_transcriptomic,
        run_proteins=args.run_proteins,
        run_masking=args.run_masking,
        run_dust=args.run_dust,
        run_trf=args.run_trf,
	    run_cpg=args.run_cpg,
	    run_eponine=args.run_eponine,
	    run_cmsearch=args.run_cmsearch,
	    run_trnascan=args.run_trnascan,
	    run_star=args.run_star,
	    run_scallop=args.run_scallop,
	    run_stringtie=args.run_stringtie,
	    run_minimap2=args.run_minimap2,
	    run_genblast=args.run_genblast,
	    run_busco=args.run_busco,
        rfam_accessions_file=args.rfam_accessions_file,
        short_read_fastq_dir=args.short_read_fastq_dir,
        long_read_fastq_dir=args.long_read_fastq_dir,
        protein_file=args.protein_file,
        busco_protein_file=args.busco_protein_file,
    )

    # Get seq_region_names
    seq_region_names = _utils.seq_region_names(genome_file)
    logger.info("Got seq_region_names")

    #################################
    # Repeat analyses
    #################################
    if analysis_flags['run_masking']:
        logger.info("Running masking via Red")
        masked_genome_file = red.run_red(red_bin = red_path,
                                         output_dir = work_dir,
                                         genome_file = genome_file)
        logger.info("Masked genome file: " + masked_genome_file)

    else:
        logger.info("Not running masking, presuming the genome file is softmasked")

    if analysis_flags['run_dust']:
        logger.info("Annotating low complexity regions")
        logger.info("run_dust genome file %s", genome_file)
        dust.run_dust(
            genome_file = genome_file,
	        dust_bin = dust_path,
	        output_dir = work_dir,
	        num_threads = num_threads
        )

    if analysis_flags['run_trf']:
        logger.info("Annotating tandem repeats")
        logger.info("run_trf genome file %s", genome_file)
        trf.run_trf(genome_file = genome_file,
                    trf_bin = trf_path,
                    output_dir = work_dir,
                    num_threads = num_threads)

    if analysis_flags['run_repeatmasker']:
        logger.info("Annotating repeats with RepeatMasker")
        logger.info("run_repeatmasker genome file %s", genome_file)
        repeatmasker.run_repeatmasker(
            genome_file = genome_file,
            repeatmasker_bin = repeatmasker_path,
            library = repeatmasker_library,
            species = species,
            output_dir = work_dir,
            num_threads = num_threads,
        )
    if not analysis_flags['run_repeatmasker'] or analysis_flags['run_trf'] or analysis_flags['run_dust']:
        logger.info("Repeat analyses were not requested")

    #################################
    # Simple feature analyses
    #################################
    logger.info("Checking simple features ")
    if analysis_flags['run_cpg']:
        logger.info("Annotating CpG islands")
        logger.info("run_cpg genome file %s", genome_file)
        cpg.run_cpg(genome_file = genome_file,
                    cpg_bin = cpg_path,
                    output_dir = work_dir,
                    num_threads = num_threads)

    if analysis_flags['run_eponine']:
        logger.info("Running Eponine to find transcription start sites")
        logger.info("run_eponine genome file %s", genome_file)
        eponine.run_eponine(
            genome_file = genome_file,
	        java_bin = java_path,
	        eponine_bin = eponine_path,
	        output_dir = work_dir,
	        num_threads = num_threads
        )
    if not (analysis_flags['run_cpg'] or analysis_flags['run_eponine']):
        logger.info("Simple features were not requested")

    #################################
    # sncRNA analyses
    #################################
    if analysis_flags['run_cmsearch']:
        logger.info("Annotating sncRNAs")
        logger.info("run_cmsearch genome file %s", genome_file)
        cmsearch.run_cmsearch(
            genome_file = genome_file,
            output_dir = work_dir,
            rfam_accession_file = rfam_accessions_file,
            num_threads = num_threads,
        )

    if analysis_flags['run_trnascan']:
        logger.info("Annotating tRNAs")
        logger.info("run_trnascan genome file %s", genome_file)
        trnascan.run_trnascan(
            genome_file = genome_file,
            output_dir = work_dir,
            num_threads = num_threads,
	        trnascan_bin = trnascan_path,
	        trnascan_filter = trnascan_filter_path

        )
    if not (analysis_flags['run_cmsearch'] or analysis_flags['run_trnascan']):
        logger.info("sncRNA analyses were not requested")

    #################################
    # Transcriptomic analyses
    #################################
    if args.trim_fastq:
        star.run_trimming(work_dir, short_read_fastq_dir, delete_pre_trim_fastq, num_threads)

    if analysis_flags['run_star']:
        logger.info("Running Star")
        logger.info("run_star genome file %s", genome_file)
        star.run_star(
            star_bin = star_path,
            trim_fastq = trim_fastq,
	        delete_pre_trim_fastq = delete_pre_trim_fastq,
            output_dir = work_dir,
            short_read_fastq_dir = short_read_fastq_dir,
            genome_file = genome_file,
            max_reads_per_sample = max_reads_per_sample,
            max_intron_length = max_intron_length,
            num_threads = num_threads,
	        samtools_bin = samtools_path,
	        trim_galore_bin = trim_galore_path,
        )

    if analysis_flags['run_scallop']:
        logger.info("Running Scallop")
        scallop.run_scallop(scallop_bin = scallop_path,
                            stringtie_bin = stringtie_path,
                            output_dir = work_dir)

    if analysis_flags['run_stringtie']:
        logger.info("Running Stringtie")
        logger.info("run_stringtie genome file %s", genome_file)
        stringtie.run_stringtie(
            stringtie_bin = stringtie_path,
	        output_dir = work_dir,
	        num_threads = num_threads
        )

    if analysis_flags['run_minimap2']:
        logger.info("Running minimap2")
        logger.info("run_minimap2 genome file %s", genome_file)
        minimap.run_minimap2(
            minimap2_bin = minimap2_path,
            paftools_bin = paftools_path,
            output_dir = work_dir,
            long_read_fastq_dir = long_read_fastq_dir,
            genome_file = genome_file,
            max_intron_length = max_intron_length,
            num_threads = num_threads,
        )

    if run_transcriptomic:
        _utils.check_transcriptomic_output(work_dir)

    #################################
    # Protein analyses
    #################################
    if analysis_flags['run_genblast']:
        logger.info("Running GenBlast")
        logger.info("run_genblast genome file %s", masked_genome_file)
        genblast.run_genblast(
	        genblast_bin = genblast_path,
	        convert2blastmask_bin = convert2blastmask_path,
	        makeblastdb_bin = makeblastdb_path,
	        output_dir = work_dir,
	        protein_dataset = protein_file,
	        masked_genome = masked_genome_file,
	        max_intron_length = max_intron_length,
	        num_threads = num_threads,
	        genblast_timeout_secs = genblast_timeout,
	        protein_set = "uniprot",
        )

    # Run GenBlast on OrthoDB or other set, gives higher priority when creating the
    # final genes in cases where transcriptomic data are missing or fragmented
    if analysis_flags['run_busco']:
        logger.info("Running GenBlast of OrthoDB or FungiDB proteins")
        logger.info("run_busco genome file %s", masked_genome_file)
        genblast.run_genblast(
	        genblast_bin = genblast_path,
	        convert2blastmask_bin = convert2blastmask_path,
	        makeblastdb_bin = makeblastdb_path,
	        output_dir = work_dir,
	        protein_dataset = busco_protein_file,
	        masked_genome = masked_genome_file,
	        max_intron_length = max_intron_length,
	        num_threads = num_threads,
	        genblast_timeout_secs = genblast_timeout,
	        protein_set="orthodb"
        )

    #################################
    # Finalisation not yet modularised
    #################################
    # Do some magic
    if analysis_flags['finalise_geneset']:
        logger.info("Finalise geneset")
        logger.info("Finalise geneset genome file %s", genome_file)
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
	        repeatmasker_analysis,
        )












if __name__ == "__main__":
    main()