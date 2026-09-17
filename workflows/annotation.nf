include { TRANSCRIPTOMICS_ANNOTATION } from '../subworkflows/transcriptomic_annotation.nf'
include { SPLIT_FASTA } from '../subworkflows/split_fasta.nf'
include { REPEATS } from '../subworkflows/repeats.nf'
include { SIMPLE_FEATURE_ANNOTATION } from '../subworkflows/simple_feature_annotation.nf'
include { SMALL_NCRNA_ANNOTATION } from '../subworkflows/small_ncRNA_annotation.nf'
include { PROTEINS } from '../subworkflows/proteins.nf'

workflow ANNOTATION {
    take:
    fasta                   // [fasta]
    short_reads             // [basename, [fastq1, optional fastq2]]
    long_reads              // [basename, [fastq]]
    proteins                // [source, [fasta]]

    fasta_slicing_params    // [param1: val, param2: val]
    transcriptomics_params  // [param1: val, param2: val]
    simple_features_params  // [param1: val, param2: val]
    repeats_params          // [param1: val, param2: val]
    protein_params          // [param1: val, param2: val]

    rfam_accession_file     // [path]
    rfam_cm_db              // [path]
    rfam_seeds_file         // [path]
    genblast_alignscore     // [path]

    main:

    sliced_fastas = SPLIT_FASTA(fasta, fasta_slicing_params)

    REPEATS(fasta, sliced_fastas, repeats_params)
    // SIMPLE_FEATURE_ANNOTATION(sliced_fastas, simple_features_params)
    // SMALL_NCRNA_ANNOTATION(fasta, 
    //                        sliced_fastas, 
    //                        rfam_accession_file,
    //                        rfam_cm_db,
    //                        rfam_seeds_file
    // )

    // TRANSCRIPTOMICS_ANNOTATION(short_reads, long_reads, fasta, transcriptomics_params)
    // PROTEINS(REPEATS.out.red_masked_genome, proteins, genblast_alignscore, protein_params)

}
