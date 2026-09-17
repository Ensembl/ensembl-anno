#!/usr/bin/env nextflow

include { ANNOTATION } from './workflows/annotation.nf'

nextflow.enable.dsl = 2

nextflow.enable.strict = true

def generate_short_reads_ch(short_read_dir) {
    // This function generates a channel of short read fastqs
    // The channel has the structure [basename, [fastq1, optional_fastq2]]

    // The logic that follows enables the pipeline to be able to handle both paired and single ended short reads
    // First strip any trailing / from the path to the directory containing the short read fastqs
    def short_read_dir_str = short_read_dir.replaceAll('/$', '')

    // Next gather together the paired end reads fastqs in a channel - this is relatively straightforwards
    // This generates a channel with the following structure [ sampleName, [sampleName_1.fq, sampleName_2.fq]]
    def short_paired_read_ch = channel.fromFilePairs(
        short_read_dir_str + '/*_{1,2,R1,R2}.{fq,fastq}{.gz,}', 
        size: 2,  // paired reads only
        checkIfExists: true
    ).map{
        it -> [[id: it[0]], it[1]]
    }

    // Now make a channel containing the single end read fastqs. This is a bit more complex.
    // We end up with a channel with the following structure [ sampleName, [sampleName.fq]]
    // First just grab all of the fastqs (both single and paired end)
    def short_single_read_ch = channel.fromPath(short_read_dir_str + '/*.fastq')
    .mix(channel.fromPath(short_read_dir_str + '/*.fq'))
    .mix(channel.fromPath(short_read_dir_str + '/*.fastq.gz'))
    .mix(channel.fromPath(short_read_dir_str + '/*.fq.gz'))
    .filter { file ->                                   // use the filtering step to remove paired end files here
            !file.baseName.split('\\.')[0].endsWith('_R1') &&
            !file.baseName.split('\\.')[0].endsWith('_R2') &&
            !file.baseName.split('\\.')[0].endsWith('_1') &&
            !file.baseName.split('\\.')[0].endsWith('_2')
    }
    .map { file ->      // this converts a channel with structure [sampleName.fq] to [sampleName, [sampleName.fq]] 
    tuple([id: file.baseName.split('\\.')[0]], [file]) 
    }

    // Mix together the single and paired end reads
    // Downstream we use the length of the list of fastqs in the channel to determine if we are handling se or pe data
    def short_read_ch = short_paired_read_ch.mix(short_single_read_ch)
    return short_read_ch
}

def generate_long_reads_ch(long_read_dir){
    // This function generates a channel of long read fastqs
    // The channel has the structure [basename, [fastq]]

    def long_reads_ch = channel.fromPath(long_read_dir + '/*.fastq')   
    .mix(channel.fromPath(long_read_dir + '/*.fq'))
    .mix(channel.fromPath(long_read_dir + '/*.fastq.gz'))
    .mix(channel.fromPath(long_read_dir + '/*.fq.gz'))
    .map { 
        file ->      // this converts a channel with structure [sampleName.fq] to [sampleName, [sampleName.fq]] 
    tuple([id: file.baseName.split('\\.')[0]], [file]) 
    }
    return long_reads_ch

}


workflow {

    // Set up reads channels for the transcriptomics pipeline
    // First short reads:
    short_read_ch = channel.empty()
    if (params.short_read_dir != null){
        short_read_ch = generate_short_reads_ch(params.short_read_dir)

    }

    // Then long reads:
    long_read_ch = channel.empty()
    if (params.long_read_dir != null){
        long_read_ch = generate_long_reads_ch(params.long_read_dir)
    }

    // Initialise a fasta channel (containing the unsliced fasta)
    fasta_ch = channel.fromPath(params.fasta)

    // Set up a channel containing protein fastas
    orthodb_ch = channel.fromPath(params.orthodb).map{
        it -> tuple('orthodb', it)
    }
    uniprot_ch = channel.fromPath(params.uniprot).map{
        it -> tuple('uniprot', it)
    }
    protein_ch = orthodb_ch.concat(uniprot_ch)



    // Now initialise param channels:

    fasta_slicing_param_ch = channel.value([
        slice_size: params.slice_size,
        min_seq_length: params.min_seq_length
    ])

    transcriptomics_param_ch = channel.value([
        max_intron_length: params.max_intron_length,
        min_total_transcriptomic_gtf_lines: params.min_total_transcriptomic_gtf_lines
    ])

    simple_features_param_ch = channel.value([
        eponine_bin: params.eponine_bin,
        eponine_threshold: params.eponine_threshold
    ])

    repeats_param_ch = channel.value([
        library: params.library,
        species: params.species,
        match_score: params.match_score,
        mismatch_score: params.mismatch_score,
        delta: params.delta,
        pm: params.pm,
        pi: params.pi,
        minscore: params.minscore,
        maxperiod: params.maxperiod
    ])


    // rfam files
    rfam_accession_file_ch = channel.fromPath(params.rfam_accession_file).collect()
    rfam_cm_db_ch = channel.fromPath(params.rfam_cm_db).collect()
    rfam_seeds_file_ch = channel.fromPath(params.rfam_seeds_file).collect()
    
    // This file is required to run genblast
    genblast_alignscore = channel.fromPath(params.genblast_alignscore).collect()


    protein_param_ch = channel.value([
        max_intron_length: params.max_intron_length
    ])

    short_read_ch.view()
    long_read_ch.view()

    // Run annotation
    ANNOTATION( fasta_ch,
                short_read_ch,
                long_read_ch,
                protein_ch,

                fasta_slicing_param_ch,
                transcriptomics_param_ch,
                simple_features_param_ch,
                repeats_param_ch,
                protein_param_ch,

                rfam_accession_file_ch,
                rfam_cm_db_ch,
                rfam_seeds_file_ch,
                genblast_alignscore)
}
