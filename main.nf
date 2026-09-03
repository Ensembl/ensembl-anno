#!/usr/bin/env nextflow

include { TRANSCRIPTOMICS_ANNOTATION } from './subworkflows/transcriptomic_annotation.nf'
include { SPLIT_FASTA } from './subworkflows/split_fasta.nf'
include { REPEATS } from './subworkflows/repeats.nf'
include { SIMPLE_FEATURE_ANNOTATION } from './subworkflows/simple_feature_annotation.nf'

nextflow.enable.dsl = 2

nextflow.enable.strict = true

def generate_short_reads_ch(short_read_dir) {
    // The logic that follows enables the pipeline to be able to handle both paired and single ended short reads
    // First strip any trailing / from the path to the directory containing the short read fastqs
    def short_read_dir_str = short_read_dir.replaceAll('/$', '')

    // Next gather together the paired end reads fastqs in a channel - this is relatively straightforwards
    // This generates a channel with the following structure [ sampleName, [sampleName_1.fq, sampleName_2.fq]]
    def short_paired_read_ch = channel.fromFilePairs(
        short_read_dir_str + '/*_{1,2,R1,R2}.{fq,fastq}{.gz,}', 
        size: 2,  // paired reads only
        checkIfExists: true
    )

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
    tuple(file.baseName.split('\\.')[0], [file]) 
    }

    // Mix together the short and long reads
    // Downstream we use the length of the list of fastqs in the channel to determine if we are handling se or pe data
    def short_read_ch = short_paired_read_ch.mix(short_single_read_ch)
    return short_read_ch
}

def generate_long_reads_ch(long_read_dir){
    def long_reads_ch = channel.fromPath(long_read_dir + '/*.fastq')   
    .mix(channel.fromPath(long_read_dir + '/*.fq'))
    .mix(channel.fromPath(long_read_dir + '/*.fastq.gz'))
    .mix(channel.fromPath(long_read_dir + '/*.fq.gz'))
    .map { 
        file ->      // this converts a channel with structure [sampleName.fq] to [sampleName, [sampleName.fq]] 
    tuple(file.baseName.split('\\.')[0], [file]) 
    }
    return long_reads_ch

}


workflow {


    // ToDo work out where all this logic should live

    short_read_ch = channel.empty()
    if (params.short_read_dir != null){
        short_read_ch = generate_short_reads_ch(params.short_read_dir)

    }

    long_read_ch = channel.empty()
    if (params.long_read_dir != null){
        long_read_ch = generate_long_reads_ch(params.long_read_dir)
    }

    fasta_ch = channel.fromPath(params.fasta)

    sliced_fastas = SPLIT_FASTA(fasta_ch)
    sliced_fastas.view()
    REPEATS(fasta_ch, sliced_fastas)
    SIMPLE_FEATURE_ANNOTATION(sliced_fastas)
    // TRANSCRIPTOMICS_ANNOTATION(short_read_ch, long_read_ch, fasta_ch)


    // TODO populate :)
}
