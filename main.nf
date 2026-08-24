#!/usr/bin/env nextflow

include { TRANSCRIPTOMICS_ANNOTATION } from './subworkflows/transcriptomic_annotation.nf'

nextflow.enable.dsl = 2

nextflow.enable.strict = true

/*
* Reverse the sequences - delete this just here for dirty testing during dev work
*/
process reverse {

    input:
    tuple val(y), path(x)

    output:
    path 'helloworld.txt'

    script:
    """
    echo ${y} > helloworld.txt
    """
}
workflow {


    // ToDo work out where all this logic should live


    // The logic that follows enables the pipeline to be able to handle both paired and single ended short reads
    // First strip any trailing / from the path to the directory containing the short read fastqs
    short_read_dir_str = params.short_read_dir.replaceAll('/$', '')

    // Next gather together the paired end reads fastqs in a channel - this is relatively straightforwards
    // This generates a channel with the following structure [ sampleName, [sampleName_1.fq, sampleName_2.fq]]
    short_paired_read_ch = channel.fromFilePairs(
        short_read_dir_str + '/*_{1,2,R1,R2}.{fq,fastq}{.gz,}', 
        size: 2,  // paired reads only
        checkIfExists: true
    )

    // Now make a channel containing the single end read fastqs. This is a bit more complex.
    // We end up with a channel with the following structure [ sampleName, [sampleName.fq]]
    // First just grab all of the fastqs (both single and paired end)
    short_single_read_ch = channel.fromPath(short_read_dir_str + '/*.fastq')
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
    short_read_ch = short_paired_read_ch.mix(short_single_read_ch)
    //short_read_ch.view()

    // ToDo update these - just take random dir for stub testing currently
    fasta_ch = channel.fromPath(params.fasta)
    long_reads_ch = channel.fromPath(params.long_reads + '/*.fastq').map { 
        file ->      // this converts a channel with structure [sampleName.fq] to [sampleName, [sampleName.fq]] 
    tuple(file.baseName.split('\\.')[0], [file]) 
    }
    long_reads_ch.view()

    TRANSCRIPTOMICS_ANNOTATION(short_read_ch, long_reads_ch, fasta_ch)


    // TODO populate :)
}
