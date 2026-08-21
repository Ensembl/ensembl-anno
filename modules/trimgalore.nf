process TRIMGALORE {

    publishDir "${params.outdir}/trim_galore",
        mode: 'copy'

    input:
    tuple val(meta), path(input_fastq)

    output:
    tuple val(meta), path("${meta}_trimmed.fq.gz"),    emit: trimmed_reads
    path "versions.yml",                               emit: versions

    script:
    paired_args = ''
    if (input_fastq.size == 2){
        paired_args = ' --paired'
    }
    """
    # TODO test this - trim_galore is funny with filenames
    # TODO decide what to do about args
    # TODO investigate gzip compression behaviour
    trim_galore --illumina --quality 20 --length 50 ${paired_args} ${input_fastq}
    trim_galore --version > versions.yml
    """

    stub:
    """
    touch ${meta}_trimmed.fq.gz
    touch versions.yml
    """
}