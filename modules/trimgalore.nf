process TRIMGALORE {
    label 'process_low'

    publishDir "${params.outdir}/trim_galore",
        mode: 'copy'

    input:
    tuple val(meta), path(input_fastq)

    output:
    tuple val(meta), path("${meta.id}_trimmed.fq.gz"),    emit: trimmed_reads
    path "versions.yml",                               emit: versions

    script:
    def fastqs = input_fastq instanceof List ? input_fastq : [input_fastq]
    def paired_args = fastqs.size() == 2 ? ' --paired' : ''
    """
    # TODO test this - trim_galore is funny with filenames
    # TODO decide what to do about args
    # TODO investigate gzip compression behaviour
    trim_galore --illumina --quality 20 --length 50 ${paired_args} ${input_fastq}
    trim_galore --version > versions.yml
    """

    stub:
    """
    touch ${meta.id}_trimmed.fq.gz
    touch versions.yml
    """
}