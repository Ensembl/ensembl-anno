process TRF {
    label 'process_medium'

    publishDir "${params.outdir}/trf",
        mode: 'copy'

    input:
    tuple val(coords), path(sliced_fasta)
 
    output:
    tuple val(coords), path('*.trf.gtf'),          emit: trf_repeats
    path "versions.yml",                     emit: versions

    script:
    """
    trf ${sliced_fasta} \
        ${params.match_score} \
        ${params.mismatch_score} \
        ${params.delta} \
        ${params.pm} \
        ${params.pi} \
        ${params.minscore} \
        ${params.maxperiod} \
        -d -h

    trf -v >> versions.yml
    """

    stub:
    """
    touch coords.trf.gtf
    touch versions.yml
    """
}