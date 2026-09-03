process TRF {
    label 'process_medium'

    publishDir "${params.outdir}/trf",
        mode: 'copy'

    input:
    tuple val(coords), path(sliced_fasta)
 
    output:
    tuple val(coords), path('*.dat'),          emit: trf_repeats
    path "versions.yml",                     emit: versions

    script:
    """
    bash -c 'trf ${sliced_fasta} \
        ${params.match_score} \
        ${params.mismatch_score} \
        ${params.delta} \
        ${params.pm} \
        ${params.pi} \
        ${params.minscore} \
        ${params.maxperiod} \
        -d -h' || echo 'processed $? TRs'

    trf -v >> versions.yml
    """

    stub:
    """
    touch coords.dat
    touch versions.yml
    """
}