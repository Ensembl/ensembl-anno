process TRF {
    label 'process_medium'

    publishDir "${params.outdir}/trf",
        mode: 'copy'

    input:
    tuple val(coords), path(sliced_fasta)
    val(repeats_params)
 
    output:
    tuple val(coords), path('*.dat'),          emit: trf_repeats
    path "versions.yml",                     emit: versions

    script:
    """
    bash -c 'trf ${sliced_fasta} \
        ${repeats_params.match_score} \
        ${repeats_params.mismatch_score} \
        ${repeats_params.delta} \
        ${repeats_params.pm} \
        ${repeats_params.pi} \
        ${repeats_params.minscore} \
        ${repeats_params.maxperiod} \
        -d -h' || echo 'processed $? TRs'

    trf -v >> versions.yml
    """

    stub:
    """
    touch coords.dat
    touch versions.yml
    """
}