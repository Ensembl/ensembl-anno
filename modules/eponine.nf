process EPONINE {
    label 'process_medium'

    publishDir "${params.outdir}/eponine",
        mode: 'copy'

    input:
    tuple val(coords), path(sliced_fastas)
    val(simple_features_params)
 
    output:
    tuple val(coords), path('*.epo'),          emit: epo
    // NB no versions file, info not available on command line

    script:
    """
    java -jar ${simple_features_params.eponine_bin} \
        -threshold ${simple_features_params.eponine_threshold} \
        -seq ${sliced_fastas} >> ${coords}.epo
    """

    stub:
    """
    touch coords.epo
    """
}