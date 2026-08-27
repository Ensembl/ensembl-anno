process PAFTOOLS {
    label 'process_medium'

    publishDir "${params.outdir}/paftools",
        mode: 'copy'

    input:
    tuple val(meta), path(sam)
 
    output:
    tuple val(meta), path("${meta}.bed"),          emit: bed
    path "versions.yml",                              emit: versions

    script:
    """
    paftools.js splice2bed ${sam} > ${meta}.bed

    echo 'paftools.js ' > versions.yml
    paftools.js version  >> versions.yml
    """

    stub:
    """
    touch ${meta}.bed

    touch versions.yml
    """
}