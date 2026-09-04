process CMSEARCH {
    label 'process_medium'

    publishDir "${params.outdir}/cmsearch",
        mode: 'copy'

    input:
    tuple val(coords), path(sliced_fastas)
    path(rfam_models)
 
    output:
    tuple val(coords), path('*.cpg'),          emit: cpgs
    //skip versions file as no version info provided by cmsearch

    script:
    """
    cmsearch --rfam \
        --cpu ${params.n_threads} \
        --nohmmonly \
        --cut_ga \
        --tblout", ${coords}.tblout \
        ${rfam_models} ${sliced_fastas}
    """

    stub:
    """
    touch coords.tblout
    """
}