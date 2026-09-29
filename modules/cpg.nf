process CPG {
    label 'process_medium'

    publishDir "${params.outdir}/cpg",
        mode: 'copy'

    input:
    tuple val(coords), path(sliced_fastas)
 
    output:
    tuple val(coords), path('*.cpg'),          emit: cpgs
    // skipping version file, no version info output by tool

    script:
    """
    cpg_lh ${sliced_fastas} > ${coords}.cpg
    """

    stub:
    """
    touch coords.cpg
    """
}