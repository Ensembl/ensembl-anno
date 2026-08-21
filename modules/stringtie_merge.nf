process STRINGTIE_MERGE {

    publishDir "${params.outdir}/merged_gtfs",
        mode: 'copy'


    input:
    tuple val(meta), path(gtf)
    val toolname
 
    output:
    path("${toolname}.gtf"),          emit: merged_gtf
    path "versions.yml",              emit: versions

    script:
    """
    # ToDo test whether need to skip emtpy gtfs
    stringtie --merge -o ${toolname}.gtf ${gtf}

    echo 'stringtie' > versions.yml
    stringtie --version >> versions.yml
    """

    stub:
    """
    touch ${toolname}.gtf

    touch versions.yml
    """
}