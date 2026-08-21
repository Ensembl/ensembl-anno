process SCALLOP {

    publishDir "${params.outdir}/scallop",
        mode: 'copy'


    input:
    tuple val(meta), path(bam)
 
    output:
    tuple val(meta), path("${meta}.scallop.gtf"),          emit: scallop_gtf
    path "versions.yml",                                   emit: versions

    script:
    """
    scallop -i ${bam} -o ${meta}.scallop.gtf --min_flank_length 10

    echo 'scallop' > versions.yml
    scallop --version >> versions.yml
    """

    stub:
    """
    touch ${meta}.scallop.gtf

    touch versions.yml
    """
}