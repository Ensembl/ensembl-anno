process SAMTOOLS {

    publishDir "${params.outdir}/samtools",
        mode: 'copy'

    input:
    tuple val(meta), path(sam)
 
    output:
    tuple val(meta), path("${meta}Aligned.out.bam"),          emit: bam
    path "versions.yml",                                      emit: versions

    script:
    """
    samtools sort -@ ${params.n_threads} -o ${params.meta}Aligned.out.bam ${sam}

    samtools --version >> versions.yml
    """

    stub:
    """
    touch ${meta}Aligned.out.bam    
    touch versions.yml
    """
}