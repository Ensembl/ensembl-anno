process SAMTOOLS {
    label 'process_high'

    publishDir "${params.outdir}/samtools",
        mode: 'copy'

    input:
    tuple val(meta), path(sam)
 
    output:
    tuple val(meta), path("${meta.id}Aligned.out.bam"),          emit: bam
    path "versions.yml",                                      emit: versions

    script:
    """
    samtools sort -@ ${params.n_threads} -o ${meta.id}Aligned.out.bam ${sam}

    samtools --version >> versions.yml
    """

    stub:
    """
    touch ${meta.id}Aligned.out.bam    
    touch versions.yml
    """
}