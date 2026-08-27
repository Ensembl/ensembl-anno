process STRINGTIE {
    label 'process_medium'

    publishDir "${params.outdir}/stringtie",
        mode: 'copy'

    input:
    tuple val(meta), path(bam)
 
    output:
    tuple val(meta), path("${meta}.stringtie.gtf"),               emit: stringtie_gtf
    path "versions.yml",                                        emit: versions

    script:
    """
    stringtie ${bam} \
    -o ${meta}.stringtie.gtf \
    -p ${params.n_threads} \
    -t -a 15   # disable trimming of predicted transcripts based on coverage + minimum anchor length for junctions
 

    echo 'stringtie' > versions.yml
    stringtie --version >> versions.yml
    """

    stub:
    """
    touch ${meta}.stringtie.gtf

    touch versions.yml
    """
}