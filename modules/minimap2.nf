process MINIMAP2 {

    publishDir "${params.outdir}/minimap2",
        mode: 'copy'


    // TODO change
    input:
    tuple val(meta), path(fastq)
    path(index)
 
    output:
    tuple val(meta), path("${meta}.sam"),          emit: sam
    path "versions.yml",                              emit: versions

    script:
    """
    minimap2 \
    -G ${params.max_intron_length} \
    -t ${params.n_threads}
    --cs \
    --secondary=no \
    -ax splice \
    -u b \
    ${index} \
    ${fastq} \
    -o ${meta}.sam

    echo 'minimap2 ' > versions.yml
    minimap2 --version >> versions.yml
    """

    stub:
    """
    touch ${meta}.sam

    touch versions.yml
    """
}