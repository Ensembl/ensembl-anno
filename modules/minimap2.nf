process MINIMAP2 {
    label 'process_high'

    publishDir "${params.outdir}/minimap2",
        mode: 'copy'

    input:
    tuple val(meta), path(fastq)
    path(index)
    val(transcriptomics_params)
 
    output:
    tuple val(meta), path("${meta.id}.sam"),          emit: sam
    path "versions.yml",                              emit: versions

    script:
    """
    minimap2 \
    -G ${transcriptomics_params.max_intron_length} \
    -t ${params.n_threads} \
    --cs \
    --secondary=no \
    -ax splice \
    -u b \
    ${index} \
    ${fastq} \
    -o ${meta.id}.sam

    echo 'minimap2 ' > versions.yml
    minimap2 --version >> versions.yml
    """

    stub:
    """
    touch ${meta.id}.sam

    touch versions.yml
    """
}