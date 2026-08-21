process MINIMAP2 {

    publishDir "${params.outdir}/minimap2",
        mode: 'copy'


    // TODO change
    input:
    tuple val(meta), path(input_file)
 
    output:
    tuple val(meta), path("${meta.id}_trimmed.fastq"),          emit: trimmed_reads
    path "versions.yml",                                        emit: versions

    script:
    """
    echo "A output" > ${meta.id}_A.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool_a: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_A.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tool_a: 1.0.0
    END_VERSIONS
    """
}