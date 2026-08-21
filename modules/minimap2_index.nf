process MINIMAP2_INDEX {

    publishDir "${params.outdir}/minimap2_index",
        mode: 'copy'


    input:
    path(genome_fasta)
 
    output:
    path("minimap2_index.mmi"),                emit: index
    path("versions.yml"),                      emit: versions

    script:
    """
    minimap2 -t ${params.n_threads} -d minimap2_index.mmi ${genome_fasta}

    echo 'minimap2 ' > versions.yml
    minimap2 --version >> versions.yml
    """

    stub:
    """
    touch minimap2_index.mmi

    touch versions.yml
    """

}