process MINIPROT_INDEX {
    label 'process_medium'

    publishDir "${params.outdir}/miniprot_index",
        mode: 'copy'

    input:
    path(masked_fasta)
 
    output:
    path("miniprot.mpi"),                              emit: index
    path("versions.yml"),                      emit: versions

    script:
    """
    /homes/ereboperezsilva/.local/bin/miniprot \
        -t ${params.n_threads} \
        -d miniprot.mpi \
        ${masked_fasta}

    /homes/ereboperezsilva/.local/bin/miniprot --version >> versions.yml
    """

    stub:
    """
    touch miniprot.mpi
    touch versions.yml
    """

}