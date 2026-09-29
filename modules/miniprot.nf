process MINIPROT {
    label 'process_medium'

    publishDir "${params.outdir}/miniprot",
        mode: 'copy'

    input:
    path(miniprot_index)
    tuple val(protein_source), path(protein_file)

    output:
    tuple val("${protein_source}"),  path("${protein_source}.miniprot.gff"),    emit: gff
    path("versions.yml"),                      emit: versions

    script:
    """
    /homes/ereboperezsilva/.local/bin/miniprot \
        -t ${params.n_threads} \
        -N 1 \
        --outs=1.0 \
        --gff \
        ${miniprot_index} \
        ${protein_file} > ${protein_source}.miniprot.gff

    /homes/ereboperezsilva/.local/bin/miniprot --version >> versions.yml
    """

    stub:
    """
    touch miniprot.mpi
    touch versions.yml
    """
}