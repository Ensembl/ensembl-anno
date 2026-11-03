process COLLAPSE_IG_ANNOTATIONS {
    label 'process_low'

    publishDir "${params.outdir}/collapsed_ig_annotation",
        mode: 'copy'

    input:
    tuple val(protein_source), path(miniprot_gff)

    output:
    path('ig_annotation.gtf'),        emit: gtf

    script:
    """
    python ${params.projectdir}/bin/src/python/ensembl/tools/anno/nextflow_utils/ig_collapse.py \
        --miniprot_gff ${miniprot_gff} \
        --collapsed_gtf ig_annotation.gtf
    """

    stub:
    """
    touch ig_annotation.gtf
    """

}