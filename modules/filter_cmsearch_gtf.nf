process FILTER_CMSEARCH_GTF {
    label 'process_low'

    publishDir "${params.outdir}/gtfs/rnafold_filtered_cmsearch_gtfs/",
        mode: 'copy'

    input:
    tuple val(coords), path(rnafold_results), path(input_gtf) 

    output:
    path('*.gtf'),        optional:true, emit: gtf

    script:

    """
    python ${params.projectdir}/bin/src/python/ensembl/tools/anno/nextflow_utils/filter_cmsearch_gtf.py \
    --input_gtf ${input_gtf} \
    --output_gtf ${coords}.gtf \
    --rnafold_predictions ${rnafold_results} 
    """

    stub:
    """
    touch filtered.gtf
    """

}