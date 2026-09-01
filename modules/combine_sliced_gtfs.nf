process COMBINE_SLICED_GTFS {
    label 'process_low'

    publishDir "${params.outdir}/annotation_gtf",
        mode: 'copy'

    input:
    val tool
    path(sliced_gtfs)
 
    output:
    path('*annotation.gtf'),       emit: gtf

    script:
    """
    python ${params.projectdir}/bin/src/python/ensembl/tools/anno/nextflow_utils/combine_gtf_slices.py \
    --sliced_gtfs ${sliced_gtfs} \
    --output_gtf ${tool}.annotation.gtf 

    """

    stub:
    """
    touch annotation.gtf
    """

}