process COMBINE_SLICED_GTFS {
    label 'process_low'

    publishDir "${params.outdir}/gtfs",
        mode: 'copy'

    input:
    val tool
    path(sliced_gtfs)
 
    output:
    path('*/*annotation.gtf'),       emit: gtf

    script:
    """
    mkdir ${tool[0]}
    python ${params.projectdir}/bin/src/python/ensembl/tools/anno/nextflow_utils/combine_gtf_slices.py \
    --sliced_gtfs ${sliced_gtfs} \
    --output_gtf ${tool[0]}/${tool[0]}_annotation.gtf \
    --tool ${tool[0]}

    """

    stub:
    """
    mkdir tool
    touch tool/annotation.gtf
    """

}