process COMBINE_SLICED_GTFS {
    label 'process_low'

    publishDir "${params.outdir}/gtfs",
        mode: 'copy'

    input:
    tuple val(tool), path(rfam_seed), path(rfam_selected_model)
    path(sliced_gtfs)
 
    output:
    path('*/*annotation.gtf'),       emit: gtf

    script:
    """
    mkdir ${tool}
    python ${params.projectdir}/bin/src/python/ensembl/tools/anno/nextflow_utils/combine_gtf_slices.py \
    --sliced_gtfs ${sliced_gtfs} \
    --output_gtf ${tool}/${tool}_annotation.gtf \
    --tool ${tool}


    """

    stub:
    """
    mkdir tool
    touch tool/annotation.gtf
    """

}