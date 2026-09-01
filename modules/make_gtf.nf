process MAKE_GTF {
    label 'process_low'

    publishDir "${params.outdir}/sliced_gtfs",
        mode: 'copy'

    input:
    val tool
    tuple val(coords), path(input_file)
 
    output:
    path('*.gtf'),       emit: gtf


    script:
    """
    python ${params.projectdir}/bin/src/python/ensembl/tools/anno/nextflow_utils/make_gtf.py \
    --input_file ${input_file} \
    --output_gtf ${tool}_${coords}.gtf \
    --region_name ${coords} \
    --${tool}


    """

    stub:
    """
    touch tool.gtf

    """

}