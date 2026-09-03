process MAKE_GTF {
    label 'process_low'

    publishDir "${params.outdir}/gtfs",
        mode: 'copy'

    input:
    val tool
    tuple val(coords), path(input_file)
 
    output:
    path('*/*.gtf'),       emit: gtf


    script:
    """
    mkdir ${tool[0]}

    python ${params.projectdir}/bin/src/python/ensembl/tools/anno/nextflow_utils/make_gtf.py \
    --input_file ${input_file.join(' ')} \
    --output_gtf ${tool[0]}/${tool[0]}_${coords}.gtf \
    --region_name ${coords} \
    --${tool[0]}


    """

    stub:
    """
    touch tool.gtf

    """

}