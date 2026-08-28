process MAKE_GTF {
    label 'process_low'

    publishDir "${params.outdir}/red_gtf",
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
    --output_gtf ${tool}.gtf \
    --coords ${coords} \
    --${tool}

    Red --version >> versions.yml
    """

    stub:
    """
    touch tool.gtf
    touch versions.yml
    """

}