process MAKE_GTF {
    label 'process_low'

    publishDir "${params.outdir}/gtfs",
        mode: 'copy'

    input:
    tuple val(tool), path(rfam_seed), path(rfam_selected_model)
    tuple val(coords), path(input_file)

 
    output:
    tuple val(coords), path('*/*.gtf'),        optional:true, emit: gtf
    tuple val(coords), path('*/*.bed'),       optional:true, emit: bed


    script:
    def region_name = coords.split(':')[0]
    cmsearch_args = ''
    if (tool == "cmsearch"){
        cmsearch_args = " --rfam_seed_descriptions ${rfam_seed} \
        --rfam_selected_models_file ${rfam_selected_model} \
        --output_bed ${tool}/${tool}_${coords}.bed"
    }
    """
    mkdir ${tool}

    python ${params.projectdir}/bin/src/python/ensembl/tools/anno/nextflow_utils/make_gtf.py \
    --input_file ${input_file.join(' ')} \
    --output_gtf ${tool}/${tool}_${coords}.gtf \
    --region_name ${region_name} \
    --${tool} ${cmsearch_args}


    """

    stub:
    """
    mkdir tool
    touch tool/tool.gtf

    """

}