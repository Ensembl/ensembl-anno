process SELECT_RFAM_MODELS {
    label 'process_low'

    publishDir "${params.outdir}/rfam_models",
        mode: 'copy'

    input:
    path(rfam_accession_file)
    path(rfam_cm_db)
 
    output:
    path('rfam_models'),       emit: rfam_models

    script:
    """
    python ${params.projectdir}/bin/src/python/ensembl/tools/anno/nextflow_utils/select_rfam_models.py \
    --rfam_accesion_file ${rfam_accession_file} \
    --rfam_cm_db ${rfam_cm_db} \
    --rfam_selected_models_file rfam_models

    """

    stub:
    """
    touch rfam_models
    """

}