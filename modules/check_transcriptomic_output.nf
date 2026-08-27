process CHECK_TRANSCRIPTOMIC_OUTPUT {
    label 'process_low'

    publishDir "${params.outdir}/check_transcriptomic_gtfs",
        mode: 'copy'

    input:
    path(gtf)
 
    output:
    path('transcript_output_checks_passed.txt'), emit: checks_passed
    path('transcript_output_checks.log'), emit: checks_log


    script:
    """
    python ${params.projectdir}/bin/src/python/ensembl/tools/anno/nextflow_utils/check_transcriptomic_output.py \
        --gtfs ${gtf} \
        --min_lines ${params.min_total_transcriptomic_gtf_lines} \
        --logfile transcript_output_checks.log \
        --checks_passed_file transcript_output_checks_passed.txt
    """

    stub:
    """
    touch transcript_output_checks.log
    touch transcript_output_checks_passed.txt
    """

}