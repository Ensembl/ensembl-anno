process CHECK_TRANSCRIPTOMIC_OUTPUT {

    input:
    path(gtf)
 
    output:
    path('transcript_output_checks_passed.txt'), emit: checks_passed
    path('transcript_output_checks.log'), emit: checks_log


    script:
    """
    python bin/src/python/ensembl/tools/anno/nextflow_utils/beds_to_gtf.py \
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