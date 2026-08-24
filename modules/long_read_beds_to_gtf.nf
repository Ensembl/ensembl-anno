process LONG_READ_BEDS_TO_GTF {

    input:
    path(bed)
 
    output:
    path('long_read_transcripts.gtf'), emit: gtf

    script:
    """
    python bin/src/python/ensembl/tools/anno/nextflow_utils/beds_to_gtf.py \
    --bedfile_list ${bed} \
    --gtf_path long_read_transcripts.gtf
    """

    stub:
    """
    touch long_read_transcripts.gtf
    """

}