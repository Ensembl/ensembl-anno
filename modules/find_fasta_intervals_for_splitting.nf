process FIND_FASTA_INTERVALS {

    publishDir "${params.outdir}/fasta_intervals",
        mode: 'copy'
    label 'process_light'

    input:
    path(genome_fasta)
 
    output:
    path('*.bed'), emit:beds

    script:
    """
    python ${params.projectdir}/bin/src/python/ensembl/tools/anno/nextflow_utils/fasta_operations.py \
    --splitFasta --genome_file ${genome_fasta} --min_seq_length 5000 --slice_size 1000000
    """

    stub:
    """
    touch chr1_1_2.bed
    """

}