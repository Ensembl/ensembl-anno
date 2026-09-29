process CALCULATE_GENOMESAINDEXNBASES {

    publishDir "${params.outdir}/star_index",
        mode: 'copy'
    label 'process_light'

    input:
    path(genome_fasta)
 
    output:
    stdout emit: genomeSAindexNbases

    script:
    """
    python ${params.projectdir}/bin/src/python/ensembl/tools/anno/nextflow_utils/fasta_operations.py \
    --calculate_genomeSAindexNbases --genome_file ${genome_fasta} --min_seq_length 0 --maximum_index_bases 14
    """

    stub:
    """
    echo 14
    """

}