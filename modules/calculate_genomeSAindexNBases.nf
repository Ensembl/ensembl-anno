process CALCULATE_GENOMESAINDEXNBASES {

    input:
    path(genome_fasta)
 
    output:
    stdout emit: genomeSAindexNbases

    script:
    """
    python bin/src/python/ensembl/tools/anno/nextflow_utils/calculate_star_index_bases.py \
    --genome_file ${genome_fasta}
    """

    stub:
    """
    echo 14
    """

}