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
    python ${params.projectdir}/bin/src/python/ensembl/tools/anno/nextflow_utils/calculate_genomeSAindexNbases.py \
    --genome_file ${genome_fasta}
    """

    stub:
    """
    echo 14
    """

}