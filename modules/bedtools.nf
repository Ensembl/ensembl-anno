process BEDTOOLS {
    label 'process_medium'

    publishDir "${params.outdir}/fasta_intervals",
        mode: 'copy'

    input:
    path(fasta)
    tuple val(coords), path(bed)
 
    output:
    tuple val(coords), path("${coords}.fa"),               emit: fasta_slice             
    path "versions.yml",                                        emit: versions

    script:
    """
    bedtools getfasta -fi ${fasta} -bed ${bed} > ${coords}.fa

    bedtools --version >> versions.yml
    """

    stub:
    """
    touch ${coords}.fa


    touch versions.yml
    """
}