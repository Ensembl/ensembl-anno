process CONVERT_TO_BLASTMASK {
    label 'process_medium'

    publishDir "${params.outdir}/convert2blastmask",
        mode: 'copy'

    input:
    path(masked_fasta)
 
    output:
    path "${masked_fasta.name}.asnb",             emit: asnb             
    path "versions.yml",                              emit: versions

    script:
    """
    convert2blastmask \
        -in ${masked_fasta} \
        -parse_seqids \
        -masking_algorithm other \
        -masking_options \"REpeatDetector, default\" \
        -outfmt maskinfo_asn1_bin \
        -out ${masked_fasta.name}.asnb

    convert2blastmask -version >> versions.yml
    """

    stub:
    """
    touch ${masked_fasta.name}.asnb

    touch versions.yml
    """
}