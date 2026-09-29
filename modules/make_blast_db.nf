process MAKE_BLAST_DB {
    label 'process_medium'

    publishDir "${params.outdir}/blast_db",
        mode: 'copy'

    input:
    path(masked_fasta)
    path(asnb)
 
    output:
    tuple val(masked_fasta.name), path("fasta_db"),             emit: fasta_db          
    path "versions.yml",         emit: versions
    //path "fasta_db/${masked_fasta.name}.ndb" //only purpose of this is to ensure makeblastdb actually runs

    script:
    """
    mkdir fasta_db
    makeblastdb \
        -in ${masked_fasta} \
        -dbtype nucl \
        -parse_seqids \
        -mask_data ${asnb} \
        -max_file_sz 4GB \
        -out fasta_db/${masked_fasta.name}
    
    cp ${masked_fasta} fasta_db/
    cp ${asnb} fasta_db/

    makeblastdb -version >> versions.yml
    """

    stub:
    """
    mkdir fasta_db
    touch fasta_db/${masked_fasta.name}.ndb

    touch versions.yml
    """
}