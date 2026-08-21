process STAR_INDEX {

    publishDir "${params.outdir}/star",
        mode: 'copy'


    input:
    path(genome_fasta)
 
    output:
    path("star"),                              emit: index
    path("versions.yml"),                      emit: versions

    script:
    """
    # TODO make a much smaller python package containing utility function to get genome size
    # and call it here or in another process

    STAR \
    --runThreadN ${params.n_threads} \
    --runMode genomeGenerate \
    --outFileNamePrefix star/ \
    --genomeDir star \
    --genomeSAindexNbases str(index_bases) \ #TODO fix this
    --genomeFastaFiles ${genome_fasta}


    echo 'STAR ' > versions.yml
    STAR --version >> versions.yml
    """

    stub:
    """
    mkdir star
    touch star/Genome
    touch star/Log.out
    touch star/SA
    touch star/SAindex
    touch star/chrLength.txt
    touch star/chrName.txt
    touch star/chrNameLength.txt
    touch star/chrStart.txt
    touch star/exonGeTrInfo.tab
    touch star/exonInfo.tab
    touch star/geneInfo.tab
    touch star/genomeParameters.txt
    touch star/sjdbInfo.txt
    touch star/sjdbList.fromGTF.out.tab
    touch star/sjdbList.out.tab
    touch star/transcriptInfo

    touch versions.yml
    """

}