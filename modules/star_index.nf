process STAR_INDEX {
    label 'process_high'

    publishDir "${params.outdir}/star_index",
        mode: 'copy'

    input:
    path(genome_fasta)
    val(genomeSAindexNbases)
 
    output:
    path("star"),                              emit: index
    path("versions.yml"),                      emit: versions

    script:
    """
    STAR \
    --runThreadN ${params.n_threads} \
    --runMode genomeGenerate \
    --outFileNamePrefix star/ \
    --genomeDir star \
    --genomeFastaFiles ${genome_fasta} \
    --genomeSAindexNbases ${genomeSAindexNbases}
    


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