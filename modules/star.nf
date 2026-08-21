process STAR {

    publishDir "${params.outdir}/star",
        mode: 'copy'

    input:
    tuple val(meta), path(fastqs)
    path(index)
 
    output:
    tuple val(meta), path("${meta}Aligned.out.sam"),          emit: sam
    tuple val(meta), path("${meta}SJ.out.tab"),               emit: junctions
    path "versions.yml",                                      emit: versions

    script:
    gzip_args = ''
    if (fastqs[0].endsWith('.gz')){
        gzip_args = '--readFilesCommand gunzip -c'
    }
    """
    STAR \
        --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --outSAMstrandField intronMotif \
        --runThreadN ${params.n_threads} \
        --twopassMode Basic \
        --runMode alignReads \
        --genomeDir ${index} \
        --readFilesIn ${fastqs.join(",")} ${gzip_args}\ 
        --outFileNamePrefix ${meta} \
        --outSAMtype SAM \
        --alignIntronMax ${params.star_max_intron_length} 

    echo 'STAR ' > versions.yml
    STAR --version >> versions.yml
    """

    stub:
    """
    touch ${meta}Aligned.out.sam
    touch ${meta}SJ.out.tab

    touch versions.yml
    """
}