process STAR {
    label 'process_high'

    publishDir "${params.outdir}/star",
        mode: 'copy'

    input:
    tuple val(meta), path(fastqs)
    path(index)
    val(transcriptomics_params)
 
    output:
    tuple val(meta), path("${meta.id}Aligned.out.sam"),          emit: sam
    tuple val(meta), path("${meta.id}SJ.out.tab"),               emit: junctions
    path "versions.yml",                                      emit: versions

    script:
    gzip_args = ''
    if (fastqs.join(",").endsWith('.gz')){
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
        --readFilesIn ${fastqs.join(",")} \
        ${gzip_args} --outFileNamePrefix ${meta.id} \
        --outSAMtype SAM \
        --alignIntronMax ${transcriptomics_params.max_intron_length} 

    echo 'STAR ' > versions.yml
    STAR --version >> versions.yml
    """

    stub:
    """
    touch ${meta.id}Aligned.out.sam
    touch ${meta.id}SJ.out.tab

    touch versions.yml
    """
}