process GENBLAST {
    label 'process_medium'

    publishDir "${params.outdir}/genblast",
        mode: 'copy'

    input:
    tuple val(masked_fasta_filename), path(fasta_db)
    tuple val(protein_source), val(slice_id), path(protein_slice)
    path(alignscores)
 
    output:
    tuple val("${protein_source}_${slice_id}"), path('*.gff'),  emit: gff

    script:
    """
    /hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/bin/genblast -p genblastg \
        -q ${protein_slice} \
        -t ${fasta_db}/${masked_fasta_filename} \
        -g T -pid -r 1 -P blast -gff -e 1e-1 -c 0.8 \
        -W 3 -softmask -scodon 50 -i 30 -x 10 -n 30 \
        -d ${params.max_intron_length} \
        -o ${protein_slice} || echo "genblast finished"

    """

    stub:
    """
    touch ${protein_slice.name}.gff
    """
}