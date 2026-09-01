process RED {
    label 'process_high'

    publishDir "${params.outdir}/red",
        mode: 'copy'

    input:
    path(fasta)
 
    output:
    path('red_masked_genome_dir/'),          emit: red_masked_genome_dir
    path('red_repeats_dir/'),                emit: red_repeats_dir
    path('red_masked_genome_dir/*.msk'),     emit: red_masked_genome_files
    tuple val('whole_genome'), path('red_repeats_dir/*.rpt'),           emit: red_repeats_files
    path "versions.yml",                     emit: versions

    script:
    """
    mkdir red_masked_genome_dir
    mkdir red_repeats_dir

    # This is a bit silly. We have to pass the directoy containing the genome fasta to Red rather
    # than the fasta itself, hence the logic below:
    mkdir input_genome_dir
    cp ${fasta} input_genome_dir/ #TODO make sure this doesn't persist
    
    Red -gnm input_genome_dir -msk red_masked_genome_dir -rpt red_repeats_dir

    Red -v > versions.yml
    """

    stub:
    """
    mkdir red_masked_genome_dir
    mkdir red_repeats_dir
    touch red_masked_genome_dir/genome.msk
    touch red_repeats_dir/genome.rpt

    touch versions.yml
    """
}