process RNAFOLD {
    label 'process_medium'

    publishDir "${params.outdir}/rnafold",
        mode: 'copy'

    input:
    tuple val(coords), path(rna_fasta)
 
    output:
    tuple val(coords), path('*rnafold_energy_scores.txt'),          emit: rnafold_energy_scores
    // skipping version file, no version info output by tool

    script:
    """
    RNAfold --infile ${rna_fasta} >> ${coords}_rnafold_energy_scores.txt
    """

    stub:
    """
    touch rnafold_energy_scores.txt
    """
}