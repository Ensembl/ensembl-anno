process REPEATMASKER {
    label 'process_high'
    label 'repeats'

    publishDir "${params.outdir}/repeatmasker",
        mode: 'copy'

    input:
    tuple val(coords), path(sliced_fasta)
    val(repeats_params)
 
    output:
    tuple val(coords), path('*/*.out'),          emit: repeatmasker_repeats
    path "versions.yml",                     emit: versions

    script:
    def args = ''
    if (repeats_params.library==null){
        def species = repeats_params.species
        if (species == null){
            species = 'homo'
        }
        args = ' --species ' + species
    } else {
        args = ' --lib ' + repeats_params.library
    }
    """
    mkdir repeatmasker_out
    RepeatMasker -nolow -engine rmblast -dir repeatmasker_out ${args} ${sliced_fasta}
    RepeatMasker -v >> versions.yml
    """

    stub:
    """
    mkdir repeatmasker_out
    touch repeatmasker_out/coords.out
    touch versions.yml
    """
}