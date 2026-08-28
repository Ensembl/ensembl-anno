process REPEATMASKER {
    label 'process_medium'
    label 'repeats'

    publishDir "${params.outdir}/repeatmasker",
        mode: 'copy'

    input:
    tuple val(coords), path(sliced_fasta)
 
    output:
    tuple val(coords), path('*.out'),          emit: repeatmasker_repeats
    path "versions.yml",                     emit: versions

    script:
    def args = ''
    if (params.library==null){
        def species = params.species
        if (species == null){
            species = 'homo'
        }
        args = ' --species ' + species
    } else {
        args = ' --lib ' + params.library
    }
    """
    mkdir repeatmasker_out
    RepeatMasker -nolow -engine rmblast -dir repeatmasker_out ${args} ${sliced_fasta}
    RepeatMasker -v >> versions.yml
    """

    stub:
    """
    touch coords.out
    touch versions.yml
    """
}