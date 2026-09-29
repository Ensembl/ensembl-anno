process DUST {
    label 'process_medium'

    publishDir "${params.outdir}/dust",
        mode: 'copy'

    input:
    tuple val(coords), path(sliced_fastas)
 
    output:
    tuple val(coords), path('*.dust'),          emit: dust_repeats
    path "versions.yml",                     emit: versions

    script:
    """
    dustmasker -in ${sliced_fastas} -out ${coords}.dust
    dustmasker -version-full >> versions.yml
    """

    stub:
    """
    touch coords.dust
    touch versions.yml
    """
}