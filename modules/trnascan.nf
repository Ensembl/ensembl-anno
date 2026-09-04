process TRNASCAN {
    label 'process_medium'

    publishDir "${params.outdir}/trnascan",
        mode: 'copy'

    input:
    tuple val(coords), path(sliced_fastas)
 
    output:
    tuple val(coords), path('*.trna'),          emit: trna
    tuple val(coords), path('*.ss'),            emit: ss
    path('tRNAscan_output/*.out'),              emit: filter_out
    path('tRNAscan_output/*.log'),              emit: log
    path('tRNAscan_output/*.ss'),               emit: ss
    path('versions.yml'),                       emit: versions

    script:
    """
    mkdir tRNAscan_output

    tRNAscan-SE ${sliced_fastas} \
        -o ${coords}.trna \
        -f ${coords}.ss \
        -H -q --detail -Q

    EukHighConfidenceFilter \
        --result ${coords}.trna \
        --ss ${coords}.ss \
        --output tRNAscan_output \
        --prefix ${coords}.filt

    tRNAscan-SE -h >> versions.yml
    """

    stub:
    """
    touch coords.trna
    touch coords.ss

    mkdir tRNAscan_output
    touch tRNAscan_output/coords.out
    touch tRNAscan_output/coords.log 
    touch tRNAscan_output/coords.ss

    touch versions.yml
    """
}