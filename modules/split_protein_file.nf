process SPLIT_PROTEIN_FILE {
    label 'process_low'

    publishDir "${params.outdir}/sliced_proteins",
        mode: 'copy'

    input:
    tuple val(source), path(proteins)
 
    output:
    tuple val(source), path("${source}/bin*/*.fa"),       emit: sliced_proteins

    script:
    """
    python ${params.projectdir}/bin/src/python/ensembl/tools/anno/nextflow_utils/slice_protein_file.py \
    --proteins ${proteins} \
    --protein_source ${source} \

    """

    stub:
    """
    mkdir ${source}
    mkdir ${source}/bin1
    mkdir ${source}/bin2
    mkdir ${source}/bin3
    mkdir ${source}/bin4
    mkdir ${source}/bin5
    mkdir ${source}/bin6
    mkdir ${source}/bin7
    mkdir ${source}/bin8
    mkdir ${source}/bin9
    mkdir ${source}/bin10
    touch ${source}/bin1/1.fa
    touch ${source}/bin5/2.fa
    """

}