include { IG_ANNOTATION } from '../subworkflows/ig_annotation.nf'


workflow MAIN_PIPELINE {
    take:
    fasta                   // [fasta]
    igmt_fasta             //  [fasta]

    main:
    IG_ANNOTATION(fasta, igmt_fasta)

}
