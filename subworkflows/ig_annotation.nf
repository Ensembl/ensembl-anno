include { MINIPROT_INDEX } from '../modules/miniprot_index.nf'
include { MINIPROT as MINIPROT_IG_ALIGNMENT } from '../modules/miniprot.nf'
include { COLLAPSE_IG_ANNOTATIONS } from '../modules/collapse_ig_annotations.nf'


workflow IG_ANNOTATION {
    take:
    masked_genome_fasta
    igmt_fasta

    main:

    MINIPROT_INDEX(masked_genome_fasta)
    MINIPROT_IG_ALIGNMENT(MINIPROT_INDEX.out.index.collect(), igmt_fasta)
    COLLAPSE_IG_ANNOTATIONS(MINIPROT_IG_ALIGNMENT.out.gff)

    emit:
    ig_annotation_gtf  = COLLAPSE_IG_ANNOTATIONS.out.gtf
    

}