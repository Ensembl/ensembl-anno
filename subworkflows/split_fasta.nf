include { FIND_FASTA_INTERVALS } from '../modules/find_fasta_intervals_for_splitting.nf'
include { BEDTOOLS } from '../modules/bedtools.nf'

workflow SPLIT_FASTA {
    take:
    fasta

    main:
    FIND_FASTA_INTERVALS(fasta)

    interval_bed_ch = FIND_FASTA_INTERVALS.out.beds.flatten().map{
        file -> tuple(file.baseName.split('\\.bed')[0], file) 
    }
    
    BEDTOOLS(fasta.collect(), interval_bed_ch)


    emit:
    BEDTOOLS.out.fasta_slice

}
