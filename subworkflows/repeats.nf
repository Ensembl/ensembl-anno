include { RED } from '../modules/red.nf'
include { DUST } from '../modules/dust.nf'
include { REPEATMASKER } from '../modules/repeatmasker.nf'
include { TRF } from '../modules/trf.nf'
include { MAKE_GTF as MAKE_RED_GTF} from '../modules/make_gtf.nf'
include { MAKE_GTF as MAKE_DUST_GTF} from '../modules/make_gtf.nf'
include { MAKE_GTF as MAKE_REPEATMASKER_GTF} from '../modules/make_gtf.nf'
include { MAKE_GTF as MAKE_TRF_GTF} from '../modules/make_gtf.nf'

workflow REPEATS {
    take:
    fasta
    sliced_fasta_dir
    sliced_fasta

    main:

    // tool channels
    red_ch = channel.of('red')
    dust_ch = channel.of('dust')
    repeatmasker_ch = channel.of('repeatmasker')
    trf_ch = channel.of('trf')

    RED(fasta)
    MAKE_RED_GTF(red_ch, RED.out.repeats)

    DUST(sliced_fasta)
    MAKE_DUST_GTF(dust_ch, DUST.out.dust_repeats)

    REPEATMASKER(sliced_fasta)
    MAKE_REPEATMASKER_GTF(repeatmasker_ch, REPEATMASKER.out.repeatmasker_repeats)

    TRF(sliced_fasta)
    MAKE_TRF_GTF(trf_ch, TRF.out.trf_repeats)

    //TODO collapse dust, repeatmasker and trf gtfs



    emit:
    red_gtf  = RED.out.gtf
    dust_gtf = MERGE_DUST_GTFS.out.merged_gtf
    trf_gtf = MERGE_TRF_GTFS.out.merged_gtf 
    repeatmasker_gtf = MERGE_REPEATMASKER_GTFS.out.merged_gtf


}