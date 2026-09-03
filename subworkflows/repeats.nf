include { RED } from '../modules/red.nf'
include { DUST } from '../modules/dust.nf'
include { REPEATMASKER } from '../modules/repeatmasker.nf'
include { TRF } from '../modules/trf.nf'
include { MAKE_GTF as MAKE_RED_GTF} from '../modules/make_gtf.nf'
include { MAKE_GTF as MAKE_DUST_GTF} from '../modules/make_gtf.nf'
include { MAKE_GTF as MAKE_REPEATMASKER_GTF} from '../modules/make_gtf.nf'
include { MAKE_GTF as MAKE_TRF_GTF} from '../modules/make_gtf.nf'
include { COMBINE_SLICED_GTFS as COMBINE_DUST_GTFS} from '../modules/combine_sliced_gtfs.nf'
include { COMBINE_SLICED_GTFS as COMBINE_REPEATMASKER_GTFS} from '../modules/combine_sliced_gtfs.nf'
include { COMBINE_SLICED_GTFS as COMBINE_TRF_GTFS} from '../modules/combine_sliced_gtfs.nf'

workflow REPEATS {
    take:
    fasta
    sliced_fasta

    main:

    // tool channels
    red_ch = channel.of('red')
    dust_ch = channel.of('dust')
    repeatmasker_ch = channel.of('repeatmasker')
    trf_ch = channel.of('trf')

    RED(fasta)
    MAKE_RED_GTF(red_ch.collect(), RED.out.red_repeats_files)

    DUST(sliced_fasta)
    MAKE_DUST_GTF(dust_ch.collect(), DUST.out.dust_repeats)
    COMBINE_DUST_GTFS(dust_ch.collect(), MAKE_DUST_GTF.out.gtf.collect())

    REPEATMASKER(sliced_fasta)
    MAKE_REPEATMASKER_GTF(repeatmasker_ch.collect(), REPEATMASKER.out.repeatmasker_repeats)
    COMBINE_REPEATMASKER_GTFS(repeatmasker_ch.collect(), MAKE_REPEATMASKER_GTF.out.gtf.collect())

    TRF(sliced_fasta)
    MAKE_TRF_GTF(trf_ch.collect(), TRF.out.trf_repeats)
    COMBINE_TRF_GTFS(trf_ch.collect(), MAKE_TRF_GTF.out.gtf.collect())

    emit:
    red_gtf  = MAKE_RED_GTF.out.gtf
    dust_gtf = COMBINE_DUST_GTFS.out.gtf
    trf_gtf = COMBINE_TRF_GTFS.out.gtf
    repeatmasker_gtf = COMBINE_REPEATMASKER_GTFS.out.gtf


}