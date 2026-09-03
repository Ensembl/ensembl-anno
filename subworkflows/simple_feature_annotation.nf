include { CPG } from '../modules/cpg.nf'
include { EPONINE } from '../modules/eponine.nf'
include { MAKE_GTF as MAKE_CPG_GTF} from '../modules/make_gtf.nf'
include { MAKE_GTF as MAKE_EPONINE_GTF} from '../modules/make_gtf.nf'
include { COMBINE_SLICED_GTFS as COMBINE_CPG_GTFS} from '../modules/combine_sliced_gtfs.nf'
include { COMBINE_SLICED_GTFS as COMBINE_EPONINE_GTFS} from '../modules/combine_sliced_gtfs.nf'

workflow SIMPLE_FEATURE_ANNOTATION {
    take:
    sliced_fasta

    main:

    // tool channels
    cpg_ch = channel.of('cpg').collect()
    eponine_ch = channel.of('eponine').collect()

    CPG(sliced_fasta)
    MAKE_CPG_GTF(cpg_ch, CPG.out.cpgs)
    COMBINE_CPG_GTFS(cpg_ch, MAKE_CPG_GTF.out.gtf.collect())

    EPONINE(sliced_fasta)
    MAKE_EPONINE_GTF(eponine_ch, EPONINE.out.epo)
    COMBINE_EPONINE_GTFS(eponine_ch, MAKE_EPONINE_GTF.out.gtf.collect())

    emit:
    cpg_gtf  = MAKE_CPG_GTF.out.gtf
    eponine_gtf = COMBINE_EPONINE_GTFS.out.gtf


}