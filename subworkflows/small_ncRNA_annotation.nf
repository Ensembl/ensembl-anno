include { SELECT_RFAM_MODELS } from '../modules/select_rfam_models.nf'
include { CMSEARCH } from '../modules/cmsearch.nf'
include { TRNASCAN } from '../modules/trnascan.nf'
include { MAKE_GTF as MAKE_CMSEARCH_GTF} from '../modules/make_gtf.nf'
include { MAKE_GTF as MAKE_TRNASCAN_GTF} from '../modules/make_gtf.nf'
include { COMBINE_SLICED_GTFS as COMBINE_CMSEARCH_GTFS} from '../modules/combine_sliced_gtfs.nf'
include { COMBINE_SLICED_GTFS as COMBINE_TRNASCAN_GTFS} from '../modules/combine_sliced_gtfs.nf'

workflow SMALL_NCRNA_ANNOTATION {
    take:
    sliced_fasta

    main:

    // tool channels
    cmsearch_ch = channel.of('cmsearch').collect()
    trnascan_ch = channel.of('trnascan').collect()

    // TODO add these to config
    rfam_accesion_file_ch = channel.fromPath(params.rfam_accesion_file)
    rfam_cm_db_ch = channel.fromPath(params.rfam_cm_db)
    
    SELECT_RFAM_MODELS(rfam_accesion_file_ch, rfam_cm_db_ch)
    CMSEARCH(sliced_fasta, SELECT_RFAM_MODELS.out.rfam_models.collect())
    MAKE_CMSEARCH_GTF(cmsearch_ch, CMSEARCH.out.something)
    COMBINE_CMSEARCH_GTFS(cmsearch_ch, MAKE_CMSEARCH_GTF.out.gtf.collect())

    TRNASCAN(sliced_fasta)
    MAKE_TRNASCAN_GTF(trnascan_ch, TRNASCAN.out.filter_out)
    COMBINE_TRNASCAN_GTFS(trnascan_ch, MAKE_TRNASCAN_GTF.out.gtf.collect())

    emit:
    cmsearch_gtf  = MAKE_CMSEARCH_GTF.out.gtf
    trnascan_gtf = COMBINE_TRNASCAN_GTFS.out.gtf


}