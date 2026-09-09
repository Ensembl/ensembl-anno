include { SELECT_RFAM_MODELS } from '../modules/select_rfam_models.nf'
include { CMSEARCH } from '../modules/cmsearch.nf'
include { TRNASCAN } from '../modules/trnascan.nf'
include { MAKE_GTF as MAKE_CMSEARCH_GTF} from '../modules/make_gtf.nf'
include { MAKE_GTF as MAKE_TRNASCAN_GTF} from '../modules/make_gtf.nf'
include { COMBINE_SLICED_GTFS as COMBINE_FILTERED_CMSEARCH_GTFS} from '../modules/combine_sliced_gtfs.nf'
include { COMBINE_SLICED_GTFS as COMBINE_TRNASCAN_GTFS} from '../modules/combine_sliced_gtfs.nf'
include { BEDTOOLS } from '../modules/bedtools.nf'
include { RNAFOLD } from '../modules/rnafold.nf'
include { FILTER_CMSEARCH_GTF } from '../modules/filter_cmsearch_gtf.nf'

workflow SMALL_NCRNA_ANNOTATION {
    take:
    fasta
    sliced_fasta

    main:

    rfam_accession_file_ch = channel.fromPath(params.rfam_accession_file)
    rfam_cm_db_ch = channel.fromPath(params.rfam_cm_db)
    
    SELECT_RFAM_MODELS(rfam_accession_file_ch, rfam_cm_db_ch)
    CMSEARCH(sliced_fasta, SELECT_RFAM_MODELS.out.rfam_models.collect())
    
    cmsearch_ch = channel.of(tuple('cmsearch', file(params.rfam_seeds_file))).combine(
        SELECT_RFAM_MODELS.out.rfam_models
    ).collect()
    MAKE_CMSEARCH_GTF(cmsearch_ch, CMSEARCH.out.tblout)

    // This could be further optimised. It is a bit inefficient that we search
    // against the entire fasta when we know the sequence is in a particular slice
    BEDTOOLS(fasta.collect(), MAKE_CMSEARCH_GTF.out.bed) 
    RNAFOLD(BEDTOOLS.out.fasta_slice)

    rnafold_cmsearch_ch = RNAFOLD.out.rnafold_energy_scores.join(
        MAKE_CMSEARCH_GTF.out.gtf
    )
    FILTER_CMSEARCH_GTF(rnafold_cmsearch_ch)
    COMBINE_FILTERED_CMSEARCH_GTFS(cmsearch_ch, FILTER_CMSEARCH_GTF.out.gtf.collect())

    // Run TRNAscan
    trnascan_ch = channel.of(tuple('trnascan'), file('optional1'), file('optional2')).collect()
    TRNASCAN(sliced_fasta)
    MAKE_TRNASCAN_GTF(trnascan_ch, TRNASCAN.out.filter_out)
    COMBINE_TRNASCAN_GTFS(trnascan_ch, MAKE_TRNASCAN_GTF.out.gtf.map{it -> it[1]}.collect())

    emit:
    cmsearch_gtf  = MAKE_CMSEARCH_GTF.out.gtf
    trnascan_gtf = COMBINE_TRNASCAN_GTFS.out.gtf

}