include { CALCULATE_GENOMESAINDEXNBASES} from '../modules/calculate_genomeSAindexNBases.nf'
include { STAR_INDEX } from '../modules/star_index.nf'
include { TRIMGALORE } from '../modules/trimgalore.nf'
include { STAR } from '../modules/star.nf'
include { SAMTOOLS } from '../modules/samtools.nf'
include { SCALLOP } from '../modules/scallop.nf'
include { STRINGTIE } from '../modules/stringtie.nf'
include { STRINGTIE_MERGE as MERGE_SCALLOP_GTFS } from '../modules/stringtie_merge.nf'
include { STRINGTIE_MERGE as MERGE_STRINGTIE_GTFS } from '../modules/stringtie_merge.nf'
include { MINIMAP2_INDEX } from '../modules/minimap2_index.nf'
include { MINIMAP2 } from '../modules/minimap2.nf'
include { PAFTOOLS } from '../modules/paftools.nf'
include { LONG_READ_BEDS_TO_GTF } from '../modules/long_read_beds_to_gtf.nf'
include { CHECK_TRANSCRIPTOMIC_OUTPUT } from '../modules/check_transcriptomic_output.nf'

workflow TRANSCRIPTOMICS_ANNOTATION {
    take:
    short_reads  // channel: [ val(meta), [path(read1), path(read2)] ]
    long_reads   // channel: [ val(meta), path(input_file) ]
    genome_fasta // channel: [fasta_path]

    main:

    // Decide whether to make star and minimap2 indexes based on whether short and long read dirs were provided:
    star_index_ch = channel.empty()
    minimap2_index_ch = channel.empty()

    if (params.short_read_dir != null){
        CALCULATE_GENOMESAINDEXNBASES(genome_fasta)
        STAR_INDEX(genome_fasta, CALCULATE_GENOMESAINDEXNBASES.out.genomeSAindexNbases)
        
        // It is neccessary to run collect() to convert a queue channel to a value channel
        // Otherwise the index is consumed after the first fastq is aligned and the rest are skipped
        star_index_ch = STAR_INDEX.out.index.collect()
    }   
    if (params.long_read_dir != null){
        MINIMAP2_INDEX(genome_fasta)

        // It is neccessary to run collect() to convert a queue channel to a value channel
        // Otherwise the index is consumed after the first fastq is aligned and the rest are skipped
        minimap2_index_ch = MINIMAP2_INDEX.out.index.collect()
    }

    // Short reads pipeline
    // This will only run if the short reads channel has been populated with fastqs
    // If short_reads is an empty channel we skip from here to the long reads pipeline
    if (params.trim_reads){
        TRIMGALORE(short_reads)
        STAR(TRIMGALORE.out.trimmed_reads, star_index_ch)
    } else {
        STAR(short_reads, star_index_ch)
    }

    SAMTOOLS(STAR.out.sam)

    SCALLOP(SAMTOOLS.out.bam)
    STRINGTIE(SAMTOOLS.out.bam)

    // Generate a channel containing a list of gtfs to combine for scallop and stringtie
    scallop_gtf_ch = SCALLOP.out.scallop_gtf.map {
        gtf -> gtf[1]
    }.collect()
    stringtie_gtf_ch = STRINGTIE.out.stringtie_gtf.map {
        gtf -> gtf[1]
    }.collect()

    // These processes merge the sample level gtfs into a single gtf of all predictions for each tool
    MERGE_SCALLOP_GTFS(scallop_gtf_ch, channel.value('scallop') )
    MERGE_STRINGTIE_GTFS(stringtie_gtf_ch, channel.value('stringtie') )

    // Long reads pipeline
    // This will only run if the short reads channel has been populated with fastqs
    // If long_reads is an empty channel we skip from here to the end
    MINIMAP2(long_reads, minimap2_index_ch)
    PAFTOOLS(MINIMAP2.out.sam)
    paftools_bed_ch = PAFTOOLS.out.bed.map {
        bed -> bed[1]
    }.collect()
    LONG_READ_BEDS_TO_GTF(paftools_bed_ch)

    gtf_ch = MERGE_SCALLOP_GTFS.out.merged_gtf
            .concat((MERGE_STRINGTIE_GTFS.out.merged_gtf)
            .concat(LONG_READ_BEDS_TO_GTF.out.gtf))
            .collect()

    CHECK_TRANSCRIPTOMIC_OUTPUT(gtf_ch)

    emit:
    scallop_gtf  = MERGE_SCALLOP_GTFS.out.merged_gtf
    stringtie_gtf = MERGE_STRINGTIE_GTFS.out.merged_gtf
    long_reads_gtf = LONG_READ_BEDS_TO_GTF.out.gtf    

}