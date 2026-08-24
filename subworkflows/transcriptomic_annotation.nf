
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

workflow TRANSCRIPTOMICS_ANNOTATION {
    take:
    short_reads  // channel: [ val(meta), [path(read1), path(read2)] ]
    long_reads   // channel: [ val(meta), path(input_file) ]
    genome_fasta // channel: [fasta_path]

    main:

    // make STAR index
    // It is neccessary to run collect() to convert a queue channel to a value channel
    // Otherwise the index is consumed after the first fastq is aligned and the rest are skipped
    STAR_INDEX(genome_fasta)
    star_index_ch = STAR_INDEX.out.index.collect()

    // ToDo split up input files (in channels)
    // ToDo add logic to handle if short and or long reads not available

    // Short reads pipeline
    // Trim and align reads in a splice aware manner
    if (params.trim_reads){
        TRIMGALORE(short_reads)
        STAR(TRIMGALORE.out.trimmed_reads, star_index_ch)
    } else {
        STAR(short_reads, star_index_ch)
    }

    SAMTOOLS(STAR.out.sam)

    // Two alternative models for building transcript models
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
    // Use minimap2 to align long reads
    // It is neccessary to run collect() to convert a queue channel to a value channel
    // Otherwise the index is consumed after the first fastq is aligned and the rest are skipped
    MINIMAP2_INDEX(genome_fasta)
    minimap2_index_ch = MINIMAP2_INDEX.out.index.collect()

    MINIMAP2(long_reads, minimap2_index_ch)
    PAFTOOLS(MINIMAP2.out.sam)

    // Combine outputs into single gtf
    // TODO: Does this work? Not sure if we actually want this
    // transcriptomics_gtf = MINIMAP2.out.gtf
    //     .concat(SCALLOP.out.gtf.concat(STRINGTIE.out.gtf))

    emit:
    // ToDo update this once python package is in place
    scallop_gtf  = MERGE_SCALLOP_GTFS.out.merged_gtf
    stringtie_gtf = MERGE_STRINGTIE_GTFS.out.merged_gtf
    long_reads_bed = PAFTOOLS.out.bed     

}