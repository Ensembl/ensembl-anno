
include { TRIMGALORE } from '../modules/trimgalore.nf'
include { STAR } from '../modules/star.nf'
include { SCALLOP } from '../modules/scallop.nf'
include { STRINGTIE } from '../modules/stringtie.nf'
include { MINIMAP2 } from '../modules/minimap2.nf'

workflow TRANSCRIPTOMICS_ANNOTATION {
    take:
    short_reads  // channel: [ val(meta), [path(read1), path(read2)] ]
    long_reads   // channel: [ val(meta), path(input_file) ]
    genome_fasta

    main:

    // ToDo split up input files (in channels)
    // ToDo add logic to handle if short and or long reads not available

    // Short reads pipeline
    // Trim and align reads in a splice aware manner
    if (params.trim_reads){
        TRIMGALORE(short_reads)
        STAR(TRIMGALORE.out.trimmed_reads)
    } else {
        STAR(short_reads)
    }
    

    // Two alternative models for building transcript models
    SCALLOP(STAR.out.bam)
    STRINGTIE(STAR.out.bam)


    // Long reads pipeline
    // Use minimap2 to align long reads
    MINIMAP2(long_reads)

    // Combine outputs into single gtf
    // TODO: Does this work?
    transcriptomics_gtf = MINIMAP2.out.gtf
        .concat(SCALLOP.out.gtf.concat(STRINGTIE.out.gtf))
    
    // TODO: Does this work?
    versions = TRIMGALORE.out.versions.concat(
                STAR.out.versions.concat(
                    SCALLOP.out.versions.concat(
                        STRINGTIE.out.versions.concat(
                            MINIMAP2.out.versions
                )
            )
        )
    )

    emit:
    transcriptomics_gtf  = transcriptomics_gtf      // channel: [ path(gtf) ]
    versions = versions                             // channel: [ path(versions) ]
}