include { CONVERT_TO_BLASTMASK } from '../modules/convert2blastmask.nf'
include { MAKE_BLAST_DB } from '../modules/make_blast_db.nf'
include { SPLIT_PROTEIN_FILE } from '../modules/split_protein_file.nf'
include { GENBLAST } from '../modules/genblast.nf'
include { MAKE_GTF as MAKE_GENBLAST_GTF} from '../modules/make_gtf.nf'
include { COMBINE_SLICED_GTFS as COMBINE_UNIPROT_GENBLAST_GTFS} from '../modules/combine_sliced_gtfs.nf'
include { COMBINE_SLICED_GTFS as COMBINE_ORTHODB_GENBLAST_GTFS} from '../modules/combine_sliced_gtfs.nf'

//include { MINIPROT } from '../modules/miniprot.nf'

workflow PROTEINS {
    take:
    masked_fasta
    proteins

    main:

    CONVERT_TO_BLASTMASK(masked_fasta)
    MAKE_BLAST_DB(masked_fasta, CONVERT_TO_BLASTMASK.out.asnb)

    SPLIT_PROTEIN_FILE(proteins)
    //SPLIT_PROTEIN_FILE.out.sliced_proteins.view()
    split_proteins_with_ids = SPLIT_PROTEIN_FILE.out.sliced_proteins.flatMap { val, files ->
        files.collect { file -> tuple(val, file.baseName, file) }
    }
    split_proteins_with_ids.view()
    
    GENBLAST(MAKE_BLAST_DB.out.fasta_db, split_proteins_with_ids)
    
    
    genblast_ch = channel.of(tuple('genblast', file('optional1'), file('optional2'))).collect()
    MAKE_GENBLAST_GTF(genblast_ch, GENBLAST.out.gff)



    genblast_branched_gtf_ch = MAKE_GENBLAST_GTF.out.gtf.branch{
        it ->
            orthodb: it[0].startsWith('orthodb')
                return it[1]

            uniprot: it[0].startsWith('uniprot')
                return it[1]

    }
    genblast_orthodb_ch = channel.of(tuple('genblast_orthodb', file('optional1'), file('optional2'))).collect()
    genblast_uniprot_ch = channel.of(tuple('genblast_uniprot', file('optional1'), file('optional2'))).collect()


    COMBINE_ORTHODB_GENBLAST_GTFS(genblast_orthodb_ch, genblast_branched_gtf_ch.orthodb.collect())
    COMBINE_UNIPROT_GENBLAST_GTFS(genblast_uniprot_ch, genblast_branched_gtf_ch.uniprot.collect())

    genblast_annot_gtf = COMBINE_ORTHODB_GENBLAST_GTFS.out.gtf.concat(
        COMBINE_UNIPROT_GENBLAST_GTFS.out.gtf
    )

    emit:
    genblast_gtf  = genblast_annot_gtf
    //miniprot_gtf =  <>


}