# Legacy finalisation module
import glob
import multiprocessing
import os
from pathlib import Path
import re
import subprocess
import logging
import logging.config
import json
from typing import Any, Dict, Union

from src.python.ensembl.tools.anno.utils._utils import (
    check_file,
    create_dir,
    check_gtf_content,
    get_seq_region_lengths,
    split_protein_file,
)

def load_json(path: Union[str, Path]) -> Dict[str, Any]:
    path = Path(path)
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)
   
# File check helpers:
def file_ok(path: Union[str, Path], min_size: int = 1) -> bool:
    path = Path(path)
    return path.exists() and path.is_file() and path.stat().st_size >= min_size


def gtf_ok(path: Union[str, Path]) -> bool:
    if not file_ok(path):
        return False
    try:
        return check_gtf_content(path, "transcript") > 0
    except Exception:
        return False


def skip_if_exists(desc: str, path: Union[str, Path], check_fn=file_ok) -> bool:
    if check_fn(path):
        logging.info("%s exists (%s), skipping", desc, path)
        return True
    return False

# Use same config path layout as monolithic script
config = load_json(Path(os.environ["ENSCODE"]) / "ensembl-anno" / "conf" /"config.json")

def run_finalise_geneset(
    main_script_dir,
    main_output_dir,
    genome_file,
    seq_region_names,
    validation_type,
    diamond_validation_db,
    num_threads,
):

    if validation_type is None:
        logging.info("Setting validation type to relaxed")
    else:
        logging.info("Setting validation type to " + validation_type)

    final_annotation_dir = create_dir(main_output_dir, "annotation_output")
    region_annotation_dir = create_dir(final_annotation_dir, "initial_region_gtfs")
    final_region_annotation_dir = create_dir(
        final_annotation_dir, "final_region_gtfs"
    )
    utr_region_annotation_dir = create_dir(final_annotation_dir, "utr_region_gtfs")
    validation_dir = create_dir(final_annotation_dir, "cds_validation")
    seq_region_lengths = get_seq_region_lengths(genome_file, 0)

    logging.info("Check if gtf file already exists")
    output_file = os.path.join(final_annotation_dir, "annotation.gtf")
    if os.path.exists(output_file):
        logging.info("final_annotation_dir exists")
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logging.info("Final gtf file exists")
            return
    else:
        logging.info("No gtf file, go on with the analysis")
    # This used to be a list of output dirs and a loop which was neat,
    # I'm coverting to a list of conditions as
    # it's more straightforward with the renaming
    # and having to merge scallop and stringtie
    protein_annotation_raw = os.path.join(
        main_output_dir, "uniprot_output", "annotation.gtf"
    )
    minimap2_annotation_raw = os.path.join(
        main_output_dir, "minimap2_output", "annotation.gtf"
    )
    stringtie_annotation_raw = os.path.join(
        main_output_dir, "stringtie_output", "annotation.gtf"
    )
    scallop_annotation_raw = os.path.join(
        main_output_dir, "scallop_output", "annotation.gtf"
    )
    busco_annotation_raw = os.path.join(main_output_dir, "orthodb_output", "annotation.gtf")

    transcript_selector_script = os.path.join(
        main_script_dir, "support_scripts_perl", "select_best_transcripts.pl"
    )
    finalise_geneset_script = os.path.join(
        main_script_dir, "support_scripts_perl", "finalise_geneset.pl"
    )
    clean_geneset_script = os.path.join(
        main_script_dir, "support_scripts_perl", "clean_geneset.pl"
    )
    clean_utrs_script = os.path.join(
        main_script_dir, "support_scripts_perl", "clean_utrs_and_lncRNAs.pl"
    )
    gtf_to_seq_script = os.path.join(
        main_script_dir, "support_scripts_perl", "gtf_to_seq.pl"
    )

    transcriptomic_annotation_raw = os.path.join(
        final_annotation_dir, "transcriptomic_raw.gtf"
    )
    if not skip_if_exists("Transcriptomic raw GTF", transcriptomic_annotation_raw):
        with open(transcriptomic_annotation_raw, "w+") as file_out:
            for transcriptomic_file in [
                minimap2_annotation_raw,
                scallop_annotation_raw,
                stringtie_annotation_raw,
            ]:
                if not os.path.exists(transcriptomic_file):
                    logging.info("Missing %s, skipping", transcriptomic_file)
                    continue

                with open(transcriptomic_file) as file_in:
                    for line in file_in:
                        file_out.write(line)

    # Copy the raw files into the annotation dir, this is not needed
    # as such, but collecting them in one place and relabelling is
    # helpful for a user
    if os.path.exists(busco_annotation_raw):
        subprocess.run(
            [
                "cp",
                busco_annotation_raw,
                os.path.join(final_annotation_dir, "busco_raw.gtf"),
            ]
        )

    if os.path.exists(protein_annotation_raw):
        subprocess.run(
            [
                "cp",
                protein_annotation_raw,
                os.path.join(final_annotation_dir, "protein_raw.gtf"),
            ]
        )

    gtf_files = ["transcriptomic_raw.gtf", "protein_raw.gtf", "busco_raw.gtf"]
    generic_select_cmd = [
        "perl",
        transcript_selector_script,
        "-genome_file",
        genome_file,
    ]
    pool = multiprocessing.Pool(int(num_threads))
    for seq_region_name in seq_region_names:
        # The selection script needs different params depending on
        # whether the seqs are from transcriptomic data or not
        region_details = (
            seq_region_name + ".rs1" + ".re" + str(seq_region_lengths[seq_region_name])
        )
        transcriptomic_region_gtf_path = os.path.join(
            region_annotation_dir, (region_details + ".trans.gtf")
        )
        busco_region_gtf_path = os.path.join(
            region_annotation_dir, (region_details + ".busco.gtf")
        )
        protein_region_gtf_path = os.path.join(
            region_annotation_dir, (region_details + ".protein.gtf")
        )

        # Transcriptomic data
        if os.path.exists(transcriptomic_annotation_raw):
            if not gtf_ok(transcriptomic_region_gtf_path):
                logging.info("Finalising transcriptomic data for: " + seq_region_name)
                transcriptomic_annotation_select = re.sub(
                    "_raw.gtf", "_sel.gtf", transcriptomic_annotation_raw
                )
                cmd = generic_select_cmd.copy()
                cmd.extend(
                    [
                        "-region_details",
                        region_details,
                        "-input_gtf_file",
                        transcriptomic_annotation_raw,
                        "-output_gtf_file",
                        transcriptomic_region_gtf_path,
                        "-cds_search",
                        "-final_biotype",
                        "transcriptomic",
                    ]
                )
                pool.apply_async(multiprocess_finalise_geneset, args=(cmd,))
            else:
                logging.info("Transcriptomic region GTF exists, skipping: %s",
                     transcriptomic_region_gtf_path)
                
        # Protein data
        if os.path.exists(busco_annotation_raw):
            if not gtf_ok(busco_region_gtf_path):
                logging.info("Finalising OrthoDB data for: " + seq_region_name)
                busco_annotation_select = re.sub("_raw.gtf", "_sel.gtf", busco_annotation_raw)
                cmd = generic_select_cmd.copy()
                cmd.extend(
                    [
                        "-region_details",
                        region_details,
                        "-input_gtf_file",
                        busco_annotation_raw,
                        "-output_gtf_file",
                        busco_region_gtf_path,
                        "-all_cds_exons",
                        "-final_biotype",
                        "busco",
                    ]
                )
                pool.apply_async(multiprocess_finalise_geneset, args=(cmd,))
            else:
                logging.info("OrthoDB/BUSCO region GTF exists, skipping: %s",
                     busco_region_gtf_path)

        if os.path.exists(protein_annotation_raw):
            if not gtf_ok(protein_region_gtf_path):
                logging.info("Finalising Uniprot data for: " + seq_region_name)
                protein_annotation_select = re.sub(
                    "_raw.gtf", "_sel.gtf", protein_annotation_raw
                )
                cmd = generic_select_cmd.copy()
                cmd.extend(
                    [
                        "-region_details",
                        region_details,
                        "-input_gtf_file",
                        protein_annotation_raw,
                        "-output_gtf_file",
                        protein_region_gtf_path,
                        "-clean_transcripts",
                        "-all_cds_exons",
                        "-final_biotype",
                        "protein",
                    ]
                )
                pool.apply_async(multiprocess_finalise_geneset, args=(cmd,))
            else:
                logging.info("UniProt region GTF exists, skipping: %s",
                            protein_region_gtf_path)

    pool.close()
    pool.join()

    # At this point we will have the region files for all the,
    merge_finalise_output_files(
        final_annotation_dir,
        region_annotation_dir,
        ".trans.gtf",
        "transcriptomic",
    )
    merge_finalise_output_files(
        final_annotation_dir, region_annotation_dir, ".busco.gtf", "busco"
    )
    merge_finalise_output_files(
        final_annotation_dir, region_annotation_dir, ".protein.gtf", "protein"
    )

    # Create a single GTF file with all the selected transcripts
    # now that they have proper ids
    fully_merged_gtf_path = os.path.join(
        final_annotation_dir, "all_selected_transcripts.gtf"
    )
    fully_merged_gtf_out = open(fully_merged_gtf_path, "w+")

    merge_gtf_cmd = ["cat"]
    merge_gtf_cmd.extend(glob.glob(str(final_annotation_dir) + "/*_sel.gtf"))
    subprocess.run(merge_gtf_cmd, stdout=fully_merged_gtf_out)
    fully_merged_gtf_out.close()

    # Now collapse the gene set
    generic_finalise_cmd = [
        "perl",
        finalise_geneset_script,
        "-genome_file",
        genome_file,
    ]

    pool = multiprocessing.Pool(int(num_threads))
    for seq_region_name in seq_region_names:
        region_details = (
            seq_region_name + ".rs1" + ".re" + str(seq_region_lengths[seq_region_name])
        )
        final_region_gtf_path = os.path.join(
            final_region_annotation_dir, (region_details + ".final.gtf")
        )
        if gtf_ok(final_region_gtf_path):
            logging.info("Final region GTF exists, skipping: %s", final_region_gtf_path)
        else:
            cmd = generic_finalise_cmd.copy()
            cmd.extend(
                [
                    "-region_details",
                    region_details,
                    "-input_gtf_file",
                    fully_merged_gtf_path,
                    "-output_gtf_file",
                    final_region_gtf_path,
                ]
            )
            pool.apply_async(multiprocess_finalise_geneset, args=(cmd,))

    pool.close()
    pool.join()

    merge_finalise_output_files(
        final_annotation_dir,
        final_region_annotation_dir,
        ".final.gtf",
        "prevalidation",
    )
    merged_gtf_file = os.path.join(final_annotation_dir, ("prevalidation_sel.gtf"))
    merged_cdna_file = os.path.join(final_annotation_dir, ("prevalidation_sel.cdna.fa"))
    merged_amino_acid_file = os.path.join(
        final_annotation_dir, ("prevalidation_sel.prot.fa")
    )

    postvalidation_gtf_file = os.path.join(final_annotation_dir, ("postvalidation.gtf"))

    # Check if already run
    if gtf_ok(postvalidation_gtf_file):
        logging.info("Post-validation GTF exists, skipping validation")
    
    else:
        updated_gtf_lines = validate_coding_transcripts(merged_cdna_file,
        merged_amino_acid_file,
        validation_dir,
        validation_type,
        diamond_validation_db,
        merged_gtf_file,
        num_threads,)
        with open(postvalidation_gtf_file, "w+") as file_out:
            file_out.writelines(updated_gtf_lines)

    cleaned_initial_gtf_file = os.path.join(final_annotation_dir, ("cleaned_pre_utr.gtf"))
    cleaned_utr_gtf_file = os.path.join(final_annotation_dir, ("annotation.gtf"))

    logging.info("Cleaning initial set")
    cleaning_cmd = [
        "perl",
        clean_geneset_script,
        "-genome_file",
        genome_file,
        "-gtf_file",
        postvalidation_gtf_file,
        "-output_gtf_file",
        cleaned_initial_gtf_file,
    ]
    logging.info(" ".join(cleaning_cmd))
    if not gtf_ok(cleaned_utr_gtf_file):
        subprocess.run(cleaning_cmd)
    else:
        logging.info("Cleaned annotation.gtf exists, skipping cleaning")

    # Clean UTRs
    generic_clean_utrs_cmd = [
        "perl",
        clean_utrs_script,
        "-genome_file",
        genome_file,
        "-input_gtf_file",
        cleaned_initial_gtf_file,
    ]
    pool = multiprocessing.Pool(int(num_threads))
    for seq_region_name in seq_region_names:
        region_details = (
            seq_region_name + ".rs1" + ".re" + str(seq_region_lengths[seq_region_name])
        )
        utr_region_gtf_path = os.path.join(
            utr_region_annotation_dir, (region_details + ".utr.gtf")
        )

        cmd = generic_clean_utrs_cmd.copy()
        cmd.extend(
            [
                "-region_details",
                region_details,
                "-input_gtf_file",
                cleaned_initial_gtf_file,
                "-output_gtf_file",
                utr_region_gtf_path,
            ]
        )
        pool.apply_async(multiprocess_generic, args=(cmd,))
    pool.close()
    pool.join()

    merge_finalise_output_files(
        final_annotation_dir,
        utr_region_annotation_dir,
        ".utr.gtf",
        "annotation",
    )
    subprocess.run(
        [
            "mv",
            os.path.join(final_annotation_dir, "annotation_sel.gtf"),
            cleaned_utr_gtf_file,
        ]
    )

    logging.info("Dumping transcript and translation sequences")
    dumping_cmd = [
        "perl",
        gtf_to_seq_script,
        "-genome_file",
        genome_file,
        "-gtf_file",
        cleaned_utr_gtf_file,
    ]
    logging.info(" ".join(dumping_cmd))
    subprocess.run(dumping_cmd)

    logging.info("Finished creating gene set")


def validate_coding_transcripts(
    cdna_file,
    amino_acid_file,
    validation_dir,
    validation_type,
    diamond_validation_db,
    gtf_file,
    num_threads,
):

    logging.info("Running CDS validation with RNAsamba")
    rnasamba_weights = config["rnasamba"]["weights"]
    rnasamba_output_path = Path(validation_dir) / "rnasamba.tsv.txt"

    if not rnasamba_output_path.exists():
        rnasamba_volume = str(validation_dir) + "/:/app:rw"
        rnasamba_cmd = [
            "singularity",
            "exec",
            "--bind",
            rnasamba_volume,
            config["rnasamba"]["software"],
            "rnasamba",
            "classify",
            str(rnasamba_output_path),
            cdna_file,
            rnasamba_weights,
        ]
        logging.info(" ".join(rnasamba_cmd))
        subprocess.run(rnasamba_cmd)

        check_file(Path(rnasamba_output_path))
        logging.info("RNAsamba run sucessfully")
    else:
        logging.info("Found RNAsamba output file. Will not run again.")

    cpc2_output_file_path = Path(validation_dir) / "cpc2.tsv.txt"
    cpc2_output_path = Path(validation_dir) / "cpc2.tsv"
    
    if not cpc2_output_file_path.exists():
        logging.info("Running CDS validation with CPC2")
        cpc2_volume = str(validation_dir) + "/:/app:rw"
        cpc2_cmd = [
            "singularity",
            "exec",
            "--bind",
            cpc2_volume,
            config["cpc2"]["software"],
            "python3",
            "/CPC2_standalone-1.0.1/bin/CPC2.py",
            "-i",
            cdna_file,
            "--ORF",
            "-o",
            cpc2_output_path,
        ]
        logging.info(" ".join(map(str, cpc2_cmd)))
        subprocess.run(cpc2_cmd)

        check_file(Path(cpc2_output_file_path))
        logging.info("CPC2 run sucessfully")
    
    else:
        logging.info("Found CPC2 output file. Will not run again.")

    logging.info("Running Diamond validation")
    diamond_results = None
    if diamond_validation_db is not None:
        diamond_output_dir = create_dir(validation_dir, "diamond_output")
        diamond_validation(
            diamond_validation_db,
            amino_acid_file,
            str(diamond_output_dir),
            num_threads,
        )
        diamond_results = read_diamond_results(diamond_output_dir)

    logging.info("Reading validation results")
    rnasamba_results = read_rnasamba_results(rnasamba_output_path)
    cpc2_results = read_cpc2_results(cpc2_output_file_path)
    combined_results = combine_results(rnasamba_results, cpc2_results, diamond_results)
    logging.info("Reading gtf genes")
    parsed_gtf_genes = read_gtf_genes(gtf_file)
    updated_gtf_lines = update_gtf_genes(
        parsed_gtf_genes, combined_results, validation_type
    )

    return updated_gtf_lines


def diamond_validation(
    diamond_validation_db, amino_acid_file, diamond_output_dir, num_threads
):
    logging.info("Starting Diamond validation")
    logging.info("Diamond DB: %s", diamond_validation_db)
    logging.info("Amino acid file: %s", amino_acid_file)
    logging.info("Diamond output dir: %s", diamond_output_dir)
    logging.info("Threads: %s", num_threads)

    batched_protein_files = split_protein_file(amino_acid_file, Path(diamond_output_dir), 100)

    logging.info(
            "Split protein file into %d batches", len(batched_protein_files)
        )

    if not batched_protein_files:
        logging.error("No batched protein files created")
        return
    
    pool = multiprocessing.Pool(int(num_threads))
    for batched_protein_file in batched_protein_files:
        logging.info(
            "Submitting Diamond job for batch: %s", batched_protein_file
        )
        pool.apply_async(
            multiprocess_diamond,
            args=(
                Path(batched_protein_file),
                Path(diamond_output_dir),
                diamond_validation_db,
            ),
        )
    pool.close()
    pool.join()

    logging.info("Diamond validation finished")


def multiprocess_diamond(
    batched_protein_file: Path,
    diamond_output_dir: Path,
    diamond_validation_db: str,
):
    try:
        logging.info("Diamond worker started for %s", batched_protein_file)

        if not batched_protein_file.exists():
            logging.error("Batch file does not exist: %s", batched_protein_file)
            return

        diamond_output_file = diamond_output_dir / (
            batched_protein_file.name + ".dmdout"
        )

        logging.info(
            "Diamond output file will be: %s", diamond_output_file
        )

        diamond_cmd = [
            "diamond",
            "blastp",
            "--query",
            str(batched_protein_file),
            "--db",
            diamond_validation_db,
            "--out",
            str(diamond_output_file),
            "--outfmt",
            "6",
            "--threads",
            "1",
        ]

        logging.info("Running Diamond command:")
        logging.info(" ".join(diamond_cmd))

        result = subprocess.run(
            diamond_cmd,
            capture_output=True,
            text=True,
        )

        logging.info(
            "Diamond return code for %s: %s",
            batched_protein_file,
            result.returncode,
        )

        if result.stdout:
            logging.debug("Diamond stdout:\n%s", result.stdout)

        if result.stderr:
            logging.warning("Diamond stderr:\n%s", result.stderr)

        if result.returncode != 0:
            logging.error(
                "Diamond failed for %s", batched_protein_file
            )
            return

        if not diamond_output_file.exists():
            logging.error(
                "Diamond produced no output file for %s",
                batched_protein_file,
            )
            return

        if diamond_output_file.stat().st_size == 0:
            logging.error(
                "Diamond output file is empty: %s",
                diamond_output_file,
            )
            return

        logging.info(
            "Diamond completed successfully for %s",
            batched_protein_file,
        )

    except Exception as e:
        logging.exception(
            "Unhandled exception in Diamond worker for %s", batched_protein_file
        )


def update_gtf_genes(parsed_gtf_genes, combined_results, validation_type):

    output_lines = []

    for gene_id in parsed_gtf_genes.keys():
        transcript_ids = parsed_gtf_genes[gene_id].keys()
        for transcript_id in transcript_ids:
            transcript_line = parsed_gtf_genes[gene_id][transcript_id]["transcript"]
            single_cds_exon_transcript = 0
            translation_match = re.search(
                r'; translation_coords "([^"]+)";', transcript_line
            )
            if translation_match:
                translation_coords = translation_match.group(1)
                translation_coords_list = translation_coords.split(":")
                # If the start exon coords of both exons are the same,
                # then it's the same exon and thus a single exon cds
                if translation_coords_list[0] == translation_coords_list[3]:
                    single_cds_exon_transcript = 1

            exon_lines = parsed_gtf_genes[gene_id][transcript_id]["exons"]
            validation_results = combined_results[transcript_id]
            rnasamba_coding_probability = float(validation_results[0])
            rnasamba_coding_potential = validation_results[1]
            cpc2_coding_probability = float(validation_results[2])
            cpc2_coding_potential = validation_results[3]
            transcript_length = int(validation_results[4])
            peptide_length = int(validation_results[5])
            diamond_e_value = None
            if len(validation_results) == 7:
                diamond_e_value = validation_results[6]

            avg_coding_probability = (
                rnasamba_coding_probability + cpc2_coding_probability
            ) / 2
            max_coding_probability = max(
                rnasamba_coding_probability, cpc2_coding_probability
            )

            match = re.search(r'; biotype "([^"]+)";', transcript_line)
            biotype = match.group(1)
            if biotype == "busco" or biotype == "protein":
                transcript_line = re.sub(
                    '; biotype "' + biotype + '";',
                    '; biotype "protein_coding";',
                    transcript_line,
                )
                output_lines.append(transcript_line)
                output_lines.extend(exon_lines)
                continue

            min_single_exon_pep_length = 100
            min_multi_exon_pep_length = 75
            min_single_source_probability = 0.8
            min_single_exon_probability = 0.9

            # Note that the below looks at validating things
            # under different levels of strictness
            # There are a few different continue statements,
            # where transcripts will be skipped resulting
            # in a smaller post validation file. It mainly
            # removes single coding exon genes with no real
            # support or for multi-exon lncRNAs that are less than 200bp long
            if single_cds_exon_transcript == 1 and validation_type == "relaxed":
                if diamond_e_value is not None:
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    rnasamba_coding_potential == "coding"
                    and cpc2_coding_potential == "coding"
                    and peptide_length >= min_single_exon_pep_length
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    (
                        rnasamba_coding_potential == "coding"
                        or cpc2_coding_potential == "coding"
                    )
                    and peptide_length >= min_single_exon_pep_length
                    and max_coding_probability >= min_single_source_probability
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                else:
                    continue
            elif single_cds_exon_transcript == 1 and validation_type == "moderate":
                if (
                    diamond_e_value is not None
                    and peptide_length >= min_single_exon_pep_length
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    (
                        rnasamba_coding_potential == "coding"
                        and cpc2_coding_potential == "coding"
                    )
                    and peptide_length >= min_single_exon_pep_length
                    and avg_coding_probability >= min_single_exon_probability
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                else:
                    continue
            else:
                if diamond_e_value is not None:
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    rnasamba_coding_potential == "coding"
                    and cpc2_coding_potential == "coding"
                    and peptide_length >= min_multi_exon_pep_length
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    (
                        rnasamba_coding_potential == "coding"
                        or cpc2_coding_potential == "coding"
                    )
                    and peptide_length >= min_multi_exon_pep_length
                    and max_coding_probability >= min_single_source_probability
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif transcript_length >= 200:
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "lncRNA";',
                        transcript_line,
                    )
                    transcript_line = re.sub(
                        ' translation_coords "[^"]+";', "", transcript_line
                    )
                else:
                    continue

            output_lines.append(transcript_line)
            output_lines.extend(exon_lines)

    return output_lines


def read_rnasamba_results(file_path):

    results = []

    file_in = open(file_path)
    line = file_in.readline()
    while line:
        line = line.rstrip()
        match = re.search(r"^sequence_name", line)
        if match:
            line = file_in.readline()
            continue

        eles = line.split("\t")
        if not len(eles) == 3:
            line = file_in.readline()
            continue

        transcript_id = eles[0]
        coding_probability = eles[1]
        coding_potential = eles[2]
        results.append([transcript_id, coding_probability, coding_potential])
        line = file_in.readline()
    file_in.close()

    return results


def read_cpc2_results(file_path):

    results = []

    file_in = open(file_path)
    line = file_in.readline()
    while line:
        line = line.rstrip()
        match = re.search(r"^#ID", line)
        if match:
            line = file_in.readline()
            continue

        eles = line.split("\t")
        if not len(eles) == 9:
            line = file_in.readline()
            continue

        transcript_id = eles[0]
        transcript_length = eles[1]
        peptide_length = eles[2]
        coding_probability = eles[7]
        coding_potential = eles[8]
        results.append(
            [
                transcript_id,
                coding_probability,
                coding_potential,
                transcript_length,
                peptide_length,
            ]
        )
        line = file_in.readline()
    file_in.close()

    return results


def read_diamond_results(diamond_output_dir):

    results = []
    diamond_files = glob.glob(str(diamond_output_dir) + "/*.dmdout")
    for file_path in diamond_files:
        file_in = open(file_path)
        line = file_in.readline()
        while line:
            line = line.rstrip()

            eles = line.split("\t")
            if not len(eles) == 12:
                line = file_in.readline()
                continue

            transcript_id = eles[0]
            e_value = eles[10]
            results.append([transcript_id, e_value])
            line = file_in.readline()
    file_in.close()

    return results


def combine_results(rnasamba_results, cpc2_results, diamond_results):

    transcript_ids = {}

    for result in rnasamba_results:
        transcript_id = result[0]
        coding_probability = result[1]
        coding_potential = result[2]

        if transcript_id not in transcript_ids:
            transcript_ids[transcript_id] = [
                coding_probability,
                coding_potential,
            ]

    for result in cpc2_results:
        transcript_id = result[0]
        coding_probability = result[1]
        coding_potential = result[2]
        transcript_length = result[3]
        peptide_length = result[4]
        transcript_ids[transcript_id].extend(
            [
                coding_probability,
                coding_potential,
                transcript_length,
                peptide_length,
            ]
        )

    if diamond_results is not None:
        for result in diamond_results:
            transcript_id = result[0]
            e_value = result[1]
            if transcript_id in transcript_ids:
                transcript_ids[transcript_id].extend([e_value])

    return transcript_ids


def read_gtf_genes(gtf_file):

    gtf_genes = {}
    gtf_in = open(gtf_file)
    line = gtf_in.readline()
    while line:
        eles = line.split("\t")
        if not len(eles) == 9:
            line = gtf_in.readline()
            continue

        match = re.search(r'gene_id "([^"]+)".+transcript_id "([^"]+)"', line)

        if not match:
            line = gtf_in.readline()
            continue

        gene_id = match.group(1)
        transcript_id = match.group(2)
        feature_type = eles[2]
        if gene_id not in gtf_genes:
            gtf_genes[gene_id] = {}
        if feature_type == "transcript":
            gtf_genes[gene_id][transcript_id] = {}
            gtf_genes[gene_id][transcript_id]["transcript"] = line
            gtf_genes[gene_id][transcript_id]["exons"] = []
        elif feature_type == "exon":
            gtf_genes[gene_id][transcript_id]["exons"].append(line)
        line = gtf_in.readline()
    gtf_in.close()

    return gtf_genes


def merge_finalise_output_files(
    final_annotation_dir, region_annotation_dir, extension, id_label
):

    gtf_files = glob.glob(str(region_annotation_dir) + "/*" + extension)

    merged_gtf_file = os.path.join(final_annotation_dir, (id_label + "_sel.gtf"))
    merged_cdna_file = os.path.join(final_annotation_dir, (id_label + "_sel.cdna.fa"))
    merged_amino_acid_file = os.path.join(
        final_annotation_dir, (id_label + "_sel.prot.fa")
    )

    # The below is not great, it's a bit messy because there might be
    # some cases where there aren't translations. So it's not as
    # straightforward as reading the records across all three files
    # in parallel. The solution is to just load the seqs into
    # memory and index them on the current header, which should
    # correspond to a transcript/gene id in the GTF. When writing the
    # results into the single merged files the ids will be updated to
    # be unique and consistent across the header,  three file types

    gene_id_counter = 0
    transcript_id_counter = 0
    gtf_out = open(merged_gtf_file, "w+")
    cdna_out = open(merged_cdna_file, "w+")
    amino_acid_out = open(merged_amino_acid_file, "w+")
    for gtf_file in gtf_files:
        logging.info("GTF file: " + gtf_file)
        cdna_seq_index = {}
        amino_acid_seq_index = {}
        cdna_file = gtf_file + ".cdna"
        amino_acid_file = gtf_file + ".prot"
        cdna_in = open(cdna_file)
        amino_acid_in = open(amino_acid_file)
        cdna_seq_index = fasta_to_dict(cdna_in.readlines())
        amino_acid_seq_index = fasta_to_dict(amino_acid_in.readlines())
        cdna_in.close()
        amino_acid_in.close()

        current_gene_id = ""
        gtf_in = open(gtf_file)
        line = gtf_in.readline()
        while line:
            if re.search(r"^#", line):
                line = gtf_in.readline()
                continue

            eles = line.split("\t")
            if not len(eles) == 9:
                line = gtf_in.readline()
                continue

            match = re.search(r'gene_id "([^"]+)".+transcript_id "([^"]+)"', line)
            if match and eles[2] == "transcript":
                transcript_id_counter += 1

            gene_id = match.group(1)
            transcript_id = match.group(2)

            if not current_gene_id:
                gene_id_counter += 1
                current_gene_id = gene_id

            if not gene_id == current_gene_id:
                gene_id_counter += 1
                current_gene_id = gene_id

            new_gene_id = id_label + "_" + str(gene_id_counter)
            new_transcript_id = id_label + "_" + str(transcript_id_counter)
            line = re.sub(
                'gene_id "' + gene_id + '"',
                ('gene_id "' + new_gene_id + '"'),
                line,
            )
            line = re.sub(
                'transcript_id "' + transcript_id + '"',
                ('transcript_id "' + new_transcript_id + '"'),
                line,
            )
            gtf_out.write(line)
            line = gtf_in.readline()

            if eles[2] == "transcript":
                new_header = ">" + new_transcript_id + "\n"
                cdna_out.write(new_header + cdna_seq_index[transcript_id])

                if transcript_id in amino_acid_seq_index:
                    amino_acid_out.write(new_header + amino_acid_seq_index[transcript_id])

    gtf_out.close()
    cdna_out.close()
    amino_acid_out.close()


def fasta_to_dict(fasta_list):

    index = {}
    it = iter(fasta_list)
    for header in it:
        match = re.search(r">(.+)\n$", header)
        header = match.group(1)
        seq = next(it)
        index[header] = seq
    return index


# def merge_gtf_files(file_paths,id_label):

#  gtf_file_path = os.path.join(output_dir,'annotation.gtf')
#  gtf_out = open(gtf_file_path,'w+')
#  for gtf_file_path in gtf_files:
#    gtf_file_name = os.path.basename(gtf_file_path)
#    match = re.search(r'\.rs(\d+)\.re(\d+)\.',gtf_file_name)
#    start_offset = int(match.group(1))
#    gtf_in = open(gtf_file_path,'r')
#    line = gtf_in.readline()
#    while line:
#      values = line.split("\t")
#      if len(values) == 9 and (values[2] in feature_types):
#        values[3] = str(int(values[3]) + (start_offset - 1))
#        values[4] = str(int(values[4]) + (start_offset - 1))
#        gtf_out.write("\t".join(values))
#        line = gtf_in.readline()
#    gtf_in.close()
#  gtf_out.close()


def multiprocess_generic(cmd):
    print(" ".join(cmd))
    subprocess.run(cmd)


def multiprocess_finalise_geneset(cmd):

    print(" ".join(cmd))
    subprocess.run(cmd)
