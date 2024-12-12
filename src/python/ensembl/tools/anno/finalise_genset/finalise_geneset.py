# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Set of functions used to collect the results of the subpipelines and generate the final output
This logic will be revised in the nextflow pipeline
"""
__all__ = ["run_finalise_geneset"]

import argparse
import logging
import logging.config
import multiprocessing
from pathlib import Path
import re
import subprocess
import shutil
from typing import List, Dict, Optional

from ensembl.tools.anno.utils._utils import (
    get_seq_region_length,
    fasta_to_dict,
    create_dir,
    check_gtf_content,
    read_gtf_genes,
    run_command,
)
from ensembl.tools.anno.finalise_genset.cpc2 import run_cpc2
from ensembl.tools.anno.finalise_genset.rnasamba import run_rnasamba
from ensembl.tools.anno.finalise_genset.diamond import run_diamond


logger = logging.getLogger(__name__)


def run_finalise_geneset(#pylint:disable=too-many-branches
    main_output_dir: Path,
    genome_file: Path,
    seq_region_names: List[str],
    diamond_validation_db: Path,
    validation_type: str = "relaxed",
    num_threads: int = 1,
    cpc2_bin: Path = Path(
        "/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/test_cpc2.sif"
    ),
    rnasamba_bin: Path = Path(
        "/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/rnasamba_latest.sif"
    ),
    rnasamba_weights: Path = Path(
        "/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/rnasamba_data/full_length_weights.hdf5"  # pylint:disable=line-too-long
    ),
    diamond_bin: Path = Path("diamond"),
    min_single_exon_pep_length: int = 100,
    min_multi_exon_pep_length: int = 75,
    min_single_source_probability: float = 0.8,
    min_single_exon_probability: float = 0.9,
):
    """Collect results

    :param main_output_dir: Path for the output dir
    :type main_output_dir: Path
    :param genome_file: Path for genome fasta file
    :type genome_file: Path
    :param seq_region_names: list of seq region names
    :type seq_region_names: List[str]
    :param diamond_validation_db: Path for diamond validation db
    :type diamond_validation_db: Path
    :param validation_type: type of validation ("relaxed" or "moderate").
    :type validation_type: str, default "relaxed"
    :param num_threads: Num of threads. Defaults to 1.
    :type num_threads: int, default 1
    :param cpc2_bin: CPC2 software path
    :type cpc2_bin: Path
    :param rnasamba_bin: RNASamba software path
    :type rnasamba_bin: Path
    :param rnasamba_weights: RNASamba weights path
    :type rnasamba_weights: Path
    :param diamond_bin: Path to the Diamond software.
    :type diamond_bin: Path
    :param min_single_exon_pep_length: Minimum peptide length for single exon.
    :type min_single_exon_pep_length: int = 100,
    :param min_multi_exon_pep_length: Minimum peptide length for multiple exons.
    :type min_multi_exon_pep_length: int, default 75
    :param min_single_source_probability: Minimum probability for single source.
    :type min_single_source_probability: float, default 0.8
    :param min_single_exon_probability: Minimum average probability for single exon.
    :type min_single_exon_probability: float, default 0.9

    """

    final_annotation_dir = create_dir(main_output_dir, "annotation_output")
    region_annotation_dir = create_dir(final_annotation_dir, "initial_region_gtfs")
    final_region_annotation_dir = create_dir(final_annotation_dir, "final_region_gtfs")
    utr_region_annotation_dir = create_dir(final_annotation_dir, "utr_region_gtfs")
    validation_dir = create_dir(final_annotation_dir, "cds_validation")
    seq_region_lengths = get_seq_region_length(genome_file, 0)
    output_file = final_annotation_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Final gtf file exists, skipping analysis")
            return

    protein_annotation_raw = main_output_dir / "genblast_output" / "annotation.gtf"
    minimap2_annotation_raw = main_output_dir / "minimap2_output" / "annotation.gtf"
    stringtie_annotation_raw = main_output_dir / "stringtie_output" / "annotation.gtf"
    scallop_annotation_raw = main_output_dir / "scallop_output" / "annotation.gtf"
    busco_annotation_raw = main_output_dir / "busco_output" / "annotation.gtf"

    transcript_selector_script = (
        Path(__file__).parents[5]
        / "perl"
        / "ensembl"
        / "tools"
        / "anno"
        / "support_scripts_perl"
        / "select_best_transcripts.pl"
    )
    finalise_geneset_script = (
        Path(__file__).parents[5]
        / "perl"
        / "ensembl"
        / "tools"
        / "anno"
        / "support_scripts_perl"
        / "finalise_geneset.pl"
    )
    clean_geneset_script = (
        Path(__file__).parents[5]
        / "perl"
        / "ensembl"
        / "tools"
        / "anno"
        / "support_scripts_perl"
        / "clean_geneset.pl"
    )
    clean_utrs_script = (
        Path(__file__).parents[5]
        / "perl"
        / "ensembl"
        / "tools"
        / "anno"
        / "support_scripts_perl"
        / "clean_utrs_and_lncRNAs.pl"
    )
    gtf_to_seq_script = (
        Path(__file__).parents[5]
        / "perl"
        / "ensembl"
        / "tools"
        / "anno"
        / "support_scripts_perl"
        / "gtf_to_seq.pl"
    )

    transcriptomic_annotation_raw = main_output_dir / "transcriptomic_raw.gtf"

    with transcriptomic_annotation_raw.open("w") as file_out:
        for transcriptomic_file in [
            minimap2_annotation_raw,
            scallop_annotation_raw,
            stringtie_annotation_raw,
        ]:
            if not transcriptomic_file.exists():
                logger.info("No annotation.gtf file found in %s, skipping", transcriptomic_file)
                continue

            with transcriptomic_file.open("r") as file_in:
                for line in file_in:
                    file_out.write(line.rstrip())

    # Copy the raw files into the annotation dir
    copy_raw_files(busco_annotation_raw, final_annotation_dir / "busco_raw.gtf")
    copy_raw_files(protein_annotation_raw, final_annotation_dir / "protein_raw.gtf")

    # select best transcript
    generic_select_cmd = [
        "perl",
        str(transcript_selector_script),
        "-genome_file",
        str(genome_file),
    ]
    pool = multiprocessing.Pool(int(num_threads))  # pylint: disable=consider-using-with
    for seq_region_name in seq_region_names:
        # The selection script needs different params depending on
        # whether the seqs are from transcriptomic data or not
        region_details = f"{seq_region_name}.rs1.re{seq_region_lengths[seq_region_name]}"
        transcriptomic_region_gtf_path = region_annotation_dir / f"{region_details}.trans.gtf"
        busco_region_gtf_path = region_annotation_dir / f"{region_details}.busco.gtf"
        protein_region_gtf_path = region_annotation_dir / f"{region_details}.protein.gtf"
        if transcriptomic_annotation_raw:
            cmd = generic_select_cmd.copy()
            cmd.extend(
                [
                    "-region_details",
                    region_details,
                    "-input_gtf_file",
                    str(transcriptomic_annotation_raw),
                    "-output_gtf_file",
                    str(transcriptomic_region_gtf_path),
                    "-cds_search",
                    "-final_biotype",
                    "transcriptomic",
                ]
            )
            pool.apply_async(run_command, args=(cmd,))
        if busco_annotation_raw:
            cmd = generic_select_cmd.copy()
            cmd.extend(
                [
                    "-region_details",
                    region_details,
                    "-input_gtf_file",
                    str(busco_annotation_raw),
                    "-output_gtf_file",
                    str(busco_region_gtf_path),
                    "-all_cds_exons",
                    "-final_biotype",
                    "busco",
                ]
            )
            pool.apply_async(run_command, args=(cmd,))
        if protein_annotation_raw:
            cmd = generic_select_cmd.copy()
            cmd.extend(
                [
                    "-region_details",
                    region_details,
                    "-input_gtf_file",
                    str(protein_annotation_raw),
                    "-output_gtf_file",
                    str(protein_region_gtf_path),
                    "-clean_transcripts",
                    "-all_cds_exons",
                    "-final_biotype",
                    "protein",
                ]
            )
            pool.apply_async(run_command, args=(cmd,))

    pool.close()
    pool.join()

    # At this point we will have the region files for each seq region and we can merge them
    merge_finalise_output_files(
        final_annotation_dir,
        region_annotation_dir,
        ".trans.gtf",
        "transcriptomic",
    )
    merge_finalise_output_files(final_annotation_dir, region_annotation_dir, ".busco.gtf", "busco")
    merge_finalise_output_files(final_annotation_dir, region_annotation_dir, ".protein.gtf", "protein")

    # Create a single GTF file with all the selected transcripts
    # now that they have proper ids
    fully_merged_gtf_path = final_annotation_dir / "all_selected_transcripts.gtf"
    # Collect all *_sel.gtf files
    gtf_files = list(final_annotation_dir.glob("*_sel.gtf"))

    # Ensure there are files to merge
    if not gtf_files:
        raise FileNotFoundError(  # pylint:disable=raising-format-tuple
            "No *_sel.gtf files found in %s", str(final_annotation_dir)
        )

    try:
        with fully_merged_gtf_path.open("w+") as fully_merged_gtf_out:
            # at thei stage we shouldn't have header duplicated so cat should be robust
            merge_gtf_cmd = ["cat"] + [str(gtf_file) for gtf_file in gtf_files]
            subprocess.run(merge_gtf_cmd, stdout=fully_merged_gtf_out, check=True)
            logger.info("Merged GTF files into %s", fully_merged_gtf_path)
    except Exception as e:  # pylint:disable=raising-format-tuple,broad-exception-caught
        print("An error occurred while merging GTF files: %s", e)

    # Now collapse the gene set
    generic_finalise_cmd = [
        "perl",
        str(finalise_geneset_script),
        "-genome_file",
        str(genome_file),
    ]

    pool = multiprocessing.Pool(int(num_threads))  # pylint:disable=consider-using-with
    for seq_region_name in seq_region_names:
        region_details = f"{seq_region_name}.rs1.re{seq_region_lengths[seq_region_name]}"
        final_region_gtf_path = final_region_annotation_dir / f"{region_details}.final.gtf"
        cmd = generic_finalise_cmd.copy()
        cmd.extend(
            [
                "-region_details",
                region_details,
                "-input_gtf_file",
                str(fully_merged_gtf_path),
                "-output_gtf_file",
                str(final_region_gtf_path),
            ]
        )
        pool.apply_async(run_command, args=(cmd,))

    pool.close()
    pool.join()

    merge_finalise_output_files(
        final_annotation_dir,
        final_region_annotation_dir,
        ".final.gtf",
        "prevalidation",
    )
    merged_gtf_file = final_annotation_dir / "prevalidation_sel.gtf"
    merged_cdna_file = final_annotation_dir / "prevalidation_sel.cdna.fa"
    merged_amino_acid_file = final_annotation_dir / "prevalidation_sel.prot.fa"
    # validation
    updated_gtf_lines = validate_coding_transcripts(
        merged_cdna_file,
        merged_amino_acid_file,
        validation_dir,
        validation_type,
        diamond_validation_db,
        merged_gtf_file,
        num_threads,
        cpc2_bin,
        rnasamba_bin,
        rnasamba_weights,
        diamond_bin,
        min_single_exon_pep_length,
        min_multi_exon_pep_length,
        min_single_source_probability,
        min_single_exon_probability,
    )

    postvalidation_gtf_file = final_annotation_dir / "postvalidation.gtf"
    with postvalidation_gtf_file.open("w+") as file_out:
        for line in updated_gtf_lines:
            file_out.write(line)
    cleaned_initial_gtf_file = final_annotation_dir / "cleaned_pre_utr.gtf"
    cleaned_utr_gtf_file = final_annotation_dir / "annotation.gtf"

    logger.info("Cleaning initial set")
    cleaning_cmd = [
        "perl",
        str(clean_geneset_script),
        "-genome_file",
        str(genome_file),
        "-gtf_file",
        str(postvalidation_gtf_file),
        "-output_gtf_file",
        str(cleaned_initial_gtf_file),
    ]
    logger.info("Running cleaning command: %s", " ".join(map(str, cleaning_cmd)))

    run_command(cleaning_cmd)

    # Clean UTRs
    generic_clean_utrs_cmd = [
        "perl",
        str(clean_utrs_script),
        "-genome_file",
        str(genome_file),
        "-input_gtf_file",
        str(cleaned_initial_gtf_file),
    ]
    pool = multiprocessing.Pool(int(num_threads))  # pylint:disable=consider-using-with
    for seq_region_name in seq_region_names:
        region_details = f"{seq_region_name}.rs1.re{seq_region_lengths[seq_region_name]}"
        utr_region_gtf_path = utr_region_annotation_dir / f"{region_details}.utr.gtf"

        cmd = generic_clean_utrs_cmd.copy()
        cmd.extend(
            [
                "-region_details",
                region_details,
                "-output_gtf_file",
                str(utr_region_gtf_path),
            ]
        )
        pool.apply_async(run_command, args=(cmd,))
    pool.close()
    pool.join()

    merge_finalise_output_files(
        final_annotation_dir,
        utr_region_annotation_dir,
        ".utr.gtf",
        "annotation",
    )
    cmd = [
        "mv",
        str(final_annotation_dir / "annotation_sel.gtf"),
        str(cleaned_utr_gtf_file),
    ]

    run_command(cmd)

    logger.info("Dumping transcript and translation sequences")
    dumping_cmd = [
        "perl",
        str(gtf_to_seq_script),
        "-genome_file",
        str(genome_file),
        "-gtf_file",
        str(cleaned_utr_gtf_file),
    ]

    logger.info("Running dumping command: %s", " ".join(map(str, dumping_cmd)))

    run_command(dumping_cmd)

    logger.info("Finished creating gene set")


def validate_coding_transcripts(
    cdna_file: Path,
    amino_acid_file: Path,
    validation_dir: Path,
    validation_type: str,
    diamond_validation_db: Path,
    gtf_file: Path,
    num_threads: int,
    cpc2_bin: Path = Path(
        "/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/test_cpc2.sif"
    ),
    rnasamba_bin: Path = Path(
        "/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/rnasamba_latest.sif"
    ),
    rnasamba_weights: Path = Path(
        "/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/rnasamba_data/full_length_weights.hdf5"  # pylint:disable=line-too-long
    ),
    diamond_bin: Path = Path("diamond"),
    min_single_exon_pep_length: int = 100,
    min_multi_exon_pep_length: int = 75,
    min_single_source_probability: float = 0.8,
    min_single_exon_probability: float = 0.9,
):
    """Run validation on protein coding transcript
    using CPC2, RNASamba, Diamond

    Args:
        cdna_file (Path): Input cdna file
        amino_acid_file (Path): Input protein file
        validation_dir (Path): Validation output directory
        validation_type (str): _description_
        diamond_validation_db (Path): _description_
        gtf_file (Path): _description_
        num_threads (int): number of threads

    Returns:
        _type_: _description_
    """

    logger.info("Running CDS validation with RNAsamba and CPC2")
    # rnasamba_weights = config["rnasamba"]["weights"]
    rnasamba_results = run_rnasamba(validation_dir, cdna_file, rnasamba_weights, rnasamba_bin)
    cpc2_results = run_cpc2(validation_dir, cdna_file, cpc2_bin)
    diamond_results = run_diamond(
        validation_dir,
        amino_acid_file,
        diamond_validation_db,
        num_threads,
        diamond_bin,
    )
    combined_results = combine_results(rnasamba_results, cpc2_results, diamond_results)
    logger.info("read gtf genes")
    parsed_gtf_genes = read_gtf_genes(gtf_file)
    updated_gtf_lines = update_gtf_genes(
        parsed_gtf_genes,
        combined_results,
        validation_type,
        min_single_exon_pep_length,
        min_multi_exon_pep_length,
        min_single_source_probability,
        min_single_exon_probability,
    )

    return updated_gtf_lines


def combine_results(
    rnasamba_results: List[List[str]],
    cpc2_results: List[List[str]],
    diamond_results: Optional[List[List[str]]] = None,
) -> Dict[str, List[str]]:
    """
    Combine results from RNA-Samba, CPC2, and DIAMOND analyses into a single dictionary.

    This function integrates coding potential predictions from RNA-Samba,
    CPC2, and optional similarity results from DIAMOND for a set of transcripts.

    Args:
        rnasamba_results (List[List[Union[str, float]]]):
            Results from RNA-Samba, where each sublist contains:
                - transcript_id (str): Transcript identifier.
                - coding_probability (float): Probability of being coding.
                - coding_potential (float): Coding potential value.

        cpc2_results (List[List[Union[str, float, int]]]):
            Results from CPC2, where each sublist contains:
                - transcript_id (str): Transcript identifier.
                - coding_probability (float): Probability of being coding.
                - coding_potential (float): Coding potential value.
                - transcript_length (int): Length of the transcript.
                - peptide_length (int): Length of the predicted peptide.

        diamond_results (Optional[List[List[Union[str, float]]]], optional):
            Results from DIAMOND (default is None). Each sublist contains:
                - transcript_id (str): Transcript identifier.
                - e_value (float): E-value for sequence similarity.

    Returns:
        Dict[str, List[Union[float, int]]]:
            A dictionary mapping `transcript_id` to a list of combined values:
            RNA-Samba predictions, CPC2 predictions, and optional DIAMOND results.
    """
    transcript_ids: Dict[str, List[str]] = {}

    # Process RNA-Samba results
    for result in rnasamba_results:
        transcript_id, coding_probability, coding_potential = result[:3]
        if transcript_id not in transcript_ids:
            transcript_ids[transcript_id] = [coding_probability, coding_potential]

    # Process CPC2 results
    for result in cpc2_results:
        transcript_id, coding_probability, coding_potential, transcript_length, peptide_length = result[:5]
        if transcript_id in transcript_ids:
            transcript_ids[transcript_id].extend(
                [coding_probability, coding_potential, transcript_length, peptide_length]
            )

    # Process DIAMOND results (if provided)
    if diamond_results:
        for result in diamond_results:
            transcript_id, e_value = result[:2]
            if transcript_id in transcript_ids:
                transcript_ids[transcript_id].append(e_value)

    return transcript_ids


def update_gtf_genes(#pylint:disable=too-many-branches
    parsed_gtf_genes: Dict[str, Dict[str, Dict[str, List[str]]]],
    combined_results: Dict[str, List[str]],
    validation_type: str,
    min_single_exon_pep_length: int = 100,
    min_multi_exon_pep_length: int = 75,
    min_single_source_probability: float = 0.8,
    min_single_exon_probability: float = 0.9,
):
    """Update GTF genes based on validation criteria.

    Args:
        parsed_gtf_genes (Dict): Parsed GTF genes data.
        combined_results (Dict): Combined validation results.
        validation_type (str): Validation strictness ("relaxed" or "moderate").
        min_single_exon_pep_length (int): Minimum peptide length for single exon.
        min_multi_exon_pep_length (int): Minimum peptide length for multiple exons.
        min_single_source_probability (float): Minimum probability for single source.
        min_single_exon_probability (float): Minimum average probability for single exon.

    Returns:
        List[str]: Updated GTF lines.
    """
    output_lines = []

    for gene_id, transcripts in parsed_gtf_genes.items():#pylint:disable=unused-variable
        for transcript_id, transcript_data in transcripts.items():
            transcript_line = "".join(transcript_data["transcript"])
            single_cds_exon_transcript = 0

            # Extract translation coordinates
            translation_match = re.search(r'; translation_coords "([^"]+)";', str(transcript_line))
            if translation_match:
                translation_coords_list = translation_match.group(1).split(":")
                # Determine if it's a single-exon CDS
                if translation_coords_list[0] == translation_coords_list[3]:
                    single_cds_exon_transcript = 1

            exon_lines = transcript_data["exons"]
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

            # Calculate coding probabilities
            avg_coding_probability = (rnasamba_coding_probability + cpc2_coding_probability) / 2
            max_coding_probability = max(rnasamba_coding_probability, cpc2_coding_probability)

            match = re.search(r'; biotype "([^"]+)";', str(transcript_line))
            if match:
                biotype = match.group(1)
                if biotype == "busco" or biotype == "protein":#pylint:disable=consider-using-in
                    transcript_line = update_biotype(str(transcript_line), biotype, "protein_coding")
                    output_lines.append(transcript_line)
                    output_lines.extend(exon_lines)
                    continue

            # Note that the below looks at validating things
            # under different levels of strictness
            # There are a few different continue statements,
            # where transcripts will be skipped resulting
            # in a smaller post validation file. It mainly
            # removes single coding exon genes with no real
            # support or for multi-exon lncRNAs that are less than 200bp long
            if single_cds_exon_transcript == 1 and validation_type == "relaxed":
                if diamond_e_value is not None:
                    transcript_line = update_biotype(transcript_line, biotype, "protein_coding")

                elif (
                    rnasamba_coding_potential == "coding"
                    and cpc2_coding_potential == "coding"
                    and peptide_length >= min_single_exon_pep_length
                ):
                    transcript_line = update_biotype(transcript_line, biotype, "protein_coding")
                elif (
                    (rnasamba_coding_potential == "coding" or cpc2_coding_potential == "coding")
                    and peptide_length >= min_single_exon_pep_length
                    and max_coding_probability >= min_single_source_probability
                ):
                    transcript_line = update_biotype(transcript_line, biotype, "protein_coding")
                else:
                    continue
            elif single_cds_exon_transcript == 1 and validation_type == "moderate":
                if diamond_e_value is not None and peptide_length >= min_single_exon_pep_length:
                    transcript_line = update_biotype(transcript_line, biotype, "protein_coding")
                elif (
                    (rnasamba_coding_potential == "coding" and cpc2_coding_potential == "coding")
                    and peptide_length >= min_single_exon_pep_length
                    and avg_coding_probability >= min_single_exon_probability
                ):
                    transcript_line = update_biotype(transcript_line, biotype, "protein_coding")
                else:
                    continue
            else:
                if diamond_e_value is not None:
                    transcript_line = update_biotype(transcript_line, biotype, "protein_coding")
                elif (
                    rnasamba_coding_potential == "coding"
                    and cpc2_coding_potential == "coding"
                    and peptide_length >= min_multi_exon_pep_length
                ):
                    transcript_line = update_biotype(transcript_line, biotype, "protein_coding")
                elif (
                    (rnasamba_coding_potential == "coding" or cpc2_coding_potential == "coding")
                    and peptide_length >= min_multi_exon_pep_length
                    and max_coding_probability >= min_single_source_probability
                ):
                    transcript_line = update_biotype(transcript_line, biotype, "protein_coding")
                elif transcript_length >= 200:
                    transcript_line = update_biotype(transcript_line, biotype, "lncRNA")
                    transcript_line = re.sub(' translation_coords "[^"]+";', "", transcript_line)
                else:
                    continue

            output_lines.append(transcript_line)
            output_lines.extend(exon_lines)

    return output_lines


def update_biotype(transcript_line: str, current_biotype: str, new_biotype: str) -> str:
    """Update the biotype in a transcript line.

    Args:
        transcript_line (str): The line containing the transcript information.
        current_biotype (str): The biotype to be replaced.
        new_biotype (str): The new biotype to set.

    Returns:
        str: The updated transcript line with the new biotype.
    """
    return re.sub(f'; biotype "{current_biotype}";', f'; biotype "{new_biotype}";', transcript_line)


def merge_finalise_output_files(
    final_annotation_dir: Path, region_annotation_dir: Path, extension: str, id_label: str
) -> None:
    """Merge the resulting gene models from all the seq regions for the specified extension

    Args:
        final_annotation_dir (Path): Output directory
        region_annotation_dir (Path): Input directory
        extension (str): File extension
        id_label (str): Type of data or process to merge
    """
    gtf_files = list(region_annotation_dir.glob(f"*{extension}"))
    merged_gtf_file = final_annotation_dir / f"{id_label}_sel.gtf"
    merged_cdna_file = final_annotation_dir / f"{id_label}_sel.cdna.fa"
    merged_amino_acid_file = final_annotation_dir / f"{id_label}_sel.prot.fa"

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

    with merged_gtf_file.open("w") as gtf_out, merged_cdna_file.open(
        "w"
    ) as cdna_out, merged_amino_acid_file.open("w") as amino_acid_out:
        for gtf_file in gtf_files:
            logger.info("GTF file: %s", gtf_file)
            cdna_seq_index = {}
            amino_acid_seq_index = {}
            cdna_file = gtf_file.with_suffix(".cdna")
            amino_acid_file = gtf_file.with_suffix(".prot")

            # Check that the files exist and are not empty
            if (not Path(cdna_file).is_file() or Path(cdna_file).stat().st_size == 0) and (
                Path(amino_acid_file).is_file() or Path(amino_acid_file).stat().st_size == 0
            ):
                logger.error(
                    "DNA file %s or amino acid file %s \
                    does not exist or is empty.",
                    cdna_file,
                    amino_acid_file,
                )
            else:
                # Open and process the files
                with open(cdna_file) as cdna_in, open(amino_acid_file) as amino_acid_in:
                    cdna_seq_index = fasta_to_dict(cdna_in.readlines())
                    amino_acid_seq_index = fasta_to_dict(amino_acid_in.readlines())

            current_gene_id = ""
            with gtf_file.open() as gtf_in:
                for line in gtf_in:
                    if re.search(r"^#", line):
                        # if line.startswith("#"):
                        continue
                    eles = line.split("\t")
                    if len(eles) != 9:
                        continue

                    match = re.search(r'gene_id "([^"]+)".+transcript_id "([^"]+)"', line)
                    gene_id = ""
                    transcript_id = ""
                    if match and eles[2] == "transcript":
                        transcript_id_counter += 1
                        gene_id = match.group(1)
                        transcript_id = match.group(2)

                    if not current_gene_id:
                        gene_id_counter += 1
                        current_gene_id = gene_id

                    if gene_id != current_gene_id:
                        gene_id_counter += 1
                        current_gene_id = gene_id
                    new_gene_id = f"{id_label}_{gene_id_counter}"
                    new_transcript_id = f"{id_label}_{transcript_id_counter}"

                    line = re.sub(f'gene_id "{gene_id}"', f'gene_id "{new_gene_id}"', line)
                    line = re.sub(
                        f'transcript_id "{transcript_id}"', f'transcript_id "{new_transcript_id}"', line
                    )
                    gtf_out.write(line)
                    if eles[2] == "transcript":
                        new_header = f">{new_transcript_id}\n"
                        cdna_out.write(new_header + cdna_seq_index[transcript_id])
                        if transcript_id in amino_acid_seq_index:
                            amino_acid_out.write(new_header + amino_acid_seq_index[transcript_id])


def copy_raw_files(raw_file: Path, destination_file: Path) -> None:
    """Copy file in a new destination

    Args:
        raw_file (Path): File to copy
        destination_file (Path): New destination
    """
    if raw_file.exists():
        try:
            shutil.copy(raw_file, destination_file)
            logger.info("Copied %s to %s", raw_file, destination_file)
        except Exception as e:  # pylint:disable=broad-exception-caught
            logger.error("Failed to copy %s to %s: %s", raw_file, destination_file, e)
    else:
        logger.info("No file found at %s, skipping", raw_file)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Genblast arguments")

    parser.add_argument("--main_output_dir", required=True, help="Output directory path")
    parser.add_argument("--genome_file", required=True, help="Path for the genome file")
    parser.add_argument("--seq_region_names", type=List[str], help="List of seq region names")
    parser.add_argument("--diamond_validation_db", required=True, help="Diamond validation db file path")
    parser.add_argument(
        "--validation_type", type=int, choices=["relaxed", "moderate"], help="Type of validation"
    )
    parser.add_argument(
        "--diamond_bin",
        help="DIAMOND executable path",
    )
    parser.add_argument(
        "--cpc2_bin",
        help="CPC2 executable path",
    )
    parser.add_argument(
        "--rnasamba_bin",
        help="RNAsamba executable path",
    )
    parser.add_argument(
        "--rnasamba_weights",
        help="Rnasamba weights path",
    )
    parser.add_argument("--num_threads", type=int, default=1, help="Number of threads")
    parser.add_argument(
        "--min_single_exon_pep_length", type=int, default=100, help="Minimum peptide length for single exon."
    )
    parser.add_argument(
        "--min_multi_exon_pep_length", type=int, default=75, help="Minimum peptide length for multiple exons."
    )
    parser.add_argument(
        "--min_single_source_probability",
        type=float,
        default=0.8,
        help="Minimum probability for single source.",
    )
    parser.add_argument(
        "--min_single_exon_probability",
        type=float,
        default=0.9,
        help="Minimum average probability for single exon.",
    )
    return parser.parse_args()


def main():
    """Finalisation's entry-point."""
    args = parse_args()

    log_file_path = create_dir(args.output_dir, "log") / "finalisation.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"

    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )

    run_finalise_geneset(
        Path(args.main_output_dir),
        Path(args.genome_file),
        args.seq_region_names,
        Path(args.diamond_validation_db),
        args.validation_type,
        args.num_threads,
        Path(args.cpc2_bin),
        Path(args.rnasamba_bin),
        Path(args.rnasamba_weights),
        Path(args.diamond_bin),
        args.min_single_exon_pep_length,
        args.min_multi_exon_pep_length,
        args.min_single_source_probability,
        args.min_single_exon_probability,
    )


if __name__ == "__main__":
    main()
