import argparse
from pathlib import Path
from os import PathLike
import re
from typing import Union, Dict, List, Any, cast
import numpy as np
from numpy.typing import NDArray


def create_red_gtf(repeat_coords_file: Path, output_file: Path):
    """
    Create Red gtf file from masked genome file

    Args:
        repeat_coords_file: Coordinates for repeats.
        output_file : GTF file with the final results.
    """
    with (
        open(repeat_coords_file, "r", encoding="utf8") as red_in,
        open(output_file, "w+", encoding="utf8") as red_out,
    ):
        for repeat_id, line in enumerate(red_in, start=1):
            result_match = re.search(r"^\>(.+)\:(\d+)\-(\d+)", line)
            if result_match:
                region_name = result_match.group(1)
                # Note that Red is 0-based, so add 1
                start = int(result_match.group(2)) + 1
                end = int(result_match.group(3)) + 1
                gtf_line = (
                    f"{region_name}\tRed\trepeat\t{start}\t"
                    f'{end}\t.\t+\t.\trepeat_id "{repeat_id}";\n'  # pylint:disable=line-too-long
                )
                red_out.write(gtf_line)


def create_dust_gtf(
    input_file: Path,
    output_gtf: Path,
    region_name: str,
) -> None:
    """
    All the genomic slices are collected in a single gtf output
    Args:
        input_file : GTF file with final results.
        output_gtf : GTF file with the results per region.
        region_name :Coordinates of genomic slice.
    """
    with (
        open(input_file, "r", encoding="utf8") as dust_in,
        open(output_gtf, "w+", encoding="utf8") as dust_out,
    ):
        repeat_count = 1
        for line in dust_in:
            result_match = re.search(r"(\d+)\ - (\d+)", line)
            if result_match:
                start = int(result_match.group(1)) + 1
                end = int(result_match.group(2)) + 1
                gtf_line = (
                    f"{region_name}\tDust\trepeat\t{start}\t"
                    f'{end}\t.\t+\t.\trepeat_id "{repeat_count}";\n'  # pylint:disable=line-too-long
                )
                dust_out.write(gtf_line)
                repeat_count += 1


# Function to find the repeat class based on the mappings
def get_repeat_type(repeat_type: str) -> str:
    """Get the repeat type based on the provided repeat_type string.

    Args:
        repeat_type (str): The repeat type string to match against the mappings.

    Returns:
        str: The corresponding repeat type description if a match is found,
        otherwise "Unknown".
    """
    mappings = {
        r"^Low_Comp": "Low complexity regions",
        r"^LINE": "Type I Transposons/LINE",
        r"^SINE": "Type I Transposons/SINE",
        r"^DNA": "Type II Transposons",
        r"^LTR": "LTRs",
        r"^Other": "Other repeats",
        r"^Satelli": "Satellite repeats",
        r"^Simple": "Simple repeats",
        r"^Tandem": "Tandem repeats",
        r"^TRF": "Tandem repeats",
        r"^Waterman": "Waterman",
        r"^Recon": "Recon",
        r"^Tet_repeat": "Tetraodon repeats",
        r"^MaskRegion": "Mask region",
        r"^dust": "Dust",
        r"^Unknown": "Unknown",
        r"RNA$": "RNA repeats",
    }
    for pattern, description in mappings.items():
        if re.match(pattern, repeat_type):
            return description
    return "Unknown"  # Default if no match is found


def create_repeatmasker_gtf(  # pylint: disable=too-many-locals
    input_file: Path,
    output_gtf: Path,
    region_name: str,
) -> None:
    """

    All the genomic slices are collected in a single gtf output with the following format:
    SW    perc perc perc query    position in query matching repeat       position in repeat
    score div. del. ins. sequence begin end (left)  repeat   class/family begin end  (left)  ID
    Args:
        input_file : GTF file with final results.
        output_gtf_path : GTF file with results per region.
        region_name : Coordinates of genomic slice.
    """
    with (
        open(input_file, "r", encoding="utf8") as repeatmasker_in,
        open(output_gtf, "w+", encoding="utf8") as repeatmasker_out,
    ):
        repeat_count = 1
        for line in repeatmasker_in:
            result_match = re.search(r"^\s*\d+\s+", line)
            if result_match:
                results = line.split()
                if results[-1] == "*":
                    results.pop()
                if len(results) != 15:
                    continue
                score = results[0]
                start = results[5]
                end = results[6]
                strand = results[8]
                repeat_name = results[9]
                repeat_class = results[10]
                repeat_type = get_repeat_type(results[10])
                if strand == "+":
                    repeat_start = results[11]
                    repeat_end = results[12]
                else:
                    repeat_start = results[13]
                    repeat_end = results[12]
                    strand = "-"
                gtf_line = (
                    f"{region_name}\tRepeatMasker\trepeat\t{start}\t{end}\t.\t"
                    f"{strand}\t.\trepeat_id {repeat_count}; "
                    f'repeat_name "{repeat_name}"; repeat_class "{repeat_class}"; '
                    f'repeat_type "{repeat_type}"; repeat_start "{repeat_start}"; '
                    f'repeat_end "{repeat_end}"; score "{score}";\n'
                )
                repeatmasker_out.write(gtf_line)
                repeat_count += 1


def create_trf_gtf(  # pylint:disable=too-many-locals, too-many-branches
    input_file: Path,
    output_gtf: Path,
    region_name: str,
) -> None:
    """

    TRF output format:
    cols 1+2:  Indices of the repeat relative to the start of the sequence
    col 3:     Period size of the repeat
    col 4:     Number of copies aligned with the consensus pattern
    col 5:     Size of consensus pattern (may differ slightly from the period size)
    col 6:     Percent of matches between adjacent copies overall
    col 7:     Percent of indels between adjacent copies overall
    col 8:     Alignment score
    cols 9-12: Percent composition for each of the four nucleotides
    col 13:    Entropy measure based on percent composition
    col 14:    Consensus sequence
    col 15:    Repeat sequence
    Args:
       input_file : GTF file with final results.
       output_gtf : GTF file with results per region.
       region_name : Coordinates of genomic slice.
    """
    with (
        open(input_file, "r", encoding="utf8") as trf_in,
        open(output_gtf, "w+", encoding="utf8") as trf_out,
    ):
        repeat_count = 1
        for line in trf_in:
            result_match = re.search(r"^\d+", line)
            if result_match:
                results = line.split()
                if len(results) != 15:
                    continue
                start = results[0]
                end = results[1]
                period = float(results[2])
                copy_number = float(results[3])
                percent_matches = float(results[5])
                score = float(results[7])
                repeat_consensus = results[13]
                if (  # pylint: disable=too-many-boolean-expressions
                    score < 50 and percent_matches >= 80 and copy_number > 2 and period < 10
                ) or (copy_number >= 2 and percent_matches >= 70 and score >= 50):
                    gtf_line = (
                        f"{region_name}\tTRF\trepeat\t{start}\t{end}\t.\t+\t.\t"
                        f'repeat_id "{repeat_count}"; score "{score}"; '
                        f'repeat_consensus "{repeat_consensus}";\n'
                    )
                    trf_out.write(gtf_line)
                    repeat_count += 1


def create_cpg_gtf(  # pylint:disable=too-many-arguments, too-many-locals, too-many-branches, too-many-positional-arguments
    input_file: Path,
    output_gtf: Path,
    region_name: str,
    cpg_min_length: int = 400,
    cpg_min_gc_content: int = 50,
    cpg_min_oe: float = 0.6,
) -> None:
    """
    Read the fasta file and save the content in gtf format
    All the genomic slices are collected in a single gtf output
    Args:
        input_file : GTF file with final results.
        output_gtf : GTF file with the results per region.
        region_name :Coordinates of genomic slice.
        cpg_dir : Output dir.
        cpg_min_length : Min length of CpG islands
        cpg_min_gc_content : Min GC frequency percentage
        cpg_min_oe :  Min ratio of the observed to expected number of CpG (CpGo/e)
    """
    with (
        open(input_file, "r", encoding="utf8") as cpg_in,
        open(output_gtf, "w+", encoding="utf8") as cpg_out,
    ):
        feature_count = 1
        for line in cpg_in:
            result_match = re.search(r"^" + region_name, line)
            if result_match:
                results = line.split()
                start = int(results[1])
                end = int(results[2])
                length = end - start + 1
                score = float(results[3])
                gc_content = float(results[6])
                oe_score_str = results[7]
                oe_score: Union[float, int]
                if oe_score_str in ("-", "inf"):
                    oe_score = 0
                else:
                    oe_score = float(oe_score_str)
                if (
                    int(length) >= int(cpg_min_length)
                    and gc_content >= int(cpg_min_gc_content)
                    and oe_score >= float(cpg_min_oe)
                ):
                    gtf_line = (
                        f"{region_name}\tCpG\tsimple_feature\t{start}\t"
                        f'{end}\t.\t+\t.\tfeature_id "{feature_count}"; score "{score}";\n'
                    )
                    cpg_out.write(gtf_line)


def create_eponine_gtf(
    input_file: Path,
    output_gtf: Path,
    region_name: str,
) -> None:
    """
    Read the fasta file and save the content in gtf format
    All the genomic slices are collected in a single gtf output
    Args:
        input_file: GTF file with final results.
        output_gtf: GTF file with the results per region.
        region_name: Coordinates of genomic slice.
    """
    with (
        open(input_file, "r", encoding="utf8") as eponine_in,
        open(output_gtf, "w+", encoding="utf8") as eponine_out,
    ):
        feature_count = 1
        for line in eponine_in:
            result_match = re.search(r"^" + region_name, line)
            if result_match:
                results = line.split()
                start = int(results[3])
                end = int(results[4])
                score = float(results[5])
                strand = results[6]
                # There's a one base offset on the reverse strand
                if strand == "-":
                    start -= 1
                    end -= 1

                gtf_line = (
                    f"{region_name}\tEponine\tsimple_feature\t"
                    f"{start}\t{end}\t.\t{strand}\t.\t"
                    f'feature_id "{feature_count}"; score "{score}";\n'
                )
                eponine_out.write(gtf_line)
                feature_count += 1


def create_trnascan_gtf(input_file: Path, output_gtf: Path, region_name: str) -> None:
    """
    Read the fasta file and save the content in gtf format
    All the genomic slices are collected in a single gtf output
    Args:
        output_gtf : GTF file with the results per region.
        filter_file : GTF file with the filtered results per region.
        region_name :Coordinates of genomic slice.

    tRNAscan-SE output format:
    col0: GtRNAdb Gene Symbol - gene ID in corresponding genome
    col1: tRNAscan-SE ID - tRNA ID in tRNAscan-SE prediction results
    col2-3: Locus - Genomic coordinates of predicted gene
    col4: Isotype (from Anticodon) - tRNA isotype determined by anticodon
    col5: Anticodon - anticodon of predicted tRNA gene
    col6-7: Intron boundaries
    col8: General tRNA Model Score - covariance model bit score from tRNAscan-SE results
    col9: Best Isotype Model - best matching (highest scoring) isotype determined
    by isotype-specific covariance model classification
    col10-11-12: Anticodon and Isotype Model Agreement - consistency between anticodon
    from predicted gene sequence and best isotype model
    col13: Features - special gene features that may include gene set categorization,
    number of introns, possible pseudogenes, possible truncation, or base-pair mismatches
    """
    with (
        open(input_file, "r", encoding="utf8") as trna_in,
        open(output_gtf, "w+", encoding="utf8") as trna_out,
    ):
        gene_counter = 1
        for line in trna_in:
            result_match = re.search(r"^" + region_name, line)
            if result_match:
                results = line.split()
                start = int(results[2])
                end = int(results[3])
                strand = "+"
                if start > end:
                    strand = "-"
                    start, end = end, start
                biotype = "tRNA" if re.search(r"high confidence set", line) else "tRNA_pseudogene"
                transcript_string = (
                    f"{region_name}\ttRNAscan\ttranscript\t{start}\t{end}\t.\t"
                    f'{strand}\t.\tgene_id "{gene_counter}"; transcript_id '
                    f'"{gene_counter}"; biotype "{biotype}";\n'
                )
                exon_string = (
                    f"{region_name}\ttRNAscan\texon\t{start}\t{end}\t.\t"
                    f'{strand}\t.\tgene_id "{gene_counter}"; transcript_id '
                    f'"{gene_counter}"; exon_number "1"; biotype "{biotype}";\n'
                )
                trna_out.write(transcript_string)
                trna_out.write(exon_string)
                trna_out.flush()
                gene_counter += 1


def get_rfam_seed_descriptions(rfam_seeds_file: PathLike) -> Dict[str, Dict[str, Any]]:
    """Get Rfam seed description

    Args:
        rfam_seeds_file (PathLike): File of Rfam seeds

    Returns:
        dict: List of Rfam seeds with description,name, type
    """
    descriptions: Dict[str, Dict[str, Any]] = {}
    rfam_seeds = []
    domain = ""
    # NOTE: for some reason the decoder breaks on the seeds file,
    # so I have made this ignore errors
    with open(rfam_seeds_file, encoding="utf-8", errors="ignore") as rfam_seeds_in:
        rfam_seeds = rfam_seeds_in.read().splitlines()

    for seed in rfam_seeds:
        matches = re.findall(r"^\#=GF (AC|DE|ID|TP)\s+(.+)", seed)

        if matches:
            key, value = matches[0]
            if key == "AC":
                domain = value
                descriptions[domain] = {}
            elif key == "DE":
                assert domain is not None, "Domain should not be None at this point."
                descriptions[domain]["description"] = value
            elif key == "ID":
                assert domain is not None, "Domain should not be None at this point."
                descriptions[domain]["name"] = value
            elif key == "TP":
                assert domain is not None, "Domain should not be None at this point."
                descriptions[domain]["type"] = value
    return descriptions


def extract_rfam_metrics(rfam_selected_models: PathLike) -> Dict[str, Dict[str, Any]]:
    """Get name, description, length, max length, threshold of each Rfam model.

    Args:
        rfam_selected_models_file : Path for Rfam models.

    Returns:
        parsed_cm_data: Rfam metrics.
    """
    with open(rfam_selected_models, "r", encoding="utf-8") as rfam_cm_in:
        rfam_models = rfam_cm_in.read().split("//\n")
        parsed_cm_data: Dict[str, Dict[str, Any]] = {}
        for model in rfam_models:
            model_name_match = re.search(r"NAME\s+(\S+)", model)
            match_infernal = re.search(r"INFERNAL", model)
            if model_name_match and match_infernal:
                model_name = model_name_match.group(1)
                parsed_cm_data[model_name] = {}
                parse_regex = {
                    r"^NAME\s+(\S+)": "-name",
                    r"^DESC\s+(\S+)": "-description",
                    r"^CLEN\s+(\d+)": "-length",
                    r"^W\s+(\d+)": "-maxlength",
                    r"^ACC\s+(\S+)": "-accession",
                    r"^GA\s+(\d+)": "-threshold",
                }
                for line in model.split("\n"):
                    for pattern, value_type in parse_regex.items():
                        match = re.search(pattern, line)
                        if match:
                            parsed_cm_data[model_name][value_type] = match.group(1)
                            continue

    return parsed_cm_data


def parse_rfam_tblout(region_tblout: Path, region_name: str) -> List[Dict[str, Any]]:
    """Parse cmsearch output
    col 0 Target Name : This is the name of the target sequence or sequence
    region that matched the query.
    col 2 Query name : This is the name of the query sequence or model that
    was used for the search.
    col 3 Accession : This usually refers to a unique identifier for the
    target sequence.
    col 5 Query Start : The position where the match starts on the query
    sequence.
    col 6 Query End : The position where the match ends on the query
    sequence.
    col 7 Target Start : The position where the match starts on the
    target sequence.
    col 8 Target End : The position where the match ends on the target
    sequence.
    col 9 Strand : Indicates the orientation of the match on the target
    sequence.
            It could be + for the forward strand or - for the reverse strand.
    col 14 Hit Score : The score assigned to this match. Higher scores
    generally indicate better matches.
    col 15 E-value : This is a statistical measure of the number of hits
    one can expect to see when searching a database of a particular size.

    Args:
        region_tblout : Cmsearch output for the region name.
        region_name : Region name.

    Returns:
        Formatted cmsearch output
    """

    with open(region_tblout, "r", encoding="utf-8") as rfam_tbl_in:
        rfam_tbl_data = rfam_tbl_in.read()

    results = []
    for line in rfam_tbl_data.splitlines():
        if not line or line.startswith("#"):
            continue
        # Must start with region_name exactly
        if not re.match(rf"^{re.escape(region_name)}\b", line):
            continue
        hit = line.split()
        if len(hit) < 16:
            continue
        results.append(
            {
                "accession": hit[3],
                "start": hit[7],
                "end": hit[8],
                "strand": 1 if hit[9] == "+" else -1,
                "query_name": hit[2],
                "score": hit[14],
            }
        )
    return results


def remove_rfam_overlap(  # pylint: disable=too-many-locals, too-many-branches
    parsed_tbl_data: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """
    Remove Rfam mdoels overlapping and with a lower score.

    Args:
        parsed_tbl_data : Cmsearch output

    Returns:
        Final Rfam models
    """
    excluded_structures = {}
    chosen_structures: List[Dict[str, Any]] = []
    for structure_x in parsed_tbl_data:
        chosen_structure = structure_x
        structure_x_start = int(structure_x["start"])
        structure_x_end = int(structure_x["end"])
        structure_x_score = float(structure_x["score"])
        structure_x_accession = structure_x["accession"]
        structure_x_string = (
            f"{structure_x_start}:{structure_x_end}:{structure_x_score}:{structure_x_accession}"
        )
        for structure_y in parsed_tbl_data:
            structure_y_start = int(structure_y["start"])
            structure_y_end = int(structure_y["end"])
            structure_y_score = float(structure_y["score"])
            structure_y_accession = structure_y["accession"]
            structure_y_string = (
                f"{structure_y_start}:{structure_y_end}:{structure_y_score}:{structure_y_accession}"
            )
            if structure_y_string in excluded_structures:
                continue
            if structure_x_start <= structure_y_end and structure_x_end >= structure_y_start:
                if structure_x_score < structure_y_score:
                    chosen_structure = structure_y
                    excluded_structures[structure_x_string] = 1
                else:
                    excluded_structures[structure_y_string] = 1
        chosen_structures.append(chosen_structure)
    return chosen_structures


def filter_rfam_results(
    unfiltered_tbl_data: List[Dict[str, Any]], cm_models: Dict[str, Dict[str, Any]]
) -> List[Dict[str, Any]]:
    """Filter Rfam models according to type and a set of thresholds.

    Args:
        unfiltered_tbl_data : Unfiltered Rfam output.
        cm_models : Rfam models.

    Returns:
        filtered_results: List of filtered models
    """
    filtered_results: List[Dict[str, Any]] = []
    thresholds = {
        "LSU_rRNA_eukarya": 1700,
        "SSU_rRNA_eukarya": 1600,
        "5_8S_rRNA": 85,
        "5S_rRNA": 75,
    }
    for structure in unfiltered_tbl_data:
        query = structure["query_name"]
        if query in ["LSU_rRNA_archaea", "LSU_rRNA_bacteria"]:
            threshold = cm_models.get(str(query), {}).get("-threshold")
        else:
            threshold = thresholds.get(str(query), cm_models.get(str(query), {}).get("-threshold"))
        if threshold is not None and float(structure["score"]) >= float(threshold):
            filtered_results.append(structure)
    return filtered_results


def create_cmsearch_gtf(  # pylint: disable=too-many-arguments, too-many-locals, too-many-positional-arguments
    filtered_results: List[Dict[str, Any]],
    cm_models: Dict[str, Dict[str, Any]],
    seed_descriptions: Dict[str, Dict[str, Any]],
    region_name: str,
    output_gtf: Path,
    output_bed: Path,
) -> None:
    """Convert RFam output per single region in gtf format

    Args:
        filtered_results : Filtered Rfam results without overlapping.
        cm_models :  Rfam database.
        seed_descriptions : Rfam seed file.
        region_name : Slice name.
        output_gtf : Rfam output file.
        genome_file : Genome file.
        rfam_dir : Output file.
        rnafold_bin: RNAfold software path.
    """
    if not filtered_results:
        return

    biotype_type_mapping = {
        r"snRNA; snoRNA; scaRNA": "scaRNA",
        r"snRNA; snoRNA": "snoRNA",
        r"snRNA": "snRNA",
        r"rRNA;": "rRNA",
        r"antisense;": "antisense",
        r"antitoxin;": "antitoxin",
        r"ribozyme;": "ribozyme",
    }
    biotype_name_mapping = {
        r"Vault": "Vault_RNA",
        r"Y_RNA": "Y_RNA",
        r"^RNaseP": "RNase_P_RNA",
        r"^RNase_M": "RNase_MRP_RNA",
    }
    with open(output_gtf, "w+", encoding="utf-8") as rfam_gtf_out:
        with open(output_bed, "w+", encoding="utf-8") as rfam_bed_out:
            gene_counter = 1
            for structure in filtered_results:
                query = structure["query_name"]
                accession = structure["accession"]
                if query in cm_models:
                    model = cm_models[query]  # pylint: disable=unused-variable
                    description = seed_descriptions.get(accession, {})
                    rfam_type = description.get("type", "misc_RNA")
                    domain = structure["query_name"]
                    # padding = model["-length"]
                    gtf_strand = structure["strand"]
                    rnafold_strand = structure["strand"]
                    if gtf_strand == 1:
                        start = structure["start"]
                        end = structure["end"]
                        gtf_strand = "+"
                    else:
                        start = structure["end"]
                        end = structure["start"]
                        # score = structure["score"]
                        gtf_strand = "-"
                        rnafold_strand = -1

                    biotype = "misc_RNA"
                    # Flag to track if a match is found
                    match_found = False
                    # Check each pattern and update biotype if a match is found
                    for pattern, mapped_biotype in biotype_type_mapping.items():
                        if re.search(pattern, str(rfam_type)):
                            biotype = mapped_biotype
                            match_found = True
                            break  # Break out of the loop once a match is found
                    # If no match is found, check matches in biotype_name_mapping
                    if not match_found:
                        # Check each pattern and update biotype if a match is found
                        for pattern, mapped_biotype in biotype_name_mapping.items():
                            if re.search(pattern, domain):
                                biotype = mapped_biotype
                                match_found = True
                                break  # Break out of the loop once a match is found

                    transcript_string = (
                        region_name
                        + "\tRfam\ttranscript\t"
                        + str(start)
                        + "\t"
                        + str(end)
                        + "\t.\t"
                        + gtf_strand
                        + "\t.\t"
                        + 'gene_id "'
                        + str(gene_counter)
                        + '"; transcript_id "'
                        + str(gene_counter)
                        + '"; biotype "'
                        + biotype  # pylint: disable=undefined-loop-variable
                        + '";\n'
                    )
                    exon_string = (
                        region_name
                        + "\tRfam\texon\t"
                        + str(start)
                        + "\t"
                        + str(end)
                        + "\t.\t"
                        + gtf_strand
                        + "\t.\t"
                        + 'gene_id "'
                        + str(gene_counter)
                        + '"; transcript_id "'
                        + str(gene_counter)
                        + '"; exon_number "1"; biotype "'
                        + biotype  # pylint: disable=undefined-loop-variable
                        + '";\n'
                    )
                    bed_string = (
                        region_name
                        + "\t"
                        + str(int(start) - 1)
                        + "\t"
                        + str(end)
                        + "\t"
                        + str(gene_counter)
                        + " \t.\t"
                        + gtf_strand
                        + "\n"
                    )

                    rfam_gtf_out.write(transcript_string)
                    rfam_gtf_out.write(exon_string)
                    rfam_bed_out.write(bed_string)
                    gene_counter += 1


def orchestrate_cmsearch_gtf(
    input_file: Path,
    output_gtf: Path,
    region_name: str,
    rfam_seed_descriptions: Path,
    rfam_selected_models_file: Path,
    output_bed: Path,
):
    seed_descriptions = get_rfam_seed_descriptions(rfam_seed_descriptions)
    cm_models = extract_rfam_metrics(rfam_selected_models_file)
    initial_table_results = parse_rfam_tblout(input_file, region_name)
    unique_table_results = remove_rfam_overlap(initial_table_results)
    filtered_table_results = filter_rfam_results(unique_table_results, cm_models)

    create_cmsearch_gtf(
        filtered_results=filtered_table_results,
        cm_models=cm_models,
        seed_descriptions=seed_descriptions,
        region_name=region_name,
        output_gtf=output_gtf,
        output_bed=output_bed,
    )


def _convert_genblast_gff_to_gtf(gff_file: Path) -> str:
    """
    Convert the content of gtf file in gff format
    gff_file: Path for the gff file
    """
    gtf_string = ""
    with open(gff_file, "r", encoding="utf8") as file_in:
        for line in file_in:
            results = line.split()
            if len(results) == 9:
                results[2] = "exon" if results[2] == "coding_exon" else results[2]
                attributes = set_genblast_attributes(str(results[8]), str(results[2]))
                results[8] = attributes
                converted_line = "\t".join(results)
                gtf_string += converted_line + "\n"
    return gtf_string


def set_genblast_attributes(attributes: str, feature_type: str) -> str:
    """
    Given the list of attributes in the genblast output,
    define the new attributes for the gtf file.
    attributes: GenBlast attribute list
    feature_type: transcript or exon
    Example genBlast output #pylint: disable=line-too-long, trailing-whitespace
    1       genBlastG       transcript      131128674       131137049       252.729 -       .       ID=259447-R1-1-A1;Name=259447;PID=84.65;Coverage=94.22;Note=PID:84.65-Cover:94.22
    1       genBlastG       coding_exon     131137031       131137049       .       -       .       ID=259447-R1-1-A1-E1;Parent=259447-R1-1-A1
    1       genBlastG       coding_exon     131136260       131136333       .       -       .       ID=259447-R1-1-A1-E2;Parent=259447-R1-1-A1
    1       genBlastG       coding_exon     131128674       131130245       .       -       .       ID=259447-R1-1-A1-E3;Parent=259447-R1-1-A1
    """
    converted_attributes = ""
    split_attributes = attributes.split(";")
    if feature_type == "transcript":
        match = re.search(r"Name\=(.+)$", split_attributes[1])
        assert match
        name = match.group(1)
        converted_attributes = f'gene_id "{name}"; transcript_id "{name}";'
    elif feature_type == "exon":
        match = re.search(r"\-E(\d+);Parent\=(.+)\-R\d+\-\d+\-", attributes)
        assert match
        exon_rank = match.group(1)
        name = match.group(2)
        converted_attributes = f'gene_id "{name}"; transcript_id "{name}"; exon_number "{exon_rank}";'  # pylint:disable=line-too-long

    return converted_attributes


def create_genblast_gtf(gff_file: Path, gtf_file: Path) -> str:
    """
    Convert the content of gtf file in gff format
    gff_file: Path for the gff file
    """
    gtf_string = ""
    with open(gff_file, "r", encoding="utf8") as file_in:
        for line in file_in:
            results = line.split()
            if len(results) == 9:
                results[2] = "exon" if results[2] == "coding_exon" else results[2]
                attributes = set_genblast_attributes(str(results[8]), str(results[2]))
                results[8] = attributes
                converted_line = "\t".join(results)
                gtf_string += converted_line + "\n"

    with open(gtf_file, "w", encoding="utf8") as file_out:
        file_out.write(gtf_string)

    return "genblast gtf created"


def create_miniprot_gtf(  # pylint: disable=too-many-locals, too-many-branches, too-many-statements
    input_file: Union[str, Path],
    output_file: Union[str, Path],
) -> None:
    """Convert Miniprot GFF output into GTF format."""

    input_file = Path(input_file)
    output_file = Path(output_file)

    with open(input_file, "r", encoding="utf-8") as input_handle:
        blocks = input_handle.read().split("\n#")

    with open(output_file, "w", encoding="utf-8") as file_out:
        for block in blocks:
            nblock_lines: List[str] = [line for line in block.split("\n") if line]

            if not nblock_lines:
                continue

            header_line = nblock_lines[0]

            nblock_list: List[List[str]] = [line.split("\t") for line in nblock_lines[1:]]

            nblock: NDArray[np.object_] = cast(
                NDArray[np.object_],
                np.array(
                    nblock_list,
                    dtype=object,
                ),
            )

            if "fs:i:" in header_line:
                match_fs = re.search(
                    r"fs:i:(\d+)",
                    header_line,
                )
                if match_fs and int(match_fs.group(1)) != 0:
                    continue

            if "st:i:" in header_line:
                match_st = re.search(
                    r"st:i:(\d+)",
                    header_line,
                )
                if match_st and int(match_st.group(1)) != 0:
                    continue

            if nblock.shape[0] == 0:
                continue

            nrows = nblock.shape[0]

            nblock[0, 2] = nblock[0, 2].replace(
                "mRNA",
                "transcript",
            )

            nblock[1:nrows, 2] = [value.replace("CDS", "exon") for value in nblock[1:nrows, 2]]

            target_info = [
                value.replace("Target=", "")
                for value in re.split(
                    r";|\s",
                    str(nblock[0, 8]),
                )
                if "Target" in value
            ][0]

            gene_transcript = f'gene_id "{target_info}"; ' f'transcript_id "{target_info}";'

            for index in range(nrows):
                if index == 0:
                    nblock[index, 8] = gene_transcript
                else:
                    nblock[index, 8] = f"{gene_transcript} " f'exon_number "{index}";'

            if nblock[nrows - 1, 2] == "stop_codon":
                if nblock[0, 6] == "-":
                    nblock[nrows - 2, 3] = nblock[
                        nrows - 1,
                        3,
                    ]
                    nblock = nblock[:-1]

                if nblock[0, 6] == "+":
                    nblock[nrows - 2, 4] = nblock[
                        nrows - 1,
                        4,
                    ]
                    nblock = nblock[:-1]

            for element in nblock:
                file_out.write("%s\n" % "\t".join(element))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Arguments for script to check contents of transcriptomic gtfs"
    )
    parser.add_argument("--input_file", help="Path to input to convert to gtf")
    parser.add_argument("--output_gtf", help="Path to output logfile recording status of each gtf")
    parser.add_argument("--region_name", default=None, help="Optional region name field")
    parser.add_argument(
        "--rfam_seed_descriptions",
        default=None,
        help="path to rfam seed description file (only required for cmsearch)",
    )
    parser.add_argument(
        "--rfam_selected_models_file",
        default=None,
        help="path to rfam selected model file (only required for cmsearch)",
    )
    parser.add_argument(
        "--output_bed", default=None, help="path to output bedfile (only required for cmsearch)"
    )
    parser.add_argument("--red", action="store_true", help="convert red output to gtf")
    parser.add_argument("--dust", action="store_true", help="convert red output to gtf")
    parser.add_argument("--repeatmasker", action="store_true", help="convert red output to gtf")
    parser.add_argument("--trf", action="store_true", help="convert trf output to gtf")
    parser.add_argument("--cpg", action="store_true", help="convert cpg output to gtf")
    parser.add_argument("--eponine", action="store_true", help="convert eponine output to gtf")
    parser.add_argument("--trnascan", action="store_true", help="convert trnascan output to gtf")
    parser.add_argument("--cmsearch", action="store_true", help="convert cmsearch/rfam output to gtf")
    parser.add_argument("--genblast", action="store_true", help="convert genblast gff to gtf")
    parser.add_argument("--miniprot", action="store_true", help="convert miniprot gff to gtf")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    if args.red:
        create_red_gtf(args.input_file, args.output_gtf)

    if args.dust:
        create_dust_gtf(args.input_file, args.output_gtf, args.region_name)

    if args.repeatmasker:
        create_repeatmasker_gtf(args.input_file, args.output_gtf, args.region_name)

    if args.trf:
        create_trf_gtf(args.input_file, args.output_gtf, args.region_name)

    if args.cpg:
        create_cpg_gtf(args.input_file, args.output_gtf, args.region_name)

    if args.eponine:
        create_eponine_gtf(args.input_file, args.output_gtf, args.region_name)

    if args.trnascan:
        create_trnascan_gtf(args.input_file, args.output_gtf, args.region_name)

    if args.cmsearch:
        orchestrate_cmsearch_gtf(
            input_file=args.input_file,
            output_gtf=args.output_gtf,
            region_name=args.region_name,
            rfam_seed_descriptions=args.rfam_seed_descriptions,
            rfam_selected_models_file=args.rfam_selected_models_file,
            output_bed=args.output_bed,
        )

    if args.genblast:
        create_genblast_gtf(gff_file=args.input_file, gtf_file=args.output_gtf)

    if args.miniprot:
        create_miniprot_gtf(input_file=args.input_file, output_file=args.output_gtf)
