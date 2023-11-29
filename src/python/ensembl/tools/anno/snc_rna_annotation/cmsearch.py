# See the NOTICE file distributed with this work for additional information #pylint: disable=missing-module-docstring
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
Infernal and its "cmsearch" tool are used for detecting sncRNAs in sequence databases.
sncRNA diversity: Small non-coding RNAs (sncRNAs) constitute a diverse group of RNA
molecules that play critical roles in various cellular processes, including gene regulation,
RNA interference, and post-transcriptional modifications. There are different types of sncRNAs,
such as microRNAs (miRNAs), small interfering RNAs (siRNAs), and small nucleolar RNAs (snoRNAs).
Despite their small size, many sncRNAs exhibit conserved structural features or sequence motifs
across species, essential to identify and study them.
Covariance models (CMs) can represent conserved RNA secondary structures as well as conserved
sequence patterns. This makes them well-suited for detecting sncRNAs in sequence databases.

Nawrocki, E. P., Kolbe, D. L., & Eddy, S. R. (2009). Infernal 1.0: inference of RNA alignments.
Bioinformatics, 25(10), 1335-1337.
"""
__all__ = ["run_cmsearch"]

import argparse
import logging
import logging.config
import multiprocessing
import os
from os import PathLike
from pathlib import Path
import re
import subprocess
import tempfile
from typing import List, Dict, Union, Any
import psutil

from ensembl.tools.anno.utils._utils import (
    check_exe,
    create_dir,
    check_gtf_content,
    get_seq_region_length,
    get_slice_id,
    slice_output_to_gtf,
    get_sequence,
)

logger = logging.getLogger(__name__)


def run_cmsearch(
    genome_file: PathLike,
    output_dir: Path,
    rfam_accession_file: Path,
    rfam_cm_db: Path = Path("/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.0/Rfam.cm"),
    rfam_seeds_file: Path = Path("/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.0/Rfam.seed"),
    cmsearch_bin: Path = Path("cmsearch"),
    rnafold_bin: Path = Path("RNAfold"),
    num_threads: int = 1,
) -> None:
    """
    Search CM(s) against a Rfam database

        :param genome_file : Genome file path.
        :type genome_file: PathLike
        :param output_dir : Working directory path.
        :type output_dir : Path
        :param rfam_accessions : list of Rfam accessions.
        :type rfam_accessions : Path
        :param rfam_cm_db : Rfam database with cm models.
        :type rfam_cm_db : Path
        :param rfam_seed : Rfam seeds file.
        :type rfam_seed : Path
        :param cmsearch_bin : cmsearch software path.
        :type cmsearch_bin : Path
        :param rnafold_bin: RNAfold software path.
        :type rnafold_bin : Path
        :param num_threads : int, number of threads.
        :type num_threads :int

        :return: None
        :rtype: None
    """
    check_exe(cmsearch_bin)
    rfam_dir = create_dir(output_dir, "rfam_output")
    os.chdir(rfam_dir)
    output_file = rfam_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Cmsearch gtf file exists, skipping analysis")
            return
    else:
        logger.info("No gtf file, go on with the analysis")

    rfam_selected_models_file = rfam_dir / "rfam_models.cm"
    with open(rfam_accession_file) as rfam_accessions_in:
        rfam_accessions = rfam_accessions_in.read().splitlines()

    with open(rfam_cm_db, "r") as rfam_cm_in:
        rfam_data = rfam_cm_in.read()

    rfam_models = rfam_data.split("//\n")
    with open(rfam_selected_models_file, "w+") as rfam_cm_out:
        for model in rfam_models:
            # The Rfam.cm file has INFERNAL and HMMR models, both are needed at this point
            # Later we just want the INFERNAL ones for looking at thresholds
            match = re.search(r"(RF\d+)", model)
            if match:
                model_accession = match.group(1)
                if model_accession in rfam_accessions:
                    rfam_cm_out.write(model + "//\n")

    seed_descriptions = _get_rfam_seed_descriptions(rfam_seeds_file)
    cm_models = _extract_rfam_metrics(rfam_selected_models_file)

    logger.info("Creating list of genomic slices")
    seq_region_to_length = get_seq_region_length(genome_file, 5000)
    slice_ids_per_region = get_slice_id(seq_region_to_length, slice_size=1000000, overlap=0, min_length=5000)

    cmsearch_cmd = [
        cmsearch_bin,
        "--rfam",
        "--cpu",
        int(num_threads),
        "--nohmmonly",
        "--cut_ga",
        "--tblout",
    ]
    logger.info("Running Rfam")

    pool = multiprocessing.Pool(int(num_threads))  # pylint: disable=consider-using-with
    for slice_id in slice_ids_per_region:
        pool.apply_async(
            _multiprocess_cmsearch,
            args=(
                cmsearch_cmd,
                slice_id,
                str(genome_file),
                str(rfam_dir),
                str(rfam_selected_models_file),
                str(cm_models),
                str(seed_descriptions),
                rnafold_bin
                # memory_limit,
            ),
        )
    pool.close()
    pool.join()
    # Retry failed runs
    _handle_failed_jobs(
        rfam_dir, cmsearch_cmd, genome_file, rfam_selected_models_file, cm_models, seed_descriptions
    )
    slice_output_to_gtf(output_dir=rfam_dir, unique_ids=True, file_extension=".rfam.gtf")


def _handle_failed_jobs(
    rfam_dir: Path,
    cmsearch_cmd: List,
    genome_file: PathLike,
    rfam_selected_models_file: Path,
    cm_models: Dict,
    seed_descriptions: Dict,
    rnafold_bin: Path = Path("RNAfold"),
) -> None:
    """Retry Rfam failed jobs using available cores and memory

    Args:
        rfam_dir (Path): Rfam output dir.
        cmsearch_cmd (List): Cmsearch command.
        genome_file (PathLike): Genome softamsked file.
        rfam_selected_models_file (Path): Rfam model file.
        cm_models (dict): Metrics of Rfam models.
        seed_descriptions (dict): Rfam seed models file.
        rnafold_bin: RNAfold software path.
    """
    exception_file_list = rfam_dir.glob("*.rfam.except")

    # Get the available system memory and CPU cores
    available_memory = psutil.virtual_memory().available
    # psutil.cpu_count(logical=False): number of physical CPU cores.
    # psutil.cpu_count(logical=True): total number of logical CPU cores (physical cores and virtual cores)
    # For computationally intensive tasks, using physical cores might be more appropriate,
    # while for tasks that can benefit from parallelism, using logical cores might be beneficial.
    available_cores = psutil.cpu_count(logical=False)

    # Calculate the optimal memory per core and the number of cores to use
    memory_per_core = available_memory // available_cores
    num_cores = min(available_cores, available_memory // memory_per_core)
    pool = multiprocessing.Pool(num_cores)  # pylint: disable=consider-using-with
    for exception_file in exception_file_list:
        # exception_file_name = os.path.basename(exception_file_path)
        logger.info("Running himem job for failed region:%s\n", exception_file)
        match = re.search(r"(.+)\.rs(\d+)\.re(\d+)\.", exception_file.name)
        if match:
            except_region = match.group(1)
            except_start = match.group(2)
            except_end = match.group(3)
            except_slice_id = [except_region, except_start, except_end]
            pool.apply_async(
                _multiprocess_cmsearch,
                args=(
                    cmsearch_cmd,
                    except_slice_id,
                    str(genome_file),
                    str(rfam_dir),
                    str(rfam_selected_models_file),
                    str(cm_models),
                    str(seed_descriptions),
                    memory_per_core,
                    rnafold_bin,
                ),
            )
    pool.close()
    pool.join()


def _get_rfam_seed_descriptions(rfam_seeds_file: PathLike) -> dict:
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

    """
    # pylint: disable=W0105
    TO DELETE IF THE ABOVE CODE WORKS #pylint: disable=pointless-string-statement
    # pylint: disable=W0105
        domain_match = re.search(
            r"^\\#=GF AC   (.+)", seed # pylint: disable=anomalous-backslash-in-string
        )
        if domain_match:
            domain = domain_match.group(1)
            descriptions[domain] = {}
            continue

        description_match = re.search(
            r"^\\#=GF DE   (.+)", seed #pylint:disable=anomalous-backslash-in-string
        )  # pylint: disable =anomalous-backslash-in-string
        if description_match:
            description = description_match.group(1)
            descriptions[domain]["description"] = description
            continue

        name_match = re.search(
            "^\\#=GF ID   (.+)", seed #pylint:disable=anomalous-backslash-in-string
        )  # pylint: disable =anomalous-backslash-in-string
        if name_match:
            name = name_match.group(1)
            descriptions[domain]["name"] = name
            continue

        type_match = re.search(
            "^\\#=GF TP   Gene; (.+)", seed #pylint:disable=anomalous-backslash-in-string
        )  # pylint: disable =anomalous-backslash-in-string
        if type_match:
            rfam_type = type_match.group(1)
            descriptions[domain]["type"] = rfam_type
            continue
            """
    return descriptions


def _extract_rfam_metrics(rfam_selected_models: PathLike) -> dict:
    """Get name, description, length, max length, threshold of each Rfam model.

    Args:
        rfam_selected_models_file (PathLike): Path for Rfam models.

    Returns:
        parsed_cm_data: disctionary with Rfam metrics.
    """
    with open(rfam_selected_models, "r") as rfam_cm_in:
        rfam_models = rfam_cm_in.read().split("//\n")

        # rfam_models = rfam_data.split("//\n")
        parsed_cm_data: Dict[str, Dict[str, Union[str, int]]] = {}
        for model in rfam_models:
            # temp = model.split("\n")
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
                """ TO DELETE IF ABOVE SOLUTION WORKS
                for line in model.split("\n"):
                    name_match = re.search("^NAME\\s+(\\S+)", line)
                    if name_match:
                        parsed_cm_data[model_name]["-name"] = name_match.group(1)
                        continue

                    description_match = re.search(r"^DESC\\s+(\\S+)", line) #pylint: disable=anomalous-backslash-in-string
                    if description_match:
                        parsed_cm_data[model_name]["-description"] = description_match.group(
                            1
                        )
                        continue

                    length_match = re.search(r"^CLEN\\s+(\\d+)", line) #pylint:disable=anomalous-backslash-in-string
                    if length_match:
                        parsed_cm_data[model_name]["-length"] = length_match.group(1)
                        continue

                    max_length_match = re.search(r"^W\\s+(\\d+)", line) #pylint:disable=anomalous-backslash-in-string
                    if max_length_match:
                        parsed_cm_data[model_name]["-maxlength"] = max_length_match.group(1)
                        continue

                    accession_match = re.search(r"^ACC\\s+(\\S+)", line) #pylint:disable=anomalous-backslash-in-string
                    if accession_match:
                        parsed_cm_data[model_name]["-accession"] = accession_match.group(1)
                        continue

                    threshold_match = re.search(r"^GA\\s+(\\d+)", line) #pylint:disable=anomalous-backslash-in-string
                    if threshold_match:
                        parsed_cm_data[model_name]["-threshold"] = threshold_match.group(1)
                        continue
                        """

    return parsed_cm_data


def _multiprocess_cmsearch(
    cmsearch_cmd: List[str],
    slice_id: List[str],
    genome_file: Path,
    rfam_dir: Path,
    rfam_selected_models_file: Path,
    cm_models: Dict,
    seed_descriptions: Dict,
    rnafold_bin: Path = Path("RNAfold"),
) -> None:
    """Run cmsearch on multiprocess on genomic slices
    Args:
        cmsearch_cmd (List): Cmsearch command to execute.
        slice_id: List of slice IDs.
        genome_file (PathLike): Genome softamsked file.
        rfam_dir (Path): Rfam output dir.
        rfam_selected_models_file (Path): Rfam model file.
        cm_models (dict): Metrics of Rfam models.
        seed_descriptions (dict): Rfam seed models file.
        rnafold_bin: RNAfold software path.
    """
    region_name, start, end = slice_id
    logger.info(
        "Processing Rfam data using cmsearch against slice: %s:%s:%s",
        region_name,
        start,
        end,
    )
    seq = get_sequence(region_name, int(start), int(end), 1, genome_file, rfam_dir)
    slice_name = f"{region_name}.rs{start}.re{end}"
    # I NEED TO SEE THE OUTPUT TO SET BIOTYPES
    # with tempfile.TemporaryDirectory(dir=rfam_dir) as tmpdirname:
    slice_file = rfam_dir / f"{slice_name}.fa"
    with open(slice_file, "w+", encoding="utf8") as region_out:
        region_out.write(f">{region_name}\n{seq}\n")
    region_tblout = rfam_dir / f"{slice_file}.tblout"
    region_results = rfam_dir / f"{slice_name}.rfam.gtf"
    exception_results = rfam_dir / f"{slice_name}.rfam.except"
    cmsearch_cmd.append(str(region_tblout))
    cmsearch_cmd.append(str(rfam_selected_models_file))
    cmsearch_cmd.append(str(slice_file))
    logger.info("cmsearch_cmd: %s", cmsearch_cmd)
    # to TEST
    # if memory_limit is not None:
    #    cmsearch_cmd = prlimit_command(cmsearch_cmd, memory_limit)

    return_value = None
    try:
        return_value = subprocess.check_output(cmsearch_cmd)
    except subprocess.CalledProcessError as ex:
        # Note that writing to file was the only option here.
        # If return_value was passed back, eventually it would clog the
        # tiny pipe that is used by the workers to send info back.
        # That would mean that it would eventually just be one
        # worked running at a time
        logger.error(
            "Issue processing the following region with cmsearch: %s %s-%s ", region_name, start, end
        )
        logger.error("Return value: %s", return_value)
        with open(exception_results, "w+") as exception_out:
            exception_out.write(f"{region_name} {start} {end}\n")
        region_results.unlink()
        region_tblout.unlink()
        raise ex  # Re-raise the exception to allow caller to handle it further if needed

    initial_table_results = _parse_rfam_tblout(region_tblout, region_name)
    unique_table_results = _remove_rfam_overlap(initial_table_results)
    filtered_table_results = _filter_rfam_results(unique_table_results, cm_models)
    _create_rfam_gtf(
        filtered_table_results,
        cm_models,
        seed_descriptions,
        region_name,
        region_results,
        genome_file,
        rfam_dir,
        rnafold_bin,
    )
    # os.remove(slice_file)
    # os.remove(region_tblout_file_path)
    # gc.collect()


def _parse_rfam_tblout(region_tblout: Path, region_name: str) -> List[Dict[str, Union[str, int]]]:
    """Parse cmsearch output
    col 0 Target Name : This is the name of the target sequence or sequence region that matched the query.
    col 2 Query name : This is the name of the query sequence or model that was used for the search.
    col 3 Accession : This usually refers to a unique identifier for the target sequence.
    col 5 Query Start : The position where the match starts on the query sequence.
    col 6 Query End : The position where the match ends on the query sequence.
    col 7 Target Start : The position where the match starts on the target sequence.
    col 8 Target End : The position where the match ends on the target sequence.
    col 9 Strand : Indicates the orientation of the match on the target sequence.
            It could be + for the forward strand or - for the reverse strand.
    col 14 Hit Score : The score assigned to this match. Higher scores generally indicate better matches.
    col 15 E-value : This is a statistical measure of the number of hits one can expect to see
            when searching a database of a particular size.

    Args:
        region_tblout (Path): Cmsearch output for the region name.
        region_name (str): Region name.

    Returns:
        List[Dict[str, Union[str, int]]]: Formatted output
    """
    with open(region_tblout, "r") as rfam_tbl_in:
        rfam_tbl_data = rfam_tbl_in.read()

    tbl_results = rfam_tbl_data.split("\n")

    all_parsed_results: List[Dict[str, Union[str, int]]] = []
    for result in tbl_results:
        parsed_tbl_data: Dict[str, Union[str, int]] = {}
        if re.compile(rf"^{region_name}").match(result):
            hit = result.split()
            parsed_tbl_data = {
                "accession": str(hit[3]),
                "start": str(hit[7]),
                "end": str(hit[8]),
                "strand": 1 if hit[9] == "+" else -1,
                "query_name": str(hit[2]),
                "score": str(hit[14]),
            }
            all_parsed_results.append(parsed_tbl_data)
    return all_parsed_results


def _remove_rfam_overlap(
    parsed_tbl_data: List[Dict[str, Union[str, int]]]
) -> List[Dict[str, Union[str, int]]]:
    """
    Remove Rfam mdoels overlapping and with a lower score.

    Args:
        parsed_tbl_data (List[Dict[str, Union[str, int]]]): Cmsearch output

    Returns:
        List[Dict[str, Union[str, int]]]: Final Rfam models
    """
    excluded_structures = {}
    chosen_structures = []
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


def _filter_rfam_results(unfiltered_tbl_data: List[Dict[str, Union[str, int]]], cm_models: Dict) -> List:
    """Filter Rfam models according to type and a set of thresholds.

    Args:
        unfiltered_tbl_data (List[Dict[str, Union[str, int]]]): Unfiltered Rfam output.
        cm_models (Dict): Rfam models.

    Returns:
        List: List of filtered models
    """

    """
    for structure in unfiltered_tbl_data:
        query = structure["query_name"]
        if query in cm_models:
            threshold = cm_models[query]["-threshold"]
            if query == "LSU_rRNA_eukarya":
                threshold = 1700
            elif query == "LSU_rRNA_archaea":
                continue
            elif query == "LSU_rRNA_bacteria":
                continue
            elif query == "SSU_rRNA_eukarya":
                threshold = 1600
            elif query == "5_8S_rRNA":
                threshold = 85
            elif query == "5S_rRNA":
                threshold = 75
            if threshold and float(structure["score"]) >= float(threshold):
                filtered_results.append(structure)
                """
    filtered_results = []
    thresholds = {"LSU_rRNA_eukarya": 1700, "SSU_rRNA_eukarya": 1600, "5_8S_rRNA": 85, "5S_rRNA": 75}
    for structure in unfiltered_tbl_data:
        query = structure["query_name"]
        if query in ["LSU_rRNA_archaea", "LSU_rRNA_bacteria"]:
            threshold = cm_models.get(str(query), {}).get("-threshold")
        else:
            threshold = thresholds.get(str(query), cm_models.get(str(query), {}).get("-threshold"))
        if threshold is not None and float(structure["score"]) >= float(threshold):
            filtered_results.append(structure)
    return filtered_results


# NOTE: The below are notes from the perl code (which has extra code)
# about possible improvements that are not implemented there
# Although not included in RefSeq filters, additional filters that
# consider sizes and score_to_size ratios can be applied
# in future work to further exclude FPs
#
# my $is_valid_size = $mapping_length > $min_length && $mapping_length < $max_length ? 1 : 0;
# my $score_size_ratio = $result->{'score'} / $mapping_length;


def _create_rfam_gtf(
    filtered_results: List,
    cm_models: Dict,
    seed_descriptions: Dict,
    region_name: str,
    region_results: Path,
    genome_file: Path,
    rfam_dir: Path,
    rnafold_bin: Path = Path("RNAfold"),
) -> None:
    """Convert RFam output per single region in gtf format

    Args:
        filtered_results : Filtered Rfam results without overlapping.
        cm_models :  Rfam database.
        seed_descriptions : Rfam seed file.
        region_name : Slice name.
        region_results : Rfam output file.
        genome_file : Genome file.
        rfam_dir : Output file.
        rnafold_bin: RNAfold software path.
    """
    if not filtered_results:
        return

    with open(region_results, "w+") as rfam_gtf_out:
        gene_counter = 1
        for structure in filtered_results:
            query = structure["query_name"]
            accession = structure["accession"]
            if query not in cm_models and accession not in seed_descriptions:
                continue
            # model = cm_models[query]
            description = seed_descriptions[accession]
            if "type" in description:
                rfam_type = description["type"]
            else:
                description = None
                rfam_type = "misc_RNA"
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
            biotype_mapping = {
                r"^snRNA;": "snRNA",
                r"^snRNA; snoRNA": "snoRNA",
                r"^snRNA; snoRNA; scaRNA;": "scaRNA",
                r"rRNA;": "rRNA",
                r"antisense;": "antisense",
                r"antitoxin;": "antitoxin",
                r"ribozyme;": "ribozyme",
                # domain: domain,
                # domain: "Vault_RNA",  # Note: This will overwrite the previous domain assignment
                domain: "Y_RNA",  # Note: This will overwrite the previous domain assignment
            }

            for pattern, biotype in biotype_mapping.items():
                if not re.match(pattern, rfam_type):
                    biotype = "misc_RNA"

            # TO DELETE IF ABOVE WORKS
            """
            biotype = "misc_RNA"
            if re.match(r"^snRNA;", rfam_type):
                biotype = "snRNA"
            if re.match(r"^snRNA; snoRNA", rfam_type):
                biotype = "snoRNA"
            if re.match(r"^snRNA; snoRNA; scaRNA;", rfam_type):
                biotype = "scaRNA"
            if re.match(r"rRNA;", rfam_type):
                biotype = "rRNA"
            if re.match(r"antisense;", rfam_type):
                biotype = "antisense"
            if re.match(r"antitoxin;", rfam_type):
                biotype = "antitoxin"
            if re.match(r"ribozyme;", rfam_type):
                biotype = "ribozyme"
                
                
            if re.match(r"" + domain, rfam_type):
                biotype = domain
            if re.match(r"" + domain, rfam_type):
                biotype = "Vault_RNA"
            if re.match(r"" + domain, rfam_type):
                biotype = "Y_RNA"
            """
            rna_seq = get_sequence(
                region_name,
                start,
                end,
                rnafold_strand,
                genome_file,
                rfam_dir,
            )
            valid_structure = check_rnafold_structure(rna_seq, rfam_dir, rnafold_bin)

            if not valid_structure:
                continue

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

            rfam_gtf_out.write(transcript_string)
            rfam_gtf_out.write(exon_string)
            gene_counter += 1


def check_rnafold_structure(seq: str, rfam_dir: Path, rnafold_bin: Path = Path("RNAfold")) -> float:
    """RNAfold reads RNA sequences, calculates their minimum free energy (mfe)
       structure and prints the mfe structure in bracket notation and its free energy.

    Args:
        seq : sequence
        rfam_dir : cmsearch working directory
        rnafold_bin : Software path. Defaults to Path("RNAfold").

    Returns:
        float: mfe structure
    """
    # Note there's some extra code in the RNAfold Perl module
    # for encoding the structure into an attrib
    # Could consider implementing this when running
    # for loading into an Ensembl db
    structure = 0
    with tempfile.NamedTemporaryFile(mode="w+t", delete=False, dir=rfam_dir) as rna_temp_in:
        rna_temp_in.writelines(">seq1\n" + seq + "\n")
        rna_in_file_path = rna_temp_in.name
        try:
            rnafold_cmd = [str(rnafold_bin), "--infile", rna_in_file_path]
            # rnafold_output = subprocess.Popen(rnafold_cmd, stdout=subprocess.PIPE)
            # for line in io.TextIOWrapper(rnafold_output.stdout, encoding="utf-8"):
            rnafold_output = subprocess.check_output(rnafold_cmd, encoding="utf-8")
            for line in rnafold_output.splitlines():
                match = re.search(r"([().]+)\s\(\s*(-*\d+.\d+)\)\n$", line)
                if match:
                    structure = int(match.group(1))
                    # free_energy_score = match.group(2)  #TOADD DAF
                    break
        except (subprocess.CalledProcessError, OSError) as e:
            logging.error("Error while running RNAfold: %s", e)
        finally:
            os.remove(rna_in_file_path)

    return structure


# this function is useful for the DnaAlignFeature and it should help when we load into the ensembl db
# NOTE some of the code above and the code commented out here is to do with creating
# a DAF. As we don't have a python concept of this I'm leaving it out for the moment
# but the code below is a reference
#    my $daf = Bio::EnsEMBL::DnaDnaAlignFeature->new(
#      -slice          => $self->queries,
#      -start          => $strand == 1 ? $start : $end,
#      -end            => $strand == 1 ? $end : $start,
#      -strand         => $strand,
#      -hstart         => $hstart,
#      -hend           => $hend,
#      -hstrand        => $strand,
#      -score          => $score,
#      -hseqname       => length($target_name) > 39 ?
# substr($target_name, 0, 39) : $target_name,,
#      -p_value  => $evalue,
#      -align_type => 'ensembl',
#      -cigar_string  => abs($hend - $hstart) . "M",
# -hcoverage    => $RNAfold->score,
#   );


# this function is useful to upload the RNAfold structure in the transcript attribute table
# should be moved into load db module
def encode_str(input_string: str) -> List:
    """Encode the number of repeated sequences reporting the sequence of characters and the number
    of occurences; when the encoding sequence reach 200 it is splitted reporting the start and the
    stop position in the transcript's structure

    Args:
        input_string : input sequence

    Returns:
        List: List of encoded sequences
    """
    codes = []
    start = 1
    count = 0
    code = ""
    last_chrom = ""
    array = []

    for chrom in input_string:
        count += 1

        if chrom == last_chrom:
            array.append(chrom)
        else:
            if code and len(code) > 200 and len(array) == 1:
                codes.append(f"{start}:{count}\t{code}")
                code = ""
                start = count + 1

            if len(array) > 1:
                code += str(len(array))
                array = []

            code += chrom
            last_chrom = chrom

    if len(array) > 1:
        code += str(len(array))

    codes.append(f"{start}:{count}\t{code}")
    return codes


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="RepeatMasker's arguments")
    parser.add_argument("--genome_file", required=True, help="Genome file path")
    parser.add_argument("--output_dir", required=True, help="Output directory path")
    parser.add_argument("--rfam_accessions", required=True, help="List of Rfam accessions.")
    parser.add_argument(
        "--rfam_cm_db",
        default="/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.0/Rfam.cm",
        help="Rfam database with cm models.",
    )
    parser.add_argument(
        "--rfam_seeds_file",
        default="/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.0/Rfam.seed",
        help="CRfam seeds file.",
    )
    parser.add_argument(
        "--cmsearch_bin",
        default="cmsearch",
        help="cmsearch software path.",
    )
    parser.add_argument(
        "--rnafold_bin",
        default="RNAfold",
        help="RNAfold software path.",
    )
    parser.add_argument(
        "--num_threads",
        type=int,
        default=1,
        help="Number of threads",
    )
    return parser.parse_args()


def main():
    """cmsearch's entry-point."""
    args = parse_args()

    log_file_path = create_dir(args.output_dir, "log") / "cmsearch.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"

    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )
    run_cmsearch(
        args.genome_file,
        args.output_dir,
        args.rfam_accession_file,
        args.rfam_cm_db,
        args.rfam_seeds_file,
        args.cmsearch_bin,
        args.rnafold_bin,
        args.num_threads,
    )


if __name__ == "__main__":
    main()
