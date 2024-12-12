# pylint: disable=missing-module-docstring
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
Loads GTF records for multiple features into an Ensembl database.
"""
__all__ = ["load_results_to_ensembl_db"]

import argparse
import logging
import logging.config
import multiprocessing
from pathlib import Path
import gc
import re
import subprocess
from typing import Generator, List, Dict
import tempfile
from ensembl.tools.anno.utils._utils import (
    create_dir,
)


def load_results_to_ensembl_db(
    genome_file: Path,
    main_output_dir: Path,
    db_details: str,
    features: Dict[str, tuple[int, str, str]],
    num_threads: int = 1,
    single_transcript_genes_loading: bool = False,
) -> None:
    """
    Loads GTF records for multiple features into an Ensembl database.

    This function processes each feature defined in `features`, batches the associated GTF records,
    and then loads them into the Ensembl database using the `generic_load_records_to_ensembl_db` function.

    :param genome_file: Path to the genome file used for the loading process.
    :type genome_file: Path
    :param main_output_dir: The directory where the feature directories and GTF files are stored.
    :type main_output_dir: Path
    :param db_details: Comma-separated database connection details \
        (db_name, db_host, db_port, db_user, db_pass).
    :type db_details: str
    :param features: A dictionary with feature names as keys, and tuples \
        with batch size, load type (e.g., 'gene'), and analysis name (e.g., 'ensembl').
    :type features: (Dict[str, Tuple[int, str, str]])
    :param num_threads: The number of threads to use for loading.
    :type num_threads: int, default 1.
    :param single_transcript_genes_loading: Flag to indicate whether to load single transcript genes.
    :type single_transcript_genes_loading: bool, default False.

    """
    db_loading_script = (
        Path(__file__).parents[5]
        / "perl"
        / "ensembl"
        / "tools"
        / "anno"
        / "support_scripts_perl"
        / "load_gtf_ensembl.pl"
    )
    db_loading_dir = create_dir(main_output_dir, "db_loading")
    # Process each selected feature
    for feature in features:
        # Retrieve settings for the current feature
        batch_size, load_type, analysis_name = features[feature]

        # Construct the path to the GTF file for the feature
        annotation_results_gtf_file = main_output_dir / feature / "annotation.gtf"

        # Check if the file exists and load the feature
        if annotation_results_gtf_file.exists():
            logging.info("Loading %s geneset to db", feature)
            gtf_records = batch_gtf_records(annotation_results_gtf_file, batch_size, load_type)
            generic_load_records_to_ensembl_db(
                single_transcript_genes_loading,
                db_loading_script,
                genome_file,
                db_details,
                db_loading_dir,
                load_type,
                analysis_name,
                gtf_records,
                num_threads,
            )
        else:
            logging.error(
                "Did not find the %s annotation file, so not loading. Path checked: %s",
                feature,
                annotation_results_gtf_file,
            )


def generic_load_records_to_ensembl_db(
    single_transcript_genes_loading: bool,
    db_loading_script: Path,
    genome_file: Path,
    db_details: str,
    db_loading_dir: Path,
    load_type: str,
    analysis_name: str,
    gtf_records: List[List[str]],
    num_threads: int,
) -> None:
    """
    Loads GTF records into the Ensembl database using multiple threads.

    This function processes the GTF records in batches and loads them asynchronously
    into the database using multiprocessing. The function delegates the actual loading
    of records to `multiprocess_load_records_to_ensembl_db` for each batch.

    Args:
        single_transcript_genes_loading (str): Option for loading transcript separately without
        defining the canonical.
        db_loading_script (Path): Path to the Perl script used for loading records into Ensembl.
        genome_file (Path): Path to the genome file.
        db_details (str): Database details.
        db_loading_dir (Path): Directory for storing temporary files related to DB loading.
        load_type (str): Type of the records being loaded (e.g., "gene", "transcript").
        analysis_name (str): Name of the analysis being performed.
        gtf_records (List[List[str]]): A list of GTF record batches to be loaded.
        num_threads (int): Number of threads to use for loading records.

    Returns:
        None: This function performs an asynchronous operation and doesn't return any value.
    """
    pool = multiprocessing.Pool(int(num_threads))#pylint:disable=consider-using-with
    for record_batch in gtf_records:
        pool.apply_async(
            multiprocess_load_records_to_ensembl_db,
            args=(
                single_transcript_genes_loading,
                db_loading_script,
                genome_file,
                db_details,
                db_loading_dir,
                load_type,
                analysis_name,
                record_batch,
            ),
        )
    try:
        pool.close()
        pool.join()
    except Exception as e:#pylint:disable=broad-exception-caught
        logging.error("Error during pool shutdown: %s", e)
    finally:
        pool.terminate()


def multiprocess_load_records_to_ensembl_db(
    single_transcript_genes_loading: bool,
    db_loading_script: str,
    genome_file: str,
    db_details: str,
    output_dir: str,
    load_type: str,
    analysis_name: str,
    record_batch: List[str],
) -> None:
    """
    Loads a batch of GTF records into an Ensembl database using a temporary GTF file.

    This function writes the provided `record_batch` to a temporary GTF file, constructs
    a command to load the records into the database using a Perl script, and executes it
    via a subprocess. After execution, the temporary GTF file is removed.

    Args:
        single_transcript_genes_loading (bool): Flag to indicate if single transcript genes
                                                should be loaded.
        db_loading_script (str): Path to the Perl script used to load records into the database.
        genome_file (str): Path to the genome file to be used by the loading script.
        db_details (str): Comma-separated string containing database connection details 
        (db_name, db_host, db_port, db_user, db_pass).
        output_dir (str): Directory where the temporary GTF file will be stored.
        load_type (str): Type of the record being loaded (e.g., "gene").
        analysis_name (str): Name of the analysis being performed (e.g., "ensembl").
        record_batch (List[str]): List of GTF records (as strings) to be written to the temporary file.

    Returns:
        None: This function performs a subprocess call and does not return any value.
    """
    with tempfile.NamedTemporaryFile(mode="w+t", delete=False, dir=output_dir) as gtf_temp_out:
        for line in record_batch:
            gtf_temp_out.writelines(line)
            # gtf_temp_file_path = gtf_temp_out.name

    (db_name, db_host, db_port, db_user, db_pass) = db_details.split(",")

    loading_cmd = [
        "perl",
        str(db_loading_script),
        "-genome_file",
        str(genome_file),
        "-dbname",
        db_name,
        "-host",
        db_host,
        "-port",
        str(db_port),
        "-user",
        db_user,
        "-pass",
        db_pass,
        "-gtf_file",
        str(gtf_temp_out),
        "-analysis_name",
        analysis_name,
        "-load_type",
        load_type,
    ]

    if load_type == "gene" and analysis_name == "ensembl":
        loading_cmd.extend(
            [
                "-protein_coding_biotype",
                "anno_protein_coding",
                "-non_coding_biotype",
                "anno_lncRNA",
            ]
        )

        if single_transcript_genes_loading:
            loading_cmd.append("-make_single_transcript_genes")

    logging.info(" ".join(loading_cmd))
    try:
        # Run the loading command as a subprocess
        subprocess.run(loading_cmd, check=True)
    except subprocess.CalledProcessError as e:
        logging.error("Error executing loading command: %s", e)
    finally:
        # Clean up temporary files
        try:
            gtf_temp_out.unlink()
            logging.info("Temporary file %s removed.", gtf_temp_out)
        except OSError as e:
            logging.error("Error removing temporary file %s, %s", gtf_temp_out, e)
        # Run garbage collection to free up memory
        gc.collect()
        logging.info("Finished processing batch from %s", gtf_temp_out)


# The following function assumes the file has unique ids for the parent features.
# It then batches them into
# sets of records based on the batch size passed in
def batch_gtf_records(input_gtf_file: Path, batch_size: int, record_type: str) -> List[List[str]]:
    """
    Processes a GTF file in batches and returns a list of record batches.

    Args:
        input_gtf_file (Path): The path to the input GTF file.
        batch_size (int): The number of records per batch.
        record_type (str): The type of records to process ('gene' or 'single_line_feature').

    Returns:
        List[List[str]]: A list of batches, each containing GTF lines as lists of strings.
    """

    def read_gtf_lines(gtf_file: Path) -> Generator[str, None, None]:
        """
        Generator that reads and yields valid GTF lines from the file.

        Args:
            gtf_file (Path): The path to the input GTF file.

        Yields:
            str: A valid GTF line.
        """
        with open(gtf_file) as gtf_in:
            for line in gtf_in:
                if re.match(r"^#", line):  # Skip comment lines
                    continue
                elems = line.split("\t")
                if len(elems) == 9:  # Ensure valid GTF line
                    yield line

    records: List[List[str]] = []
    current_record_batch: List[str] = []
    record_counter: int = 0
    current_gene_id: str = ""

    # Iterate over the GTF lines
    for line in read_gtf_lines(input_gtf_file):
        if record_type == "gene":
            match = re.search(r'gene_id "([^"]+)"', line)
            if match:
                gene_id = match.group(1)

                if not current_gene_id:
                    record_counter += 1
                    current_gene_id = gene_id

                if gene_id != current_gene_id:
                    record_counter += 1
                    if record_counter % batch_size == 0:
                        records.append(current_record_batch)
                        current_record_batch = []
                    current_gene_id = gene_id

                current_record_batch.append(line)

        elif record_type == "single_line_feature":
            record_counter += 1
            if record_counter % batch_size == 0:
                records.append(current_record_batch)
                current_record_batch = []

            current_record_batch.append(line)

    # Append the last batch of records
    if current_record_batch:
        records.append(current_record_batch)

    return records


def build_feature_settings(
    args: argparse.Namespace, repeatmasker_analysis: str = "repeatmask_repbase_human"
) -> Dict[str, tuple[int, str, str]]:
    """
    Builds the feature settings dictionary based on the enabled flags from the argument object.

    This function uses the provided argument flags (`args`) to determine which features
    to include in the returned settings dictionary. The dictionary contains the configuration
    for each feature, including batch size, record type, and analysis name.

    Args:
        args (object): An object containing flags (e.g., `args.annotation`, `args.rfam`)
                       that determine which features should be included.
        repeatmasker_analysis (str, optional): The analysis name for RepeatMasker output.
                                                Defaults to "repeatmask_repbase_human".

    Returns:
        Dict[str, tuple[int, str, str]]: A dictionary where keys are feature names (e.g., "annotation_output")
                                          and values are tuples containing batch size (int),
                                          record type (str), and analysis name (str).
    """
    # Base configurations for all features
    all_features = {
        "annotation_output": (200, "gene", "ensembl"),
        "rfam_output": (500, "gene", "ncrna"),
        "trnascan_output": (500, "gene", "ncrna"),
        "dust_output": (500, "single_line_feature", "dust"),
        "red_output": (500, "single_line_feature", "repeatdetector"),
        "trf_output": (500, "single_line_feature", "trf"),
        "repeatmasker_output": (500, "single_line_feature", repeatmasker_analysis),
        "cpg_output": (500, "single_line_feature", "cpg"),
        "eponine_output": (500, "single_line_feature", "eponine"),
    }

    # Start with an empty dictionary
    feature_settings = {}

    # Add features based on enabled flags
    if args.annotation:
        feature_settings["annotation_output"] = all_features["annotation_output"]
    if args.rfam:
        feature_settings["rfam_output"] = all_features["rfam_output"]
    if args.trnascan:
        feature_settings["trnascan_output"] = all_features["trnascan_output"]
    if args.dust:
        feature_settings["dust_output"] = all_features["dust_output"]
    if args.red:
        feature_settings["red_output"] = all_features["red_output"]
    if args.trf:
        feature_settings["trf_output"] = all_features["trf_output"]
    if args.repeatmasker:
        feature_settings["repeatmasker_output"] = all_features["repeatmasker_output"]
    if args.cpg:
        feature_settings["cpg_output"] = all_features["cpg_output"]
    if args.eponine:
        feature_settings["eponine_output"] = all_features["eponine_output"]

    return feature_settings


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Load multiple features into the Ensembl database.")
    parser.add_argument("--output_dir", required=True, help="Output directory path")
    parser.add_argument("--genome_file", required=True, help="Path for the genome fasta file")
    parser.add_argument("--annotation", action="store_true", help="Enable annotation output feature")
    parser.add_argument("--rfam", action="store_true", help="Enable RFAM feature")
    parser.add_argument("--trnascan", action="store_true", help="Enable tRNA scan feature")
    parser.add_argument("--dust", action="store_true", help="Enable dust feature")
    parser.add_argument("--red", action="store_true", help="Enable RED feature")
    parser.add_argument("--trf", action="store_true", help="Enable TRF feature")
    parser.add_argument("--repeatmasker", action="store_true", help="Enable RepeatMasker feature")
    parser.add_argument("--cpg", action="store_true", help="Enable CpG feature")
    parser.add_argument("--eponine", action="store_true", help="Enable Eponine feature")
    parser.add_argument(
        "--db_details",
        required=True,
        type=str,
        help="Comma-separated string containing database details <db_name,db_host,db_port,db_user,db_pass>",
    )
    parser.add_argument("--num_threads", type=int, default=1, help="Number of threads")

    parser.add_argument(
        "--single_transcript_genes_loading",
        action="store_true",
        help="Flag to indicate whether to load single transcript genes (default is False).",
    )
    return parser.parse_args()


def main():
    """Loading into core database's entry-point."""
    args = parse_args()

    # Build the feature settings dictionary dynamically
    feature_setting = build_feature_settings(args, repeatmasker_analysis="repeatmask_repbase_human")

    load_results_to_ensembl_db(
        args.genome_file,
        args.output_dir,
        args.db_details,
        feature_setting,
        args.num_threads,
        args.single_transcript_genes_loading,
    )
