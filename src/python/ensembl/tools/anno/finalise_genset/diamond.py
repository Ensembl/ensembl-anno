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
DIAMOND is a sequence aligner for protein and translated DNA searches,
designed for high performance analysis of big sequence data..

Buchfink B, Reuter K, Drost HG, "Sensitive protein alignments at tree-of-life scale
using DIAMOND", Nature Methods 18, 366–368 (2021). doi:10.1038/s41592-021-01101-x
"""
__all__ = ["run_diamond"]

import argparse
import logging
import logging.config
import multiprocessing
from pathlib import Path

from typing import List, Optional, Union

from ensembl.tools.anno.utils._utils import (
    create_dir,
    run_command,
    split_protein_file,
)

logger = logging.getLogger(__name__)


def run_diamond(  # pylint:disable=dangerous-default-value
    validation_dir: Path,
    amino_acid_file: Path,
    diamond_validation_db: Path,
    num_threads: int = 1,
    diamond_bin: Path = Path("diamond"),
) -> List[List[str]]:
    """
    Reads the output file from CPC2 and extracts coding potential information.

        :param validation_dir: Path to the directory where the CPC2 results will be stored.
        :type validation_dir: Path
        :param diamond_validation_db: Path to the DIAMOND database.
        :type diamond_validation_db: Path
        :param amino_acid_file: Path to the input FASTA file containing amino acid sequences.
        :type amino_acid_file: Path
        :param num_threads: Number of threads.
        :type num_threads: int, default 1
        :param diamond_bin: Path to the Diamond software.
        :type diamond_bin: Path

        :return: DIAMOND results, where each sublist contains:
            - `transcript_id` (str): Identifier of the protein sequence.
            - `e_value` (str): E-value of the alignment.
        :rtype: List[List[str]]

    """

    diamond_results = []
    if diamond_validation_db is not None:
        diamond_output_dir = create_dir(validation_dir, "diamond_output")
        diamond_validation(
            diamond_validation_db,
            amino_acid_file,
            diamond_output_dir,
            num_threads,
            diamond_bin,
        )
        diamond_results = read_diamond_results(diamond_output_dir)

    return diamond_results


def diamond_validation(
    diamond_validation_db: Path,
    amino_acid_file: Path,
    diamond_output_dir: Path,
    num_threads: int = 1,
    diamond_bin: Path = Path("diamond"),
) -> None:
    """
    Perform DIAMOND validation on batched protein sequences.

    Args:
        diamond_validation_db (Path): Path to the DIAMOND database.
        amino_acid_file (Path): Path to the input FASTA file containing amino acid sequences.
        diamond_output_dir (Path): Directory to store DIAMOND output files.
        num_threads (int, optional): Number of threads to use for multiprocessing. Defaults to 1.
        diamond_bin (Path, optional): Path to the DIAMOND binary. Defaults to "diamond".

    Returns:
        None
    """
    batched_protein_files = split_protein_file(amino_acid_file, diamond_output_dir, 100)

    pool = multiprocessing.Pool(int(num_threads))  # pylint: disable=consider-using-with
    for batched_protein_file in batched_protein_files:
        pool.apply_async(
            multiprocess_diamond,
            args=(
                batched_protein_file,
                diamond_output_dir,
                diamond_validation_db,
                diamond_bin,
            ),
        )
    pool.close()
    pool.join()


def multiprocess_diamond(
    batched_protein_file: Path,
    diamond_output_dir: Path,
    diamond_validation_db: Path,
    diamond_bin: Path = Path("diamond"),
) -> None:
    """
    Run DIAMOND on a single batch of protein sequences.

    Args:
        batched_protein_file (Path): Path to the batched protein sequence file.
        diamond_output_dir (Path): Directory to store DIAMOND output files.
        diamond_validation_db (Path): Path to the DIAMOND database.
        diamond_bin (Path, optional): Path to the DIAMOND binary. Defaults to "diamond".

    Returns:
        None
    """
    # batch_num = os.path.splitext(batched_protein_file)[0]
    # batch_dir = os.path.dirname(batched_protein_file)
    diamond_output_file = batched_protein_file.with_suffix("dmdout")
    logger.info("Running diamond on %s :", batched_protein_file)

    diamond_cmd = [
        str(diamond_bin),
        "blastp",
        "--query",
        str(batched_protein_file),
        "--db",
        str(diamond_validation_db),
        "--out",
        str(diamond_output_file),
    ]

    logger.info("Running diamond command: %s", " ".join(map(str, diamond_cmd)))
    run_command(diamond_cmd)
    cmd = ["mv", diamond_output_file, diamond_output_dir]
    run_command(cmd)


def read_diamond_results(diamond_output_dir: Path) -> List[List[str]]:
    """
    Parse DIAMOND output files and extract results.

    Args:
        diamond_output_dir (Path): Directory containing DIAMOND output files.

    Returns:
        List[List[str]]: Parsed DIAMOND results, where each sublist contains:
            - `transcript_id` (str): Identifier of the protein sequence.
            - `e_value` (str): E-value of the alignment.
    """
    results = []
    diamond_files = list(diamond_output_dir.glob("*.dmdout"))
    for file_path in diamond_files:
        with file_path.open("r") as file_in:
            for line in file_in:
                line = line.rstrip()

                eles = line.split("\t")
                if len(eles) != 12:
                    continue

                transcript_id = eles[0]
                e_value = eles[10]
                results.append([transcript_id, e_value])
    return results


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="CPC2 arguments")
    parser.add_argument("--validation_dir", required=True, help="Validation directory path")
    parser.add_argument(
        "--amino_acid_file",
        required=True,
        help="Path to the input FASTA file containing amino acid sequences.",
    )
    parser.add_argument("--diamond_validation_db", required=True,help="Path to the DIAMOND database.")
    parser.add_argument(
        "--diamond_bin",
        help="DIAMOND executable path",
    )
    parser.add_argument("--num_threads", type=int, default=1, help="Number of threads")

    return parser.parse_args()


def main():
    """DIAMOND's entry-point."""
    args = parse_args()

    log_file_path = create_dir(args.output_dir, "log") / "diamond.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"

    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )

    run_diamond(
        Path(args.validation_dir),
        Path(args.amino_acid_file),
        Path(args.diamond_validation_db),
        int(args.num_threads),
        Path(args.diamond_bin),
    )


if __name__ == "__main__":
    main()
