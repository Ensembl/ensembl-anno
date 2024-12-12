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
CPC2 (Coding Potential Calculator 2) analyzes RNA sequences to determine whether
they are likely coding or non-coding based on sequence features and machine learning.

CPC2: a fast and accurate coding potential calculator based on sequence intrinsic features 
Yu-Jian Kang,   De-Chang Yang,   Lei Kong,   Mei Hou,   Yu-Qi Meng,   Liping Wei, Ge Gao
Nucleic Acids Research, Volume 45, Issue W1, 3 July 2017, Pages W12–W16, https://doi.org/10.1093/nar/gkx428
"""
__all__ = ["run_cpc2"]

import argparse
import logging
import logging.config
from pathlib import Path
import re
from typing import List

from ensembl.tools.anno.utils._utils import (
    create_dir,
    check_file,
    run_command,
)

logger = logging.getLogger(__name__)


def run_cpc2(  # pylint:disable=dangerous-default-value
    validation_dir: Path,
    cdna_file: Path,
    cpc2_bin: Path = Path(
        "/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/rnasamba_latest.sif"
    ),
) -> List[List[str]]:
    """
    Reads the output file from CPC2 and extracts coding potential information.

        :param validation_dir: Path to the directory where the CPC2 results will be stored.
        :type validation_dir: Path
        :param cdna_file: Path to the input cDNA file containing transcript sequences in FASTA format.
        :type cdna_file: Path
        :param cpc2_bin: Path to the CPC2 software.
        :type cpc2_bin: Path

        :return: A list of parsed results, where each sublist contains:
            - `transcript_id` (str): Identifier of the RNA transcript.
            - `coding_probability` (str): Probability of the transcript being coding.
            - `coding_potential` (str): Predicted coding potential (e.g., coding/non-coding).
            - `transcript_length` (str): Length of the transcript in nucleotides.
            - `peptide_length` (str): Length of the predicted peptide, if applicable.
        :rtype: List[List[str]]

    """

    cpc2_output_path = validation_dir / "cpc2.tsv"
    check_file(cpc2_output_path)
    cpc2_volume = f"{validation_dir}/:/app:rw"
    cpc2_cmd = [
        "singularity",
        "exec",
        "--bind",
        str(cpc2_volume),
        str(cpc2_bin),
        "python3",
        "/CPC2_standalone-1.0.1/bin/CPC2.py",
        "-i",
        str(cdna_file),
        "--ORF",
        "-o",
        str(cpc2_output_path),
    ]
    logger.info("Running CPC2 command: %s", " ".join(map(str, cpc2_cmd)))
    run_command(cpc2_cmd)
    cpc2_output_path = cpc2_output_path.with_suffix(".txt")
    cpc2_results = read_cpc2_results(cpc2_output_path)
    return cpc2_results


def read_cpc2_results(file_path: Path) -> List[List[str]]:
    """
    Reads CPC2 results from a file and returns a list of coding potential data.

    Args:
        file_path (Path): Path to the CPC2 results file.

    Returns:
        List[List[str]]: A list of records, each containing transcript ID,
        coding probability, coding potential, transcript length, and peptide length.
    """
    results = []

    with file_path.open("r") as file_in:
        for line in file_in:
            line = line.rstrip()
            if re.search(r"^#ID", line):
                continue  # Skip header lines

            eles = line.split("\t")
            if len(eles) != 9:
                continue  # Skip malformed lines

            (
                transcript_id,
                transcript_length,
                peptide_length,
                _,
                _,
                _,
                _,
                coding_probability,
                coding_potential,
            ) = eles  # Unpack all three elements
            results.append(
                [
                    transcript_id,
                    coding_probability,
                    coding_potential,
                    transcript_length,
                    peptide_length,
                ]
            )

    return results


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="CPC2 arguments")
    parser.add_argument("--validation_dir", required=True, help="Validation directory path")
    parser.add_argument("--cdna_file", required=True, help="Path for the cDNA fasta file")
    parser.add_argument(
        "--cpc2_bin",
        help="CPC2 executable path",
    )
    return parser.parse_args()


def main():
    """CPC2's entry-point."""
    args = parse_args()

    log_file_path = create_dir(args.output_dir, "log") / "cpc2.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"

    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )

    run_cpc2(
        Path(args.validation_dir),
        Path(args.cdna_file),
        Path(args.cpc2_bin),
    )


if __name__ == "__main__":
    main()
