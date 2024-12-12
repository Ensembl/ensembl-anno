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
Run RNA-Samba to predict the coding potential of transcripts.

RNA-Samba is a machine learning-based tool for distinguishing coding 
and non-coding RNA transcripts based on sequence features. 

Camargo, A. P., Sourkov, V., Pereira, G. A. G. & Carazzolle, M. F.. 
"RNAsamba: neural network-based assessment of the protein-coding potential 
of RNA sequences" NAR Genomics and Bioinformatics 2, lqz024 (2020).
"""
__all__ = ["run_rnasamba"]

import argparse
import logging
import logging.config
from pathlib import Path
import re
from typing import List, Union

from ensembl.tools.anno.utils._utils import (
    create_dir,
    check_file,
    run_command,
)

logger = logging.getLogger(__name__)


def run_rnasamba(  # pylint:disable=dangerous-default-value
    validation_dir: Path,
    cdna_file: Path,
    rnasamba_weights: Path = Path(
        "/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/rnasamba_data/full_length_weights.hdf5"  # pylint:disable=line-too-long
    ),
    rnasamba_bin: Path = Path(
        "/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/rnasamba_latest.sif"
    ),
) -> List[List[str]]:
    """
    Runs the RNA-Samba classification tool within a Singularity container to
    evaluate the coding potential of cDNA sequences.

    This function uses RNA-Samba to classify sequences in the provided cDNA file
    and outputs results in the specified validation directory.

        :param validation_dir: Path to the directory where the RNA-Samba results will be stored.
        :type validation_dir: Path
        :param cdna_file: Path to the input cDNA file containing transcript sequences in FASTA format.
        :type cdna_file: Path
        :param rnasamba_weights: Path to the pre-trained RNA-Samba weights file.
        :type rnasamba_weights: Path
        :param rnasamba_bin: Path to the RNA-Samba software.
        :type rnasamba_bin: Path

        :return: A list of RNA-Samba results, where each entry contains
                - Transcript ID
                - Coding probability
                - Coding potential
        :rtype: List[List[str]]

    """

    rnasamba_output_path = validation_dir / "rnasamba.tsv.txt"
    check_file(rnasamba_output_path)
    rnasamba_output_path = validation_dir / "rnasamba.tsv.txt"
    rnasamba_volume = f"{validation_dir}/:/app:rw"
    rnasamba_cmd = [
        "singularity",
        "exec",
        "--bind",
        str(rnasamba_volume),
        str(rnasamba_bin),
        "rnasamba",
        "classify",
        str(rnasamba_output_path),
        str(cdna_file),
        str(rnasamba_weights),
    ]

    logger.info("Running Rnasamba command: %s", " ".join(map(str, rnasamba_cmd)))
    run_command(rnasamba_cmd)
    rnasamba_results = read_rnasamba_results(rnasamba_output_path)
    return rnasamba_results


def read_rnasamba_results(file_path: Path) -> List[List[str]]:
    """
    Reads RNA-Samba results from a file and returns a list of coding potential data.

    Args:
        file_path (Path): Path to the RNA-Samba results file.

    Returns:
        List[List[str]]: A list of records, each containing transcript ID,
        coding probability, and coding potential.
    """
    results = []

    with file_path.open("r") as file_in:
        for line in file_in:
            line = line.rstrip()
            if re.search(r"^sequence_name", line):
                continue  # Skip header lines

            eles = line.split("\t")
            if len(eles) != 3:
                continue  # Skip malformed lines

            transcript_id, coding_probability, coding_potential = eles  # Unpack all three elements
            results.append([transcript_id, coding_probability, coding_potential])

    return results


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Rnasamba arguments")
    parser.add_argument("--validation_dir", required=True, help="Validation directory path")
    parser.add_argument("--cdna_file", required=True, help="Path for the cDNA fasta file")

    parser.add_argument(
        "--rnasamba_weights",
        help="Rnasamba weight path",
    )
    parser.add_argument(
        "--rnasamba_bin",
        help="Rnasamba executable path",
    )
    return parser.parse_args()


def main():
    """Rnasamba's entry-point."""
    args = parse_args()

    log_file_path = create_dir(args.output_dir, "log") / "rnasamba.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"

    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )

    run_rnasamba(
        Path(args.validation_dir), Path(args.cdna_file), Path(args.rnasamba_weights), Path(args.rnasamba_bin)
    )


if __name__ == "__main__":
    main()
