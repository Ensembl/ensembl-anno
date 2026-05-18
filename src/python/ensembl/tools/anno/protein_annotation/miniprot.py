
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
Miniprot aligns protein sequences to genomic DNA to recover spliced,
protein-to-genome alignments useful for gene annotation.
"""

__all__ = ["run_miniprot"]

import argparse
import logging
import logging.config
import os
import re
import subprocess

from pathlib import Path
from typing import Union

import numpy as np

from src.ensembl.tools.anno.utils._utils import (
    check_exe,
    check_gtf_content,
    create_dir,
)

logger = logging.getLogger(__name__)


def run_miniprot(  # pylint: disable=too-many-arguments
    masked_genome: Path,
    output_dir: Path,
    protein_dataset: Path,
    miniprot_bin: Path = Path("miniprot"),
    num_threads: int = 1,
    top_n: int = 1,
    outs: float = 1.0,
    protein_set: str = "uniprot",
) -> None:
    """Run Miniprot protein-to-genome alignment workflow."""

    check_exe(miniprot_bin)

    if protein_set not in {"uniprot", "orthodb"}:
        raise ValueError(
            "protein_set must be either 'uniprot' or 'orthodb'"
        )

    if protein_set == "uniprot":
        miniprot_dir = create_dir(
            output_dir,
            "miniprot_output/uniprot_output",
        )
    else:
        miniprot_dir = create_dir(
            output_dir,
            "miniprot_output/orthodb_output",
        )

    initial_output_file = miniprot_dir / "initial_annotation.gff"
    output_file = miniprot_dir / "annotation.gtf"

    if output_file.exists():
        transcript_count = check_gtf_content(
            output_file,
            "transcript",
        )
        if transcript_count > 0:
            logger.info(
                "Miniprot GTF file exists, skipping analysis"
            )
            return

    if not masked_genome.exists():
        raise IOError(
            f"Masked genome file does not exist: {masked_genome}"
        )

    if not protein_dataset.exists():
        raise IOError(
            f"Protein file does not exist: {protein_dataset}"
        )

    miniprot_index_file = Path(f"{masked_genome}.mpi")

    if not miniprot_index_file.exists():
        run_miniprot_index(
            miniprot_bin,
            masked_genome,
            miniprot_index_file,
            num_threads,
        )
    else:
        logger.info(
            "Found an existing miniprot index, skipping indexing"
        )

    miniprot_cmd = [
        str(miniprot_bin),
        "-t",
        str(num_threads),
        "-N",
        str(top_n),
        f"--outs={outs}",
        "--gff",
        str(miniprot_index_file),
        str(protein_dataset),
    ]

    logger.info("Running miniprot mapping")
    logger.info(" ".join(miniprot_cmd))

    with open(initial_output_file, "w", encoding="utf-8") as output_handle:
        subprocess.run(
            miniprot_cmd,
            stdout=output_handle,
            check=True,
        )

    logger.info("Completed running miniprot")
    logger.info("Creating standardised GFF")
    generate_miniprot_gtf(miniprot_dir)


def generate_miniprot_gtf(miniprot_dir: str) -> None:
    """Generate a standardised GTF file from Miniprot output."""

    logger.info("generate_miniprot_gff")

    file_out_name = os.path.join(
        miniprot_dir,
        "annotation.gtf",
    )

    for root, _, files in os.walk(miniprot_dir):
        for miniprot_file in files:
            miniprot_path = os.path.join(
                root,
                miniprot_file,
            )

            if miniprot_path.endswith(".gff"):
                convert_miniprot_gff_to_gtf(
                    input_file=miniprot_path,
                    output_file=file_out_name,
                )


def convert_miniprot_gff_to_gtf(
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
            nblock_lines = [
                line for line in block.split("\n")
                if line
            ]

            if not nblock_lines:
                continue

            header_line = nblock_lines[0]

            nblock = [
                line.split("\t")
                for line in nblock_lines[1:]
            ]

            nblock = np.array(
                nblock,
                dtype=object,
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

            nblock[1:nrows, 2] = [
                value.replace("CDS", "exon")
                for value in nblock[1:nrows, 2]
            ]

            target_info = [
                value.replace("Target=", "")
                for value in re.split(
                    ";|\\s",
                    nblock[0, 8],
                )
                if "Target" in value
            ][0]

            gene_transcript = (
                f'gene_id "{target_info}"; '
                f'transcript_id "{target_info}";'
            )

            for index in range(nrows):
                if index == 0:
                    nblock[index, 8] = gene_transcript
                else:
                    nblock[index, 8] = (
                        f'{gene_transcript} '
                        f'exon_number "{index}";'
                    )

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
                file_out.write(
                    "%s\n" % "\t".join(element)
                )


def run_miniprot_index(
    miniprot_bin: Union[str, Path],
    masked_genome: Union[str, Path],
    miniprot_index_file: Union[str, Path],
    num_threads: int,
) -> None:
    """Execute Miniprot indexing on a masked genome."""

    miniprot_bin = Path(miniprot_bin)
    masked_genome = Path(masked_genome)
    miniprot_index_file = Path(miniprot_index_file)

    miniprot_index_cmd = [
        str(miniprot_bin),
        f"-t{num_threads}",
        "-d",
        str(miniprot_index_file),
        str(masked_genome),
    ]

    logger.info(" ".join(miniprot_index_cmd))

    subprocess.run(
        miniprot_index_cmd,
        check=True,
    )

    logger.info(
        "Completed running miniprot indexing"
    )


def parse_args():
    """Parse command line arguments."""

    parser = argparse.ArgumentParser(
        description="Miniprot arguments"
    )

    parser.add_argument(
        "--masked_genome_file",
        required=True,
        help="Masked genome file path",
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory path",
    )

    parser.add_argument(
        "--protein_file",
        required=True,
        help="Path for the protein dataset",
    )

    parser.add_argument(
        "--miniprot_bin",
        default="miniprot",
        help="Miniprot executable path",
    )

    parser.add_argument(
        "--num_threads",
        type=int,
        default=1,
        help="Number of threads",
    )

    parser.add_argument(
        "--protein_set",
        required=True,
        choices=["uniprot", "orthodb"],
        help="Protein set [uniprot, orthodb]",
    )

    return parser.parse_args()


def main() -> None:
    """Miniprot entry-point."""

    args = parse_args()

    log_file_path = (
        create_dir(args.output_dir, "log")
        / "miniprot.log"
    )

    loginipath = (
        Path(__file__).parents[6]
        / "conf"
        / "logging.conf"
    )

    logging.config.fileConfig(
        loginipath,
        defaults={
            "logfilename": str(log_file_path),
        },
        disable_existing_loggers=False,
    )

    run_miniprot(
        Path(args.masked_genome_file),
        Path(args.output_dir),
        Path(args.protein_file),
        Path(args.miniprot_bin),
        args.num_threads,
        protein_set=args.protein_set,
    )


if __name__ == "__main__":
    main()
