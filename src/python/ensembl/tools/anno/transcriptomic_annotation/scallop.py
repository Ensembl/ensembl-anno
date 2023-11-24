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
Scallop is a high-performance tool designed for the accurate and efficient quantification 
of transcriptome assembly. 
It's capable of handling large-scale transcriptomic data while providing precise estimates 
of transcript abundances.
Scallop's algorithmic approach allows it to efficiently reconstruct transcript structures 
and quantify their expression levels, making it a valuable resource for studying gene 
expression and transcriptome analysis.

Shao M, Kingsford C. Accurate assembly of transcripts through phase-preserving graph 
decomposition. Nat Biotechnol.
2017 Dec;35(12):1167-1169. doi: 10.1038/nbt.4020. Epub 2017 Nov 13. PMID: 29131147; PMCID: PMC5722698.
"""

__all__ = ["run_scallop"]

import argparse
import logging
import logging.config
from pathlib import Path
import re
import subprocess

from ensembl.tools.anno.utils._utils import (
    check_exe,
    create_dir,
    check_gtf_content,
)


def run_scallop(
    output_dir: Path,
    scallop_bin: Path = Path("scallop"),
    prlimit_bin: Path = Path("prlimit"),
    stringtie_bin: Path = Path("stringtie"),
    memory_limit: int = 40 * 1024**3,
) -> None:
    """
    Run Scallop assembler on short read data after STAR alignment.

        :param output_dir: Working directory path.
        :type output_dir: Path
        :param scallop_bin: Software path.
        :type scallop_bin: Path, default scallop
        :param prlimit_bin: Software path.
        :type prlimit_bin: Path, default prlimit
        :param stringtie_bin: Software path.
        :type stringtie_bin: Path, default stringtie
        :param memory_limit: Memory limit Scallop command Defaults to 40*1024**3.
        :type memory_limit: int

        :return: None
        :rtype: None
    """
    check_exe(scallop_bin)
    check_exe(stringtie_bin)
    scallop_dir = create_dir(output_dir, "scallop_output")
    logging.info("Skip analysis if the gtf file already exists")
    output_file = scallop_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logging.info("Scallop gtf file exists, skipping analysis")
            return

    star_dir = Path(f"{output_dir}/star_output")

    if star_dir.exists() and len(list(star_dir.glob("*.bam"))) != 0:
        for sorted_bam_file in star_dir.glob("*.bam"):
            transcript_file_name = re.sub(".bam", ".scallop.gtf", sorted_bam_file.name)
            transcript_file = scallop_dir / transcript_file_name
            if transcript_file.exists():
                logging.info(
                    "Found an existing stringtie gtf file, will not overwrite. \
                        File found: %s",
                    transcript_file,
                )
            else:
                logging.info("Running Scallop on: %s", sorted_bam_file.name)
                try:
                    scallop_cmd = [
                        scallop_bin,
                        "-i",
                        sorted_bam_file,
                        "-o",
                        transcript_file,
                        "--min_flank_length",
                        "10",
                    ]
                    if memory_limit is not None:
                        scallop_cmd = _prlimit_command(prlimit_bin, scallop_cmd, memory_limit)
                    subprocess.check_output(scallop_cmd, stderr=subprocess.STDOUT, universal_newlines=True)
                    # This combines the standard output and error streams into a single
                    # string and ensures that the output is in text mode

                except subprocess.CalledProcessError as ex:
                    logging.error("Error occurred while running Scallop:")
                    logging.error("Command: %s\n", " ".join(scallop_cmd))
                    logging.error("Return code: %s\n", str(ex.returncode))
                    logging.error("Output and error messages: %s\n", ex.output)
    else:
        raise IndexError(f"The list of sorted bam files is empty, Star output dir: {star_dir}")

    # Now need to merge
    logging.info("Merge Scaalop's output.")
    _scallop_merge(scallop_dir, stringtie_bin)


def _scallop_merge(scallop_dir: Path, stringtie_bin: Path = Path("stringtie")) -> None:
    """
    Merge Scallop result in a single gtf file

    scallop_dir : Input directory's path.
    stringtie_bin : Software path.
    """
    scallop_input_to_file = scallop_dir / "scallop_assemblies.txt"
    scallop_merge_output_file = scallop_dir / "annotation.gtf"
    with open(scallop_input_to_file, "w+", encoding="utf8") as gtf_list_out:
        for gtf_file in scallop_dir.glob("*.scallop.gtf"):
            transcript_count = check_gtf_content(gtf_file, "transcript")
            if transcript_count > 0:
                gtf_list_out.write(str(gtf_file) + "\n")
            else:
                logging.warning("Warning, skipping file with no transcripts. Path:%s\n", gtf_file)

    try:
        subprocess.check_output(
            [
                stringtie_bin,
                "--merge",
                "-o",
                scallop_merge_output_file,
                scallop_input_to_file,
            ],
            stderr=subprocess.STDOUT,
            text=True,
        )

    except subprocess.CalledProcessError as e:
        print("StringTie execution failed with an error:%s", e.output)


def _prlimit_command(prlimit_bin, command_list, virtual_memory_limit) -> list:
    """
    Prepend memory limiting arguments to a command list to be run with subprocess.

    This method uses the `prlimit` program to set the memory limit.

    The `virtual_memory_limit` size is in bytes.

    prlimit arguments:
    -v, --as[=limits]
           Address space limit.
    """
    return [str(prlimit_bin), f"-v{virtual_memory_limit}"] + command_list


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Scallop's arguments")
    parser.add_argument("--output_dir", required=True, help="Output directory path")
    parser.add_argument("--scallop_bin", default="scallop", help="Scallop software path")
    parser.add_argument("--stringtie_bin", default="stringtie", help="StringTie software path")
    parser.add_argument("--prlimit_bin", default="prlimit", help="Prlimit software path")
    parser.add_argument(
        "--memory_limit", type=int, default=40 * 1024**3, help="Memory's limit for Scallop command"
    )
    return parser.parse_args()


def main():
    """Scallop's entry-point."""
    args = parse_args()

    log_file_path = create_dir(args.output_dir, "log") / "scallop.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"

    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )

    run_scallop(args.output_dir, args.scallop_bin, args.prlimit_bin, args.stringtie_bin, args.memory_limit)


if __name__ == "__main__":
    main()
