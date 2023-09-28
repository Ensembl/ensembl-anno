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
StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts.
It uses a novel network flow algorithm as well as an optional de novo assembly step to assemble and
quantitate full-length transcripts representing multiple splice variants for each gene locus.
Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT & Salzberg SL. StringTie enables improved 
reconstruction of a transcriptome from RNA-seq reads Nature Biotechnology 2015, doi:10.1038/nbt.3122
"""

__all__ = ["run_stringtie"]
import logging
import logging.config
from pathlib import Path
import re
import subprocess
import argschema

from ensembl.tools.anno.utils._utils import (
    check_exe,
    create_dir,
    check_gtf_content,
)


def run_stringtie(
    output_dir: Path,
    stringtie_bin: Path = Path("stringtie"),
    num_threads: int = 1,
) -> None:
    """
    StringTie assembler of short read data.
        :param output_dir: Working directory path.
        :type output_dir: Path
        :param stringtie_bin: Software path.
        :type stringtie_bin: Path, default stringtie
        :param num_threads: Number of available threads.
        :type num_threads: int, default 1
                        
        :return: None
        :rtype: None
    """
    check_exe(stringtie_bin)
    stringtie_dir = create_dir(output_dir, "stringtie_output")
    logging.info("Skip analysis if the gtf file already exists")
    output_file = stringtie_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logging.info("Stringtie gtf file exists, skipping analysis")
            return

    stringtie_merge_input_file = stringtie_dir / "stringtie_assemblies.txt"
    stringtie_merge_output_file = stringtie_dir / "annotation.gtf"
    star_dir = output_dir / "star_output"

    if star_dir.exists() and len(list(star_dir.glob("*.bam"))) != 0:
        for sorted_bam_file in star_dir.glob("*.bam"):
            transcript_file_name = re.sub(".bam", ".stringtie.gtf", sorted_bam_file.name)
            transcript_file = stringtie_dir / transcript_file_name
            if transcript_file.exists():
                logging.info(
                    "Found an existing stringtie gtf file, will not overwrite. \
                        File found: %s",
                    transcript_file,
                )
            else:
                logging.info("Running Stringtie on: %s", sorted_bam_file.name)
                try:
                    subprocess.check_output(  # pylint:disable=subprocess-run-check
                        [
                            stringtie_bin,
                            sorted_bam_file,
                            "-o",
                            transcript_file,
                            "-p",
                            str(num_threads),
                            "-t",  # disable trimming of predicted transcripts based on coverage
                            "-a",  # minimum anchor length for junctions
                            "15",
                        ]
                    )
                except subprocess.CalledProcessError as e:
                    logging.error("Error running Stringtie command: %s", e)
                    logging.error("Return code: %s", str(e.returncode))
                    logging.error("Output and error messages:%s\n", e.output)
    else:
        raise IndexError(f"The list of sorted bam files is empty, Star output dir: {star_dir}")

    logging.info("Creating Stringtie merge input file: %s", stringtie_merge_input_file)
    with open(stringtie_merge_input_file, "w+", encoding="utf8") as gtf_list_out:
        for gtf_file in stringtie_dir.glob("*.stringtie.gtf"):
            transcript_count = check_gtf_content(gtf_file, "transcript")
            if transcript_count > 0:
                gtf_list_out.write(f"{gtf_file}\n")
            else:
                logging.warning("Warning, skipping file with no transcripts. Path:%s", gtf_file)
    logging.info("Merging Stringtie results.")
    try:
        subprocess.run(  # pylint:disable=subprocess-run-check
            [
                stringtie_bin,
                "--merge",
                "-o",
                stringtie_merge_output_file,
                stringtie_merge_input_file,
            ]
        )
    except subprocess.CalledProcessError as e:
        logging.error("Error running Stringtie merging command: %s", e)


class InputSchema(argschema.ArgSchema):
    """Input arguments expected to run StringTie software."""

    output_dir = argschema.fields.OutputDir(required=True, description="Output directory path")
    stringtie_bin = argschema.fields.String(
        required=False,
        default="stringtie",
        description="StringTie software path",
    )
    num_threads = argschema.fields.Integer(required=False, default=1, description="Number of threads")


def main() -> None:
    """StringTie's entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    log_file_path = create_dir(mod.args["output_dir"], "log") / "stringtie.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"
    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )
    run_stringtie(
        mod.args["output_dir"],
        mod.args["stringtie_bin"],
        mod.args["num_threads"],
    )
