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
Red is the first repeat-detection tool capable of labeling its training data
and training itself automatically on an entire genome.
Girgis, H.Z. Red: an intelligent, rapid, accurate tool for detecting repeats
de-novo on the genomic scale. BMC Bioinformatics 16, 227 (2015).
https://doi.org/10.1186/s12859-015-0654-5
"""
__all__ = ["run_red"]

import logging
import logging.config
from os import PathLike
from pathlib import Path
import re
import subprocess
import argschema

from ensembl.tools.anno.utils._utils import (
    check_exe,
    create_dir,
)

logger = logging.getLogger(__name__)


def run_red(genome_file: Path, output_dir: Path, red_bin: Path = Path("Red"),) -> str:
    """
    Run Red on genome file
        :param genome_file: Genome file path.
        :type genome_file: Path
        :param output_dir: Working directory path.
        :type output_dir: Path
        :param red_bin: Red software path.
        :type red_bin: Path, default Red
        
        :return: Masked genome file
        :rtype: str
    """
    check_exe(red_bin)
    red_dir = create_dir(output_dir, "red_output")
    red_mask_dir = create_dir(red_dir, "mask_output")
    red_repeat_dir = create_dir(red_dir, "repeat_output")
    red_genome_dir = create_dir(red_dir, "genome_dir")

    sym_link_genome_cmd = "ln -s " + str(genome_file) + " " + str(red_genome_dir)
    genome_file_name = genome_file.name
    red_genome_file = red_genome_dir / genome_file_name
    genome_file_stem = genome_file.stem
    masked_genome_file = red_mask_dir / f"{genome_file_stem}.msk"
    repeat_coords_file = red_repeat_dir / f"{genome_file_stem}.rpt"
    output_file = red_dir / "annotation.gtf"

    if masked_genome_file.exists():
        logger.warning(
            "Masked Genome file already found on the path to the Red mask output dir. \
            Will not create a new file"
        )
        # _create_red_gtf(repeat_coords_file, output_file)
        return str(masked_genome_file)
    if red_genome_file.exists():
        logger.warning(
            "Unmasked genome file already found on the path to the Red genome dir, \
            will not create a sym link"
        )

    else:
        logger.info(
            "Preparing to sym link the genome file to the Red genome dir. Cmd\n %s",
            sym_link_genome_cmd,
        )
        # subprocess.run(["ln", "-s", genome_file, red_genome_dir])
        red_genome_file.symlink_to(genome_file)
    try:
        if red_genome_file.exists():
         logger.info("Running Red")
         subprocess.run(
            [
                red_bin,
                "-gnm",
                red_genome_dir,
                "-msk",
                red_mask_dir,
                "-rpt",
                red_repeat_dir,
            ],
            check=True,
        )
    except:
        logger.error(
            "Could not find the genome file in the Red genome dir or sym link \
            to the original file. Path expected:\n%s",
            genome_file,
        )
    _create_red_gtf(repeat_coords_file, output_file)
    return str(masked_genome_file)


def _create_red_gtf(repeat_coords_file: Path, output_file: Path):
    """
    Create Red gtf file from masked genome file

    Args:
        repeat_coords_file: Coordinates for repeats.
        output_file : GTF file with the final results.
    """
    with open(repeat_coords_file, "r", encoding="utf8") as red_in, open(
        output_file, "w+", encoding="utf8"
    ) as red_out:
        for repeat_id, line in enumerate(red_in, start=1):
            result_match = re.search(r"^\>(.+)\:(\d+)\-(\d+)", line)
            if result_match:
                region_name = result_match.group(1)
                # Note that Red is 0-based, so add 1
                start = int(result_match.group(2)) + 1
                end = int(result_match.group(3)) + 1
                gtf_line = (
                    f"{region_name}\tRed\trepeat\t{start}\t"
                    f'{end}\t.\t+\t.\trepeat_id "{repeat_id}";\n'
                )
                red_out.write(gtf_line)


class InputSchema(argschema.ArgSchema):
    """Input arguments expected to run Red."""

    genome_file = argschema.fields.InputFile(
        required=True, description="Genome file path"
    )
    output_dir = argschema.fields.OutputDir(
        required=True, description="Output directory path"
    )
    red_bin = argschema.fields.String(
        required=False, default="Red", description="Red executable path",
    )


def main() -> None:
    """Red's entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    log_file_path = create_dir(mod.args["output_dir"], "log") / "red.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"
    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )
    run_red(
        Path(mod.args["genome_file"]), mod.args["output_dir"], mod.args["red_bin"],
    )


if __name__ == "__main__":
    main()
