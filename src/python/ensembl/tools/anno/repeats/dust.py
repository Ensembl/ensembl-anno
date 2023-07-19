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
DustMasker is a program that identifies and masks out low complexity
parts of a genome using a new and improved DUST algorithm.

Morgulis A, Gertz EM, Schaffer AA, Agarwala R. A Fast and Symmetric
DUST Implementation to Mask Low-Complexity DNA Sequences.
"""
__all__ = ["run_dust"]

import logging
import logging.config
import multiprocessing
import os
from os import PathLike
from pathlib import Path
import re
import subprocess
import tempfile
from typing import List
import argschema


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


def run_dust(
    genome_file: PathLike,
    output_dir: Path,
    dust_bin: Path = Path("dustmasker"),
    num_threads: int = 1,
) -> None:
    """
    Run Dust on genomic slices with mutiprocessing
    Args:
        genome_file : Genome file path.
        output_dir : Working directory path.
        dust_bin : Dust software path.
        num_threads: Number of threads.
    """

    check_exe(dust_bin)
    dust_dir = create_dir(output_dir, "dust_output")
    os.chdir(str(dust_dir))
    output_file = dust_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "repeat")
        if transcript_count > 0:
            logger.info("Dust gtf file exists, skipping analysis")
            return
    logger.info("Creating list of genomic slices")
    seq_region_to_length = get_seq_region_length(genome_file, 5000)
    slice_ids_per_region = get_slice_id(
        seq_region_to_length, slice_size=1000000, overlap=0, min_length=5000
    )
    dust_cmd = [dust_bin, "-in"]
    pool = multiprocessing.Pool(num_threads)  # pylint: disable=consider-using-with
    for slice_id in slice_ids_per_region:
        pool.apply_async(
            _multiprocess_dust,
            args=(
                dust_cmd,
                slice_id,
                dust_dir,
                genome_file,
            ),
        )
    pool.close()
    pool.join()
    slice_output_to_gtf(dust_dir, "repeat_id", "dust", True, ".dust.gtf")
    for gtf_file in dust_dir.glob("*.dust.gtf"):
        gtf_file.unlink()


def _multiprocess_dust(  # pylint: disable=too-many-locals
    dust_cmd: List[str],
    slice_id: List[str],
    dust_dir: Path,
    genome_file: Path,
) -> None:
    """
    Run Dust on multiprocess on genomic slices
    Args:
        dust_cmd: Dust command to execute.
        slice_id: List of slice IDs.
        dust_dir : Dust output directory path.
        genome_file : Genome file.
    """
    region_name, start, end = slice_id
    logger.info(
        "Processing slice to find low complexity regions with Dust: %s:%s:%s",
        region_name,
        start,
        end,
    )
    seq = get_sequence(region_name, int(start), int(end), 1, genome_file, dust_dir)
    slice_name = f"{region_name}.rs{start}.re{end}"
    with tempfile.TemporaryDirectory(dir=dust_dir) as tmpdirname:
        slice_file = dust_dir / tmpdirname / f"{slice_name}.fa"
        with open(slice_file, "w+", encoding="utf8") as region_out:
            region_out.write(f">{region_name}\n{seq}\n")
        region_results = dust_dir / f"{slice_name}.dust.gtf"
        output_file = Path(f"{slice_file}.dust")
        dust_cmd = dust_cmd.copy()
        logger.info("dust_cmdPRIMA: %s", dust_cmd)
        dust_cmd.append(str(slice_file))
        logger.info("dust_cmd: %s", dust_cmd)
        with open(output_file, "w+", encoding="utf8") as dust_out:
            subprocess.run(dust_cmd, stdout=dust_out, check=True)
        _create_dust_gtf(output_file, region_results, region_name)
        slice_file.unlink()
        region_results.unlink()


def _create_dust_gtf(
    output_file: Path,
    region_results: Path,
    region_name: str,
) -> None:
    """
    Read the fasta file and save the content in gtf format
    All the genomic slices are collected in a single gtf output
    Args:
        output_file : GTF file with final results.
        region_results : GTF file with the results per region.
        region_name :Coordinates of genomic slice.
    """
    with open(output_file, "r", encoding="utf8") as dust_in, open(
        region_results, "w+", encoding="utf8"
    ) as dust_out:
        repeat_count = 1
        for line in dust_in:
            result_match = re.search(r"(\d+)\ - (\d+)", line)
            if result_match:
                start = int(result_match.group(1)) + 1
                end = int(result_match.group(2)) + 1
                gtf_line = (
                    f"{region_name}\tDust\trepeat\t{start}\t"
                    f'{end}\t.\t+\t.\trepeat_id "{repeat_count}";\n'
                )
                dust_out.write(gtf_line)
                repeat_count += 1


class InputSchema(argschema.ArgSchema):
    """Input arguments expected to run DustMasker."""

    genome_file = argschema.fields.InputFile(
        required=True, description="Genome file path"
    )
    output_dir = argschema.fields.OutputDir(
        required=True, description="Output directory path"
    )
    dust_bin = argschema.fields.String(
        required=False,
        default="dustmasker",
        description="Dust executable path",
    )
    num_threads = argschema.fields.Integer(
        required=False, default=1, description="Number of threads"
    )


def main() -> None:
    """Dust's entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    log_file_path = create_dir(mod.args["output_dir"], "log") / "dust.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"
    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )
    run_dust(
        mod.args["genome_file"],
        mod.args["output_dir"],
        mod.args["dust_bin"],
        mod.args["num_threads"],
    )


if __name__ == "__main__":
    main()
