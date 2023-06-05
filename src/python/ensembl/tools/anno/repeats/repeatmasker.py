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
    RepeatMasker is a program that screens DNA sequences for interspersed
    repeats and low complexity DNA sequences.
    Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-4.0
"""

__all__ = ["run_repeatmasker"]

import json
import logging
import logging.config
import multiprocessing
import os
from os import PathLike
from pathlib import Path
import re
import subprocess
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
logger = logging.getLogger('__name__')


def run_repeatmasker(
    genome_file: PathLike,
    output_dir: Path,
    repeatmasker_bin: Path = Path("RepeatMasker"),
    library: str = "",
    repeatmasker_engine: str = "rmblast",
    species: str = "",
    num_threads: int = 1,
) -> None:

    """
    Executes RepeatMasker on the genome slices and stores the final annotation.gtf in repeatmasker_output
    Args:
        genome_file : Genome file path.
        repeatmasker_path : RepeatMasker executable path.
        library : Custom repeat library.
        species :Species name.
        output_dir : Output directory path.
        num_threads: Number of threads.

    """
    check_exe(repeatmasker_bin)
    repeatmasker_dir = create_dir(output_dir, "repeatmasker_output")

    output_file = repeatmasker_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "repeat")
        if transcript_count > 0:
            logger.info("Repeatmasker gtf file exists")
            return
    logger.info("Creating list of genomic slices")
    seq_region_to_length = get_seq_region_length(genome_file, 5000)
    slice_ids_per_region = get_slice_id(
        seq_region_to_length, slice_size=1000000, overlap=0, min_length=5000
    )
    repeatmasker_cmd = [
        str(repeatmasker_bin),
        "-nolow",
        "-engine",
        repeatmasker_engine,
        "-dir",
        str(repeatmasker_dir),
    ]
    if not library:
        if not species:
            species = "homo"
        repeatmasker_cmd.extend(["-species", species])
    else:
        repeatmasker_cmd.extend(["-lib", library])
    logger.info(f"Running RepeatMasker {repeatmasker_cmd}")
    pool = multiprocessing.Pool(num_threads)  # pylint: disable=consider-using-with
    for slice_id in slice_ids_per_region:
        pool.apply_async(
            _multiprocess_repeatmasker,
            args=(
                repeatmasker_cmd,
                slice_id,
                genome_file,
                repeatmasker_dir,
            ),
        )
    pool.close()
    pool.join()
    slice_output_to_gtf(repeatmasker_dir, "repeat_id", "repeatmask", True, ".rm.gtf")
    for gtf_file in repeatmasker_dir.glob("*.rm.gtf"):
        gtf_file.unlink()

def _multiprocess_repeatmasker(  # pylint: disable=too-many-locals
    repeatmasker_cmd: List[str],
    slice_id: List[str],
    genome_file: Path,
    repeatmasker_output_dir: Path,
) -> None:
    """
    Run Repeatmasker on genomic slice

    Args:
        generic_repeatmasker_cmd: RepeatMasker command to execute.
        slice_id: Slice ID to run RepeatMasker on.
        genome_file : Genome file path.
        repeatmasker_output_dir : RepeatMasker output directory path.
    """

    region_name, start, end = slice_id
    logger.info(
        "Processing slice to find repeats with RepeatMasker: %s:%s:%s",
        region_name,
        start,
        end,
    )
    seq = get_sequence(
        region_name, int(start), int(end), 1, genome_file, repeatmasker_output_dir
    )
    slice_file_name = f"{region_name}.rs{start}.re{end}"
    region_file = repeatmasker_output_dir / f"{slice_file_name}.fa"
    with open(region_file, "w+", encoding="utf8") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")
    region_results_file = Path(f"{region_file}.rm.gtf")
    output_file = Path(f"{region_file}.out")
    masked_file = Path(f"{region_file}.masked")
    tbl_file = Path(f"{region_file}.tbl")
    log_file = Path(f"{region_file}.log")
    cat_file = Path(f"{region_file}.cat")
    repeatmasker_cmd = repeatmasker_cmd.copy()
    repeatmasker_cmd.append(str(region_file))
    logger.info(repeatmasker_cmd)
    subprocess.run(repeatmasker_cmd, check=True)
    _create_repeatmasker_gtf(output_file, region_results_file, region_name)
    output_file.unlink()
    region_file.unlink()
    masked_file.unlink(missing_ok=True)
    tbl_file.unlink(missing_ok=True)
    log_file.unlink(missing_ok=True)
    cat_file.unlink(missing_ok=True)


def _create_repeatmasker_gtf(  # pylint: disable=too-many-locals
    output_file: Path,
    region_results_file: Path,
    region_name: str,
) -> None:
    """
    Read the fasta file and save the content in gtf format

    All the genomic slices are collected in a single gtf output with the following format:
    SW    perc perc perc query    position in query matching repeat       position in repeat
    score div. del. ins. sequence begin end (left)  repeat   class/family begin end  (left)  ID
    Args:
        output_file : GTF file with final results.
        region_results_file_path : GTF file with results per region.
        region_name : Coordinates of genomic slice.
    """
    with open(output_file, "r", encoding="utf8") as repeatmasker_in, open(
        region_results_file, "w+", encoding="utf8"
    ) as repeatmasker_out:
        repeat_count = 1
        for line in repeatmasker_in:
            result_match = re.search(r"^\s*\d+\s+", line)
            if result_match:
                results = line.split()
                if results[-1] == "*":
                    results.pop()
                if len(results) != 15:
                    continue
                score = results[0]
                start = results[5]
                end = results[6]
                strand = results[8]
                repeat_name = results[9]
                repeat_class = results[10]
                if strand == "+":
                    repeat_start = results[11]
                    repeat_end = results[12]
                else:
                    repeat_start = results[13]
                    repeat_end = results[12]
                    strand = "-"
                gtf_line = (
                    f"{region_name}\tRepeatMasker\trepeat\t{start}\t{end}\t.\t"
                    f"{strand}\t.\trepeat_id{repeat_count}; "
                    f'repeat_name "{repeat_name}"; repeat_class "{repeat_class}"; '
                    f'repeat_start "{repeat_start}"; '
                    f'repeat_end "{repeat_end}"; score "{score}";\n'
                )
                repeatmasker_out.write(gtf_line)
                repeat_count += 1


class InputSchema(argschema.ArgSchema):
    """Input arguments expected to run RepeatMasker."""

    genome_file = argschema.fields.InputFile(
        required= True, description= "Genome file path"
    )
    output_dir = argschema.fields.OutputDir(
        required= True, description= "Output directory path"
    )
    repeatmasker_bin = argschema.fields.String(
        required= False, default= "RepeatMasker",
            description = "RepeatMasker executable path",

    )
    library = argschema.fields.String(
            required= False, default= "", description= "Custom repeat library"
    )
    repeatmasker_engine = argschema.fields.String(
            required= False, default= "rmblast", description= "RepeatMasker engine"
    )
    species = argschema.fields.String(
            required= False,
            default="homo",
            description="Species name (used if no library is provided)"
    )
    num_threads = argschema.fields.Integer(
            required= False, default= 1, description= "Number of threads"
    )


def main() -> None:
    """RepeatMasker's entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    log_file_path = create_dir(mod.args["output_dir"], "log") /"repeatmasking.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"
    logging.config.fileConfig(loginipath, defaults={"logfilename": str(log_file_path)}, disable_existing_loggers=False,)
    run_repeatmasker(
        mod.args["genome_file"],
        mod.args["output_dir"],
        mod.args["repeatmasker_bin"],
        mod.args["library"],
        mod.args["repeatmasker_engine"],
        mod.args["species"],
        mod.args["num_threads"],
        )

if __name__ == "__main__":
    main()

