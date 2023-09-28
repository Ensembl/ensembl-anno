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
    Tandem Repeats Finder is a program to locate and display tandem repeats in DNA sequences.
    Benson G. Tandem repeats finder: a program to analyze DNA sequences.
    Nucleic Acids Res. 1999; 27(2):573–580. doi:10.1093/nar/27.2.573
"""
__all__ = ["run_trf"]

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


def run_trf(
    genome_file: PathLike,
    output_dir: Path,
    num_threads: int = 1,
    trf_bin: Path = Path("trf"),
    match_score: int = 2,
    mismatch_score: int = 5,
    delta: int = 7,
    pm: int = 80,
    pi: int = 10,
    minscore: int = 40,
    maxperiod: int = 500,
) -> None:
    """
    Executes TRF on genomic slices
            :param genome_file: Genome file path.
            :type genome_file: PathLike
            :param output_dir:  Working directory path.
            :type output_dir: Path
            :param num_threads: int, number of threads.
            :type num_threads: int, default 1
            :param trf_bin: TRF software path.
            :type trf_bin: Path, default trf
            :param match_score: Matching weight.
            :type match_score: int, default 2
            :param mismatch_score: Mismatching penalty.
            :type mismatch_score: int, default 5
            :param delta: Indel penalty.
            :type delta: int, default 7
            :param pm: Match probability (whole number).
            :type pm: int, default 80
            :param pi: Indel probability (whole number).
            :type pi: int, default 10
            :param minscore: Minimum alignment score to report.
            :type minscore: int, default 40
            :param maxperiod: Maximum period size to report.
            :type maxperiod: int, default 500
                    
            :return: None
            :rtype: None
    """
    check_exe(trf_bin)
    trf_dir = create_dir(output_dir, "trf_output")
    os.chdir(str(trf_dir))
    output_file = trf_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "repeat")
        if transcript_count > 0:
            logger.info("Trf gtf file exists, skipping analysis")
            return
    logger.info("Creating list of genomic slices")
    seq_region_to_length = get_seq_region_length(genome_file, 5000)
    slice_ids_per_region = get_slice_id(
        seq_region_to_length, slice_size=1000000, overlap=0, min_length=5000
    )
    trf_output_extension = (
        f".{match_score}.{mismatch_score}.{delta}."
        f"{pm}.{pi}.{minscore}.{maxperiod}.dat"
    )
    trf_cmd = [
        trf_bin,
        None,
        str(match_score),
        str(mismatch_score),
        str(delta),
        str(pm),
        str(pi),
        str(minscore),
        str(maxperiod),
        "-d",
        "-h",
    ]
    logger.info("Running TRF")
    pool = multiprocessing.Pool(num_threads)#pylint:disable=consider-using-with
    for slice_id in slice_ids_per_region:
        pool.apply_async(
            _multiprocess_trf,
            args=(
                trf_cmd,
                slice_id,
                trf_dir,
                trf_output_extension,
                genome_file,
            ),
        )
    pool.close()
    pool.join()
    slice_output_to_gtf(trf_dir, "repeat_id", "trf", True, ".trf.gtf")
    for gtf_file in trf_dir.glob("*.trf.gtf"):
        gtf_file.unlink()


def _multiprocess_trf(
    trf_cmd: List[str],
    slice_id: List[str],
    trf_dir: Path,
    trf_output_extension: Path,
    genome_file:Path,
) -> None:
    """
    Run TRF on multiprocess on genomic slices
    Args:
        trf_cmd: TRF command to execute.
        slice_id: Slice Id to run TRF on.
        trf_dir : TRF output dir.
        trf_output_extension: TRF file output extension.
        genome_file : Genome file.
    """
    region_name, start, end = slice_id
    logger.info(
        "Processing slice to find tandem repeats with TRF:%s:%s:%s",
        region_name,
        start,
        end,
    )
    seq = get_sequence(region_name, int(start), int(end), 1, genome_file, trf_dir)
    slice_name = f"{region_name}.rs{start}.re{end}"
    with tempfile.TemporaryDirectory(dir=trf_dir) as tmpdirname:
        slice_file = trf_dir / tmpdirname / f"{slice_name}.fa"
        with open(slice_file, "w+", encoding="utf8") as region_out:
            region_out.write(f">{region_name}\n{seq}\n")
        region_results = trf_dir / f"{slice_name}.trf.gtf"
        # TRF writes to the current dir, so swtich to the output dir for it
        # os.chdir(str(trf_output_dir))
        output_file = Path(f"{slice_file}{trf_output_extension}")
        trf_cmd = trf_cmd.copy()
        trf_cmd[1] = str(slice_file)
        logger.info("trf_cmd: %s", trf_cmd)
        # with open(trf_output_file_path, "w+") as trf_out:
        subprocess.run(trf_cmd, cwd=trf_dir / tmpdirname)#pylint:disable=subprocess-run-check
        _create_trf_gtf(output_file, region_results, region_name)
        slice_file.unlink()
        output_file.unlink()


def _create_trf_gtf(
    output_file: Path,
    region_results: Path,
    region_name: str,
) -> None:
    """
    Read the fasta file and save the content in gtf format

    TRF output format:
    cols 1+2:  Indices of the repeat relative to the start of the sequence
    col 3:     Period size of the repeat
    col 4:     Number of copies aligned with the consensus pattern
    col 5:     Size of consensus pattern (may differ slightly from the period size)
    col 6:     Percent of matches between adjacent copies overall
    col 7:     Percent of indels between adjacent copies overall
    col 8:     Alignment score
    cols 9-12: Percent composition for each of the four nucleotides
    col 13:    Entropy measure based on percent composition
    col 14:    Consensus sequence
    col 15:    Repeat sequence
    Args:
       output_file : GTF file with final results.
       region_results : GTF file with results per region.
       region_name : Coordinates of genomic slice.
    """
    with open(output_file, "r", encoding="utf8") as trf_in, open(
        region_results, "w+", encoding="utf8"
    ) as trf_out:
        repeat_count = 1
        for line in trf_in:
            result_match = re.search(r"^\d+", line)
            if result_match:
                results = line.split()
                if len(results) != 15:
                    continue
                start = results[0]
                end = results[1]
                period = float(results[2])
                copy_number = float(results[3])
                percent_matches = float(results[5])
                score = float(results[7])
                repeat_consensus = results[13]
                if (  # pylint: disable=too-many-boolean-expressions
                    score < 50
                    and percent_matches >= 80
                    and copy_number > 2
                    and period < 10
                ) or (copy_number >= 2 and percent_matches >= 70 and score >= 50):
                    gtf_line = (
                        f"{region_name}\tTRF\trepeat\t{start}\t{end}\t.\t+\t.\t"
                        f'repeat_id "{repeat_count}"; score "{score}"; '
                        f'repeat_consensus "{repeat_consensus}";\n'
                    )
                    trf_out.write(gtf_line)
                    repeat_count += 1


class InputSchema(argschema.ArgSchema):
    """Input arguments expected to run TRF."""

    genome_file = argschema.fields.InputFile(
        required=True, description="Genome file path"
    )
    output_dir = argschema.fields.OutputDir(
        required=True, description="Output directory path"
    )
    trf_bin = argschema.fields.String(
        required=False,
        default="trf",
        description="TRF executable path",
    )
    match_score = argschema.fields.Integer(
        required=False, default=2, description="Matching weight"
    )
    mismatch_score = argschema.fields.Integer(
        required=False, default=5, description="Mismatching penalty"
    )
    delta = argschema.fields.Integer(
        required=False, default=7, description="Indel penalty"
    )
    pm = argschema.fields.Integer(
        required=False, default=80, description="Match probability"
    )
    pi = argschema.fields.Integer(
        required=False, default=10, description="Indel probability"
    )
    minscore = argschema.fields.Integer(
        required=False, default=40, description="Minimum alignment score to report"
    )
    maxperiod = argschema.fields.Integer(
        required=False, default=500, description="Maximum period size to report"
    )
    num_threads = argschema.fields.Integer(
        required=False, default=1, description="Number of threads"
    )


def main() -> None:
    """TRF's entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    log_file_path = create_dir(mod.args["output_dir"], "log") / "trf.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"
    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )
    run_trf(
        mod.args["genome_file"],
        mod.args["output_dir"],
        mod.args["num_threads"],
        mod.args["trf_bin"],
        mod.args["match_score"],
        mod.args["mismatch_score"],
        mod.args["delta"],
        mod.args["pm"],
        mod.args["pi"],
        mod.args["minscore"],
        mod.args["maxperiod"],
    )


if __name__ == "__main__":
    main()
