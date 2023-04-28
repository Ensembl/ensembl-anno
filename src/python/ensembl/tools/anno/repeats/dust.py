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
import json
import logging
import multiprocessing
import os
import pathlib
import re
import subprocess
import tempfile
import typing

__all__ = ["run_red", "run_dust_regions", "run_trf_repeats", "run_repeatmasker_regions"]

from utils._utils import (
    check_exe,
    create_dir,
    check_gtf_content,
    get_seq_region_lengths,
    create_slice_ids,
    slice_output_to_gtf,
    get_sequence,
)

logger = logging.getLogger(__name__)
with open(os.environ["ENSCODE"] + "/ensembl-anno/config.json", "r") as f:
    config = json.load(f)

__all__ = ["run_repeatmasker_regions"]

def run_dust_regions(
    genome_file: typing.Union[pathlib.Path, str],
    dust_path: str,
    main_output_dir: str,
    num_threads: int,
):
    """
    Run Dust on genomic slices

    Args:
        genome_file : pathlib.Path
        dust_path : str
        main_output_dir : pathlib.Path
        num_threads: int

    Return:
        gtfs with the masked sequence for each genome slice

    """
    if not dust_path:
        dust_path = config["dust"]["software"]

    check_exe(dust_path)
    dust_output_dir = pathlib.Path(create_dir(main_output_dir, "dust_output"))
    os.chdir(str(dust_output_dir))
    output_file = dust_output_dir / "annotation.gtf"
    logger.info("dust output %s", output_file)
    if output_file.is_file():
        transcript_count = check_gtf_content(output_file, "repeat")
        if transcript_count > 0:
            logger.info("Dust gtf file exists")
            return 0

    logger.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, 5000)
    slice_ids = create_slice_ids(
        seq_region_lengths, slice_size=1000000, overlap=0, min_length=5000
    )
    generic_dust_cmd = [dust_path, "-in"]
    logger.info("Running Dust")
    pool = multiprocessing.Pool(int(num_threads))
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_dust,
            args=(
                generic_dust_cmd,
                slice_id,
                genome_file,
                dust_output_dir,
            ),
        )

    pool.close()
    pool.join()
    slice_output_to_gtf(str(dust_output_dir), ".dust.gtf", 1, "repeat_id", "dust")
    for gtf_file in pathlib.Path(dust_output_dir).glob("*.dust.gtf"):
        gtf_file.unlink()


def multiprocess_dust(  # pylint: disable=too-many-locals
    generic_dust_cmd: list,
    slice_id: str,
    genome_file: pathlib.Path,
    dust_output_dir: pathlib.Path,
):
    """
    Run Dust on multiprocess on genomic slices
    Args:
        generic_dust_cmd:list
        slice_id:str
        genome_file : pathlib.Path
        dust_output_dir : pathlib.Path
    """
    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing slice to find low complexity regions with Dust: %s:%s:%s",
        region_name,
        start,
        end,
    )
    seq = get_sequence(region_name, start, end, 1, genome_file, str(dust_output_dir))

    slice_file_name = f"{region_name}.rs{start}.re{end}"
    with tempfile.TemporaryDirectory(dir=dust_output_dir) as tmpdirname:
        region_fasta_file_path = dust_output_dir / tmpdirname / f"{slice_file_name}.fa"
        with open(region_fasta_file_path, "w+") as region_fasta_out:
            region_fasta_out.write(f">{region_name}\n{seq}\n")

        region_results_file_path = dust_output_dir / f"{slice_file_name}.dust.gtf"

        dust_output_file_path = pathlib.Path(f"{region_fasta_file_path}.dust")
        dust_cmd = generic_dust_cmd.copy()
        dust_cmd.append(region_fasta_file_path)
        logger.info(dust_cmd)
        with open(dust_output_file_path, "w+") as dust_out:
            subprocess.run(dust_cmd, stdout=dust_out, check=True)

        create_dust_gtf(dust_output_file_path, region_results_file_path, region_name)

def create_dust_gtf(
    dust_output_file_path: pathlib.Path,
    region_results_file_path: pathlib.Path,
    region_name: str,
):
    """
    Read the fasta file and save the content in gtf format

    All the genomic slices are collected in a single gtf output
    Args:
        dust_output_file_path : pathlib.Path
        region_results_file_path : pathlib.Path
        region_name :str
    """
    with open(dust_output_file_path, "r") as dust_in, open(
        region_results_file_path, "w+"
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
