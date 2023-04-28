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

__all__ = ["Dust"]

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


class Dust:
    """
    DustMasker is a program that identifies and masks out low complexity
    parts of a genome using a new and improved DUST algorithm.

    Morgulis A, Gertz EM, Schaffer AA, Agarwala R. A Fast and Symmetric
    DUST Implementation to Mask Low-Complexity DNA Sequences.

    Args:
            genome_file : pathlib.Path, genome file path.
            dust_path : str, software path.
            main_output_dir : pathlib.Path, working directory path.
            num_threads: int, number of threads.

    Return:
            gtfs with the masked sequence for each genome slice

    """

    def __init__(  # pylint: disable=too-many-arguments
        self,
        genome_file: typing.Union[pathlib.Path, str],
        dust_path: str,
        main_output_dir: str,
        num_threads: int,
    ) -> None:
        self.genome_file = genome_file
        self.dust_path = dust_path
        self.main_output_dir = main_output_dir
        self.num_threads = num_threads

    def run_dust_regions(self) -> int:
        """Executes Dust on genomic slices"""
        if not self.dust_path:
            self.dust_path = config["dust"]["software"]
        check_exe(self.dust_path)
        dust_output_dir = pathlib.Path(create_dir(self.main_output_dir, "dust_output"))
        os.chdir(str(dust_output_dir))
        output_file = dust_output_dir / "annotation.gtf"
        logger.info("dust output %s", output_file)
        if output_file.is_file():
            transcript_count = check_gtf_content(output_file, "repeat")
            if transcript_count > 0:
                logger.info("Dust gtf file exists")
                return 0
        logger.info("Creating list of genomic slices")
        seq_region_lengths = get_seq_region_lengths(self.genome_file, 5000)
        slice_ids = create_slice_ids(
            seq_region_lengths, slice_size=1000000, overlap=0, min_length=5000
        )
        generic_dust_cmd = [self.dust_path, "-in"]
        logger.info("Running Dust")
        pool = multiprocessing.Pool(
            int(self.num_threads)
        )  # pylint: disable=consider-using-with
        for slice_id in slice_ids:
            pool.apply_async(
                self._multiprocess_dust,
                args=(
                    generic_dust_cmd,
                    slice_id,
                    dust_output_dir,
                ),
            )
        pool.close()
        pool.join()
        slice_output_to_gtf(str(dust_output_dir), ".dust.gtf", 1, "repeat_id", "dust")
        for gtf_file in pathlib.Path(dust_output_dir).glob("*.dust.gtf"):
            gtf_file.unlink()
        return 0

    def _multiprocess_dust(  # pylint: disable=too-many-locals
        self,
        generic_dust_cmd: list,
        slice_id: list,
        dust_output_dir: pathlib.Path,
    ) -> None:
        """
        Run Dust on multiprocess on genomic slices
        Args:
            generic_dust_cmd:list
            slice_id:list
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
        seq = get_sequence(
            region_name, start, end, 1, self.genome_file, str(dust_output_dir)
        )
        slice_file_name = f"{region_name}.rs{start}.re{end}"
        with tempfile.TemporaryDirectory(dir=dust_output_dir) as tmpdirname:
            region_fasta_file_path = (
                dust_output_dir / tmpdirname / f"{slice_file_name}.fa"
            )
            with open(region_fasta_file_path, "w+", encoding="utf8") as region_fasta_out:
                region_fasta_out.write(f">{region_name}\n{seq}\n")
            region_results_file_path = dust_output_dir / f"{slice_file_name}.dust.gtf"
            dust_output_file_path = pathlib.Path(f"{region_fasta_file_path}.dust")
            dust_cmd = generic_dust_cmd.copy()
            dust_cmd.append(region_fasta_file_path)
            logger.info(dust_cmd)
            with open(dust_output_file_path, "w+", encoding="utf8") as dust_out:
                subprocess.run(dust_cmd, stdout=dust_out, check=True)
            self._create_dust_gtf(
                dust_output_file_path, region_results_file_path, region_name
            )

    def _create_dust_gtf(
        self,
        dust_output_file_path: pathlib.Path,
        region_results_file_path: pathlib.Path,
        region_name: str,
    ) -> None:
        """
        Read the fasta file and save the content in gtf format

        All the genomic slices are collected in a single gtf output
        Args:
            dust_output_file_path : pathlib.Path
            region_results_file_path : pathlib.Path
            region_name :str
        """
        with open(dust_output_file_path, "r", encoding="utf8") as dust_in, open(
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
