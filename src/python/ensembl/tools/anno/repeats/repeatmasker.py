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
import typing

__all__ = ["RepeatMasker"]

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
with open(os.environ["ENSCODE"] + "/ensembl-anno/config.json", "r", encoding="utf8") as f:
    config = json.load(f)

__all__ = ["run_repeatmasker_regions"]


class RepeatMasker: # pylint: disable=too-few-public-methods
    """
    RepeatMasker is a program that screens DNA sequences for interspersed repeats and low complexity DNA sequences.
    Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-4.0

    Args:
        genome_file : pathlib.Path, genome file path.
        repeatmasker_path : str, RepeatMasker executable path.
        library : str, custom repeat library.
        species :str, species name.
        main_output_dir : pathlib.Path, working directory path.
        num_threads: int, number of threads.

    Return:
        A GTF file with the repeatmasked sequence for each genome slice

    """

    def __init__( #pylint: disable=too-many-arguments
        self,
        genome_file: typing.Union[pathlib.Path, str],
        repeatmasker_path: str,
        library: str,
        species: str,
        main_output_dir: str,
        num_threads: int,
    ) -> None:
        self.genome_file = genome_file
        self.repeatmasker_path = repeatmasker_path
        self.library = library
        self.species = species
        self.main_output_dir = main_output_dir
        self.num_threads = num_threads

    def run_repeatmasker_regions(self) -> None:
        """Executes RepeatMasker on the genome slices"""
        if not self.repeatmasker_path:
            self.repeatmasker_path = config["repeatmasker"]["software"]

        check_exe(self.repeatmasker_path)
        repeatmasker_output_dir = pathlib.Path(
            create_dir(self.main_output_dir, "repeatmasker_output")
        )

        output_file = repeatmasker_output_dir / "annotation.gtf"
        if output_file.exists():
            transcript_count = check_gtf_content(output_file, "repeat")
            if transcript_count > 0:
                logger.info("Repeatmasker gtf file exists")
                return

        logger.info("Creating list of genomic slices")
        seq_region_lengths = get_seq_region_lengths(self.genome_file, 5000)
        slice_ids = create_slice_ids(
            seq_region_lengths, slice_size=1000000, overlap=0, min_length=5000
        )
        generic_repeatmasker_cmd = [
            self.repeatmasker_path,
            "-nolow",
            "-engine",
            config["repeatmasker"]["engine"],
            "-dir",
            repeatmasker_output_dir,
        ]
        if not self.library:
            if not self.species:
                self.species = "homo"
                generic_repeatmasker_cmd.extend(["-species", self.species])

            else:
                generic_repeatmasker_cmd.extend(["-species", self.species])
        else:
            generic_repeatmasker_cmd.extend(["-lib", self.library])
        logger.info("Running RepeatMasker")
        pool = multiprocessing.Pool(self.num_threads) # pylint: disable=consider-using-with
        for slice_id in slice_ids:
            pool.apply_async(
                self._multiprocess_repeatmasker,
                args=(
                    generic_repeatmasker_cmd,
                    slice_id,
                    repeatmasker_output_dir,
                ),
            )

        pool.close()
        pool.join()
        slice_output_to_gtf(
            str(repeatmasker_output_dir), ".rm.gtf", 1, "repeat_id", "repeatmask"
        )
        for gtf_file in pathlib.Path(repeatmasker_output_dir).glob("*.rm.gtf"):
            gtf_file.unlink()

    def _multiprocess_repeatmasker(  # pylint: disable=too-many-locals
        self,
        generic_repeatmasker_cmd: list,
        slice_id: list,
        repeatmasker_output_dir: pathlib.Path,
    )->None:
        """
        Run Repeatmasker on multiprocess on genomic slices

        Args:
            generic_repeatmasker_cmd: list
            slice_id: list
            genome_file : pathlib.Path
            repeatmasker_output_dir : pathlib.Path
        """

        region_name = slice_id[0]
        start = slice_id[1]
        end = slice_id[2]
        logger.info(
            "Processing slice to find repeats with RepeatMasker: %s:%s:%s",
            region_name,
            start,
            end,
        )
        seq = get_sequence(
            region_name, start, end, 1, self.genome_file, str(repeatmasker_output_dir)
        )
        slice_file_name = f"{region_name}.rs{start}.re{end}"
        region_fasta_file_path = repeatmasker_output_dir / f"{slice_file_name}.fa"
        with open(region_fasta_file_path, "w+", encoding="utf8") as region_fasta_out:
            region_fasta_out.write(f">{region_name}\n{seq}\n")
        region_results_file_path = pathlib.Path(f"{region_fasta_file_path}.rm.gtf")
        repeatmasker_output_file_path = pathlib.Path(f"{region_fasta_file_path}.out")
        repeatmasker_masked_file_path = pathlib.Path(f"{region_fasta_file_path}.masked")
        repeatmasker_tbl_file_path = pathlib.Path(f"{region_fasta_file_path}.tbl")
        repeatmasker_log_file_path = pathlib.Path(f"{region_fasta_file_path}.log")
        repeatmasker_cat_file_path = pathlib.Path(f"{region_fasta_file_path}.cat")
        repeatmasker_cmd = generic_repeatmasker_cmd.copy()
        repeatmasker_cmd.append(region_fasta_file_path)
        logger.info(repeatmasker_cmd)
        subprocess.run(repeatmasker_cmd, check=True)
        self._create_repeatmasker_gtf(
            repeatmasker_output_file_path, region_results_file_path, region_name
        )
        repeatmasker_output_file_path.unlink()
        region_fasta_file_path.unlink()
        # if region_results_file_path.exists():
        #    region_results_file_path.unlink()
        if repeatmasker_masked_file_path.exists():
            repeatmasker_masked_file_path.unlink()
        if repeatmasker_tbl_file_path.exists():
            repeatmasker_tbl_file_path.unlink()
        if repeatmasker_log_file_path.exists():
            repeatmasker_log_file_path.unlink()
        if repeatmasker_cat_file_path.exists():
            repeatmasker_cat_file_path.unlink()
        reveal_locals()
    def _create_repeatmasker_gtf(  # pylint: disable=too-many-locals
        self,
        repeatmasker_output_file_path: pathlib.Path,
        region_results_file_path: pathlib.Path,
        region_name: str,
    )->None:
        """
        Read the fasta file and save the content in gtf format

        All the genomic slices are collected in a single gtf output
        Args:
            repeatmasker_output_dir : pathlib.Path
            region_results_file_path : pathlib.Path
            region_name :str

        region_results_file_path format
        SW    perc perc perc query    position in query matching repeat       position in repeat
        score div. del. ins. sequence begin end (left)  repeat   class/family begin end  (left)  ID
        """
        with open(repeatmasker_output_file_path, "r", encoding="utf8") as repeatmasker_in, open(
            region_results_file_path, "w+", encoding="utf8"
        ) as repeatmasker_out:
            repeat_count = 1
            for line in repeatmasker_in:
                result_match = re.search(r"^\s*\d+\s+", line)
                if result_match:
                    results = line.split()
                    if results[-1] == "*":
                        results.pop()
                    if not len(results) == 15:
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
