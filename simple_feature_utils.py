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

import utils


_REPO_ROOT = pathlib.Path(__file__).parent


logger = logging.getLogger(__name__)
config_file = _REPO_ROOT / "config.json"
with config_file.open("r") as f:
    config = json.load(f)


def run_eponine_regions(  # pylint: disable=too-many-locals
    genome_file: typing.Union[pathlib.Path, str],
    java_path: str,
    eponine_path: str,
    main_output_dir: str,
    num_threads: int,
):
    """
    Run Eponine on genomic slices
    Args:
        genome_file : pathlib.Path
        java_path : str
        eponine_path : str
        main_output_dir : pathlib.Path
        num_threads: int
    Return:
        gtfs with the simple feature evidence for each genome slice
    """
    if not java_path:
        java_path = config["java"]["path"]

    if not eponine_path:
        eponine_path = config["eponine"]["software"]

    utils.check_file(pathlib.Path(eponine_path))
    utils.check_exe(java_path)

    eponine_output_dir = pathlib.Path(utils.create_dir(main_output_dir, "eponine_output"))

    logger.info("Skip analysis if the gtf file already exists")
    output_file = eponine_output_dir / "annotation.gtf"
    if output_file.is_file():
        transcript_count = utils.check_gtf_content(output_file, "simple_feature")
        if transcript_count > 0:
            logger.info("Eponine gtf file exists")
            return 0

    logger.info("Creating list of genomic slices")
    seq_region_lengths = utils.get_seq_region_lengths(genome_file, 5000)
    logger.info("Creating list of genomic slices %s", seq_region_lengths)
    slice_ids = utils.create_slice_ids(
        seq_region_lengths, slice_size=1000000, overlap=0, min_length=5000
    )

    threshold = config["eponine"]["threshold"]
    generic_eponine_cmd = [
        java_path,
        "-jar",
        eponine_path,
        "-threshold",
        threshold,
        "-seq",
    ]
    logger.info("Running Eponine")
    pool = multiprocessing.Pool(int(num_threads))
    # tasks = []
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_eponine,
            args=(
                generic_eponine_cmd,
                slice_id,
                genome_file,
                eponine_output_dir,
            ),
        )

    pool.close()
    pool.join()
    utils.slice_output_to_gtf(
        str(eponine_output_dir), ".epo.gtf", 1, "feature_id", "eponine"
    )
    for gtf_file in pathlib.Path(eponine_output_dir).glob("*.epo.gtf"):
        gtf_file.unlink()
    return 0


def multiprocess_eponine(  # pylint: disable=too-many-locals
    generic_eponine_cmd: list,
    slice_id: str,
    genome_file: pathlib.Path,
    eponine_output_dir: pathlib.Path,
):
    """
    Run Eponine on multiprocess on genomic slices
    Args:
        generic_eponine_cmd:list
        slice_id: list
        genome_file : pathlib.Path
        eponine_output_dir : pathlib.Path
    """
    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing slice to find transcription start sites with Eponine: %s:%s:%s",
        region_name,
        start,
        end,
    )
    seq = utils.get_sequence(
        region_name, start, end, 1, str(genome_file), str(eponine_output_dir)
    )

    slice_file_name = f"{region_name}.rs{start}.re{end}"
    with tempfile.TemporaryDirectory(dir=eponine_output_dir) as tmpdirname:
        region_fasta_file_path = eponine_output_dir / tmpdirname / f"{slice_file_name}.fa"
        with open(region_fasta_file_path, "w+") as region_fasta_out:
            region_fasta_out.write(f">{region_name}\n{seq}\n")

        region_results_file_path = eponine_output_dir / f"{slice_file_name}.epo.gtf"

        eponine_output_file_path = f"{region_fasta_file_path}.epo"
        eponine_cmd = generic_eponine_cmd.copy()
        eponine_cmd.append(str(region_fasta_file_path))
        with open(eponine_output_file_path, "w+") as eponine_out:
            subprocess.run(eponine_cmd, stdout=eponine_out, check=True)
        create_eponine_gtf(
            eponine_output_file_path, region_results_file_path, region_name
        )


def create_eponine_gtf(
    eponine_output_file_path: pathlib.Path,
    region_results_file_path: pathlib.Path,
    region_name: str,
):
    """
    Read the fasta file and save the content in gtf format
    All the genomic slices are collected in a single gtf output
    Args:
        eponine_output_file_path : pathlib.Path
        region_results_file_path : pathlib.Path
        region_name :str
    """
    with open(eponine_output_file_path, "r") as eponine_in, open(
        region_results_file_path, "w+"
    ) as eponine_out:
        feature_count = 1
        for line in eponine_in:
            result_match = re.search(r"^" + region_name, line)
            if result_match:
                results = line.split()
                start = int(results[3])
                end = int(results[4])
                score = float(results[5])
                strand = results[6]

                # There's a one base offset on the reverse strand
                if strand == "-":
                    start -= 1
                    end -= 1

                gtf_line = (
                    f"{region_name}\tEponine\tsimple_feature\t"
                    f"{start}\t{end}\t.\t{strand}\t.\t"
                    f'feature_id "{feature_count}"; score "{score}";\n'
                )
                eponine_out.write(gtf_line)
                feature_count += 1


def run_cpg_regions(
    genome_file: typing.Union[pathlib.Path, str],
    cpg_path: str,
    main_output_dir: str,
    num_threads: int,
):
    """
    Run CpG islands on genomic slices
    Args:
        genome_file : pathlib.Path
        cpg_path : str
        main_output_dir : pathlib.Path
        num_threads: int
    Return:
        gtfs with the cpg evidence for each genome slice
    """
    if not cpg_path:
        cpg_path = config["cpg"]["software"]

    utils.check_exe(cpg_path)
    cpg_output_dir = pathlib.Path(utils.create_dir(main_output_dir, "cpg_output"))

    output_file = cpg_output_dir / "annotation.gtf"
    if output_file.is_file():
        transcript_count = utils.check_gtf_content(output_file, "simple_feature")
        if transcript_count > 0:
            logger.info("Cpg gtf file exists")
            return 0

    logger.info("Creating list of genomic slices")
    seq_region_lengths = utils.get_seq_region_lengths(genome_file, 5000)
    slice_ids = utils.create_slice_ids(
        seq_region_lengths, slice_size=1000000, overlap=0, min_length=5000
    )

    logger.info("Running CpG")
    pool = multiprocessing.Pool(int(num_threads))
    # tasks = []
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_cpg,
            args=(
                cpg_path,
                slice_id,
                genome_file,
                cpg_output_dir,
            ),
        )

    pool.close()
    pool.join()
    utils.slice_output_to_gtf(str(cpg_output_dir), ".cpg.gtf", 1, "feature_id", "cpg")
    for gtf_file in pathlib.Path(cpg_output_dir).glob("*.cpg.gtf"):
        gtf_file.unlink()
    return 0


def multiprocess_cpg(  # pylint: disable=too-many-locals
    cpg_path: str, slice_id: str, genome_file: pathlib.Path, cpg_output_dir: pathlib.Path
):
    """
    Annotation of CpG islands on multiprocess on genomic slices
    Args:
        cpg_path:str
        slice_id:str
        genome_file : pathlib.Path
        cpg_output_dir : pathlib.Path
    """
    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing slice to find CpG islands with cpg_lh: %s:%s:%s",
        region_name,
        start,
        end,
    )
    seq = utils.get_sequence(
        region_name, start, end, 1, str(genome_file), str(cpg_output_dir)
    )

    slice_file_name = f"{region_name}.rs{start}.re{end}"
    with tempfile.TemporaryDirectory(dir=cpg_output_dir) as tmpdirname:
        region_fasta_file_path = cpg_output_dir / tmpdirname / f"{slice_file_name}.fa"

        with open(region_fasta_file_path, "w+") as region_fasta_out:
            region_fasta_out.write(f">{region_name}\n{seq}\n")

        region_results_file_path = cpg_output_dir / f"{slice_file_name}.cpg.gtf"

        cpg_output_file_path = f"{region_fasta_file_path}.cpg"
        cpg_cmd = [cpg_path, str(region_fasta_file_path)]

        with open(cpg_output_file_path, "w+") as cpg_out:
            subprocess.run(cpg_cmd, stdout=cpg_out, check=True)

        create_cpg_gtf(cpg_output_file_path, region_results_file_path, region_name)


def create_cpg_gtf(  # pylint: disable=too-many-locals
    cpg_output_file_path: pathlib.Path,
    region_results_file_path: pathlib.Path,
    region_name: str,
):
    """
    Read the fasta file and save the content in gtf format
    All the genomic slices are collected in a single gtf output
    Args:
        cpg_output_file_path : pathlib.Path
        region_results_file_path : pathlib.Path
        region_name :str
    """
    cpg_min_length = config["cpg"]["cpg_min_length"]
    cpg_min_gc_content = config["cpg"]["cpg_min_gc_content"]
    cpg_min_oe = config["cpg"]["cpg_min_oe"]

    with open(cpg_output_file_path, "r") as cpg_in, open(
        region_results_file_path, "w+"
    ) as cpg_out:
        feature_count = 1
        for line in cpg_in:
            result_match = re.search(r"^" + region_name, line)
            if result_match:
                results = line.split()
                start = int(results[1])
                end = int(results[2])
                length = end - start + 1
                score = float(results[3])
                gc_content = float(results[6])
                oe_score = results[7]
                if oe_score in ("-", "inf"):
                    oe_score = 0
                else:
                    oe_score = float(oe_score)
                if (
                    int(length) >= int(cpg_min_length)
                    and gc_content >= int(cpg_min_gc_content)
                    and oe_score >= float(cpg_min_oe)
                ):
                    gtf_line = (
                        f"{region_name}\tCpG\tsimple_feature\t{start}\t"
                        f'{end}\t.\t+\t.\tfeature_id "{feature_count}"; score "{score}";\n'
                    )
                    cpg_out.write(gtf_line)
                    feature_count += 1
