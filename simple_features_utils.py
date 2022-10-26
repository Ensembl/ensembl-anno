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

# standard library
import argparse
import gc
import glob
import io
import math
import multiprocessing
import os
import pathlib
import random
import re
import shutil
import signal
import subprocess
import tempfile
import sys
import errno
import logging

from typing import List, Union

def run_eponine_regions(
    genome_file: Union[pathlib.Path, str],
    java_path,
    eponine_path,
    main_output_dir,
    num_threads: int,
):
    if not java_path:
        java_path = "java"

    if not eponine_path:
        eponine_path = "/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/opt/eponine/libexec/eponine-scan.jar"

    check_file(eponine_path)
    check_exe(java_path)

    eponine_output_dir = create_dir(main_output_dir, "eponine_output")

    output_file = eponine_output_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "simple_feature")
        if transcript_count > 0:
            logger.info("Eponine gtf file already exists, skipping analysis")
            return

    logger.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, min_seq_length=5000)
    slice_ids = create_slice_ids(
        seq_region_lengths, slice_size=1_000_000, overlap=0, min_length=5000
    )

    threshold = "0.999"
    generic_eponine_cmd = [
        java_path,
        "-jar",
        eponine_path,
        "-threshold",
        threshold,
        "-seq",
    ]
    logger.info("Running Eponine processes")
    pool = multiprocessing.Pool(num_threads)
    tasks = []
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
    slice_output_to_gtf(eponine_output_dir, ".epo.gtf", 1, "feature_id", "eponine")


def multiprocess_eponine(
    generic_eponine_cmd,
    slice_id,
    genome_file: Union[pathlib.Path, str],
    eponine_output_dir: Union[pathlib.Path, str],
):
    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing slice to find repeats with Eponine: %s:%s:%s"
        % (region_name, start, end)
    )
    seq = get_sequence(
        seq_region=region_name,
        start=start,
        end=end,
        strand=1,
        fasta_file=genome_file,
        output_dir=eponine_output_dir,
    )

    slice_file_name = f"{region_name}.rs{start}.re{end}"
    region_fasta_file_name = f"{slice_file_name}.fa"
    region_fasta_file_path = eponine_output_dir / region_fasta_file_name

    with open(region_fasta_file_path, "w+") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")

    region_results_file_name = f"{slice_file_name}.epo.gtf"
    region_results_file_path = eponine_output_dir / region_results_file_name

    eponine_output_file_path = f"{region_fasta_file_path}.epo"
    eponine_out = open(eponine_output_file_path, "w+")

    eponine_cmd = generic_eponine_cmd.copy()
    eponine_cmd.append(region_fasta_file_path)

    logger.info("eponine_cmd: %s" % " ".join(eponine_cmd))
    subprocess.run(eponine_cmd, stdout=eponine_out)
    eponine_out.close()

    create_eponine_gtf(eponine_output_file_path, region_results_file_path, region_name)
    os.remove(eponine_output_file_path)
    os.remove(region_fasta_file_path)


def create_eponine_gtf(eponine_output_file_path, region_results_file_path, region_name):
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

                gtf_line = f'{region_name}\tEponine\tsimple_feature\t{start}\t{end}\t.\t{strand}\t.\tfeature_id "{feature_count}"; score "{score}";\n'
                eponine_out.write(gtf_line)
                feature_count += 1


def run_cpg_regions(
    genome_file: Union[pathlib.Path, str], cpg_path, main_output_dir, num_threads: int
):
    if not cpg_path:
        cpg_path = "cpg_lh"

    check_exe(cpg_path)
    cpg_output_dir = create_dir(main_output_dir, "cpg_output")

    output_file = cpg_output_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "simple_feature")
        if transcript_count > 0:
            logger.info("Cpg gtf file already exists, skipping analysis")
            return

    logger.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, min_seq_length=5000)
    slice_ids = create_slice_ids(
        seq_region_lengths, slice_size=1_000_000, overlap=0, min_length=5000
    )

    logger.info("Running CpG processes")
    pool = multiprocessing.Pool(num_threads)
    tasks = []
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
    slice_output_to_gtf(cpg_output_dir, ".cpg.gtf", 1, "feature_id", "cpg")


def multiprocess_cpg(
    cpg_path, slice_id, genome_file, cpg_output_dir: Union[pathlib.Path, str]
):
    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing slice to find CpG islands with cpg_lh: %s:%s:%s"
        % (region_name, start, end)
    )
    seq = get_sequence(
        seq_region=region_name,
        start=start,
        end=end,
        strand=1,
        fasta_file=genome_file,
        output_dir=cpg_output_dir,
    )

    slice_file_name = f"{region_name}.rs{start}.re{end}"
    region_fasta_file_name = f"{slice_file_name}.fa"
    region_fasta_file_path = cpg_output_dir / region_fasta_file_name

    with open(region_fasta_file_path, "w+") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")

    region_results_file_name = f"{slice_file_name}.cpg.gtf"
    region_results_file_path = cpg_output_dir / region_results_file_name

    cpg_output_file_path = f"{region_fasta_file_path}.cpg"
    cpg_out = open(cpg_output_file_path, "w+")

    cpg_cmd = [cpg_path, region_fasta_file_path]
    logger.info("cpg_cmd: %s" % " ".join(cpg_cmd))
    subprocess.run(cpg_cmd, stdout=cpg_out)
    cpg_out.close()

    create_cpg_gtf(cpg_output_file_path, region_results_file_path, region_name)
    os.remove(cpg_output_file_path)
    os.remove(region_fasta_file_path)


def create_cpg_gtf(cpg_output_file_path, region_results_file_path, region_name):
    cpg_min_length = 400
    cpg_min_gc_content = 50
    cpg_min_oe = 0.6

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
                oe = results[7]

                if oe == "-" or oe == "inf":
                    oe = 0
                else:
                    oe = float(oe)

                if (
                    length >= cpg_min_length
                    and gc_content >= cpg_min_gc_content
                    and oe >= cpg_min_oe
                ):
                    gtf_line = f'{region_name}\tCpG\tsimple_feature\t{start}\t{end}\t.\t+\t.\tfeature_id "{feature_count}"; score "{score}";\n'
                    cpg_out.write(gtf_line)
                    feature_count += 1


