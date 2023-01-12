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

logger = logging.getLogger(__name__)
with open(os.environ["ENSCODE"] + "/ensembl-anno/config.json", "r") as f:
    config = json.load(f)


def run_eponine_regions(
    genome_file, java_path, eponine_path, main_output_dir, num_threads
):

    if not java_path:
        java_path = config["java"]["path"]

    if not eponine_path:
        eponine_path = config["eponine"]["software"]

    check_file(eponine_path)
    check_exe(java_path)

    eponine_output_dir = utils.create_dir(main_output_dir, "eponine_output")

    logger.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(eponine_output_dir, "annotation.gtf")
    if os.path.exists(output_file):
        transcript_count = utils.check_gtf_content(output_file, "simple_feature")
        if transcript_count > 0:
            logger.info("Eponine gtf file exists")
            return
    else:
        logger.info("No gtf file, go on with the analysis")

    logger.info("Creating list of genomic slices")
    seq_region_lengths = utils.get_seq_region_lengths(genome_file, 5000)
    slice_ids = utils.create_slice_ids(seq_region_lengths, 1000000, 0, 5000)

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
    utils.slice_output_to_gtf(eponine_output_dir, ".epo.gtf", 1, "feature_id", "eponine")


def multiprocess_eponine(generic_eponine_cmd, slice_id, genome_file, eponine_output_dir):

    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing slice to find transcription start sites with Eponine: "
        + region_name
        + ":"
        + str(start)
        + ":"
        + str(end)
    )
    seq = utils.get_sequence(region_name, start, end, 1, genome_file, eponine_output_dir)

    slice_file_name = region_name + ".rs" + str(start) + ".re" + str(end)
    region_fasta_file_name = slice_file_name + ".fa"
    region_fasta_file_path = os.path.join(eponine_output_dir, region_fasta_file_name)

    region_fasta_out = open(region_fasta_file_path, "w+")
    region_fasta_out.write(">" + region_name + "\n" + seq + "\n")
    region_fasta_out.close()

    region_results_file_name = slice_file_name + ".epo.gtf"
    region_results_file_path = os.path.join(eponine_output_dir, region_results_file_name)

    eponine_output_file_path = region_fasta_file_path + ".epo"
    eponine_out = open(eponine_output_file_path, "w+")

    eponine_cmd = generic_eponine_cmd.copy()
    eponine_cmd.append(region_fasta_file_path)

    logger.info(" ".join(eponine_cmd))
    subprocess.run(eponine_cmd, stdout=eponine_out)
    eponine_out.close()

    create_eponine_gtf(eponine_output_file_path, region_results_file_path, region_name)
    os.remove(eponine_output_file_path)
    os.remove(region_fasta_file_path)


def create_eponine_gtf(eponine_output_file_path, region_results_file_path, region_name):

    eponine_in = open(eponine_output_file_path, "r")
    eponine_out = open(region_results_file_path, "w+")

    line = eponine_in.readline()
    feature_count = 1
    while line:
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
                region_name
                + "\tEponine\tsimple_feature\t"
                + str(start)
                + "\t"
                + str(end)
                + "\t.\t"
                + strand
                + "\t.\t"
                + 'feature_id "'
                + str(feature_count)
                + '"; score "'
                + str(score)
                + '";\n'
            )
            eponine_out.write(gtf_line)
            feature_count += 1
        line = eponine_in.readline()
    eponine_in.close()
    eponine_out.close()


def run_cpg_regions(genome_file, cpg_path, main_output_dir, num_threads):

    if not cpg_path:
        cpg_path = config["cpg"]["software"]

    utils.check_exe(cpg_path)
    cpg_output_dir = utils.create_dir(main_output_dir, "cpg_output")

    logger.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(cpg_output_dir, "annotation.gtf")
    if os.path.exists(output_file):
        transcript_count = utils.check_gtf_content(output_file, "simple_feature")
        if transcript_count > 0:
            logger.info("Cpg gtf file exists")
            return
    else:
        logger.info("No gtf file, go on with the analysis")

    logger.info("Creating list of genomic slices")
    seq_region_lengths = utils.get_seq_region_lengths(genome_file, 5000)
    slice_ids = utils.create_slice_ids(seq_region_lengths, 1000000, 0, 5000)

    logger.info("Running CpG")
    pool = multiprocessing.Pool(int(num_threads))
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
    utils.slice_output_to_gtf(cpg_output_dir, ".cpg.gtf", 1, "feature_id", "cpg")


def multiprocess_cpg(cpg_path, slice_id, genome_file, cpg_output_dir):

    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing slice to find CpG islands with cpg_lh: "
        + region_name
        + ":"
        + str(start)
        + ":"
        + str(end)
    )
    seq = utils.get_sequence(region_name, start, end, 1, genome_file, cpg_output_dir)

    slice_file_name = region_name + ".rs" + str(start) + ".re" + str(end)
    region_fasta_file_name = slice_file_name + ".fa"
    region_fasta_file_path = os.path.join(cpg_output_dir, region_fasta_file_name)

    region_fasta_out = open(region_fasta_file_path, "w+")
    region_fasta_out.write(">" + region_name + "\n" + seq + "\n")
    region_fasta_out.close()

    region_results_file_name = slice_file_name + ".cpg.gtf"
    region_results_file_path = os.path.join(cpg_output_dir, region_results_file_name)

    cpg_output_file_path = region_fasta_file_path + ".cpg"
    cpg_out = open(cpg_output_file_path, "w+")

    cpg_cmd = [cpg_path, region_fasta_file_path]
    logger.info(" ".join(cpg_cmd))
    subprocess.run(cpg_cmd, stdout=cpg_out)
    cpg_out.close()

    create_cpg_gtf(cpg_output_file_path, region_results_file_path, region_name)
    os.remove(cpg_output_file_path)
    os.remove(region_fasta_file_path)


def create_cpg_gtf(cpg_output_file_path, region_results_file_path, region_name):

    cpg_min_length = config["cpg"]["cpg_min_length"]
    cpg_min_gc_content = config["cpg"]["cpg_min_gc_content"]
    cpg_min_oe = config["cpg"]["cpg_min_oe"]

    cpg_in = open(cpg_output_file_path, "r")
    cpg_out = open(region_results_file_path, "w+")
    line = cpg_in.readline()
    feature_count = 1
    while line:
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
                gtf_line = (
                    region_name
                    + "\tCpG\tsimple_feature\t"
                    + str(start)
                    + "\t"
                    + str(end)
                    + "\t.\t+\t.\t"
                    + 'feature_id "'
                    + str(feature_count)
                    + '"; score "'
                    + str(score)
                    + '";\n'
                )
                cpg_out.write(gtf_line)
                feature_count += 1
        line = cpg_in.readline()
    cpg_in.close()
    cpg_out.close()
