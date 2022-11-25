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

import multiprocessing
import os
import pathlib
import re
import subprocess

from pathlib import Path

from typing import Union

def run_repeatmasker_regions(
    genome_file: Union[pathlib.Path, str],
    repeatmasker_path: str,
    library: str,
    species: str,
    main_output_dir: Union[pathlib.Path, str],
    num_threads: int,
):
    """
    Run Repeatmasker on genomic slices

    engine = crossmatch
    Args:
        genome_file : pathlib.Path
        repeatmasker_path : str
        library : str
        species :str
        main_output_dir : pathlib.Path
        num_threads: int

    Return:
        gtfs with the repeatmasked sequence for each genome slice

    """
    if not repeatmasker_path:
        repeatmasker_path = "RepeatMasker"

    if not library:
        library = "homo"

    check_exe(repeatmasker_path)
    repeatmasker_output_dir = Path(create_dir(main_output_dir, "repeatmasker_output"))
    os.chdir(repeatmasker_output_dir)

    output_file = repeatmasker_output_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "repeat")
        if transcript_count > 0:
            logger.info("Repeatmasker gtf file exists")
            return

    logger.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, 5000)
    slice_ids = create_slice_ids(seq_region_lengths, 1000000, 0, 5000)

    if not library:
        if not species:
            species = "homo"
        generic_repeatmasker_cmd = [
            repeatmasker_path,
            "-nolow",
            "-species",
            species,
            "-engine",
            "crossmatch",
            "-dir",
            repeatmasker_output_dir,
        ]

    else:
        generic_repeatmasker_cmd = [
            repeatmasker_path,
            "-nolow",
            "-lib",
            library,
            "-engine",
            "crossmatch",
            "-dir",
            repeatmasker_output_dir,
        ]

    logger.info("Running RepeatMasker")
    pool = multiprocessing.Pool(num_threads)
    tasks = []
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_repeatmasker,
            args=(
                generic_repeatmasker_cmd,
                slice_id,
                genome_file,
                repeatmasker_output_dir,
            ),
        )

    pool.close()
    pool.join()
    slice_output_to_gtf(repeatmasker_output_dir, ".rm.gtf", 1, "repeat_id", "repeatmask")


def multiprocess_repeatmasker(
    generic_repeatmasker_cmd, slice_id, genome_file, repeatmasker_output_dir
):
    """
    Run Repeatmasker on multiprocess on genomic slices

    engine = crossmatch
    Args:
        generic_repeatmasker_cmd: list
        slice_id: str
        genome_file : pathlib.Path
        repeatmasker_output_dir : pathlib.Path
    """

    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]
    logger.info(
        "Processing slice to find repeats with RepeatMasker: {region_name}:{start}:{end}"
    )
    seq = get_sequence(region_name, start, end, 1, genome_file, repeatmasker_output_dir)

    slice_file_name = f"{region_name}.rs{start}.re{end}"
    region_fasta_file_path = repeatmasker_output_dir / f"{slice_file_name}.fa"
    with open(region_fasta_file_path, "w+") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")

    region_results_file_path = f"{region_fasta_file_path}.rm.gtf"
    repeatmasker_output_file_path = f"{region_fasta_file_path}.out"
    repeatmasker_masked_file_path = f"{region_fasta_file_path}.masked"
    repeatmasker_tbl_file_path = f"{region_fasta_file_path}.tbl"
    repeatmasker_log_file_path = f"{region_fasta_file_path}.log"
    repeatmasker_cat_file_path = f"{region_fasta_file_path}.cat"

    repeatmasker_cmd = generic_repeatmasker_cmd.copy()
    repeatmasker_cmd.append(region_fasta_file_path)
    logger.info(" ".join(repeatmasker_cmd))
    subprocess.run(repeatmasker_cmd)

    create_repeatmasker_gtf(
        repeatmasker_output_file_path, region_results_file_path, region_name
    )

    region_fasta_file_path.unlink()
    if region_results_file_path.exists():
        region_results_file_path.unlink()
    if repeatmasker_masked_file_path.exists():
        repeatmasker_masked_file_path.unlink()
    if repeatmasker_tbl_file_path.exists():
        repeatmasker_tbl_file_path.unlink()
    if repeatmasker_log_file_path.exists():
        repeatmasker_log_file_path.unlink()
    if repeatmasker_cat_file_path.exists():
        repeatmasker_cat_file_path.unlink()


def create_repeatmasker_gtf(
    repeatmasker_output_file_path, region_results_file_path, region_name
):
    """
    Read the fasta file and save the content in gtf format

    All the genomic slices are collected in a single gtf output
    Args:
        repeatmasker_output_dir : pathlib.Path
        region_results_file_path : pathlib.Path
        region_name :str
    """
    with open(repeatmasker_output_file_path, "r") as repeatmasker_in, open(
        region_results_file_path, "w+"
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

                gtf_line = f'{region_name}\tRepeatMasker\trepeat\t{start}\t{end}\t.\t{strand}\t.\trepeat_id{repeat_count}; repeat_name "{repeat_name}"; repeat_class "{repeat_class}"; repeat_start "{repeat_start}"; repeat_end "{repeat_end}"; score "{score}";\n'
                repeatmasker_out.write(gtf_line)
                repeat_count += 1


def run_dust_regions(
    genome_file: Union[pathlib.Path, str],
    dust_path: Union[pathlib.Path, str],
    main_output_dir,
    num_threads,
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
        dust_path = "dustmasker"

    check_exe(dust_path)
    dust_output_dir = Path(create_dir(main_output_dir, "dust_output"))

    output_file = dust_output_dir / "annotation.gtf"
    logger.info(f"dust output {output_file}")
    if output_file.is_file():
        transcript_count = check_gtf_content(output_file, "repeat")
        if transcript_count > 0:
            logger.info("Dust gtf file exists")
            return 0

    logger.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, 5000)
    slice_ids = create_slice_ids(seq_region_lengths, 1000000, 0, 5000)

    generic_dust_cmd = [dust_path, "-in"]
    logger.info("Running Dust")
    pool = multiprocessing.Pool(int(num_threads))
    tasks = []
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
    slice_output_to_gtf(dust_output_dir, ".dust.gtf", 1, "repeat_id", "dust")


def multiprocess_dust(generic_dust_cmd, slice_id, genome_file, dust_output_dir):
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
        f"Processing slice to find low complexity regions with Dust: {region_name}:{start}:{end}"
    )
    seq = get_sequence(region_name, start, end, 1, genome_file, dust_output_dir)

    slice_file_name = f"{region_name}.rs{start}.re{end}"
    region_fasta_file_path = dust_output_dir / f"{slice_file_name}.fa"
    with open(region_fasta_file_path, "w+") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")

    region_results_file_path = dust_output_dir / f"{slice_file_name}.dust.gtf"

    dust_output_file_path = region_fasta_file_path + ".dust"
    dust_out = open(dust_output_file_path, "w+")
    dust_cmd = generic_dust_cmd.copy()
    dust_cmd.append(region_fasta_file_path)
    logger.info(" ".join(dust_cmd))
    subprocess.run(dust_cmd, stdout=dust_out)
    dust_out.close()

    create_dust_gtf(dust_output_file_path, region_results_file_path, region_name)
    dust_output_file_path.unlink()
    region_fasta_file_path.unlink()


def create_dust_gtf(dust_output_file_path, region_results_file_path, region_name):
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
                gtf_line = f'{region_name}\tDust\trepeat\t{start}\t{end}\t.\t+\t.\trepeat_id "{repeat_count}";\n'
                dust_out.write(gtf_line)
                repeat_count += 1


def run_trf_repeats(genome_file, trf_path, main_output_dir, num_threads):

    if not trf_path:
        trf_path = "trf"

    check_exe(trf_path)
    trf_output_dir = create_dir(main_output_dir, "trf_output")

    logger.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(trf_output_dir, "annotation.gtf")
    if os.path.exists(output_file):
        transcript_count = check_gtf_content(output_file, "repeat")
        if transcript_count > 0:
            logger.info("Trf gtf file exists")
            return
    else:
        logger.info("No gtf file, go on with the analysis")

    logger.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, 5000)
    slice_ids = create_slice_ids(seq_region_lengths, 1000000, 0, 5000)

    match_score = 2
    mismatch_score = 5
    delta = 7
    pm = 80
    pi = 10
    minscore = 40
    maxperiod = 500

    trf_output_extension = (
        "."
        + str(match_score)
        + "."
        + str(mismatch_score)
        + "."
        + str(delta)
        + "."
        + str(pm)
        + "."
        + str(pi)
        + "."
        + str(minscore)
        + "."
        + str(maxperiod)
        + ".dat"
    )

    generic_trf_cmd = [
        trf_path,
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
    pool = multiprocessing.Pool(int(num_threads))
    tasks = []
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_trf,
            args=(
                generic_trf_cmd,
                slice_id,
                genome_file,
                trf_output_dir,
                trf_output_extension,
            ),
        )

    pool.close()
    pool.join()
    slice_output_to_gtf(trf_output_dir, ".trf.gtf", 1, "repeat_id", "trf")


def multiprocess_trf(
    generic_trf_cmd,
    slice_id,
    genome_file,
    trf_output_dir,
    trf_output_extension,
):

    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing slice to find tandem repeats with TRF: "
        + region_name
        + ":"
        + str(start)
        + ":"
        + str(end)
    )
    seq = get_sequence(region_name, start, end, 1, genome_file, trf_output_dir)

    slice_file_name = region_name + ".rs" + str(start) + ".re" + str(end)
    region_fasta_file_name = slice_file_name + ".fa"
    region_fasta_file_path = os.path.join(trf_output_dir, region_fasta_file_name)
    region_fasta_out = open(region_fasta_file_path, "w+")
    region_fasta_out.write(">" + region_name + "\n" + seq + "\n")
    region_fasta_out.close()

    region_results_file_name = slice_file_name + ".trf.gtf"
    region_results_file_path = os.path.join(trf_output_dir, region_results_file_name)

    # TRF writes to the current dir, so swtich to the output dir for it
    os.chdir(trf_output_dir)
    trf_output_file_path = region_fasta_file_path + trf_output_extension
    trf_cmd = generic_trf_cmd.copy()
    trf_cmd[1] = region_fasta_file_path
    logger.info(" ".join(trf_cmd))
    subprocess.run(trf_cmd)
    create_trf_gtf(trf_output_file_path, region_results_file_path, region_name)
    os.remove(trf_output_file_path)
    os.remove(region_fasta_file_path)


def create_trf_gtf(trf_output_file_path, region_results_file_path, region_name):

    trf_in = open(trf_output_file_path, "r")
    trf_out = open(region_results_file_path, "w+")
    line = trf_in.readline()
    repeat_count = 1
    while line:
        result_match = re.search(r"^\d+", line)
        if result_match:
            results = line.split()
            if not len(results) == 15:
                continue
            start = results[0]
            end = results[1]
            period = float(results[2])
            copy_number = float(results[3])
            percent_matches = float(results[5])
            score = float(results[7])
            repeat_consensus = results[13]
            if (
                score < 50 and percent_matches >= 80 and copy_number > 2 and period < 10
            ) or (copy_number >= 2 and percent_matches >= 70 and score >= 50):
                gtf_line = (
                    region_name
                    + "\tTRF\trepeat\t"
                    + str(start)
                    + "\t"
                    + str(end)
                    + "\t.\t+\t.\t"
                    + 'repeat_id "'
                    + str(repeat_count)
                    + '"; score "'
                    + str(score)
                    + '"; repeat_consensus "'
                    + repeat_consensus
                    + '";\n'
                )
                trf_out.write(gtf_line)
                repeat_count += 1
        line = trf_in.readline()
    trf_in.close()
    trf_out.close()


def run_red(red_path, main_output_dir, genome_file):

    if not red_path:
        red_path = "Red"

    check_exe(red_path)
    red_dir = create_dir(main_output_dir, "red_output")
    red_mask_dir = create_dir(red_dir, "mask_output")
    red_repeat_dir = create_dir(red_dir, "repeat_output")
    red_genome_dir = create_dir(red_dir, "genome_dir")

    sym_link_genome_cmd = "ln -s " + genome_file + " " + red_genome_dir

    genome_file_name = os.path.basename(genome_file)
    red_genome_file = os.path.join(red_genome_dir, genome_file_name)
    masked_genome_file = os.path.join(
        red_mask_dir, os.path.splitext(genome_file_name)[0] + ".msk"
    )
    repeat_coords_file = os.path.join(
        red_repeat_dir, os.path.splitext(genome_file_name)[0] + ".rpt"
    )
    gtf_output_file_path = os.path.join(red_dir, "annotation.gtf")

    if os.path.exists(masked_genome_file):
        logger.warning(
            "Masked Genome file already found on the path to the Red mask output dir. Will not create a new file"
        )
        create_red_gtf(repeat_coords_file, gtf_output_file_path)
        return masked_genome_file

    if os.path.exists(red_genome_file):
        logger.warning(
            "Unmasked genome file already found on the path to the Red genome dir, will not create a sym link"
        )

    else:
        logger.info(
            "Preparing to sym link the genome file to the Red genome dir. Cmd\n%s"
            % sym_link_genome_cmd
        )
        subprocess.run(["ln", "-s", genome_file, red_genome_dir])

    if not os.path.exists(os.path.join(red_genome_dir, genome_file_name)):
        logger.error(
            "Could not find the genome file in the Red genome dir or sym link to the original file. Path expected:\n%s"
            % red_genome_file
        )

    logger.info("Running Red, this may take some time depending on the genome size")
    subprocess.run(
        [
            red_path,
            "-gnm",
            red_genome_dir,
            "-msk",
            red_mask_dir,
            "-rpt",
            red_repeat_dir,
        ]
    )

    logger.info("Completed running Red")

    create_red_gtf(repeat_coords_file, gtf_output_file_path)

    return masked_genome_file


def create_red_gtf(repeat_coords_file, gtf_output_file_path):

    red_in = open(repeat_coords_file, "r")
    red_out = open(gtf_output_file_path, "w+")
    line = red_in.readline()
    repeat_id = 1
    while line:
        result_match = re.search(r"^\>(.+)\:(\d+)\-(\d+)", line)
        if result_match:
            region_name = result_match.group(1)
            # Note that Red is 0-based, so add 1
            start = int(result_match.group(2)) + 1
            end = int(result_match.group(3)) + 1
            gtf_line = (
                region_name
                + "\tRed\trepeat\t"
                + str(start)
                + "\t"
                + str(end)
                + "\t.\t+\t.\t"
                + 'repeat_id "'
                + str(repeat_id)
                + '";\n'
            )
            red_out.write(gtf_line)
            repeat_id += 1
        line = red_in.readline()
    red_in.close()
    red_out.close()
