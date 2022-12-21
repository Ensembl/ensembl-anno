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

import utils

logger = logging.getLogger(__name__)
with open("./config.json", "r") as f:
    config = json.load(f)


def run_repeatmasker_regions(  # pylint: disable=too-many-arguments
    genome_file: typing.Union[pathlib.Path, str],
    repeatmasker_path: str,
    library: str,
    species: str,
    main_output_dir: str,
    num_threads: int,
):
    """
    Run Repeatmasker on genomic slices using the crossmatch engine.
    Args:
        genome_file : pathlib.Path
        repeatmasker_path : str path to the RepeatMasker executable
        library : str
        species :str
        main_output_dir : pathlib.Path
        num_threads: int

    Return:
        A GTF file with the repeatmasked sequence for each genome slice

    """
    if not repeatmasker_path:
        repeatmasker_path = config["repeatmasker"]["software"]

    utils.check_exe(repeatmasker_path)
    repeatmasker_output_dir = pathlib.Path(
        utils.create_dir(main_output_dir, "repeatmasker_output")
    )

    output_file = repeatmasker_output_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = utils.check_gtf_content(output_file, "repeat")
        if transcript_count > 0:
            logger.info("Repeatmasker gtf file exists")
            return

    logger.info("Creating list of genomic slices")
    seq_region_lengths = utils.get_seq_region_lengths(genome_file, 5000)
    slice_ids = utils.create_slice_ids(
        seq_region_lengths, slice_size=1000000, overlap=0, min_length=5000
    )
    generic_repeatmasker_cmd = [
        repeatmasker_path,
        "-nolow",
        "-engine",
        config["repeatmasker"]["engine"],
        "-dir",
        repeatmasker_output_dir,
    ]
    if not library:
        if not species:
            species = "homo"
            generic_repeatmasker_cmd.extend(["-species", species])

        else:
            generic_repeatmasker_cmd.extend(["-species", species])
    else:
        generic_repeatmasker_cmd.extend(["-lib", library])
    logger.info("Running RepeatMasker")
    pool = multiprocessing.Pool(num_threads)
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
    utils.slice_output_to_gtf(
        str(repeatmasker_output_dir), ".rm.gtf", 1, "repeat_id", "repeatmask"
    )


def multiprocess_repeatmasker(  # pylint: disable=too-many-locals
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
    seq = utils.get_sequence(
        region_name, start, end, 1, genome_file, str(repeatmasker_output_dir)
    )

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
    subprocess.run(repeatmasker_cmd, check=True)

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


def create_repeatmasker_gtf(  # pylint: disable=too-many-locals
    repeatmasker_output_file_path, region_results_file_path, region_name
):
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

                gtf_line = (
                    f"{region_name}\tRepeatMasker\trepeat\t{start}\t{end}\t.\t"
                    f"{strand}\t.\trepeat_id{repeat_count}; "
                    f'repeat_name "{repeat_name}"; repeat_class "{repeat_class}"; '
                    f'repeat_start "{repeat_start}"; '
                    f'repeat_end "{repeat_end}"; score "{score}";\n'
                )
                repeatmasker_out.write(gtf_line)
                repeat_count += 1


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

    utils.check_exe(dust_path)
    dust_output_dir = pathlib.Path(utils.create_dir(main_output_dir, "dust_output"))
    os.chdir(str(dust_output_dir))
    output_file = dust_output_dir / "annotation.gtf"
    logger.info("dust output %s", output_file)
    if output_file.is_file():
        transcript_count = utils.check_gtf_content(output_file, "repeat")
        if transcript_count > 0:
            logger.info("Dust gtf file exists")
            return 0

    logger.info("Creating list of genomic slices")
    seq_region_lengths = utils.get_seq_region_lengths(genome_file, 5000)
    slice_ids = utils.create_slice_ids(
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
    utils.slice_output_to_gtf(str(dust_output_dir), ".dust.gtf", 1, "repeat_id", "dust")
    return 0


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
        "Processing slice to find low complexity regions with Dust: %s:%s:%s",
        region_name,
        start,
        end,
    )
    seq = utils.get_sequence(
        region_name, start, end, 1, genome_file, str(dust_output_dir)
    )

    slice_file_name = f"{region_name}.rs{start}.re{end}"
    region_fasta_file_path = dust_output_dir / f"{slice_file_name}.fa"
    with open(region_fasta_file_path, "w+") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")

    region_results_file_path = dust_output_dir / f"{slice_file_name}.dust.gtf"

    dust_output_file_path = f"{region_fasta_file_path}.dust"
    dust_out = open(dust_output_file_path, "w+")
    dust_cmd = generic_dust_cmd.copy()
    dust_cmd.append(region_fasta_file_path)
    logger.info(dust_cmd)
    subprocess.run(dust_cmd, stdout=dust_out, check=True)
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
                gtf_line = (
                    f"{region_name}\tDust\trepeat\t{start}\t"
                    f'{end}\t.\t+\t.\trepeat_id "{repeat_count}";\n'
                )
                dust_out.write(gtf_line)
                repeat_count += 1


def run_trf_repeats(  # pylint: disable=too-many-locals
    genome_file: typing.Union[pathlib.Path, str],
    trf_path,
    main_output_dir,
    num_threads: int,
):
    """
    Run trf on genomic slices

    Args:
        genome_file : pathlib.Path
        trf_path : str
        main_output_dir : pathlib.Path
        num_threads: int

    Return:
        gtfs with the masked sequence for each genomic slice
    """

    if not trf_path:
        trf_path = config["trf"]["software"]

    utils.check_exe(trf_path)
    trf_output_dir = pathlib.Path(utils.create_dir(main_output_dir, "trf_output"))
    os.chdir(str(trf_output_dir))
    output_file = trf_output_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = utils.check_gtf_content(output_file, "repeat")
        if transcript_count > 0:
            logger.info("Trf gtf file exists, skipping analysis")
            return

    logger.info("Creating list of genomic slices")
    seq_region_lengths = utils.get_seq_region_lengths(genome_file, 5000)
    slice_ids = utils.create_slice_ids(
        seq_region_lengths, slice_size=1000000, overlap=0, min_length=5000
    )

    match_score = config["trf"]["match_score"]
    mismatch_score = config["trf"]["mismatch_score"]
    delta = config["trf"]["delta"]
    pm = config["trf"]["pm"]  # pylint: disable=invalid-name
    pi = config["trf"]["pi"]  # pylint: disable=invalid-name
    minscore = config["trf"]["minscore"]
    maxperiod = config["trf"]["maxperiod"]

    trf_output_extension = (
        f".{match_score}.{mismatch_score}.{delta}.{pm}.{pi}.{minscore}.{maxperiod}.dat"
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
    utils.slice_output_to_gtf(str(trf_output_dir), ".trf.gtf", 1, "repeat_id", "trf")


def multiprocess_trf(
    generic_trf_cmd,
    slice_id,
    genome_file,
    trf_output_dir,
    trf_output_extension,
):
    """
    Run Dust on multiprocess on genomic slices
    Args:
        generic_trf_cmd:list
        slice_id:str
        genome_file : pathlib.Path
        trf_output_dir : pathlib.Path
        trf_output_extension: str
    """
    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing slice to find tandem repeats with TRF:%s:%s:%s",
        region_name,
        start,
        end,
    )
    seq = utils.get_sequence(region_name, start, end, 1, genome_file, str(trf_output_dir))

    slice_file_name = f"{region_name}.rs{start}.re{end}"
    region_fasta_file_path = trf_output_dir / f"{slice_file_name}.fa"
    with open(region_fasta_file_path, "w+") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")

    region_results_file_path = trf_output_dir / f"{slice_file_name}.trf.gtf"

    # TRF writes to the current dir, so swtich to the output dir for it
    os.chdir(str(trf_output_dir))
    trf_output_file_path = f"{region_fasta_file_path}{trf_output_extension}"
    trf_cmd = generic_trf_cmd.copy()
    trf_cmd[1] = str(region_fasta_file_path)
    logger.info("trf_cmd: %s", trf_cmd)
    subprocess.run(trf_cmd)  # pylint: disable=subprocess-run-check
    create_trf_gtf(trf_output_file_path, region_results_file_path, region_name)
    trf_output_file_path.unlink()
    region_fasta_file_path.unlink()


def create_trf_gtf(
    trf_output_file_path, region_results_file_path, region_name
):  # pylint: disable=too-many-locals
    """
        Read the fasta file and save the content in gtf format

        All the genomic slices are collected in a single gtf output
        Args:
            trf_output_file_path : pathlib.Path
            region_results_file_path : pathlib.Path
            region_name :str

        trf_output_file_path is txt file space delimited where the colummns are
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
    """
    with open(trf_output_file_path, "r") as trf_in, open(
        region_results_file_path, "w+"
    ) as trf_out:
        repeat_count = 1
        for line in trf_in:
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


def run_red(red_path, main_output_dir, genome_file: typing.Union[pathlib.Path, str]):
    """
    Run Red on genome file

    Args:
        red_path : str
        main_output_dir : pathlib.Path
        genome_file : pathlib.Path

    Return:
        masked genome file
    """
    if not red_path:
        red_path = config["red"]["software"]
    genome_file = pathlib.Path(genome_file)
    utils.check_exe(red_path)
    red_dir = pathlib.Path(utils.create_dir(main_output_dir, "red_output"))
    red_mask_dir = pathlib.Path(utils.create_dir(red_dir, "mask_output"))
    red_repeat_dir = pathlib.Path(utils.create_dir(red_dir, "repeat_output"))
    red_genome_dir = pathlib.Path(utils.create_dir(red_dir, "genome_dir"))

    sym_link_genome_cmd = "ln -s " + str(genome_file) + " " + str(red_genome_dir)
    logger.info(genome_file)
    genome_file_name = genome_file.name
    red_genome_file = red_genome_dir / genome_file_name
    genome_file_stem = genome_file.stem
    masked_genome_file = red_mask_dir / f"{genome_file_stem}.msk"
    repeat_coords_file = red_repeat_dir / f"{genome_file_stem}.rpt"
    gtf_output_file_path = red_dir / "annotation.gtf"

    if masked_genome_file.exists():
        logger.warning(
            "Masked Genome file already found on the path to the Red mask output dir. \
            Will not create a new file"
        )
        create_red_gtf(repeat_coords_file, gtf_output_file_path)
        return str(masked_genome_file)

    if red_genome_file.exists():
        logger.warning(
            "Unmasked genome file already found on the path to the Red genome dir, \
            will not create a sym link"
        )

    else:
        logger.info(
            "Preparing to sym link the genome file to the Red genome dir. Cmd\n %s",
            sym_link_genome_cmd,
        )
        # subprocess.run(["ln", "-s", genome_file, red_genome_dir])
        red_genome_file.symlink_to(genome_file)
    if not red_genome_file.exists():
        logger.error(
            "Could not find the genome file in the Red genome dir or sym link \
            to the original file. Path expected:\n%s",
            red_genome_file,
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
        ],
        check=True,
    )

    logger.info("Completed running Red")

    create_red_gtf(repeat_coords_file, gtf_output_file_path)

    return str(masked_genome_file)


def create_red_gtf(repeat_coords_file, gtf_output_file_path):
    """
    Create Red gtf file  from masked  genome file

    Args:
        repeat_coords_file: pathlib.Path
        gtf_output_file_path : pathlib.Path
    """
    with open(repeat_coords_file, "r") as red_in, open(
        gtf_output_file_path, "w+"
    ) as red_out:
        for repeat_id, line in enumerate(red_in, start=1):
            result_match = re.search(r"^\>(.+)\:(\d+)\-(\d+)", line)
            if result_match:
                region_name = result_match.group(1)
                # Note that Red is 0-based, so add 1
                start = int(result_match.group(2)) + 1
                end = int(result_match.group(3)) + 1
                gtf_line = (
                    f"{region_name}\tRed\trepeat\t{start}\t"
                    f'{end}\t.\t+\t.\trepeat_id "{repeat_id}";\n'
                )
                red_out.write(gtf_line)
