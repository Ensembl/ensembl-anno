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
import math
import multiprocessing
import os
import pathlib
import re
import shutil
import subprocess

from typing import Union

# project imports
from utils import (
    check_exe,
    create_dir,
    create_paired_paths,
    get_seq_region_lengths,
    logger,
)


def run_star_align(
    star_path,
    trim_fastq,
    subsample_script_path,
    main_output_dir: pathlib.Path,
    short_read_fastq_dir,
    genome_file: Union[pathlib.Path, str],
    max_reads_per_sample,
    max_total_reads,
    max_intron_length: int,
    num_threads: int,
):
    # TODO
    # !!! Need to add in samtools path above instead of just using 'samtools' in command

    # TODO
    # use and pass Path objects directly as arguments to the function
    short_read_fastq_dir = pathlib.Path(short_read_fastq_dir)

    if not star_path:
        star_path = "STAR"

    check_exe(star_path)

    # If trimming has been enabled then switch the path for short_read_fastq_dir from the original location to the trimmed fastq dir
    if trim_fastq:
        short_read_fastq_dir = main_output_dir / "trim_galore_output"

    #  if not os.path.exists(subsample_script_path):
    subsample_script_path = "subsample_fastq.py"

    star_dir = create_dir(main_output_dir, "star_output")

    log_final_out = star_dir / "Log.final.out"
    log_out = star_dir / "Log.out"
    log_progress_out = star_dir / "Log.progress.out"
    if log_final_out.is_file() and log_out.is_file() and log_progress_out.is_file():
        logger.info("Skipping analysis, log files already exist")
        return

    star_tmp_dir = star_dir / "tmp"
    # delete directory if already exists
    shutil.rmtree(star_tmp_dir, ignore_errors=True)

    star_index_file = star_dir / "SAindex"

    fastq_file_list = []
    file_types = ["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"]
    for file_type in file_types:
        fastq_file_list.extend(short_read_fastq_dir.glob(file_type))

    # This works out if the files are paired or not
    fastq_file_list = create_paired_paths(fastq_file_list)

    # Subsamples in parallel if there's a value set
    if max_reads_per_sample:
        pool = multiprocessing.Pool(num_threads)
        for fastq_files in fastq_file_list:
            fastq_file = fastq_files[0]
            fastq_file_pair = ""
            if len(fastq_files) == 2:
                fastq_file_pair = fastq_files[1]

            if (
                fastq_file_pair
                and os.path.exists(f"{fastq_file}.sub")
                and os.path.exists(f"{fastq_file_pair}.sub")
            ):
                logger.info(
                    "Found existing .sub files on the fastq path for both members of the pair, will use those instead of subsampling again:\n%s\n%s"
                    % f"{fastq_file}.sub",
                    f"{fastq_file_pair}.sub",
                )
            elif fastq_file_pair:
                pool.apply_async(
                    run_subsample_script,
                    args=(
                        fastq_file,
                        fastq_file_pair,
                        subsample_script_path,
                    ),
                )
            elif os.path.exists(f"{fastq_file}.sub"):
                logger.info(
                    "Found an existing .sub file on the fastq path, will use that instead:\n%s"
                    % f"{fastq_file}.sub"
                )
            else:
                pool.apply_async(
                    run_subsample_script,
                    args=(
                        fastq_file,
                        fastq_file_pair,
                        subsample_script_path,
                    ),
                )

        pool.close()
        pool.join()

    fastq_file_list = check_for_fastq_subsamples(fastq_file_list)

    if not fastq_file_list:
        raise IndexError(f"Empty fastq files list. Fastq dir:\n{short_read_fastq_dir}")

    if not star_index_file.is_file():
        logger.info("STAR index file does not exist, generating new one.")
        seq_region_lengths = get_seq_region_lengths(genome_file, 0)
        genome_size = sum(seq_region_lengths.values())
        index_bases = min(14, math.floor((math.log(genome_size, 2) / 2) - 1))
        subprocess.run(
            [
                star_path,
                "--runThreadN",
                str(num_threads),
                "--runMode",
                "genomeGenerate",
                "--outFileNamePrefix",
                f"{star_dir}/",
                "--genomeDir",
                star_dir,
                "--genomeSAindexNbases",
                str(index_bases),
                "--genomeFastaFiles",
                genome_file,
            ]
        )
    if not star_index_file.is_file():
        raise IOError(f"STAR index file failed to be generated at:\n{star_index_file}")

    logger.info("Running STAR on the files in the fastq dir")
    for fastq_file_path in fastq_file_list:
        logger.info("fastq_file_path: %s" % fastq_file_path)
        fastq_file_name = os.path.basename(fastq_file_path)
        check_compression = re.search(r".gz$", fastq_file_name)

        # If there's a tmp dir already, the most likely cause is that STAR failed on the previous input file(s)
        # In this case STAR would effectively break for all the rest of the files, as it won't run if the tmp
        # dir exists. So clean it up and put out a warning. Another approach would be to just name the tmp dir
        # uniquely. Also there should be code checking the return on STAR anyway. So later when this is being
        # cleaned up to have proper tests, need to decide on the best implementation
        if star_tmp_dir.exists():
            logger.error(
                "Found existing tmp dir, implies potential failure on previous file, deleting"
            )
            shutil.rmtree(star_tmp_dir)

        sam_file_path = star_dir / f"{fastq_file_name}.sam"
        sam_temp_file_path = star_dir / f"{sam_file_path.name}.tmp"
        bam_sort_file_path = sam_file_path.with_suffix(".bam")

        if bam_sort_file_path.is_file() and bam_sort_file_path.stat().st_size > 0:
            logger.info(
                "Existing bam file for fastq file found, skipping processing file"
            )
            continue

        logger.info("Processing %s" % fastq_file_path)
        #    star_command = [star_path,'--outFilterIntronMotifs','RemoveNoncanonicalUnannotated','--outSAMstrandField','intronMotif','--runThreadN',str(num_threads),'--twopassMode','Basic','--runMode','alignReads','--genomeDir',star_dir,'--readFilesIn',fastq_file_path,'--outFileNamePrefix',(star_dir + '/'),'--outTmpDir',star_tmp_dir,'--outSAMtype','SAM','--alignIntronMax',str(max_intron_length),'--outSJfilterIntronMaxVsReadN','5000','10000','25000','40000','50000','50000','50000','50000','50000','100000']

        star_command = [
            star_path,
            "--outFilterIntronMotifs",
            "RemoveNoncanonicalUnannotated",
            "--outSAMstrandField",
            "intronMotif",
            "--runThreadN",
            str(num_threads),
            "--twopassMode",
            "Basic",
            "--runMode",
            "alignReads",
            "--genomeDir",
            star_dir,
            "--readFilesIn",
            fastq_file_path,
            "--outFileNamePrefix",
            f"{star_dir}/",
            "--outTmpDir",
            star_tmp_dir,
            "--outSAMtype",
            "SAM",
            "--alignIntronMax",
            str(max_intron_length),
        ]

        if check_compression:
            star_command.append("--readFilesCommand")
            star_command.append("gunzip")
            star_command.append("-c")

        subprocess.run(star_command)
        (star_dir / "Aligned.out.sam").rename(sam_file_path)
        junctions_file_path = star_dir / f"{fastq_file_name}.sj.tab"
        (star_dir / "SJ.out.tab").rename(junctions_file_path)

        logger.info("Converting samfile into sorted bam file:\n%s" % bam_sort_file_path)
        subprocess.run(
            [
                "samtools",
                "sort",
                "-@",
                str(num_threads),
                "-T",
                sam_temp_file_path,
                "-o",
                bam_sort_file_path,
                sam_file_path,
            ]
        )

        os.remove(sam_file_path)

    logger.info("Completed running STAR")


def check_for_fastq_subsamples(fastq_file_list):
    # This should probably removed at some point as it is needlessly complicated
    # Would be better to just build into the previous step
    # Mainly just about making sure that if you have subsamples they're used and if you have pairs they're paired
    for idx, fastq_files in enumerate(fastq_file_list):
        fastq_file = fastq_files[0]
        subsample_file = f"{fastq_file}.sub"

        fastq_file_pair = ""
        subsample_file_pair = ""
        if len(fastq_files) == 2:
            fastq_file_pair = fastq_files[1]
            subsample_file_pair = f"{fastq_file_pair}.sub"

        # This bit will replace the list entry with a string, don't need a list after this function for each pair/file
        if os.path.exists(subsample_file):
            logger.info(
                "Found a subsampled file extension, will use that instead of the original file:\n%s"
                % subsample_file
            )
            fastq_file_list[idx] = subsample_file
        else:
            fastq_file_list[idx] = fastq_file

        # This bit just concats the paired file (or subsampled paired file) if it exists
        if os.path.exists(subsample_file_pair):
            logger.info(
                "Found a subsampled paired file extension, will use that instead of the original file:\n%s"
                % subsample_file_pair
            )
            fastq_file_list[idx] = f"{subsample_file},{subsample_file_pair}"
        elif fastq_file_pair:
            fastq_file_list[idx] = f"{fastq_file},{fastq_file_pair}"

        logger.info("Entry at current index:\n%s" % fastq_file_list[idx])

    return fastq_file_list


def run_subsample_script(fastq_file, fastq_file_pair, subsample_script_path):
    if fastq_file_pair:
        subprocess.run(
            [
                "python3",
                subsample_script_path,
                "--fastq_file",
                fastq_file,
                "--fastq_file_pair",
                fastq_file_pair,
            ]
        )
    else:
        subprocess.run(["python3", subsample_script_path, "--fastq_file", fastq_file])
