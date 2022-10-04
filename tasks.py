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
import glob
import math
import multiprocessing
import os
import re
import shutil
import subprocess

# project imports
from utils import create_dir, logger


def run_star_align(
    star_path,
    trim_fastq,
    subsample_script_path,
    main_output_dir,
    short_read_fastq_dir,
    genome_file,
    max_reads_per_sample,
    max_total_reads,
    max_short_read_intron_length,
    num_threads,
):
    # !!! Need to add in samtools path above instead of just using 'samtools' in command

    if not star_path:
        star_path = "STAR"

    check_exe(star_path)

    # If trimming has been enabled then switch the path for short_read_fastq_dir from the original location to the trimmed fastq dir
    if trim_fastq:
        short_read_fastq_dir = os.path.join(main_output_dir, "trim_galore_output")

    #  if not os.path.exists(subsample_script_path):
    subsample_script_path = "subsample_fastq.py"

    star_dir = create_dir(main_output_dir, "star_output")

    logger.info("Skip analysis if the final log file already exists")
    log_final_out = Path(os.path.join(star_dir, "Log.final.out"))
    log_out = Path(os.path.join(star_dir, "Log.out"))
    log_progress_out = Path(os.path.join(star_dir, "Log.progress.out"))
    if log_final_out.is_file() and log_out.is_file() and log_progress_out.is_file():
        logger.info("Star gtf file exists")
        return
    else:
        logger.info("No log files, go on with the analysis")

    star_tmp_dir = os.path.join(star_dir, "tmp")
    if os.path.exists(star_tmp_dir):
        subprocess.run(["rm", "-rf", star_tmp_dir])

    star_index_file = os.path.join(star_dir, "SAindex")

    fastq_file_list = []
    file_types = ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz")
    for file_type in file_types:
        fastq_file_list.extend(glob.glob(os.path.join(short_read_fastq_dir, file_type)))

    # This works out if the files are paired or not
    fastq_file_list = create_paired_paths(fastq_file_list)

    # Subsamples in parallel if there's a value set
    if max_reads_per_sample:
        pool = multiprocessing.Pool(int(num_threads))
        for fastq_files in fastq_file_list:
            fastq_file = fastq_files[0]
            fastq_file_pair = ""
            if len(fastq_files) == 2:
                fastq_file_pair = fastq_files[1]

            if (
                fastq_file_pair
                and os.path.exists(fastq_file + ".sub")
                and os.path.exists(fastq_file_pair + ".sub")
            ):
                logger.info(
                    "Found an existing .sub files on the fastq path for both members of the pair, will use those instead of subsampling again. Files:"
                )
                logger.info(fastq_file + ".sub")
                logger.info(fastq_file_pair + ".sub")
            elif fastq_file_pair:
                pool.apply_async(
                    run_subsample_script,
                    args=(
                        fastq_file,
                        fastq_file_pair,
                        subsample_script_path,
                    ),
                )
            elif os.path.exists(fastq_file + ".sub"):
                logger.info(
                    "Found an existing .sub file on the fastq path, will use that instead. File:"
                )
                logger.info(fastq_file + ".sub")
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
        raise IndexError(
            "The list of fastq files is empty. Fastq dir:\n%s" % short_read_fastq_dir
        )

    if not os.path.exists(star_index_file):
        logger.info("Did not find an index file for Star. Will create now")
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
                (star_dir + "/"),
                "--genomeDir",
                star_dir,
                "--genomeSAindexNbases",
                str(index_bases),
                "--genomeFastaFiles",
                genome_file,
            ]
        )

    if not star_index_file:
        raise IOError(
            "The index file does not exist. Expected path:\n%s" % star_index_file
        )

    logger.info("Running Star on the files in the fastq dir")
    for fastq_file_path in fastq_file_list:
        logger.info(fastq_file_path)
        fastq_file_name = os.path.basename(fastq_file_path)
        check_compression = re.search(r".gz$", fastq_file_name)

        # If there's a tmp dir already, the most likely cause is that STAR failed on the previous input file(s)
        # In this case STAR would effectively break for all the rest of the files, as it won't run if the tmp
        # dir exists. So clean it up and put out a warning. Another approach would be to just name the tmp dir
        # uniquely. Also there should be code checking the return on STAR anyway. So later when this is being
        # cleaned up to have proper tests, need to decide on the best implementation
        if os.path.exists(star_tmp_dir):
            logger.error(
                "Found an existing tmp dir, implies potential failure on previous file. Removing tmp dir"
            )
            try:
                shutil.rmtree(star_tmp_dir)
            except OSError as e:
                logger.error("Error: %s - %s." % (e.filename, e.strerror))

        sam_file_path = os.path.join(star_dir, (fastq_file_name + ".sam"))
        junctions_file_path = os.path.join(star_dir, (fastq_file_name + ".sj.tab"))
        sam_file_name = os.path.basename(sam_file_path)
        sam_temp_file_path = os.path.join(star_dir, (sam_file_name + ".tmp"))
        bam_sort_file_path = os.path.join(
            star_dir, re.sub(".sam", ".bam", sam_file_name)
        )

        if (
            os.path.isfile(bam_sort_file_path)
            and not os.stat(bam_sort_file_path).st_size == 0
        ):
            logger.info(
                "Found an existing bam file for the fastq file, presuming the file has been processed, will skip"
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
            (star_dir + "/"),
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
        subprocess.run(["mv", os.path.join(star_dir, "Aligned.out.sam"), sam_file_path])
        subprocess.run(
            ["mv", os.path.join(star_dir, "SJ.out.tab"), junctions_file_path]
        )

        logger.info("Converting samfile into sorted bam file. Bam file:")
        logger.info(bam_sort_file_path)
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

        logger.info("Removing sam file")
        subprocess.run(["rm", sam_file_path])

    logger.info("Completed running STAR")
