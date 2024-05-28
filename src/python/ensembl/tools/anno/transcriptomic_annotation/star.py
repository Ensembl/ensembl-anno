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
"""The STAR (Spliced Transcripts Alignment to a Reference) alignment tool is widely used
in genomics research for aligning RNA-seq data to a reference genome.
Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner.
Bioinformatics. 2013;29(1):15-21. doi:10.1093/bioinformatics/bts635
"""
__all__ = ["run_star", "subsample_transcriptomic_data", "run_trimming"]

import argparse
import logging
import logging.config
import gzip
import math
import multiprocessing
from pathlib import Path
import random
import re
import shutil
import subprocess
from typing import Dict, List

from ensembl.tools.anno.utils._utils import (
    check_exe,
    create_dir,
    check_gtf_content,
    get_seq_region_length,
)


def run_star(  # pylint:disable=too-many-branches
    genome_file: Path,
    output_dir: Path,
    short_read_fastq_dir: Path,
    delete_pre_trim_fastq: bool = False,
    trim_fastq: bool = False,
    max_reads_per_sample: int = 0,
    max_intron_length: int = 100000,
    subsample_read_limit: int = 100000000,
    subsample_percentage: float = 0.25,
    sampling_via_read_limit: bool = False,
    sampling_via_percentage: bool = False,
    sampling_via_read_limit_percentage: bool = False,
    num_threads: int = 1,
    star_bin: Path = Path("star"),
    samtools_bin: Path = Path("samtools"),
    trim_galore_bin: Path = Path("trim_galore"),
) -> None:
    """
    Run STAR alignment on list of short read data.
        :param genome_file: Genome file path.
        :type genome_file: Path
        :param output_dir: Working directory path.
        :type output_dir: Path
        :param short_read_fastq_dir: Short read directory path.
        :type short_read_fastq_dir: Path
        :param delete_pre_trim_fastq: Delete the original fastq files after trimming. Defaults to False.
        :type delete_pre_trim_fastq: boolean, default False
        :param trim_fastq: Trim short read files using TrimGalore. Defaults to False.
        :type trim_fastq: boolean, default False
        :param max_reads_per_sample: Max number of reads per sample. Defaults to 0 (unlimited).
        :type max_reads_per_sample: int, default 0
        :param max_intron_length: The maximum intron size for alignments. Defaults to 100000.
        :type max_intron_length: int, default 100000
        :param subsample_read_limit: Maximum number of reads to subsample. Defaults to 100000000.
        :type subsample_read_limit:int, default 100000000,
        :param subsample_percentage: Maximun percentage of reads to subsample.
        :type subsample_percentage: int, default 0.25,
        :param sampling_via_read_limit: subsample fastq files using --subsample_read_limit.
        :type sampling_via_read_limit : boolean, False,
        :param sampling_via_percentage: subsample fastq files using --subsample_percentage.
        :type sampling_via_percentage : boolean, False,
        :param sampling_via_read_limit_percentage: use max read limit and  percentage value.
        :type sampling_via_read_limit_percentage : boolean, False,
        :param num_threads: Number of available threads.
        :type num_threads: int, default 1
        :param star_bin: Software path.
        :type star_bin: Path, default star
        :param samtools_bin: Software path.
        :type samtools_bin: Path,default samtools
        :param trim_galore_bin: Software path.
        :type trim_galore_bin: Path, default trim_galore

        :return: None
        :rtype: None
    """
    check_exe(star_bin)
    # If trimming has been enabled then switch the path for
    # short_read_fastq_dir from the original location to the trimmed fastq dir
    if trim_fastq:
        run_trimming(output_dir, short_read_fastq_dir, delete_pre_trim_fastq, num_threads, trim_galore_bin)
        short_read_fastq_dir = output_dir / "trim_galore_output"

    #  if not os.path.exists(subsample_script_path):
    # subsample_script_path = "subsample_fastq.py"

    star_dir = create_dir(output_dir, "star_output")

    for output_file in [
        Path(f"{output_dir}/stringtie_output/annotation.gtf"),
        Path(f"{output_dir}/scallop_output/annotation.gtf"),
    ]:
        if output_file.exists():
            transcript_count = check_gtf_content(output_file, "transcript")  # check a gtf
            if transcript_count > 0:
                logging.info("Transcriptomic alignment exists")
                return

    star_index_file = star_dir / "SAindex"
    fastq_file_list: List[Path] = []
    file_types = ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz")
    fastq_file_list = [
        path for file_type in file_types for path in Path(short_read_fastq_dir).rglob(file_type)
    ]
    if len(fastq_file_list) == 0:
        raise IndexError(f"The list of fastq files is empty. Fastq dir:\n{short_read_fastq_dir}")

    # for file_type in file_types:
    #    fastq_file_list.extend(glob.glob(os.path.join(short_read_fastq_dir, file_type)))

    # Get list of paired paths
    paired_fastq_file_list = _create_paired_paths(fastq_file_list)
    # Subsamples in parallel if there's a value set
    if max_reads_per_sample:
        paired_fastq_file_list_sub: List[List[Path]] = []
        for paired_fastq_files in paired_fastq_file_list:
            paired_fastq_files_sub = subsample_transcriptomic_data(
                paired_fastq_files,
                subsample_read_limit,
                subsample_percentage,
                sampling_via_read_limit,
                sampling_via_percentage,
                sampling_via_read_limit_percentage,
                num_threads,
            )
            paired_fastq_file_list_sub.append(paired_fastq_files_sub)
            paired_fastq_file_list = paired_fastq_file_list_sub
        # Get the list of the new subsampled files
        # fastq_file_list = [
        #    path for file_type in file_types for path in Path(short_read_fastq_dir).rglob(file_type)
        # ]
    # I don't think is needed
    # fastq_file_list = check_for_fastq_subsamples(fastq_file_list)

    if not star_index_file.exists():
        logging.info("Did not find an index file for Star. Will create now")
        seq_region_to_length = get_seq_region_length(genome_file, 0)
        genome_size = sum(seq_region_to_length.values())
        # This calculates the base-2 logarithm of the genome_size. The logarithm of the genome size is
        # a measure of how many bits are needed to represent the genome size in binary.
        #
        # The choice of 14 as the maximum value is likely based on empirical observations and optimization
        # considerations. Too large of a seed length can lead to increased memory usage and potentially
        # slower indexing, while a seed length that is too small might affect alignment accuracy.
        index_bases = min(14, math.floor((math.log(genome_size, 2) / 2) - 1))
        try:
            subprocess.run(  # pylint:disable=subprocess-run-check
                [
                    str(star_bin),
                    "--runThreadN",
                    str(num_threads),
                    "--runMode",
                    "genomeGenerate",
                    "--outFileNamePrefix",
                    f"{star_dir}/",
                    "--genomeDir",
                    str(star_dir),
                    "--genomeSAindexNbases",
                    str(index_bases),
                    "--genomeFastaFiles",
                    str(genome_file),
                ]
            )
        except Exception as e:  # pylint:disable=broad-exception-caught
            logging.error("An error occurred while creating star index: %s", e)

    logging.info("Running Star on the files in the fastq dir")
    for paired_files in paired_fastq_file_list:
        first_file_name = paired_files[0].name  # Get the name of the first file
        match = re.search(r"(.+)_\d+\.(fastq|fq)", first_file_name)  # Search for pattern
        if match:
            first_part_of_name = match.group(1)
        # logger.info(fastq_file_path)
        # fastq_file_name = os.path.basename(fastq_file_path)
        star_tmp_dir = star_dir / "tmp"
        if star_tmp_dir.exists():
            shutil.rmtree(star_tmp_dir)
        sam_file = Path(f"{star_dir}/{first_part_of_name}.sam")
        junctions_file = Path(f"{star_dir}/{first_part_of_name}.sj.tab")
        sam_file_name = sam_file.name
        sam_temp_file = Path(f"{star_dir}/{sam_file_name}.tmp")
        bam_file = re.sub(".sam", ".bam", sam_file_name)
        bam_sort_file = Path(f"{star_dir}/{bam_file}")
        log_out_file = Path(f"{star_dir}/{first_part_of_name}.Log.final.out")
        if log_out_file.exists() and bam_sort_file.exists() and bam_sort_file.stat().st_size != 0:
            logging.info(
                "Found an existing bam file for the fastq file, \
                presuming the file has been processed, will skip"
            )
            continue

        read_files_in = ",".join(str(fastq_file) for fastq_file in paired_files)

        logging.info("Processing %s", list(paired_files))
        star_command = [
            str(star_bin),
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
            str(star_dir),
            "--readFilesIn",
            read_files_in,
            "--outFileNamePrefix",
            f"{star_dir}/",
            "--outTmpDir",
            str(star_tmp_dir),
            "--outSAMtype",
            "SAM",
            "--alignIntronMax",
            str(max_intron_length),
        ]
        #'--outSJfilterIntronMaxVsReadN','5000','10000','25000','40000',
        #'50000','50000','50000','50000','50000','100000']
        # check_compression = re.search(r".gz$", fastq_file)
        if Path(paired_files[0].name).suffix.endswith(".gz"):
            star_command.append("--readFilesCommand")
            star_command.append("gunzip")
            star_command.append("-c")
        subprocess.run(star_command)  # pylint:disable=subprocess-run-check
        shutil.move(Path(f"{star_dir}/Aligned.out.sam"), sam_file)
        shutil.move(Path(f"{star_dir}/SJ.out.tab"), junctions_file)
        logging.info("Converting samfile into sorted bam file. Bam file: %s", bam_file)
        subprocess.run(  # pylint:disable=subprocess-run-check
            [
                str(samtools_bin),
                "sort",
                "-@",
                str(num_threads),
                "-T",
                str(sam_temp_file),
                "-o",
                str(bam_sort_file),
                str(sam_file),
            ]
        )
        shutil.move(star_dir / "Log.final.out", log_out_file)
        sam_file.unlink()
    logging.info("Completed running STAR")


def _create_paired_paths(fastq_file_paths: List[Path]) -> List[List[Path]]:
    """
    Create list of paired transcriptomic fastq files

    Args:
        fastq_file_paths (List): List of transcriptomic file paths.

    Returns:
        List[List[Path]]: List of paired transcriptomic files
    """
    path_dict: Dict[str, List[Path]] = {}
    final_list: List[List[Path]] = []
    for fastq_file in fastq_file_paths:
        paired_name = re.search(r"(.+)_\d+\.(fastq|fq)", str(fastq_file))
        if not paired_name:
            logging.exception(
                "Could not find _1 or _2 at the end of the prefix \
                for file. Assuming file is not paired: %s",
                str(fastq_file),
            )
            final_list.append([fastq_file])
            continue
        run_accession = paired_name.group(1)
        if run_accession in path_dict:
            path_dict[run_accession].append(fastq_file)
        else:
            path_dict[run_accession] = [fastq_file]
    for pair in path_dict:  # pylint:disable=consider-using-dict-items
        final_list.append(path_dict[pair])
    return final_list


# pylint:disable=pointless-string-statement
"""
For an advanced and optimised subsampling we could use 
https://github.com/lh3/seqtk 
"""


def subsample_transcriptomic_data(
    fastq_file_list: List[Path],
    subsample_read_limit: int = 100000000,
    subsample_percentage: float = 0.25,
    sampling_via_read_limit: bool = False,
    sampling_via_percentage: bool = False,
    sampling_via_read_limit_percentage: bool = True,
    num_threads: int = 2,
) -> List[Path]:
    """
    Subsample list of paired files.
    Args:
        fastq_file_list : Subsample paired fastq files.
        subsample_read_limit : Maximum number of reads to subsample, default to 100000000.
        subsample_percentage : Maximun percentage of reads to subsample, default to 0.25.
        sampling_via_read_limit : If True will subsample an input dataset of fastq files \
        using --subsample_read_limit value.
        sampling_via_percentage : If True will subsample an input dataset of fastq files \
        using --subsample_percentage value.
        sampling_via_read_limit_percentage : If True will subsample an input dataset of \
        fastq files using --subsample_read_limit and --subsample_percentage value; \
        the lowest number of reads is taken.
        num_threads : number of threads
    Returns:
        List[Path]: List of subsampled paired transcriptomic files
    """
    subsampled_fastq_files: List[Path] = []
    # for fastq_files in fastq_file_list:
    # fastq_file_1, fastq_file_2 = fastq_files
    # fastq_file_pair = ""
    # if len(fastq_files) == 2:
    #    fastq_file_pair = fastq_files[1]

    if len(fastq_file_list) == 1:
        fastq_file_1 = fastq_file_list[0]
        if Path(f"{fastq_file_1}.sub").exists():
            logging.info(
                "Found an existing .sub file on the fastq path, will use that \
                    instead. File:%s.sub",
                fastq_file_1,
            )
        else:
            _subsample_paired_fastq_files(
                fastq_file_list,
                subsample_read_limit,
                subsample_percentage,
                sampling_via_read_limit,
                sampling_via_percentage,
                sampling_via_read_limit_percentage,
                num_threads,
            )
        subsampled_fastq_files = [Path(f"{fastq_file_1}.sub")]
    if len(fastq_file_list) == 2:
        fastq_file_1, fastq_file_2 = fastq_file_list
        if Path(f"{fastq_file_1}.sub").exists() and Path(f"{fastq_file_2}.sub").exists():
            logging.info(
                "Found an existing .sub files on the fastq path for both members of the pair, will use \
                those instead of subsampling again. Files: %s.sub,%s.sub",
                fastq_file_1,
                fastq_file_2,
            )
        else:
            _subsample_paired_fastq_files(
                fastq_file_list,
                subsample_read_limit,
                subsample_percentage,
                sampling_via_read_limit,
                sampling_via_percentage,
                sampling_via_read_limit_percentage,
                num_threads,
            )
        subsampled_fastq_files = [Path(f"{fastq_file_1}.sub"), Path(f"{fastq_file_2}.sub")]
    return subsampled_fastq_files


def _subsample_paired_fastq_files(  # pylint:disable=too-many-branches
    fastq_files: List[Path],
    subsample_read_limit: int = 100000000,
    subsample_percentage: float = 0.25,
    sampling_via_read_limit: bool = False,
    sampling_via_percentage: bool = False,
    sampling_via_read_limit_percentage: bool = True,
    num_threads: int = 2,
) -> None:
    """
    Perform subsampling on two paired FastQ files in parallel using multiple threads.

    Args:
        fastq_files : Path for paired fastq files.
        subsample_read_limit : Maximum number of reads to subsample, default to 100000000.
        subsample_percentage : Maximun percentage of reads to subsample, default to 0.25.
        sampling_via_read_limit : If True will subsample an input dataset of fastq files \
        using --subsample_read_limit value.
        sampling_via_percentage : If True will subsample an input dataset of fastq files \
        using --subsample_percentage value.
        sampling_via_read_limit_percentage : If True will subsample an input dataset of \
        fastq files using --subsample_read_limit and --subsample_percentage value; \
        the lowest number of reads is taken.
        num_threads : Number of threads, defaults to 2.
        compressed : file compressed, defaults to False.
    """

    if len(fastq_files) == 2:
        fastq_file_1, fastq_file_2 = fastq_files
        output_file_1, output_file_2 = [Path(f"{fastq_file_1}.sub"), Path(f"{fastq_file_2}.sub")]
    elif len(fastq_files) == 1:
        fastq_file_1 = fastq_files[0]
        output_file_1 = Path(f"{fastq_file_1}.sub")
    else:
        raise FileNotFoundError("No fastq file found")

    if fastq_file_1.suffix.endswith(".gz$"):
        compressed = True
        num_lines = sum(1 for line in gzip.open(fastq_file_1))  # pylint:disable=consider-using-with
    else:
        compressed = False
        num_lines = sum(1 for line in open(fastq_file_1))  # pylint:disable=consider-using-with

    range_limit = int(num_lines / 4)
    sampling_size = 0
    if sampling_via_read_limit and subsample_read_limit:
        sampling_size = subsample_read_limit
    elif sampling_via_percentage and subsample_percentage:
        sampling_size = round(range_limit * subsample_percentage)
    elif sampling_via_read_limit_percentage and subsample_percentage and subsample_read_limit:
        sampling_size = min(subsample_read_limit, round(range_limit * subsample_percentage))

    if range_limit <= sampling_size:
        logging.info(
            "Number of reads (%s is less than the max allowed read count (%s), \
            no need to subsample",
            str(range_limit),
            str(sampling_size),
        )
        return

    rand_list = random.sample(range(0, range_limit - 1), sampling_size)
    random_indices = {idx * 4: 1 for idx in rand_list}
    logging.info("Processing paired files in parallel")
    if num_threads >= 2:
        pool = multiprocessing.Pool(int(num_threads))  # pylint:disable=consider-using-with
        pool.apply_async(
            _subsample_fastq_subset,
            args=(
                fastq_file_1,
                output_file_1,
                random_indices,
                compressed,
            ),
        )
        if len(fastq_files) == 2:
            pool.apply_async(
                _subsample_fastq_subset,
                args=(
                    fastq_file_2,
                    output_file_2,
                    random_indices,
                    compressed,
                ),
            )

        pool.close()
        pool.join()
    else:
        _subsample_fastq_subset(
            fastq_file_1,
            output_file_1,
            random_indices,
            compressed,
        )
        if len(fastq_files) == 2:
            _subsample_fastq_subset(
                fastq_file_2,
                output_file_2,
                random_indices,
                compressed,
            )


def _subsample_fastq_subset(
    fastq_file: Path, output_file: Path, random_indices: dict, compressed: bool
) -> None:
    """
    Selecting specific sets of four lines from an input FastQ file and writing them to an output file.
    Args:
        fastq_file : Path for the fastq file.
        output_file : Path for the output file.
        random_indices : set of random indices.
        compressed : the files is compressed
    """
    line_index = 0

    with gzip.open(fastq_file, "rt") if compressed else open(fastq_file) as file_in, open(
        output_file, "w+"
    ) as file_out:
        lines = [file_in.readline() for _ in range(4)]
        while lines[3]:  # This ensures that the loop continues until the end of the input file.
            if line_index in random_indices:
                file_out.writelines(lines)
            line_index += 4
            lines = [file_in.readline() for _ in range(4)]


def run_trimming(
    output_dir: Path,
    short_read_fastq_dir: Path,
    delete_pre_trim_fastq: bool = False,
    num_threads: int = 1,
    trim_galore_bin="trim_galore",
) -> None:
    """
    Trim list of short read fastq files.
    Args:
        output_dir : Working directory path.
        short_read_fastq_dir : Short read directory path.
        delete_pre_trim_fastq : Removing original fastq file post trimming. Defaults to False.
        num_threads : Number of threads.
        trim_galore_bin : Software path.
    """
    check_exe(trim_galore_bin)
    trim_dir = create_dir(output_dir, "trim_galore_output")

    fastq_file_list: List[Path] = []
    file_types = ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz")
    fastq_file_list = [
        path for file_type in file_types for path in Path(short_read_fastq_dir).rglob(file_type)
    ]
    paired_fastq_file_list = _create_paired_paths(fastq_file_list)

    trim_galore_cmd = [
        str(trim_galore_bin),
        "--illumina",
        "--quality",
        "20",
        "--length",
        "50",
        "--output_dir",
        str(trim_dir),
    ]

    pool = multiprocessing.Pool(int(num_threads))  # pylint:disable=consider-using-with
    for fastq_file in paired_fastq_file_list:
        file1, file2 = fastq_file
        if delete_pre_trim_fastq:
            file1.unlink()
            file2.unlink()
        pool.apply_async(
            multiprocess_trim_galore,
            args=(
                trim_galore_cmd,
                fastq_file,
                trim_dir,
            ),
        )

    pool.close()
    pool.join()

    trimmed_fastq_list = trim_dir.glob("*.fq.gz")

    for trimmed_fastq_path in trimmed_fastq_list:
        logging.info("Trimmed file path: %s", str(trimmed_fastq_path))
        sub_patterns = re.compile(r"|".join(("_val_1.fq", "_val_2.fq", "_trimmed.fq")))
        updated_file_path_name = sub_patterns.sub(".fq", trimmed_fastq_path.name)
        updated_file_path = short_read_fastq_dir / updated_file_path_name
        logging.info("Updated file path: %s", str(updated_file_path))
        trimmed_fastq_path.rename(updated_file_path)

    files_to_delete_list: List[Path] = []
    for file_type in file_types:
        files_to_delete_list.extend(short_read_fastq_dir.glob(file_type))

    for file_to_delete in files_to_delete_list:
        file_to_delete.unlink()


def multiprocess_trim_galore(trim_galore_cmd: List, fastq_paired_files: List[Path]) -> None:
    """
    Trim short paired or single short read fastq file.
    Args:
        trim_galore_cmd : Generic command.
        fastq_paired_files : List of single or paired fastq files.
    """

    fastq_file = fastq_paired_files[0]
    fastq_file_pair = None

    if len(fastq_paired_files) == 2:
        fastq_file, fastq_file_pair = fastq_paired_files
        trim_galore_cmd.append("--paired")
        trim_galore_cmd.append(fastq_file)
        trim_galore_cmd.append(fastq_file_pair)
    elif len(fastq_paired_files) == 1:
        trim_galore_cmd.append(fastq_paired_files)

    logging.info("Running Trim Galore with the following command: %s", {" ".join(trim_galore_cmd)})
    subprocess.run(trim_galore_cmd, check=True)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="STAR's arguments")
    parser.add_argument("--genome_file", required=True, help="Genome file path")
    parser.add_argument("--output_dir", required=True, help="Output directory path")
    parser.add_argument("--short_read_fastq_dir", required=True, help="Short read directory path")
    parser.add_argument(
        "--delete_pre_trim_fastq",
        action="store_true",
        default=False,
        help="Delete the original fastq files after trimming",
    )
    parser.add_argument(
        "--trim_fastq", action="store_true", default=False, help="Trim the short read files using Trim Galore"
    )
    parser.add_argument(
        "--max_reads_per_sample", type=int, default=0, help="The maximum number of reads to use per sample"
    )
    parser.add_argument(
        "--max_intron_length", type=int, default=100000, help="The maximum intron size for alignments"
    )
    parser.add_argument("--num_threads", type=int, default=1, help="Number of threads")
    parser.add_argument("--star_bin", default="STAR", help="Star software path")
    parser.add_argument("--samtools_bin", default="samtools", help="Samtools software path")
    parser.add_argument("--trim_galore_bin", default="trim_galore", help="Trim Galore software path")
    parser.add_argument(
        "--subsample_read_limit",
        type=int,
        default=100000000,
        help="Maximum number of reads to subsample. Default 1 hundred million reads",
        required=False,
    )
    parser.add_argument(
        "--subsample_percentage",
        type=float,
        default=0.25,
        help="Maximun percentage of reads to subsample (0 to 1)",
        required=False,
    )
    parser.add_argument(
        "--sampling_via_read_limit",
        type=bool,
        default=False,
        help="If True will subsample an input dataset of fastq files using --subsample_read_limit value.",
        required=False,
    )
    parser.add_argument(
        "--sampling_via_percentage",
        type=bool,
        default=False,
        help="If True will subsample an input dataset of fastq files using --subsample_percentage value.",
        required=False,
    )
    parser.add_argument(
        "--sampling_via_read_limit_percentage",
        type=bool,
        default=True,
        help="If True will subsample an input dataset of fastq files using --subsample_read_limit \
            and --subsample_percentage value; the lowest number of reads is taken.",
        required=False,
    )
    parser.add_argument(
        "--paired_file_1",
        help="Path for single or paired fastq file; used when --run_subsampling \
            or --run_trimming are enabled.",
        required=False,
    )
    parser.add_argument(
        "--paired_file_2",
        help="Path for single or paired fastq file; used when --run_subsampling \
            or --run_trimming are enabled.",
        required=False,
    )
    parser.add_argument(
        "--run_star",
        type=bool,
        default=True,
        help="If True will run STAR alignment given an input dataset of fastq files.",
        required=False,
    )
    parser.add_argument(
        "--run_subsampling",
        type=bool,
        default=False,
        help="If True will subsample an input dataset of fastq files.",
        required=False,
    )
    parser.add_argument(
        "--run_trimming",
        type=bool,
        default=False,
        help="If True will trim input dataset of fastq files using TrimGalore suite.",
        required=False,
    )
    return parser.parse_args()


def main():
    """STAR's entry-point."""
    args = parse_args()

    log_file_path = create_dir(args.output_dir, "log") / "star.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"

    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )
    if args.run_star:
        run_star(
            args.genome_file,
            args.output_dir,
            args.short_read_fastq_dir,
            args.delete_pre_trim_fastq,
            args.trim_fastq,
            args.max_reads_per_sample,
            args.max_intron_length,
            args.subsample_read_limit,
            args.subsample_percentage,
            args.sampling_via_read_limit,
            args.sampling_via_percentage,
            args.sampling_via_read_limit_percentage,
            args.num_threads,
            args.star_bin,
            args.samtools_bin,
            args.trim_galore_bin,
        )
    if args.run_subsampling:
        fastq_file_list = [Path(args.paired_file_1), Path(args.paired_file_2)]
        subsample_transcriptomic_data(
            fastq_file_list,
            args.subsample_read_limit,
            args.subsample_percentage,
            args.sampling_via_read_limit,
            args.sampling_via_percentage,
            args.num_threads,
        )
    if args.run_trimming:
        run_trimming(
            args.output_dir,
            args.short_read_fastq_dir,
            args.delete_pre_trim_fastq,
            args.num_threads,
            args.trim_galore_bin,
        )


if __name__ == "__main__":
    main()


# pylint:disable=pointless-string-statement
"""
def model_builder(work_dir):

    star_output_dir = os.path.join(work_dir, "star_output")

    all_junctions_file = os.path.join(star_output_dir, "all_junctions.sj")
    sjf_out = open(all_junctions_file, "w+")

    for sj_tab_file in glob.glob(input_dir + "/*.sj.tab"):
        sjf_in = open(sj_tab_file)
        sjf_lines = sjf_in.readlines()
        for line in sjf_lines:
            elements = line.split("\t")
            strand = "+"

            #    my $slice_name = $eles[0];
            #    my $start = $eles[1];
            #    my $end = $eles[2];
            #    my $strand = $eles[3];

            # If the strand is undefined then skip, Augustus expects a strand
            if elements[3] == "0":
                continue
            elif elements[3] == "2":
                strand = "-"

            junction_length = int(elements[2]) - int(elements[1]) + 1
            if junction_length < 100:
                continue

            if not elements[4] and elements[7] < 10:
                continue

            # For the moment treat multimapping and single
            # mapping things as a combined score
            score = float(elements[6]) + float(elements[7])
            score = str(score)
            output_line = [
                elements[0],
                "RNASEQ",
                "intron",
                elements[1],
                elements[2],
                score,
                strand,
                ".",
                ("src=W;mul=" + score + ";"),
            ]
            sjf_out.write("\t".join(output_line) + "\n")

    sjf_out.close()
"""
