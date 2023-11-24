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

"""
The STAR (Spliced Transcripts Alignment to a Reference) alignment tool is widely used
in genomics research for aligning RNA-seq data to a reference genome.
Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner.
Bioinformatics. 2013;29(1):15-21. doi:10.1093/bioinformatics/bts635
"""

__all__ = ["run_star"]

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
from typing import List

from ensembl.tools.anno.utils._utils import (
    check_exe,
    create_dir,
    check_gtf_content,
    get_seq_region_length,
)


def run_star(
    genome_file: Path,
    output_dir: Path,
    short_read_fastq_dir: Path,
    delete_pre_trim_fastq: bool = False,
    trim_fastq: bool = False,
    max_reads_per_sample: int = 0,
    max_intron_length: int = 100000,
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
    fastq_file_list = []
    file_types = ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz")
    fastq_file_list = [
        path for file_type in file_types for path in Path(short_read_fastq_dir).rglob(file_type)
    ]
    if len(fastq_file_list) == 0:
        raise IndexError(f"The list of fastq files is empty. Fastq dir:\n{short_read_fastq_dir}")

    # for file_type in file_types:
    #    fastq_file_list.extend(glob.glob(os.path.join(short_read_fastq_dir, file_type)))

    # Get list of paired paths
    fastq_file_list = _create_paired_paths(fastq_file_list)
    # Subsamples in parallel if there's a value set
    if max_reads_per_sample:
        subsample_transcriptomic_data(fastq_file_list)
        # Get the list of the new subsampled files
        fastq_file_list = [
            path for file_type in file_types for path in Path(short_read_fastq_dir).rglob(file_type)
        ]
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
        except Exception as e:#pylint:disable=broad-exception-caught
            logging.error("An error occurred while creating star index: %s", e)

    logging.info("Running Star on the files in the fastq dir")
    for fastq_file in fastq_file_list:
        # logger.info(fastq_file_path)
        # fastq_file_name = os.path.basename(fastq_file_path)
        star_tmp_dir = star_dir / "tmp"
        if star_tmp_dir.exists():
            shutil.rmtree(star_tmp_dir)
        sam_file = Path(f"{star_dir}/{fastq_file.name}.sam")
        junctions_file = Path(f"{star_dir}/{fastq_file.name}.sj.tab")
        sam_file_name = sam_file.name
        sam_temp_file = Path(f"{star_dir}/{sam_file_name}.tmp")
        bam_file = re.sub(".sam", ".bam", sam_file_name)
        bam_sort_file = Path(f"{star_dir}/{bam_file}")
        log_out_file = Path(f"{star_dir}/{fastq_file.name}.Log.final.out")
        if log_out_file.exists() and bam_sort_file.exists() and bam_sort_file.stat().st_size != 0:
            logging.info(
                "Found an existing bam file for the fastq file, \
                presuming the file has been processed, will skip"
            )
            continue

        logging.info("Processing %s", fastq_file)
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
            str(fastq_file),
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
        if fastq_file.suffix.endswith(".gz"):
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


def _create_paired_paths(fastq_file_paths: List) -> List[Path]:
    """
    Create list of paired transcriptomic fastq files

    Args:
        fastq_file_paths (List): List of transcriptomic file paths.

    Returns:
        List: List of paired transcriptomic files
    """
    path_dict = {}
    # final_list = []
    for fastq_file in fastq_file_paths:
        paired_name = re.search(r"(.+)_\d+\.(fastq|fq)", str(fastq_file))
        if not paired_name:
            logging.exception(
                "Could not find _1 or _2 at the end of the prefix \
                for file. Assuming file is not paired: %s",
                fastq_file,
            )
            # final_list.append([fastq_file])
            path_dict[fastq_file] = [fastq_file]
            continue
        run_accession = paired_name.group(1)
        if run_accession in path_dict:
            path_dict[run_accession].append(fastq_file)
        else:
            path_dict[run_accession] = [fastq_file]
    # for pair in path_dict:
    #    final_list.append(path_dict[pair])
    logging.info([value for values_list in path_dict.values() for value in values_list])
    return [value for values_list in path_dict.values() for value in values_list]


# pylint:disable=pointless-string-statement
"""
For an advanced and optimised subsampling we could use 
https://github.com/lh3/seqtk 
"""


def _subsample_paired_fastq_files(
    fastq_files: List[Path],
    subsample_read_limit: int = 100000000,
    num_threads: int = 2,
    compressed: bool = False,
) -> None:
    """
    Perform subsampling on two paired FastQ files in parallel using multiple threads.

    Args:
        fastq_files : Path for paired fastq files.
        output_files : Path for the output file.
        subsample_read_limit : Subsample size, defaults to 100000000.
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
        num_lines = sum(1 for line in open(fastq_file_1))  # pylint:disable=consider-using-with

    range_limit = int(num_lines / 4)
    if range_limit <= subsample_read_limit:
        logging.info(
            "Number of reads (%s is less than the max allowed read count (%s), \
            no need to subsample",
            str(range_limit),
            str(subsample_read_limit),
        )
        return

    rand_list = random.sample(range(0, range_limit - 1), subsample_read_limit)
    random_indices = {idx * 4: 1 for idx in rand_list}
    logging.info("Processing paired files in parallel")
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


def subsample_transcriptomic_data(fastq_file_list: List, num_threads: int = 2) -> None:
    """
    Subsample paired fastq files.

    Args:
        fastq_file_list : List of fastq file path to process.
        num_threads : number of threads
    """
    for fastq_files in fastq_file_list:
        fastq_file_1, fastq_file_2 = fastq_files
        # fastq_file_pair = ""
        # if len(fastq_files) == 2:
        #    fastq_file_pair = fastq_files[1]

        if len(fastq_files) == 1:
            fastq_file_1 = fastq_files[0]
            if Path(f"{fastq_file_1}.sub").exists():
                logging.info(
                    "Found an existing .sub file on the fastq path, will use that instead. File:%s.sub",
                    fastq_file_1,
                )
            else:
                _subsample_paired_fastq_files(fastq_files, compressed=True, num_threads=num_threads)

        elif len(fastq_files) == 2:
            fastq_file_1, fastq_file_2 = fastq_files
            if Path(f"{fastq_file_1}.sub").exists() and Path(f"{fastq_file_2}.sub").exists():
                logging.info(
                    "Found an existing .sub files on the fastq path for both members of the pair, will use \
                    those instead of subsampling again. Files: %s.sub,%s.sub",
                    fastq_file_1,
                    fastq_file_2,
                )
            elif Path(f"{fastq_file_2}.sub").exists():
                _subsample_paired_fastq_files(fastq_files, compressed=True, num_threads=num_threads)


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

    fastq_file_list = []
    file_types = ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz")
    fastq_file_list = [
        path for file_type in file_types for path in Path(short_read_fastq_dir).rglob(file_type)
    ]
    fastq_file_list = _create_paired_paths(fastq_file_list)

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
    for fastq_file in fastq_file_list:
        if delete_pre_trim_fastq:
            fastq_file.unlink()
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

    run_star(
        args.genome_file,
        args.output_dir,
        args.short_read_fastq_dir,
        args.delete_pre_trim_fastq,
        args.trim_fastq,
        args.max_reads_per_sample,
        args.max_intron_length,
        args.num_threads,
        args.star_bin,
        args.samtools_bin,
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