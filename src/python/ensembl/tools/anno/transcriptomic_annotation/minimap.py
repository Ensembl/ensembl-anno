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
"""
Minimap2 is a pairwise sequence alignment algorithm designed for efficiently comparing nucleotide sequences.
The algorithm uses a versatile indexing strategy to quickly find approximate matches between sequences, 
allowing it to efficiently align long sequences against reference genomes or other sequences.

Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094-3100.
"""

__all__ = ["run_minimap2"]

import argparse
import logging
import logging.config
from pathlib import Path
import subprocess
from typing import List

from ensembl.tools.anno.utils._utils import (
    check_exe,
    create_dir,
    check_gtf_content,
)


def run_minimap2(
    output_dir: Path,
    long_read_fastq_dir: Path,
    genome_file: Path,
    minimap2_bin: Path = Path("minimap2"),
    paftools_bin: Path = Path("paftools.js"),
    max_intron_length: int = 100000,
    num_threads: int = 1,
) -> None:
    """
    Run Minimap2 to align long read data against genome file.
    Default Minimap set for PacBio data.

        :param output_dir: Working directory path.
        :type output_dir: Path
        :param long_read_fastq_dir: Long read directory path.
        :type long_read_fastq_dir: Path
        :param genome_file: Genome file path.
        :type genome_file: Path
        :param minimap2_bin: Software path.
        :type minimap2_bin: Path, default minimap2
        :param paftools_bin: Js path.
        :type paftools_bin: Path, default paftools.js
        :param max_intron_length: The maximum intron size for alignments. Defaults to 100000.
        :type max_intron_length: int, default 100000
        :param num_threads: Number of available threads.
        :type num_threads: int, default 1

        :return: None
        :rtype: None
    """
    check_exe(minimap2_bin)
    check_exe(paftools_bin)
    minimap2_dir = create_dir(output_dir, "minimap2_output")

    logging.info("Skip analysis if the gtf file already exists")
    output_file = minimap2_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logging.info("Minimap2 gtf file exists, skipping analysis")
            return
    minimap2_index_file = minimap2_dir / f"{Path(genome_file).name}.mmi"
    # minimap2_hints_file = minimap2_dir /"minimap2_hints.gff"
    file_types = ("*.fastq", "*.fq")
    fastq_file_list = [
        path for file_type in file_types for path in Path(long_read_fastq_dir).rglob(file_type)
    ]
    if len(fastq_file_list) == 0:
        raise IndexError(f"The list of fastq files is empty. Fastq dir:\n{long_read_fastq_dir}")

    if not minimap2_index_file.exists():
        logging.info("Did not find an index file for minimap2. Will create now")
        try:
            subprocess.run(  # pylint:disable=subprocess-run-check
                [
                    minimap2_bin,
                    "-t",
                    str(num_threads),
                    "-d",
                    str(minimap2_index_file),
                    genome_file,
                ]
            )
        except subprocess.CalledProcessError as e:
            logging.error("An error occurred while creating minimap2 index: %s", e)
        except OSError as e:
            logging.error("An OS error occurred: %s", e)

    logging.info("Running minimap2 on the files in the long read fastq dir")
    for fastq_file in fastq_file_list:
        sam_file = minimap2_dir / f"{fastq_file.name}.sam"
        bed_file = minimap2_dir / f"{fastq_file.name}.bed"
        logging.info("Processing %s", fastq_file)
        with open(bed_file, "w+", encoding="utf8") as bed_file_out:
            subprocess.run(  # pylint:disable=subprocess-run-check
                [
                    minimap2_bin,
                    "-G",
                    str(max_intron_length),
                    "-t",
                    str(num_threads),
                    "--cs",
                    "--secondary=no",
                    "-ax",
                    "splice",
                    "-u",
                    "b",
                    minimap2_index_file,
                    fastq_file,
                    "-o",
                    sam_file,
                ]
            )
            logging.info("Creating bed file from SAM")
            subprocess.run(
                [paftools_bin, "splice2bed", sam_file], stdout=bed_file_out
            )  # pylint:disable=subprocess-run-check

    _bed_to_gtf(minimap2_dir)

    logging.info("Completed running minimap2")


def _bed_to_gtf(output_dir: Path) -> None:
    """
    Convert bed file into gtf file
    Args:
        output_dir : Working directory path.
    """
    gtf_file_path = output_dir / "annotation.gtf"
    with open(gtf_file_path, "w+", encoding="utf8") as gtf_out:
        gene_id = 1
        for bed_file in output_dir.glob("*.bed"):
            logging.info("Converting bed to GTF: %s", str(bed_file))
            with open(bed_file, "r", encoding="utf8") as bed_in:
                for line in bed_in:
                    elements = line.rstrip().split("\t")
                    seq_region_name = elements[0]
                    offset = int(elements[1])
                    strand = elements[5]
                    # sizes of individual block of exons
                    block_sizes = [size for size in elements[10].split(",") if size]
                    block_starts = [size for size in elements[11].split(",") if size]
                    exons = _bed_block_to_exons(block_sizes, block_starts, offset)
                    transcript_start = None
                    transcript_end = None
                    exon_records = []
                    for i, exon_coords in enumerate(exons):
                        if transcript_start is None or exon_coords[0] < transcript_start:
                            transcript_start = exon_coords[0]

                        if transcript_end is None or exon_coords[1] > transcript_end:
                            transcript_end = exon_coords[1]

                        exon_line = (
                            f"{seq_region_name}\tminimap\texon\t{exon_coords[0]}\t"
                            f"{exon_coords[1]}\t.\t{strand}\t.\t"
                            f'gene_id "minimap_{gene_id}"; transcript_id "minimap_{gene_id}"; '
                            f'exon_number "{i+ 1}";\n'
                        )
                        exon_records.append(exon_line)
                    transcript_line = (
                        f"{seq_region_name}\tminimap\ttranscript\t{transcript_start}\t"
                        f"{transcript_end}\t.\t{strand}\t.\t"
                        f'gene_id "minimap_{gene_id}"; transcript_id "minimap_{gene_id}"\n'
                    )
                    gtf_out.write(transcript_line)
                    for exon_line in exon_records:
                        gtf_out.write(exon_line)
                    gene_id += 1


def _bed_block_to_exons(block_sizes: List, block_starts: List, offset: int) -> List:
    """
    Extract exon size and start from exon feature block
    Args:
        block_sizes : Block feature size.
        block_starts : Block feature starts.
        offset : Feature offset.

    Returns:
        List of exon coordinates
    """
    exons = []
    for i, _ in enumerate(block_sizes):
        block_start = offset + int(block_starts[i]) + 1
        block_end = block_start + int(block_sizes[i]) - 1
        if block_end < block_start:
            logging.warning("Warning: block end is less than block start, skipping exon")
            continue
        exon_coords = [str(block_start), str(block_end)]
        exons.append(exon_coords)
    return exons


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Minimap2's arguments")
    parser.add_argument("--output_dir", required=True, help="Output directory path")
    parser.add_argument("--long_read_fastq_dir", required=True, help="Long read directory path")
    parser.add_argument("--genome_file", required=True, help="Genome file path")
    parser.add_argument("--minimap2_bin", default="minimap2", help="Minimap2 software path")
    parser.add_argument("--paftools_bin", default="paftools.js", help="Paftools js path")
    parser.add_argument("--max_intron_length", type=int, default=100000, help="The maximum intron length.")
    parser.add_argument("--num_threads", type=int, default=1, help="Number of threads")
    return parser.parse_args()


def main():
    """Minimap2's entry-point."""
    args = parse_args()

    log_file_path = create_dir(args.output_dir, "log") / "minimap.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"

    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )

    run_minimap2(
        args.output_dir,
        args.long_read_fastq_dir,
        args.genome_file,
        args.minimap2_bin,
        args.paftools_bin,
        args.max_intron_length,
        args.num_threads,
    )


if __name__ == "__main__":
    main()
