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
Minimap2 is a pairwise sequence alignment algorithm designed for efficiently comparing nucleotide sequences.
The algorithm uses a versatile indexing strategy to quickly find approximate matches between sequences, 
allowing it to efficiently align long sequences against reference genomes or other sequences.

Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094-3100.
"""

__all__ = ["run_minimap2"]
import logging
import logging.config
from pathlib import Path
import subprocess
from typing import List
import argschema


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
    Args:
        output_dir : Working directory path.
        long_read_fastq_dir : Long read directory path.
        genome_file : Genome file path.
        minimap2_bin : Software path.
        paftools_bin : Software path.
        max_intron_length : The maximum intron size for alignments. Defaults to 100000.
        num_threads : Number of available threads.
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
    minimap2_index_file = minimap2_dir / f"{genome_file.name}.mmi"
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

    bed_to_gtf(minimap2_dir)

    logging.info("Completed running minimap2")


def bed_to_gtf(output_dir: Path) -> None:
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
                    exons = bed_block_to_exons(block_sizes, block_starts, offset)
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


def bed_block_to_exons(block_sizes: List, block_starts: List, offset: int) -> List:
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


class InputSchema(argschema.ArgSchema):
    """Input arguments expected to run Minimap2 software."""

    output_dir = argschema.fields.OutputDir(required=True, description="Output directory path")
    long_read_fastq_dir = argschema.fields.String(
        required=True,
        description="Long read directory path",
    )
    genome_file = argschema.fields.InputFile(required=True, description="Genome file path")
    minimap2_bin = argschema.fields.String(
        required=False,
        default="minimap2",
        description="Minimap2 software path",
    )
    paftools_bin = argschema.fields.String(
        required=False,
        default="paftools.js",
        description="Paftools software path",
    )
    max_intron_length = argschema.fields.Integer(
        required=False,
        default="100000",
        description="The maximum intron length.",
    )
    max_intron_length = argschema.fields.Integer(
        required=False,
        default="100000",
        description="The maximum intron size for alignments.",
    )
    num_threads = argschema.fields.Integer(required=False, default=1, description="Number of threads")


def main() -> None:
    """Minimap2's entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    log_file_path = create_dir(mod.args["output_dir"], "log") / "minimap.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"
    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )
    run_minimap2(
        mod.args["output_dir"],
        mod.args["long_read_fastq_dir"],
        mod.args["genome_file"],
        mod.args["minimap2_bin"],
        mod.args["paftools_bin"],
        mod.args["max_intron_length"],
        mod.args["num_threads"],
    )
