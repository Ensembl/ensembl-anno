# pylint: disable=missing-module-docstring
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

import errno
import io
import logging
import os
from os import PathLike
from pathlib import Path
import random
import re
import subprocess
import shutil
import tempfile
from typing import Dict, List

logger = logging.getLogger("__name__." + __name__)

# with open(Path(__file__).parents[6] / "conf/config.json", "r", encoding="utf8") as f:
#    config = json.load(f)


def create_dir(input_dir: Path, dir_name: str | None) -> Path:
    """
    Create directory or subdirectory and log operations.
    Args:
        main_output_dir: str main output directory path
        dir_name: str optional subdirectory to be created
    Returns:
        str Path to the created directory
    """
    if dir_name:
        target_dir = Path(input_dir) / str(dir_name)
    else:
        target_dir = input_dir

    if os.path.exists(target_dir):
        logger.warning("Directory already exists, will not create again")
        return target_dir

    logger.info("Attempting to create target dir: %s", target_dir)

    try:
        os.mkdir(target_dir)
    except OSError:
        logger.error("Creation of the dir failed, path used: %s", target_dir)
    else:
        logger.info("Successfully created the dir on the following path: %s", target_dir)
    return target_dir


def check_exe(exe_bin: PathLike) -> None:
    """
    Check executable path
    Args:
        exe_bin: Executable path

    Raises:
        OSError: Executable path does not exist
    """
    if not shutil.which(exe_bin):
        raise OSError(f"Exe does not exist. Path checked: {exe_bin}")


def check_gtf_content(gtf_file: PathLike, content_obj: str) -> int:
    """
    Check number of transcript lines in the GTF

    Arg:
        gtf_file: GTF file path
        content_obj: Object to detect and count in the gtf i.e transcript, repeat

    Return: Number of occurences
    """
    obj_count = 0
    with open(gtf_file, "r", encoding="utf8") as gtf_in:
        for line in gtf_in:
            gtf_raw = line.split("\t")
            if not len(gtf_raw) == 9:
                continue
            if gtf_raw[2] == content_obj:
                obj_count += 1
    logger.info("%s; Number of %s detected: %d", gtf_file, content_obj, int(obj_count))
    return obj_count


def get_seq_region_length(genome_file: PathLike, min_seq_length: int = 0) -> Dict:
    """
    Split the genome file according to the header and store in a dictionary
    all the sequences whose length is greater than min_seq_length.
    Args:
        genome_file: Genome file path.
        min_seq_length: Minimum slice length.
    Return: Dictionary of sequence headers and the corresponding sequence length
    """
    current_header = ""
    current_seq = ""

    seq_region_to_length = {}
    with open(genome_file, "r", encoding="utf8") as file_in:
        for line in file_in:
            match = re.search(r">(.+)$", line)
            if match and current_header:
                if len(current_seq) > min_seq_length:
                    seq_region_to_length[current_header] = len(current_seq)

                current_seq = ""
                current_header = match.group(1)
            elif match:
                current_header = match.group(1)
            else:
                current_seq += line.rstrip()

        if len(current_seq) > min_seq_length:
            seq_region_to_length[current_header] = len(current_seq)
    return seq_region_to_length


def get_slice_id(
    seq_region_to_length: Dict,
    slice_size: int = 1000000,
    overlap: int = 0,
    min_length: int = 0,
) -> List:
    """
    Get list of ids for a genomic slice
    Arg:
    seq_region_to_length: Dictionary with the sequence headers
    as keys and the sequence lengths as values
    slice_size: Size of the slice
    overlap: Overlap length between two slices
    min_length: Min length of the slice
    Return: List of IDs for each genomic slice
    """

    slice_ids_per_region = []
    for seq_region in seq_region_to_length:
        seq_region_length = int(seq_region_to_length[seq_region])
        if seq_region_length < min_length:
            continue

        if seq_region_length <= slice_size:
            slice_ids_per_region.append([seq_region, 1, seq_region_length])
            continue

        start = 1
        end = start + slice_size - 1
        while end < seq_region_length:
            start = start - overlap
            start = max(start, 1)
            end = start + slice_size - 1
            end = min(end, seq_region_length)
            if (end - start + 1) >= min_length:
                slice_ids_per_region.append([seq_region, start, end])
            start = end + 1
    return slice_ids_per_region


def update_string(text, lookup_text, new_text) -> str:
    """
    Update substring
    Args:
        text : Text.
        lookup_text : String to look for.
        new_text : String to substitute.
    Return
        Text with string substitution
    """
    text = re.sub(
        lookup_text,
        new_text,
        text,
    )
    return text


#   slice_output_to_gtf(repeatmasker_dir, "repeat_id",
# "repeatmask", True, ".rm.gtf")
def slice_output_to_gtf(  # pylint: disable=too-many-branches, too-many-statements, too-many-locals
    output_dir: Path,
    feature_id_label: str = "",
    new_id_prefix: str = "",
    unique_ids: bool = True,
    file_extension: str = ".gtf",
) -> None:
    """
    Collect all the gtf files per file extension and
    merge them in a final gtf file  assigning unique ids.

    This holds keys of the current slice details with
    the gene id to form unique keys.
    Each time a new key is added the overall gene
    counter is incremented
    and the value of the key is set to the new gene id.
    Any subsequent
    lines with the same region/gene id key will then just get
    the new id without incrementing the counter.

    Args:
    output_dir : Output directory.
    feature_id_label : Feature identifier.
    new_id_prefix : New feature identifier.
    unique_ids : If True assign unique ids for the same
    feature type.
    file_extension : Input file extension.
    """
    feature_types = ["exon", "transcript", "repeat", "simple_feature"]
    # Initialise gene and feature counter
    gene_counter = 1
    feature_counter = 1
    # Initialise dictionaries that will store the list of
    # gene and transcript indexes
    gene_id_collection = {}
    gene_transcript_id_collection = {}
    transcript_id_count_gene_id = {}  # one gene can have multiple transcripts
    input_files = output_dir.glob(str("*" + file_extension))
    with open(output_dir / "annotation.gtf", "w+", encoding="utf8") as output_file:
        for input_file in input_files:
            if os.stat(input_file).st_size == 0:
                logger.info("File is empty, will skip %s", input_file.name)
                continue
            match = re.search(r"\.rs(\d+)\.re(\d+)\.", input_file.name)
            assert match is not None
            start_offset = int(match.group(1))
            with open(input_file, "r", encoding="utf8") as gtf_in:
                for line in gtf_in:
                    values = line.split("\t")
                    if len(values) == 9 and (values[2] in feature_types) and unique_ids:
                        # each slice start from 0 so we need to add the
                        # offset to get the real coordinates
                        values[3] = str(int(values[3]) + (start_offset - 1))
                        values[4] = str(int(values[4]) + (start_offset - 1))
                        # Unique id based on gene and transcript ids
                        attribs = values[8]
                        # Get the slice details, gene id and transcript id.
                        match_gene_type = re.search(
                            r'(gene_id +"([^"]+)").+(transcript_id +"([^"]+)")',
                            line,
                        )
                        # gene_id "1"; transcript_id "1"; biotype "tRNA_pseudogene";
                        if match_gene_type:
                            full_gene_id_string = match_gene_type.group(1)
                            current_gene_id = match_gene_type.group(2)
                            full_transcript_id_string = match_gene_type.group(3)
                            current_transcript_id = match_gene_type.group(4)
                            # Example key KS8000.rs1.re1000000.1
                            gene_id_slice = input_file.name + "." + str(current_gene_id)
                            # Example key KS8000.rs1.re1000000.1.transcript.1
                            transcript_id_slice = gene_id_slice + "." + str(current_transcript_id)
                            # If there is no existing entry, the gene key will be added
                            # and the gene counter is incremented.
                            # gene_id "gene1"; transcript_id "gene1.t1"
                            if gene_id_slice not in gene_id_collection:
                                new_gene_id = f"gene{gene_counter}"
                                gene_id_collection[gene_id_slice] = new_gene_id
                                transcript_id_count_gene_id[gene_id_slice] = 1
                                gene_counter += 1
                            else:
                                # If there is a key then the gene
                                # id will be replaced with the new gene id
                                new_gene_id = gene_id_collection[gene_id_slice]
                            attribs = re.sub(
                                full_gene_id_string,
                                'gene_id "' + new_gene_id + '"',
                                attribs,
                            )
                            # If there is no existing entry, the transcript key will be added
                            # and the transcript counter is incremented.
                            if transcript_id_slice not in gene_transcript_id_collection:
                                new_transcript_id = (
                                    gene_id_collection[gene_id_slice]
                                    + ".t"
                                    + str(
                                        transcript_id_count_gene_id[gene_id_slice]
                                    )  # pylint:disable=line-too-long
                                )
                                gene_transcript_id_collection[transcript_id_slice] = (
                                    new_transcript_id  # pylint:disable=line-too-long
                                )
                                transcript_id_count_gene_id[gene_id_slice] += 1
                            else:
                                # If a transcript of the same set is already present,
                                # the new id with the incremented counter is added
                                new_transcript_id = gene_transcript_id_collection[
                                    transcript_id_slice
                                ]  # pylint:disable=line-too-long
                            attribs = re.sub(
                                full_transcript_id_string,
                                'transcript_id "' + new_transcript_id + '"',
                                attribs,
                            )
                            values[8] = attribs
                            logger.info("FINAL GTF %s", values)
                            output_file.write("\t".join(values))
                            # Unique id based on the feature type
                        else:
                            if new_id_prefix == "repeatmask":
                                match_feature_type = re.search(
                                    r"("
                                    + feature_id_label
                                    + "\d+)",  # pylint:disable=line-too-long, anomalous-backslash-in-string
                                    line,
                                )
                            else:
                                match_feature_type = re.search(
                                    r"("
                                    + feature_id_label
                                    + ' +"([^"]+)")',  # "\d+)",# pylint:disable=line-too-long
                                    line,
                                )
                            if match_feature_type:
                                full_feature_id = match_feature_type.group(1)
                                new_feature_id = new_id_prefix + str(feature_counter)
                                attribs = re.sub(
                                    full_feature_id,
                                    feature_id_label + ' "' + new_feature_id + '"',
                                    attribs,
                                )
                                feature_counter += 1
                                values[8] = attribs
                                output_file.write("\t".join(values))
                    else:
                        logger.info(
                            "Feature type not recognised, will skip. Feature type: %s",
                            values[2],
                        )


def get_sequence(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    seq_region: str,
    start: int,
    end: int,
    strand: int,
    fasta_file: Path,
    output_dir: Path,
) -> str:
    """
    Creates a tempfile and writes the bed info to it based
    on whatever information
    has been passed in about the sequence. Then runs
    bedtools getfasta. The fasta file
    should have a faidx. This can be created with the
    create_faidx static method prior
    to fetching sequence

    Arg:
    seq_region: str region name
    start: int region start
    end: int region end
    strand: int strand of the sequence
    fasta_file: genome FASTA file with indexing
    output_dir: str working dir
    Return: str sequence
    """
    start -= 1
    logger.info(
        "get_sequence %s",
        f"{seq_region}\t{start}\t{end}\t{strand}\t{fasta_file}\t{output_dir}",  # pylint:disable=line-too-long
    )
    with tempfile.NamedTemporaryFile(
        mode="w+t", delete=False, dir=str(output_dir)
    ) as bed_temp_file:  # pylint:disable=line-too-long
        bed_temp_file.writelines(f"{seq_region}\t{start}\t{end}")
        bed_temp_file.close()
    bedtools_command = [
        "bedtools",
        "getfasta",
        "-fi",
        str(fasta_file),
        "-bed",
        str(bed_temp_file.name),
    ]
    bedtools_output = subprocess.Popen(  # pylint:disable=consider-using-with
        bedtools_command, stdout=subprocess.PIPE
    )
    if bedtools_output.stdout is None:
        raise RuntimeError("bedtools process did not provide stdout")

    for idx, line in enumerate(
        io.TextIOWrapper(bedtools_output.stdout, encoding="utf-8")
    ):  # pylint:disable=line-too-long, consider-using-with
        if idx == 1:
            if strand == 1:
                sequence = line.rstrip()
            else:
                sequence = reverse_complement(line.rstrip())

    os.remove(bed_temp_file.name)
    return sequence  # pylint:disable=possibly-used-before-assignment


def reverse_complement(sequence: str) -> str:
    """
    Get the reverse complement of a nucleotide sequence.
    Args:
        sequence: The nucleotide sequence
    Returns:
        The reverse complement of the sequence
    """
    rev_matrix = str.maketrans("atgcATGC", "tacgTACG")
    return sequence.translate(rev_matrix)[::-1]


def check_file(file_path: Path) -> None:
    """
    Raise an error when the file doesn't exist
    Args:
        file_path: File path
    """
    if not file_path.is_file():
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file_path)


def split_protein_file(protein_dataset: Path, output_dir: Path, batch_size: int = 20) -> List:
    """
    The protein dataset file is split by a number of sequence
    equals to the batch_size
    in batch files stored in 10 output directories.
    protein_dataset : Path for the protein dataset.
    output_dir : Output directory path.
    batch_size : Size of the batch, it needs to be equals to the
    number of threads
    to parallelise the sequence processing for each file.
    """
    batched_protein_files = []

    for i in range(0, 10):
        create_dir(output_dir, (f"bin_{i}"))
    with open(protein_dataset, "r", encoding="utf8") as file_in:
        seq_count = 0
        batch_count = 0
        current_record = ""
        initial_seq = True
        for line in file_in:
            match = re.search(r">(.+)$", line)
            # match header and is not first sequence, if the number
            # of stored sequences in each file equals
            # the number of batch_size, a new file will be created
            # and the current_record reset
            if match and not initial_seq and seq_count % batch_size == 0:
                bin_num = random.randint(0, 9)
                batch_file = output_dir / f"bin_{bin_num}" / f"{batch_count}.fa"
                with batch_file.open("w+") as file_out:
                    file_out.write(current_record)
                batch_count += 1
                seq_count += 1
                current_record = line
                batched_protein_files.append(batch_file)
            # match header and is the first sequence
            elif match:
                current_record += line
                initial_seq = False
                seq_count += 1
            # other lines
            else:
                current_record += line

        if current_record:
            bin_num = random.randint(0, 9)
            batch_file = output_dir / f"bin_{bin_num}" / f"{batch_count}.fa"
            with batch_file.open("w+") as file_out:
                file_out.write(current_record)
            batched_protein_files.append(batch_file)
    return batched_protein_files


def get_seq_region_lengths(
    genome_file: Path,
    min_seq_length: int,
) -> dict[str, int]:
    """
    Extract sequence lengths from a FASTA file.

    Only sequences longer than ``min_seq_length`` are returned.

    Args:
        genome_file: Path to the FASTA file.
        min_seq_length: Minimum sequence length to retain.

    Returns:
        Dictionary mapping sequence headers to sequence lengths.
    """
    seq_regions: dict[str, int] = {}

    current_header: str | None = None
    current_seq = ""

    with open(genome_file, "r", encoding="utf-8") as file_in:
        for line in file_in:
            match = re.search(r">(.+)$", line)
            if match and current_header:
                if len(current_seq) > min_seq_length:
                    seq_regions[current_header] = len(current_seq)

                current_seq = ""
                current_header = match.group(1)
            elif match:
                current_header = match.group(1)
            else:
                current_seq += line.rstrip()

        if current_header and len(current_seq) > min_seq_length:
            seq_regions[current_header] = len(current_seq)

    return seq_regions


def seq_region_names(genome_file: Path) -> list[str]:
    """
    Extract sequence region names from a FASTA file.

    The region named "MT" is skipped.

    Args:
        genome_file: Path to the FASTA file.

    Returns:
        List of sequence region names.
    """
    region_list: list[str] = []

    with open(genome_file, "r", encoding="utf-8") as file_in:
        line = file_in.readline()
        while line:
            match = re.search(r">([^\s]+)", line)
            if match:
                region_name = match.group(1)
                if region_name == "MT":  # pylint: disable=no-else-continue
                    logger.info("Skipping region named MT")
                    line = file_in.readline()
                    continue
                else:
                    region_list.append(match.group(1))
            line = file_in.readline()

    return region_list


def check_transcriptomic_output(main_output_dir: Path) -> None:
    """
    Validate transcriptomic annotation outputs.

    Args:
        main_output_dir: Main output directory.

    Raises:
        OSError: If no transcriptomic annotations are found or the total
            number of lines is below the expected threshold.
    """
    transcriptomic_dirs = (
        "scallop_output",
        "stringtie_output",
        "minimap2_output",
    )
    total_lines = 0
    min_lines = 100000
    for transcriptomic_dir in transcriptomic_dirs:
        annotation_file = main_output_dir / transcriptomic_dir / "annotation.gtf"
        if not annotation_file.exists():
            logger.warning(
                "No annotation.gtf found in %s. This may be expected "
                "(e.g. no long-read data were provided).",
                transcriptomic_dir,
            )
            continue
        with annotation_file.open(encoding="utf-8") as file_handle:
            num_lines = sum(1 for _ in file_handle)

        total_lines += num_lines

        logger.info(
            "Found %s lines in %s/annotation.gtf",
            num_lines,
            transcriptomic_dir,
        )
    if total_lines == 0:
        raise OSError(
            "Transcriptomic mode was enabled, but all transcriptomic annotation files are empty."
        )  # pylint:disable=line-too-long

    if total_lines <= min_lines:
        raise OSError(
            f"Transcriptomic mode was enabled, but the total number of "
            f"lines across annotation files is below the expected "
            f"minimum. Found: {total_lines}. Minimum allowed: {min_lines}."
        )
    logger.info(
        "Found %s total lines across transcriptomic annotation files. Checks passed.",
        total_lines,
    )
