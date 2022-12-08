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


import io
import logging
import os
import pathlib
import re
import subprocess
import shutil
import sys
import glob
import tempfile
from typing import Union

logger = logging.getLogger(__name__)


def add_log_file_handler(
    log_file_path: Union[pathlib.Path, str],
    logger_formatter: logging.Formatter,
    main_logger: logging.Logger,
):
    # logger: logger.Logger,
    """
    Create file handler and add to logger.
    """
    file_handler = main_logger.FileHandler(log_file_path, mode="a+")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(logger_formatter)
    return main_logger.addHandler(file_handler)


def add_log_console_handler(
    logger_formatter: logging.Formatter,
    main_logger: logging.Logger,
):
    """
    Create console handler and add to logger.
    """
    console_handler = main_logger.StreamHandler(sys.stderr)
    console_handler.setLevel(logging.DEBUG)
    console_handler.setFormatter(logger_formatter)
    return main_logger.addHandler(console_handler)


def create_dir(main_output_dir, dir_name):
    """
    Create directory or subdirectory and log operations.
    Args:
        main_output_dir: str main output directory path
        dir_name: str optional subdirectory to be created
    Returns:
        str Path to the created directory
    """
    if dir_name:
        target_dir = os.path.join(main_output_dir, dir_name)
    else:
        target_dir = main_output_dir

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


def check_exe(exe_path):
    """
    Check executable path
    Args:
        exe_path: str
            Path to the executable file

    Raises:
        OSError: If the executable file does not exist at the given path
    """
    if not shutil.which(exe_path):
        raise OSError("Exe does not exist. Path checked: %s" % exe_path)


def check_gtf_content(gtf_file, content_obj):
    """
    Check number of transcript lines in the GTF

    Arg:
      gtf_file: str path for the GTF file
      content_obj: str object to check in the gtf i.e gene_id, repeat

    Return: number of transcript lines
    """
    transcript_count = 0
    with open(gtf_file) as gtf_in:
        for line in gtf_in:
            eles = line.split("\t")
            if not len(eles) == 9:
                continue
            if eles[2] == content_obj:
                transcript_count += 1
    logger.info("%s GTF transcript count: %d", gtf_file, int(transcript_count))
    return transcript_count


def get_seq_region_lengths(genome_file, min_seq_length):
    """
    Split the genomic sequence in slices defined by  min_seq_length
    Args:
        genome_file: str path for the genome file
        min_seq_length: int slice length
    Return: Dict Dictionary with the sequence headers as keys and the sequence lengths as values
    """
    current_header = ""
    current_seq = ""

    seq_regions = {}
    with open(genome_file) as file_in:
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

        if len(current_seq) > min_seq_length:
            seq_regions[current_header] = len(current_seq)

    return seq_regions


def create_slice_ids(seq_region_lengths, slice_size, overlap, min_length):
    """
    Get list of ids for a genomic slice
    Arg:
    seq_region_lengths: dict
    Dictionary with the sequence headers as keys and the sequence lengths as values
    slice_size: int size of the slice
    overlap: int size of the overlap between two slices
    min_length: int min length of the slice
    Return: list List of IDs for the genomic slices
    """
    if not slice_size:
        slice_size = 1000000

    if not overlap:
        overlap = 0

    if not min_length:
        min_length = 0

    slice_ids = []

    for region in seq_region_lengths:
        region_length = int(seq_region_lengths[region])
        if region_length < min_length:
            continue

        if region_length <= slice_size:
            slice_ids.append([region, 1, region_length])
            continue

        start = 1
        end = start + slice_size - 1
        while end < region_length:
            start = start - overlap
            if start < 1:
                start = 1

            end = start + slice_size - 1
            if end > region_length:
                end = region_length
            if (end - start + 1) >= min_length:
                slice_ids.append([region, start, end])
            start = end + 1

    return slice_ids


def slice_output_to_gtf(  # pylint: disable=too-many-locals, too-many-branches, too-many-statements
    output_dir, extension, unique_ids, feature_id_label, new_id_prefix
):
    """
    Note that this does not make unique ids at the moment
    In many cases this is fine because the ids are unique by seq region,
    but in cases like batching it can cause problems
    So will add in a helper method to make ids unique

    This holds keys of the current slice details with the gene id to form unique keys.
    Each time a new key is added the overall gene counter is incremented
    and the value of the key is set to the new gene id. Any subsequent
    lines with the same region/gene id key will then just get
    the new id without incrementing the counter
    """
    gene_id_index = {}
    gene_transcript_id_index = {}
    gene_counter = 1

    # Similar to the gene id index, this will have a key that is based on
    # the slice details, gene id and transcript id. If there
    # is no existing entry, the transcript key will be added and the transcript
    # counter is incremented. If there is a key then
    # the transcript id will be replaced with the new transcript id
    # (which is based on the new gene id and transcript counter)
    # Example key KS8000.rs1.re1000000.gene_1.transcript_1

    transcript_id_count_index = {}

    feature_counter = 1

    feature_types = ["exon", "transcript", "repeat", "simple_feature"]
    if not extension:
        extension = ".gtf"
    gtf_files = glob.glob(output_dir + "/*" + extension)
    gtf_file_path = os.path.join(output_dir, "annotation.gtf")
    gtf_out = open(gtf_file_path, "w+")
    for gtf_file_path in gtf_files:  # pylint: disable=too-many-nested-blocks
        if os.stat(gtf_file_path).st_size == 0:
            logger.info("File is empty, will skip %s", gtf_file_path)
            continue

        gtf_file_name = os.path.basename(gtf_file_path)
        match = re.search(r"\.rs(\d+)\.re(\d+)\.", gtf_file_name)
        start_offset = int(match.group(1))
        gtf_in = open(gtf_file_path, "r")
        line = gtf_in.readline()
        while line:
            values = line.split("\t")
            if len(values) == 9 and (values[2] in feature_types):
                values[3] = str(int(values[3]) + (start_offset - 1))
                values[4] = str(int(values[4]) + (start_offset - 1))
                if unique_ids:
                    # Maybe make a unique id based on the feature type
                    # Basically region/feature id should be unique at this point,
                    # so could use region_id and current_id is key,
                    # value is the unique id that is incremented
                    attribs = values[8]

                    # This bit assigns unique gene/transcript ids if the line
                    # contains gene_id/transcript_id
                    match_gene_type = re.search(
                        r'(gene_id +"([^"]+)").+(transcript_id +"([^"]+)")',
                        line,
                    )
                    if match_gene_type:
                        full_gene_id_string = match_gene_type.group(1)
                        current_gene_id = match_gene_type.group(2)
                        full_transcript_id_string = match_gene_type.group(3)
                        current_transcript_id = match_gene_type.group(4)
                        gene_id_key = gtf_file_name + "." + str(current_gene_id)
                        transcript_id_key = gene_id_key + "." + str(current_transcript_id)
                        if gene_id_key not in gene_id_index:
                            new_gene_id = "gene" + str(gene_counter)
                            gene_id_index[gene_id_key] = new_gene_id
                            attribs = re.sub(
                                full_gene_id_string,
                                'gene_id "' + new_gene_id + '"',
                                attribs,
                            )
                            transcript_id_count_index[gene_id_key] = 1
                            gene_counter += 1
                        else:
                            new_gene_id = gene_id_index[gene_id_key]
                            attribs = re.sub(
                                full_gene_id_string,
                                'gene_id "' + new_gene_id + '"',
                                attribs,
                            )
                        if transcript_id_key not in gene_transcript_id_index:
                            new_transcript_id = (
                                gene_id_index[gene_id_key]
                                + ".t"
                                + str(transcript_id_count_index[gene_id_key])
                            )
                            gene_transcript_id_index[
                                transcript_id_key
                            ] = new_transcript_id
                            attribs = re.sub(
                                full_transcript_id_string,
                                'transcript_id "' + new_transcript_id + '"',
                                attribs,
                            )
                            transcript_id_count_index[gene_id_key] += 1
                        else:
                            new_transcript_id = gene_transcript_id_index[
                                transcript_id_key
                            ]
                            attribs = re.sub(
                                full_transcript_id_string,
                                'transcript_id "' + new_transcript_id + '"',
                                attribs,
                            )
                        values[8] = attribs

                    # If you don't match a gene line, try a feature line
                    else:
                        match_feature_type = re.search(
                            r"(" + feature_id_label + ' +"([^"]+)")',
                            line,
                        )
                        if match_feature_type:
                            full_feature_id_string = match_feature_type.group(1)
                            # current_feature_id = (
                            #    match_feature_type.group(2)
                            # )
                            new_feature_id = new_id_prefix + str(feature_counter)
                            attribs = re.sub(
                                full_feature_id_string,
                                feature_id_label + ' "' + new_feature_id + '"',
                                attribs,
                            )
                            feature_counter += 1
                            values[8] = attribs

                gtf_out.write("\t".join(values))
                line = gtf_in.readline()
            else:
                logger.info(
                    "Feature type not recognised, will skip. Feature type: %s", values[2]
                )
                line = gtf_in.readline()
        gtf_in.close()
    gtf_out.close()


def get_sequence(  # pylint: disable=too-many-arguments
    seq_region, start, end, strand, fasta_file, output_dir
):
    """
    Creates a tempfile and writes the bed info to it based on whatever information
    has been passed in about the sequence. Then runs bedtools getfasta. The fasta file
    should have a faidx. This can be created with the create_faidx static method prior
    to fetching sequence

    Arg:
    seq_region: str region name
    start: int region start
    end: int region end
    strand: int strand of the sequence
    fasta_file: str genome FASTA file
    output_dir: str working dir
    Return: str sequence
    """
    start = int(start)
    end = int(end)
    strand = int(strand)
    start -= 1
    bedtools_path = "bedtools"
    logger.info(
        "get_sequence %s",
        f"{seq_region}\t{start}\t{end}\t{strand}\t{fasta_file}\t{output_dir}",
    )
    with tempfile.NamedTemporaryFile(
        mode="w+t", delete=False, dir=output_dir
    ) as bed_temp_file:
        bed_temp_file.writelines(f"{seq_region}\t{start}\t{end}")
        bed_temp_file.close()
    bedtools_command = [
        bedtools_path,
        "getfasta",
        "-fi",
        fasta_file,
        "-bed",
        bed_temp_file.name,
    ]
    bedtools_output = subprocess.Popen(bedtools_command, stdout=subprocess.PIPE)
    for idx, line in enumerate(
        io.TextIOWrapper(bedtools_output.stdout, encoding="utf-8")
    ):
        if idx == 1:
            if strand == 1:
                sequence = line.rstrip()
            else:
                sequence = reverse_complement(line.rstrip())

    os.remove(bed_temp_file.name)
    return sequence


def reverse_complement(sequence):
    """
    Get the reverse complement of a nucleotide sequence.
    Args:
        sequence: str
            The nucleotide sequence
    Returns:
        str
            The reverse complement of the sequence
    """
    rev_matrix = str.maketrans("atgcATGC", "tacgTACG")
    return sequence.translate(rev_matrix)[::-1]
