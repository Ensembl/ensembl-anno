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
import json
import logging
import os
from os import PathLike
from pathlib import Path
import random
import re
import subprocess
import shutil
import glob
import tempfile
from typing import Dict, List

logger = logging.getLogger("__name__." + __name__)

# with open(Path(__file__).parents[6] / "conf/config.json", "r", encoding="utf8") as f:
#    config = json.load(f)


def create_dir(input_dir: Path, dir_name: str) -> Path:
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
    with open(gtf_file) as gtf_in:
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
    with open(genome_file) as file_in:
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
    #logging.info("length of the seq region %s", seq_region_to_length)
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
    seq_region_to_length: Dictionary with the sequence headers as keys and the sequence lengths as values
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
    #logging.info("list of slice ids %s", str(slice_ids_per_region))
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


#   slice_output_to_gtf(repeatmasker_dir, "repeat_id", "repeatmask", True, ".rm.gtf")
def slice_output_to_gtf(  # pylint: disable=too-many-branches
    output_dir: Path,
    feature_id_label: str = "",
    new_id_prefix: str = "",
    unique_ids: bool = True,
    file_extension: str = ".gtf",
) -> None:
    """
    Collect all the gtf files per file extension and merge them in a final gtf file  assigning unique ids.

    This holds keys of the current slice details with the gene id to form unique keys.
    Each time a new key is added the overall gene counter is incremented
    and the value of the key is set to the new gene id. Any subsequent
    lines with the same region/gene id key will then just get
    the new id without incrementing the counter.

    Args:
    output_dir : Output directory.
    feature_id_label : Feature identifier.
    new_id_prefix : New feature identifier.
    unique_ids : If True assign unique ids for the same feature type.
    file_extension : Input file extension.
    """
    feature_types = ["exon", "transcript", "repeat", "simple_feature"]
    # Initialise gene and feature counter
    gene_counter = 1
    feature_counter = 1
    # Initialise dictionaries that will store the list of gene and transcript indexes
    gene_id_collection = {}
    gene_transcript_id_collection = {}
    transcript_id_count_gene_id = {}  # one gene can have multiple transcripts
    input_files = output_dir.glob(str("*" + file_extension))
    with open(output_dir / "annotation.gtf", "w+") as output_file:
        for input_file in input_files:
            if os.stat(input_file).st_size == 0:
                logger.info("File is empty, will skip %s", input_file.name)
                continue
            match = re.search(r"\.rs(\d+)\.re(\d+)\.", input_file.name)
            assert match is not None
            start_offset = int(match.group(1))
            with open(input_file, "r") as gtf_in:
                for line in gtf_in:
                    values = line.split("\t")
                    if len(values) == 9 and (values[2] in feature_types) and unique_ids:
                        # each slice start from 0 so we need to add the offset to get the real coordinates
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
                            transcript_id_slice = (
                                gene_id_slice + "." + str(current_transcript_id)
                            )
                            # If there is no existing entry, the gene key will be added
                            # and the gene counter is incremented.
                            # gene_id "gene1"; transcript_id "gene1.t1"
                            if gene_id_slice not in gene_id_collection:
                                new_gene_id = f"gene{gene_counter}"
                                gene_id_collection[gene_id_slice] = new_gene_id
                                transcript_id_count_gene_id[gene_id_slice] = 1
                                gene_counter += 1
                            else:
                                # If there is a key then the gene id will be replaced with the new gene id
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
                                    + str(transcript_id_count_gene_id[gene_id_slice])
                                )
                                gene_transcript_id_collection[
                                    transcript_id_slice
                                ] = new_transcript_id
                                transcript_id_count_gene_id[gene_id_slice] += 1
                            else:
                                # If a transcript of the same set is already present,
                                # the new id with the incremented counter is added
                                new_transcript_id = gene_transcript_id_collection[
                                    transcript_id_slice
                                ]
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
                                    r"(" + feature_id_label + "\d+)",
                                    line,
                                )
                            else:
                                match_feature_type = re.search(  # repeat_id\d+  r"(" + feature_id_label + ' +"([^"]+)")',
                                    r"(" + feature_id_label + ' +"([^"]+)")',  # "\d+)",
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
def get_sequence(
    seq_region: str,
    start: int,
    end: int,
    strand: int,
    fasta_file: Path,
    output_dir: Path,
) -> str:
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
    fasta_file: genome FASTA file with indexing
    output_dir: str working dir
    Return: str sequence
    """
    start -= 1
    logger.info(
        "get_sequence %s",
        f"{seq_region}\t{start}\t{end}\t{strand}\t{fasta_file}\t{output_dir}",
    )
    with tempfile.NamedTemporaryFile(
        mode="w+t", delete=False, dir=str(output_dir)
    ) as bed_temp_file:
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
    bedtools_output = subprocess.Popen(bedtools_command, stdout=subprocess.PIPE)
    for idx, line in enumerate(
        io.TextIOWrapper(bedtools_output.stdout, encoding="utf-8")  # type: ignore
    ):
        if idx == 1:
            if strand == 1:
                sequence = line.rstrip()
            else:
                sequence = reverse_complement(line.rstrip())

    os.remove(bed_temp_file.name)
    # logger.info(f"sequence : {sequence}")
    return sequence


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

def fasta_to_dict(fasta_lines: List[str]) -> Dict[str, str]:
    """
    Save fastq sequences in a dictionary <header, sequence>
    fasta_lines: List[str]: List of sequences
    """
    index = {}
    header = None
    seq = []

    for line in fasta_lines:
        line = line.strip()
        if line.startswith(">"):
            if header:
                index[header] = ''.join(seq)
            header_match = re.match(r">(.+)", line)
            if header_match:
                header = header_match.group(1)
                seq = []
            else:
                raise ValueError(f"Invalid FASTA header line: {line}")
        else:
            seq.append(line)

    if header:
        index[header] = ''.join(seq)

    return index
"""
def fasta_to_dict(fasta_lines:List[str]):
    index = {}
    for line in fasta_lines:
    #it = iter(fasta_lines)
    #for header in it:
        match = re.search(r">(.+)\n$", line)
        header = match.group(1)
        seq = next(line)
        index[header] = seq
    return index
"""

def split_protein_file(
    protein_dataset: Path, output_dir: Path, batch_size: int = 20
) -> List:
    """
    The protein dataset file is splitted by a number of sequence equals to the batch_size
    in batch files stored in 10 output directories.
    protein_dataset : Path for the protein dataset.
    output_dir : Output directory path.
    batch_size : Size of the batch, it needs to be equals to the number of threads
    to parallelise the sequence processing for each file.
    """
    batched_protein_files = []

    for i in range(0, 10):
        create_dir(output_dir, (f"bin_{i}"))
    with open(protein_dataset,"r", encoding="utf8") as file_in:
        seq_count = 0
        batch_count = 0
        current_record = ""
        initial_seq = True
        for line in file_in:
            match = re.search(r">(.+)$", line)
            # match header and is not first sequence, if the number of stored sequences in each file equals
            # the number of batch_size, a new file will be created and the current_record reset
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

def load_results_to_ensembl_db(
    main_script_dir,
    load_to_ensembl_db,
    genome_file,
    main_output_dir,
    db_details,
    num_threads,
):
    db_loading_script = os.path.join(
        main_script_dir, "support_scripts_perl", "load_gtf_ensembl.pl"
    )
    db_loading_dir = create_dir(main_output_dir, "db_loading")

    # Should collapse this into a function
    annotation_results_gtf_file = os.path.join(
        main_output_dir, "annotation_output", "annotation.gtf"
    )
    if os.path.exists(annotation_results_gtf_file):
        logger.info("Loading main geneset to db")
        batch_size = 200
        load_type = "gene"
        analysis_name = "ensembl"
        gtf_records = batch_gtf_records(
            annotation_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logger.error(
            "Did not find the main gene annotation file, so not loading. Path checked:\n"
            + annotation_results_gtf_file
        )

    rfam_results_gtf_file = os.path.join(main_output_dir, "rfam_output", "annotation.gtf")
    if os.path.exists(rfam_results_gtf_file):
        logger.info("Loading Rfam-based sncRNA genes to db")
        batch_size = 500
        load_type = "gene"
        analysis_name = "ncrna"
        gtf_records = batch_gtf_records(
            rfam_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logger.error(
            "Did not find an Rfam annotation file, so not loading. Path checked:\n"
            + rfam_results_gtf_file
        )

    trnascan_results_gtf_file = os.path.join(
        main_output_dir, "trnascan_output", "annotation.gtf"
    )
    if os.path.exists(trnascan_results_gtf_file):
        logger.info("Loading tRNAScan-SE tRNA genes to db")
        batch_size = 500
        load_type = "gene"
        analysis_name = "ncrna"
        gtf_records = batch_gtf_records(
            trnascan_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logger.error(
            "Did not find an tRNAScan-SE annotation file, so not loading. Path checked:\n"
            + trnascan_results_gtf_file
        )

    dust_results_gtf_file = os.path.join(main_output_dir, "dust_output", "annotation.gtf")
    if os.path.exists(dust_results_gtf_file):
        logger.info("Loading Dust repeats to db")
        batch_size = 500
        load_type = "single_line_feature"
        analysis_name = "dust"
        gtf_records = batch_gtf_records(
            dust_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logger.error(
            "Did not find a Dust annotation file, so not loading. Path checked:\n"
            + dust_results_gtf_file
        )

    red_results_gtf_file = os.path.join(main_output_dir, "red_output", "annotation.gtf")
    if os.path.exists(red_results_gtf_file):
        logger.info("Loading Red repeats to db")
        batch_size = 500
        load_type = "single_line_feature"
        analysis_name = "repeatdetector"
        gtf_records = batch_gtf_records(
            red_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logger.error(
            "Did not find a Red annotation file, so not loading. Path checked:\n"
            + red_results_gtf_file
        )

    trf_results_gtf_file = os.path.join(main_output_dir, "trf_output", "annotation.gtf")
    if os.path.exists(trf_results_gtf_file):
        logger.info("Loading TRF repeats to db")
        batch_size = 500
        load_type = "single_line_feature"
        analysis_name = "trf"
        gtf_records = batch_gtf_records(
            trf_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logger.error(
            "Did not find a TRF annotation file, so not loading. Path checked:\n"
            + trf_results_gtf_file
        )

    cpg_results_gtf_file = os.path.join(main_output_dir, "cpg_output", "annotation.gtf")
    if os.path.exists(cpg_results_gtf_file):
        logger.info("Loading CpG islands to db")
        batch_size = 500
        load_type = "single_line_feature"
        analysis_name = "cpg"
        gtf_records = batch_gtf_records(
            cpg_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logger.error(
            "Did not find a CpG annotation file, so not loading. Path checked:\n"
            + cpg_results_gtf_file
        )

    eponine_results_gtf_file = os.path.join(
        main_output_dir, "eponine_output", "annotation.gtf"
    )
    if os.path.exists(eponine_results_gtf_file):
        logger.info("Loading Eponine repeats to db")
        batch_size = 500
        load_type = "single_line_feature"
        analysis_name = "eponine"
        gtf_records = batch_gtf_records(
            eponine_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logger.error(
            "Did not find an Eponine annotation file, so not loading. Path checked:\n"
            + eponine_results_gtf_file
        )

    logger.info("Finished loading records to db")


def generic_load_records_to_ensembl_db(
    load_to_ensembl_db,
    db_loading_script,
    genome_file,
    db_details,
    db_loading_dir,
    load_type,
    analysis_name,
    gtf_records,
    num_threads,
):
    pool = multiprocessing.Pool(int(num_threads))
    for record_batch in gtf_records:
        pool.apply_async(
            multiprocess_load_records_to_ensembl_db,
            args=(
                load_to_ensembl_db,
                db_loading_script,
                genome_file,
                db_details,
                db_loading_dir,
                load_type,
                analysis_name,
                record_batch,
            ),
        )

    pool.close()
    pool.join()


def multiprocess_load_records_to_ensembl_db(
    load_to_ensembl_db,
    db_loading_script,
    genome_file,
    db_details,
    output_dir,
    load_type,
    analysis_name,
    record_batch,
):
    with tempfile.NamedTemporaryFile(
        mode="w+t", delete=False, dir=output_dir
    ) as gtf_temp_out:
        for line in record_batch:
            gtf_temp_out.writelines(line)
            gtf_temp_file_path = gtf_temp_out.name

    (db_name, db_host, db_port, db_user, db_pass) = db_details.split(",")

    loading_cmd = [
        "perl",
        db_loading_script,
        "-genome_file",
        genome_file,
        "-dbname",
        db_name,
        "-host",
        db_host,
        "-port",
        str(db_port),
        "-user",
        db_user,
        "-pass",
        db_pass,
        "-gtf_file",
        gtf_temp_file_path,
        "-analysis_name",
        analysis_name,
        "-load_type",
        load_type,
    ]

    if load_type == "gene" and analysis_name == "ensembl":
        loading_cmd.extend(
            [
                "-protein_coding_biotype",
                "anno_protein_coding",
                "-non_coding_biotype",
                "anno_lncRNA",
            ]
        )

        if load_to_ensembl_db == "single_transcript_genes":
            loading_cmd.append("-make_single_transcript_genes")

    logger.info(" ".join(loading_cmd))
    subprocess.run(loading_cmd)
    gtf_temp_out.close()
    os.remove(gtf_temp_file_path)  # doesn't seem to be working
    logger.info("Finished: " + gtf_temp_file_path)
    gc.collect()


def batch_gtf_records(input_gtf_file, batch_size, output_dir, record_type):
    records = []
    gtf_in = open(input_gtf_file)
    if record_type == "gene":
        # NOTE that the neverending variations on GTF reading/writing/merging
        # is becoming very messy
        # need to create a set of utility methods outside of this script
        # This one assumes the file has unique ids for the parent features.
        # It then batches them into
        # sets of records based on the batch size passed in
        record_counter = 0
        current_record_batch = []
        current_gene_id = ""
        line = gtf_in.readline()
        while line:
            if re.search(r"^#", line):
                line = gtf_in.readline()
                continue

            eles = line.split("\t")
            if not len(eles) == 9:
                line = gtf_in.readline()
                continue

            match = re.search(r'gene_id "([^"]+)"', line)
            gene_id = match.group(1)

            if not current_gene_id:
                record_counter += 1
                current_gene_id = gene_id

            if not gene_id == current_gene_id:
                record_counter += 1
                if record_counter % batch_size == 0:
                    records.append(current_record_batch)
                    current_record_batch = []
                current_gene_id = gene_id

            current_record_batch.append(line)
            line = gtf_in.readline()

        records.append(current_record_batch)

    elif record_type == "single_line_feature":
        record_counter = 0
        current_record_batch = []
        current_gene_id = ""
        line = gtf_in.readline()
        while line:
            if re.search(r"^#", line):
                line = gtf_in.readline()
                continue

            eles = line.split("\t")
            if not len(eles) == 9:
                line = gtf_in.readline()
                continue

            record_counter += 1

            if record_counter % batch_size == 0:
                records.append(current_record_batch)
                current_record_batch = []

            current_record_batch.append(line)
            line = gtf_in.readline()

        records.append(current_record_batch)

    gtf_in.close()

    return records


def run_find_orfs(genome_file, main_output_dir):
    min_orf_length = 600

    orf_output_dir = create_dir(main_output_dir, "orf_output")
    seq_region_lengths = get_seq_region_lengths(genome_file, 5000)
    for region_name in seq_region_lengths:
        region_length = seq_region_lengths[region_name]
        seq = get_sequence(region_name, 1, region_length, 1, genome_file, orf_output_dir)
        for phase in range(0, 6):
            find_orf_phased_region(
                region_name, seq, phase, min_orf_length, orf_output_dir
            )


def find_orf_phased_region(region_name, seq, phase, min_orf_length, orf_output_dir):
    current_index = phase
    orf_counter = 1
    if phase > 2:
        seq = reverse_complement(seq)
        current_index = current_index % 3

    orf_file_path = os.path.join(
        orf_output_dir, (region_name + ".phase" + str(phase) + ".orf.fa")
    )
    orf_out = open(orf_file_path, "w+")

    while current_index < len(seq):
        codon = seq[current_index : current_index + 3]
        if codon == "ATG":
            orf_seq = codon
            for j in range(current_index + 3, len(seq), 3):
                next_codon = seq[j : j + 3]
                if next_codon == "TAA" or next_codon == "TAG" or next_codon == "TGA":
                    orf_seq += next_codon
                    if len(orf_seq) >= min_orf_length:
                        orf_out.write(
                            ">"
                            + region_name
                            + "_phase"
                            + str(phase)
                            + "_orf"
                            + str(orf_counter)
                            + "\n"
                        )
                        orf_out.write(orf_seq + "\n")
                        orf_counter += 1
                        orf_seq = ""
                        break

                # If there's another met in phase, then put i to the start of
                # the codon after j so that only the longest ORF is found
                if next_codon == "ATG":
                    current_index = j + 3
                orf_seq += next_codon
        current_index += 3
    orf_out.close()

def bed_to_gff(input_dir, hints_file):
    gff_out = open(hints_file, "w+")
    exons_dict = {}
    for bed_file in glob.glob(input_dir + "/*.bed"):
        logger.info("Processing file for hints:")
        logger.info(bed_file)
        bed_in = open(bed_file)
        bed_lines = bed_in.readlines()
        for line in bed_lines:
            line = line.rstrip()
            elements = line.split("\t")
            seq_region_name = elements[0]
            offset = int(elements[1])
            hit_name = elements[3]
            strand = elements[5]
            block_sizes = elements[10].split(",")
            block_sizes = list(filter(None, block_sizes))
            block_starts = elements[11].split(",")
            block_starts = list(filter(None, block_starts))
            exons = bed_block_to_exons(block_sizes, block_starts, offset)
            for i, element in enumerate(exons):
                exon_coords = exons[i]
                exon_key = (
                    seq_region_name
                    + ":"
                    + exon_coords[0]
                    + ":"
                    + exon_coords[1]
                    + ":"
                    + strand
                )
                if exon_key in exons_dict:
                    exons_dict[exon_key][5] += 1
                else:
                    gff_list = [
                        seq_region_name,
                        "CDNA",
                        "exon",
                        exon_coords[0],
                        exon_coords[1],
                        1.0,
                        strand,
                        ".",
                    ]
                    exons_dict[exon_key] = gff_list

    for exon_key, gff_list in exons_dict.items():
        gff_list[5] = str(gff_list[5])
        gff_line = "\t".join(gff_list) + "\tsrc=W;mul=" + gff_list[5] + ";\n"
        gff_out.write(gff_line)

    gff_out.close()

    sorted_hints_out = open((hints_file + ".srt"), "w+")
    subprocess.run(
        ["sort", "-k1,1", "-k7,7", "-k4,4", "-k5,5", hints_file],
        stdout=sorted_hints_out,
    )
    sorted_hints_out.close()

def splice_junction_to_gff(input_dir, hints_file):
    sjf_out = open(hints_file, "w+")

    for sj_tab_file in glob.glob(input_dir + "/*.sj.tab"):
        sjf_in = open(sj_tab_file)
        sjf_lines = sjf_in.readlines()
        for line in sjf_lines:
            elements = line.split("\t")
            strand = "+"
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

            # For the moment treat multimapping and single mapping things as a combined score
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


def split_genome(genome_file, target_dir, min_seq_length):
    # This is the lazy initial way of just splitting into a dir
    # of files based on the toplevel sequence with a min sequence length filter
    # There are a couple of obvious improvements:
    # 1) Instead of making files for all seqs, just process N seqs
    #    parallel, where N = num_threads. Then you could clean up the seq file
    #    after each seq finishes, thus avoiding potentially
    #    having thousands of file in a dir
    # 2) Split the seq into even slices and process these in parallel
    #    (which the same cleanup as in 1). For sequences smaller than the
    #    target slice size, bundle them up together into a single file. Vastly
    #    more complex, partially implemented in the splice_genome method
    #    Allows for more consistency with parallelisation (since there should
    #    be no large outliers). But require a mapping strategy for the
    #    coords and sequence names and all the hints will need to be adjusted
    current_header = ""
    current_seq = ""

    file_in = open(genome_file)
    line = file_in.readline()
    while line:
        match = re.search(r">(.+)$", line)
        if match and current_header:
            if len(current_seq) > min_seq_length:
                file_out_name = os.path.join(target_dir, (current_header + ".split.fa"))
                if not os.path.exists(file_out_name):
                    file_out = open(file_out_name, "w+")
                    file_out.write(">" + current_header + "\n" + current_seq + "\n")
                    file_out.close()

                else:
                    print(
                        "Found an existing split file, so will not overwrite. File found:"
                    )
                    print(file_out_name)

            current_seq = ""
            current_header = match.group(1)
        elif match:
            current_header = match.group(1)
        else:
            current_seq += line.rstrip()

        line = file_in.readline()

    if len(current_seq) > min_seq_length:
        file_out_name = os.path.join(target_dir, (current_header + ".split.fa"))
        if not os.path.exists(file_out_name):
            file_out = open(file_out_name, "w+")
            file_out.write(">" + current_header + "\n" + current_seq + "\n")
            file_out.close()

        else:
            logger.info(
                "Found an existing split file, so will not overwrite. File found:"
            )
            logger.info(file_out_name)

    file_in.close()


def update_gtf_genes(parsed_gtf_genes, combined_results, validation_type):
    output_lines = []

    for gene_id in parsed_gtf_genes.keys():
        transcript_ids = parsed_gtf_genes[gene_id].keys()
        for transcript_id in transcript_ids:
            transcript_line = parsed_gtf_genes[gene_id][transcript_id]["transcript"]
            single_cds_exon_transcript = 0
            translation_match = re.search(
                r'; translation_coords "([^"]+)";', transcript_line
            )
            if translation_match:
                translation_coords = translation_match.group(1)
                translation_coords_list = translation_coords.split(":")
                # If the start exon coords of both exons are the same,
                # then it's the same exon and thus a single exon cds
                if translation_coords_list[0] == translation_coords_list[3]:
                    single_cds_exon_transcript = 1

            exon_lines = parsed_gtf_genes[gene_id][transcript_id]["exons"]
            validation_results = combined_results[transcript_id]
            rnasamba_coding_probability = float(validation_results[0])
            rnasamba_coding_potential = validation_results[1]
            cpc2_coding_probability = float(validation_results[2])
            cpc2_coding_potential = validation_results[3]
            transcript_length = int(validation_results[4])
            peptide_length = int(validation_results[5])
            diamond_e_value = None
            if len(validation_results) == 7:
                diamond_e_value = validation_results[6]

            avg_coding_probability = (
                rnasamba_coding_probability + cpc2_coding_probability
            ) / 2
            max_coding_probability = max(
                rnasamba_coding_probability, cpc2_coding_probability
            )

            match = re.search(r'; biotype "([^"]+)";', transcript_line)
            biotype = match.group(1)
            if biotype == "busco" or biotype == "protein":
                transcript_line = re.sub(
                    '; biotype "' + biotype + '";',
                    '; biotype "protein_coding";',
                    transcript_line,
                )
                output_lines.append(transcript_line)
                output_lines.extend(exon_lines)
                continue

            min_single_exon_pep_length = 100
            min_multi_exon_pep_length = 75
            min_single_source_probability = 0.8
            min_single_exon_probability = 0.9

            # Note that the below looks at validating things
            # under different levels of strictness
            # There are a few different continue statements,
            # where transcripts will be skipped resulting
            # in a smaller post validation file. It mainly
            # removes single coding exon genes with no real
            # support or for multi-exon lncRNAs that are less than 200bp long
            if single_cds_exon_transcript == 1 and validation_type == "relaxed":
                if diamond_e_value is not None:
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    rnasamba_coding_potential == "coding"
                    and cpc2_coding_potential == "coding"
                    and peptide_length >= min_single_exon_pep_length
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    (
                        rnasamba_coding_potential == "coding"
                        or cpc2_coding_potential == "coding"
                    )
                    and peptide_length >= min_single_exon_pep_length
                    and max_coding_probability >= min_single_source_probability
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                else:
                    continue
            elif single_cds_exon_transcript == 1 and validation_type == "moderate":
                if (
                    diamond_e_value is not None
                    and peptide_length >= min_single_exon_pep_length
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    (
                        rnasamba_coding_potential == "coding"
                        and cpc2_coding_potential == "coding"
                    )
                    and peptide_length >= min_single_exon_pep_length
                    and avg_coding_probability >= min_single_exon_probability
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                else:
                    continue
            else:
                if diamond_e_value is not None:
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    rnasamba_coding_potential == "coding"
                    and cpc2_coding_potential == "coding"
                    and peptide_length >= min_multi_exon_pep_length
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    (
                        rnasamba_coding_potential == "coding"
                        or cpc2_coding_potential == "coding"
                    )
                    and peptide_length >= min_multi_exon_pep_length
                    and max_coding_probability >= min_single_source_probability
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif transcript_length >= 200:
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "lncRNA";',
                        transcript_line,
                    )
                    transcript_line = re.sub(
                        ' translation_coords "[^"]+";', "", transcript_line
                    )
                else:
                    continue

            output_lines.append(transcript_line)
            output_lines.extend(exon_lines)

    return output_lines


def read_gtf_genes(gtf_file):
    gtf_genes = {}
    gtf_in = open(gtf_file)
    line = gtf_in.readline()
    while line:
        eles = line.split("\t")
        if not len(eles) == 9:
            line = gtf_in.readline()
            continue

        match = re.search(r'gene_id "([^"]+)".+transcript_id "([^"]+)"', line)

        if not match:
            line = gtf_in.readline()
            continue

        gene_id = match.group(1)
        transcript_id = match.group(2)
        feature_type = eles[2]
        if gene_id not in gtf_genes:
            gtf_genes[gene_id] = {}
        if feature_type == "transcript":
            gtf_genes[gene_id][transcript_id] = {}
            gtf_genes[gene_id][transcript_id]["transcript"] = line
            gtf_genes[gene_id][transcript_id]["exons"] = []
        elif feature_type == "exon":
            gtf_genes[gene_id][transcript_id]["exons"].append(line)
        line = gtf_in.readline()
    gtf_in.close()

    return gtf_genes





# def merge_gtf_files(file_paths,id_label):

#  gtf_file_path = os.path.join(output_dir,'annotation.gtf')
#  gtf_out = open(gtf_file_path,'w+')
#  for gtf_file_path in gtf_files:
#    gtf_file_name = os.path.basename(gtf_file_path)
#    match = re.search(r'\.rs(\d+)\.re(\d+)\.',gtf_file_name)
#    start_offset = int(match.group(1))
#    gtf_in = open(gtf_file_path,'r')
#    line = gtf_in.readline()
#    while line:
#      values = line.split("\t")
#      if len(values) == 9 and (values[2] in feature_types):
#        values[3] = str(int(values[3]) + (start_offset - 1))
#        values[4] = str(int(values[4]) + (start_offset - 1))
#        gtf_out.write("\t".join(values))
#        line = gtf_in.readline()
#    gtf_in.close()
#  gtf_out.close()


def multiprocess_generic(cmd):
    print(" ".join(cmd))
    subprocess.run(cmd)


def seq_region_names(genome_file):
    region_list = []

    file_in = open(genome_file)
    line = file_in.readline()
    while line:
        match = re.search(r">([^\s]+)", line)
        if match:
            region_name = match.group(1)
            if region_name == "MT":
                logger.info("Skipping region named MT")
                line = file_in.readline()
                continue
            else:
                region_list.append(match.group(1))
        line = file_in.readline()

    return region_list


def slice_genome(genome_file, target_dir, target_slice_size):
    # The below is sort of tested
    # Without the
    target_seq_length = 50000000
    min_seq_length = 1000
    current_header = ""
    current_seq = ""
    seq_dict = {}
    for line in seq:  # pylint:disable=used-before-assignment
        match = re.search(r">(.+)$", line)
        if match and current_header:
            seq_dict[current_header] = current_seq
            current_seq = ""
            current_header = match.group(1)
        elif match:
            current_header = match.group(1)
        else:
            current_seq += line.rstrip()

    seq_dict[current_header] = current_seq

    seq_buffer = 0
    file_number = 0
    file_name = "genome_file_" + str(file_number)

    for header in seq_dict:
        seq_iterator = 0
        seq = seq_dict[header]

        while len(seq) > target_seq_length:
            file_out = open(os.path.join(target_dir, file_name), "w+")
            subseq = seq[0:target_seq_length]
            file_out.write(
                ">" + header + "_sli" + str(seq_iterator) + "\n" + subseq + "\n"
            )
            file_out.close()
            seq = seq[target_seq_length:]
            seq_iterator += 1
            file_number += 1
            file_name = "genome_file_" + str(file_number)

        if len(seq) >= min_seq_length:
            file_name = "genome_file_" + str(file_number)
            file_out = open(os.path.join(file_name), "w+")
            file_out.write(">" + header + "_sli" + str(seq_iterator) + "\n" + seq + "\n")
            file_out.close()
            file_number += 1
            file_name = "genome_file_" + str(file_number)


def coallate_results(main_output_dir):
    results_dir = create_dir(main_output_dir, "results")
    output_dirs = [
        "augustus_output",
        "cpg_output",
        "dust_output",
        "eponine_output",
        "red_output",
        "rfam_output",
        "trf_output",
        "trnascan_output",
    ]
    for output_dir in output_dirs:
        match = re.search(r"(.+)_output", output_dir)
        result_type = match.group(1)
        results_path = os.path.join(main_output_dir, output_dir, "annotation.gtf")
        copy_path = os.path.join(results_dir, (result_type + ".gtf"))
        if os.path.exists(results_path):
            cpy_cmd = ["cp", results_path, copy_path]
            subprocess.run(cpy_cmd)



def check_transcriptomic_output(main_output_dir):

    # This will check across the various transcriptomic
    # dirs and check there's actually some data
    transcriptomic_dirs = [
        "scallop_output",
        "stringtie_output",
        "minimap2_output",
    ]
    total_lines = 0
    min_lines = 100000
    for transcriptomic_dir in transcriptomic_dirs:
        full_file_path = os.path.join(
            main_output_dir, transcriptomic_dir, "annotation.gtf"
        )
        if not os.path.exists(full_file_path):
            logger.warning(
                "Warning, no annotation.gtf found for "
                + transcriptomic_dir
                + ". This might be fine, e.g. no long read data were provided"
            )
            continue
        num_lines = sum(1 for line in open(full_file_path))
        total_lines = total_lines + num_lines
        logger.info(
            "For "
            + transcriptomic_dir
            + " found a total of "
            + str(num_lines)
            + " in the annotation.gtf file"
        )
    if total_lines == 0:
        raise IOError(
            "Anno was run with transcriptomic mode enabled,\
            but the transcriptomic annotation files are empty"
        )
    elif total_lines <= min_lines:
        raise IOError(
            "Anno was run with transcriptomic mode enabled, \
            but the total number of lines in the output \
            files were less than the min expected value"
            + "\n"
            "Found: " + str(total_lines) + "\n"
            "Min allowed: " + str(min_lines)
        )

    else:
        logger.info(
            "Found "
            + str(total_lines)
            + " total lines across the transcriptomic files. Checks passed"
        )
