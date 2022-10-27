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
import errno
import logging
import os
import pathlib
import re
import shutil
import sys

from typing import Union


# logging formats
logging_formatter_time_message = logging.Formatter(
    fmt="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

# set up base logger
logger = logging.getLogger("main_logger")
logger.setLevel(logging.DEBUG)
logger.propagate = False
# create console handler and add to logger
console_handler = logging.StreamHandler(sys.stderr)
console_handler.setLevel(logging.DEBUG)
console_handler.setFormatter(logging_formatter_time_message)
logger.addHandler(console_handler)


def add_log_file_handler(
    logger: logging.Logger,
    log_file_path: Union[pathlib.Path, str],
    logging_formatter: logging.Formatter = logging_formatter_time_message,
):
    """
    Create file handler and add to logger.
    """
    file_handler = logging.FileHandler(log_file_path, mode="a+")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(logging_formatter)
    logger.addHandler(file_handler)


def check_exe(exe_path):
    if not shutil.which(exe_path):
        raise OSError('Executable file not found at "%s"' % exe_path)


def check_file(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file_path)


def check_gtf_content(gtf_file: Union[pathlib.Path, str], content_obj):
    """
    This just checks how many transcript lines are in a GTF
    """
    transcript_count = 0
    with open(gtf_file) as gtf_in:
        for line in gtf_in:
            eles = line.split("\t")
            if not len(eles) == 9:
                continue
            if eles[2] == content_obj:
                transcript_count += 1
    logger.info("%s GTF transcript count: %s" % (gtf_file, transcript_count))
    return transcript_count


def create_dir(main_output_dir: Union[pathlib.Path, str], dir_name: str = None):
    """
    Create directory or subdirectory and log operations.
    Args:
        main_output_dir: main output directory path
        dir_name: optional subdirectory to be created
    Returns:
        created directory Path object
    """
    main_output_dir = pathlib.Path(main_output_dir)

    if dir_name:
        target_dir = main_output_dir / dir_name
    else:
        target_dir = main_output_dir

    try:
        target_dir.mkdir()
    except FileExistsError:
        logger.warning('Directory "%s" already exists' % target_dir)
    except OSError:
        logger.error('Failed to create directory "%s"' % target_dir)
        sys.exit()
    else:
        logger.info('Successfully created directory "%s"' % target_dir)

    return target_dir



def get_seq_region_lengths(genome_file: Union[pathlib.Path, str], min_seq_length: int):
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


def prlimit_command(command_list: list, virtual_memory_limit: int):
    """
    Uses the `prlimit` program to set a memory limit for a command list to be run with subprocess.
    prlimit - get and set process resource limits
    -v, --as[=limits]
        Address space limit.
    Args:
        command_list: original subprocess command list
        virtual_memory_limit: virtual memory limit in bytes
    Returns:
        memory limited subprocess command list
    """
    return ["prlimit", f"-v{virtual_memory_limit}"] + command_list

def load_results_to_ensembl_db(
    main_script_dir: pathlib.Path,
    load_to_ensembl_db,
    genome_file: Union[pathlib.Path, str],
    main_output_dir: pathlib.Path,
    db_details,
    num_threads: int,
):
    db_loading_script = main_script_dir / "support_scripts_perl" / "load_gtf_ensembl.pl"
    db_loading_dir = create_dir(main_output_dir, "db_loading")

    # Should collapse this into a function
    annotation_results_gtf_file = (
        main_output_dir / "annotation_output" / "annotation.gtf"
    )
    if annotation_results_gtf_file.exists():
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
            "Main gene annotation file not found, can't load:\n%s"
            % annotation_results_gtf_file
        )

    rfam_results_gtf_file = main_output_dir / "rfam_output" / "annotation.gtf"
    if rfam_results_gtf_file.exists():
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
            "Rfam annotation file not found, can't load:\n%s" % rfam_results_gtf_file
        )

    trnascan_results_gtf_file = main_output_dir / "trnascan_output" / "annotation.gtf"
    if trnascan_results_gtf_file.exists():
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
            "tRNAScan-SE annotation file not found, can't load:\n%s"
            % trnascan_results_gtf_file
        )

    dust_results_gtf_file = main_output_dir / "dust_output" / "annotation.gtf"
    if dust_results_gtf_file.exists():
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
            "Dust annotation file not found, can't load:\n%s" % dust_results_gtf_file
        )

    red_results_gtf_file = main_output_dir / "red_output" / "annotation.gtf"
    if red_results_gtf_file.exists():
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
            "Red annotation file not found, can't load:\n%s" % red_results_gtf_file
        )

    trf_results_gtf_file = main_output_dir / "trf_output" / "annotation.gtf"
    if trf_results_gtf_file.exists():
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
            "TRF annotation file not found, can't load:\n%s" % trf_results_gtf_file
        )

    cpg_results_gtf_file = main_output_dir / "cpg_output" / "annotation.gtf"
    if cpg_results_gtf_file.exists():
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
            "CpG annotation file not found, not loading:\n%s" % cpg_results_gtf_file
        )

    eponine_results_gtf_file = main_output_dir / "eponine_output" / "annotation.gtf"
    if eponine_results_gtf_file.exists():
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
            "Eponine annotation file not found, can't load:\n%s"
            % eponine_results_gtf_file
        )

    logger.info("Finished loading records to db")


def generic_load_records_to_ensembl_db(
    load_to_ensembl_db,
    db_loading_script,
    genome_file: Union[pathlib.Path, str],
    db_details,
    db_loading_dir,
    load_type,
    analysis_name,
    gtf_records,
    num_threads: int,
):
    pool = multiprocessing.Pool(num_threads)
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
    genome_file: Union[pathlib.Path, str],
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
            gtf_temp_out.write(line)
            gtf_temp_file_path = gtf_temp_out.name

    db_name, db_host, db_port, db_user, db_pass = db_details.split(",")

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

    logger.info("loading_cmd: %s" % " ".join(loading_cmd))
    subprocess.run(loading_cmd)
    gtf_temp_out.close()
    os.remove(gtf_temp_file_path)  # NOTE: doesn't seem to be working
    logger.info("Finished: %s" % gtf_temp_file_path)
    gc.collect()


def batch_gtf_records(input_gtf_file, batch_size, output_dir, record_type):
    records = []
    with open(input_gtf_file) as gtf_in:
        if record_type == "gene":
            # NOTE that the neverending variations on GTF reading/writing/merging is becoming very messy
            # need to create a set of utility methods outside of this script
            # This one assumes the file has unique ids for the parent features. It then batches them into
            # sets of records based on the batch size passed in
            record_counter = 0
            current_record_batch = []
            current_gene_id = ""
            for line in gtf_in:
                if re.search(r"^#", line):
                    continue

                eles = line.split("\t")
                if not len(eles) == 9:
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

            records.append(current_record_batch)

        elif record_type == "single_line_feature":
            record_counter = 0
            current_record_batch = []
            current_gene_id = ""
            for line in gtf_in:
                if re.search(r"^#", line):
                    continue

                eles = line.split("\t")
                if not len(eles) == 9:
                    continue

                record_counter += 1

                if record_counter % batch_size == 0:
                    records.append(current_record_batch)
                    current_record_batch = []

                current_record_batch.append(line)

            records.append(current_record_batch)

    return records   

def slice_output_to_gtf(
    output_dir: pathlib.Path,
    extension,
    unique_ids,
    feature_id_label,
    new_id_prefix,
):
    if not extension:
        extension = ".gtf"

    # Note that this does not make unique ids at the moment
    # In many cases this is fine because the ids are unique by seq region, but in cases like batching it can cause problems
    # So will add in a helper method to make ids unique

    # This holds keys of the current slice details with the gene id to form unique keys. Each time a new key is added
    # the overall gene counter is incremented and the value of the key is set to the new gene id. Any subsequent
    # lines with the same region/gene id key will then just get the new id without incrementing the counter
    gene_id_index = {}
    gene_transcript_id_index = {}
    gene_counter = 1

    # Similar to the gene id index, this will have a key that is based on the slice details, gene id and transcript id. If there
    # is no existing entry, the transcript key will be added and the transcript counter is incremented. If there is a key then
    # the transcript id will be replaced with the new transcript id (which is based on the new gene id and transcript counter)
    # Example key KS8000.rs1.re1000000.gene_1.transcript_1 =
    transcript_id_count_index = {}

    feature_counter = 1

    feature_types = ["exon", "transcript", "repeat", "simple_feature"]
    gtf_output_file_path = output_dir / "annotation.gtf"
    with open(gtf_output_file_path, "w+") as gtf_out:
        for gtf_input_file_path in output_dir.glob("*{extension}"):
            if gtf_input_file_path.stat().st_size == 0:
                logger.info("File is empty, will skip:\n%s" % gtf_input_file_path)
                continue

            gtf_file_name = gtf_input_file_path.name
            match = re.search(r"\.rs(\d+)\.re(\d+)\.", gtf_file_name)
            start_offset = int(match.group(1))
            with open(gtf_input_file_path, "r") as gtf_in:
                for line in gtf_in:
                    values = line.split("\t")
                    if len(values) == 9 and (values[2] in feature_types):
                        values[3] = str(int(values[3]) + (start_offset - 1))
                        values[4] = str(int(values[4]) + (start_offset - 1))
                        if unique_ids:
                            # Maybe make a unique id based on the feature type
                            # Basically region/feature id should be unique at this point, so could use region_id and current_id is key, value is the unique id that is incremented
                            attribs = values[8]

                            # This bit assigns unique gene/transcript ids if the line contains gene_id/transcript_id
                            match_gene_type = re.search(
                                r'(gene_id +"([^"]+)").+(transcript_id +"([^"]+)")',
                                line,
                            )
                            if match_gene_type:
                                full_gene_id_string = match_gene_type.group(1)
                                current_gene_id = match_gene_type.group(2)
                                full_transcript_id_string = match_gene_type.group(3)
                                current_transcript_id = match_gene_type.group(4)
                                gene_id_key = f"{gtf_file_name}.{current_gene_id}"
                                transcript_id_key = (
                                    f"{gene_id_key}.{current_transcript_id}"
                                )
                                if gene_id_key not in gene_id_index:
                                    new_gene_id = f"gene{gene_counter}"
                                    gene_id_index[gene_id_key] = new_gene_id
                                    attribs = re.sub(
                                        full_gene_id_string,
                                        f'gene_id "{new_gene_id}"',
                                        attribs,
                                    )
                                    transcript_id_count_index[gene_id_key] = 1
                                    gene_counter += 1
                                else:
                                    new_gene_id = gene_id_index[gene_id_key]
                                    attribs = re.sub(
                                        full_gene_id_string,
                                        f'gene_id "{new_gene_id}"',
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
                                        f'transcript_id "{new_transcript_id}"',
                                        attribs,
                                    )
                                    transcript_id_count_index[gene_id_key] += 1
                                else:
                                    new_transcript_id = gene_transcript_id_index[
                                        transcript_id_key
                                    ]
                                    attribs = re.sub(
                                        full_transcript_id_string,
                                        f'transcript_id "{new_transcript_id}"',
                                        attribs,
                                    )
                                values[8] = attribs

                            # If you don't match a gene line, try a feature line
                            else:
                                match_feature_type = re.search(
                                    r"(" + feature_id_label + ' +"([^"]+)")', line
                                )
                                if match_feature_type:
                                    full_feature_id_string = match_feature_type.group(1)
                                    current_feature_id = match_feature_type.group(2)
                                    new_feature_id = f"{new_id_prefix}{feature_counter}"
                                    attribs = re.sub(
                                        full_feature_id_string,
                                        f'{feature_id_label} "{new_feature_id}"',
                                        attribs,
                                    )
                                    feature_counter += 1
                                    values[8] = attribs

                        gtf_out.write("\t".join(values))
                    else:
                        logger.info(
                            "Feature type not recognised, will skip: %s" % values[2]
                        )

   
def convert_gff_to_gtf(gff_file):
    gtf_string = ""
    with open(gff_file) as file_in:
        for line in file_in:
            # match = re.search(r"genBlastG",line)
            # if match:
            results = line.split()
            if not len(results) == 9:
                continue
            if results[2] == "coding_exon":
                results[2] = "exon"
            attributes = set_attributes(results[8], results[2])
            results[8] = attributes
            converted_line = "\t".join(results)
            gtf_string += f"{converted_line}\n"
    
    return gtf_string   
def set_attributes(attributes, feature_type):
    converted_attributes = ""
    split_attributes = attributes.split(";")
    if feature_type == "transcript":
        match = re.search(r"Name\=(.+)$", split_attributes[1])
        name = match.group(1)
        converted_attributes = f'gene_id "{name}"; transcript_id "{name}";'
    elif feature_type == "exon":
        match = re.search(r"\-E(\d+);Parent\=(.+)\-R\d+\-\d+\-", attributes)
        exon_rank = match.group(1)
        name = match.group(2)
        converted_attributes = (
            f'gene_id "{name}"; transcript_id "{name}"; exon_number "{exon_rank}";'
        )

    return converted_attributes


# Example genBlast output
# 1       genBlastG       transcript      131128674       131137049       252.729 -       .       ID=259447-R1-1-A1;Name=259447;PID=84.65;Coverage=94.22;Note=PID:84.65-Cover:94.22
# 1       genBlastG       coding_exon     131137031       131137049       .       -       .       ID=259447-R1-1-A1-E1;Parent=259447-R1-1-A1
# 1       genBlastG       coding_exon     131136260       131136333       .       -       .       ID=259447-R1-1-A1-E2;Parent=259447-R1-1-A1
# 1       genBlastG       coding_exon     131128674       131130245       .       -       .       ID=259447-R1-1-A1-E3;Parent=259447-R1-1-A1
##sequence-region       1_group1        1       4534
# 1       genBlastG       transcript      161503457       161503804       30.94   +       .       ID=259453-R1-1-A1;Name=259453;PID=39.46;Coverage=64.97;Note=PID:39.46-Cover:64.97
# 1       genBlastG       coding_exon     161503457       161503804       .       +       .       ID=259453-R1-1-A1-E1;Parent=259453-R1-1-A1
##sequence-region       5_group1        1       4684
# 5       genBlastG       transcript      69461063        69461741        86.16   +       .       ID=259454-R1-1-A1;Name=259454;PID=82.02;Coverage=91.67;Note=PID:82.02-Cover:91.67
# 5       genBlastG       coding_exon     69461063        69461081        .       +       .       ID=259454-R1-1-A1-E1;Parent=259454-R1-1-A1
# 5       genBlastG       coding_exon     69461131        69461741        .       +       .       ID=259454-R1-1-A1-E2;Parent=259454-R1-1-A1


def create_slice_ids(
    seq_region_lengths,
    slice_size: int = 1_000_000,
    overlap: int = 0,
    min_length: int = 0,
):
    slice_ids = []
    for region, region_length in seq_region_lengths.items():
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
                # If the start exon coords of both exons are the same, then it's the same exon and thus a single exon cds
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
                    '; biotype "{biotype}";',
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

            # Note that the below looks at validating things under different levels of strictness
            # There are a few different continue statements, where transcripts will be skipped resulting
            # in a smaller post validation file. It mainly removes single coding exon genes with no real
            # support or for multi-exon lncRNAs that are less than 200bp long
            if single_cds_exon_transcript == 1 and validation_type == "relaxed":
                if diamond_e_value is not None:
                    transcript_line = re.sub(
                        '; biotype "{biotype}";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    rnasamba_coding_potential == "coding"
                    and cpc2_coding_potential == "coding"
                    and peptide_length >= min_single_exon_pep_length
                ):
                    transcript_line = re.sub(
                        '; biotype "{biotype}";',
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
                        '; biotype "{biotype}";',
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
                        '; biotype "{biotype}";',
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
                        '; biotype "{biotype}";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                else:
                    continue
            else:
                if diamond_e_value is not None:
                    transcript_line = re.sub(
                        '; biotype "{biotype}";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    rnasamba_coding_potential == "coding"
                    and cpc2_coding_potential == "coding"
                    and peptide_length >= min_multi_exon_pep_length
                ):
                    transcript_line = re.sub(
                        '; biotype "{biotype}";',
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
                        '; biotype "{biotype}";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif transcript_length >= 200:
                    transcript_line = re.sub(
                        '; biotype "{biotype}";', '; biotype "lncRNA";', transcript_line
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
    with open(gtf_file) as gtf_in:
        for line in gtf_in:
            eles = line.split("\t")
            if not len(eles) == 9:
                continue

            match = re.search(r'gene_id "([^"]+)".+transcript_id "([^"]+)"', line)

            if not match:
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

    return gtf_genes


def fasta_to_dict(fasta_list):
    index = {}
    it = iter(fasta_list)
    for header in it:
        match = re.search(r">(.+)\n$", header)
        header = match.group(1)
        seq = next(it)
        index[header] = seq
    return index


def subprocess_run_and_log(command):
    logger.info("subprocess_run_and_log command: %s" % " ".join(command))
    subprocess.run(command)


def get_sequence(
    seq_region,
    start: int,
    end: int,
    strand: int,
    fasta_file,
    output_dir: Union[pathlib.Path, str],
):
    start -= 1
    bedtools_path = "bedtools"

    # This creates a tempfile and writes the bed info to it based on whatever information
    # has been passed in about the sequence. Then runs bedtools getfasta. The fasta file
    # should have a faidx. This can be created with the create_faidx static method prior
    # to fetching sequence
    with tempfile.NamedTemporaryFile(
        mode="w+t", delete=False, dir=output_dir
    ) as bed_temp_file:
        bed_temp_file.write(f"{seq_region}\t{start}\t{end}")
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


def reverse_complement(sequence: str) -> str:
    rev_matrix = str.maketrans("atgcATGC", "tacgTACG")
    return sequence.translate(rev_matrix)[::-1]


def get_seq_region_names(genome_file: Union[pathlib.Path, str]):
    region_list = []
    with open(genome_file) as file_in:
        for line in file_in:
            match = re.search(r">([^\s]+)", line)
            if match:
                region_name = match.group(1)
                if region_name == "MT":
                    logger.info('Skipping region named "MT"')
                    continue
                else:
                    region_list.append(match.group(1))

    return region_list


def run_find_orfs(genome_file, main_output_dir):

    min_orf_length = 600

    orf_output_dir = create_dir(main_output_dir, "orf_output")
    seq_region_lengths = get_seq_region_lengths(genome_file, 5000)
    for region_name in seq_region_lengths:
        region_length = seq_region_lengths[region_name]
        seq = get_sequence(
            region_name, 1, region_length, 1, genome_file, orf_output_dir
        )
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

                # If there's another met in phase, then put i to the start of the codon after j so that only the longest ORF is found
                if next_codon == "ATG":
                    current_index = j + 3
                orf_seq += next_codon
        current_index += 3
    orf_out.close()
 

def list_to_string(original_list: List, separator: str = " ") -> str:
    """
    Create a string with the original list elements string representations concatenated
    with spaces between them.
    Args:
        original_list: original list
        separator: character (or string) to separate the elements string representations
        in the resulting string
    Returns:
        generated list element string
    """
    return str.join(separator, [str(element) for element in original_list])

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
    # This is the lazy initial way of just splitting into a dir of files based on the toplevel sequence with a min sequence length filter
    # There are a couple of obvious improvements:
    # 1) Instead of making files for all seqs, just process N seqs parallel, where N = num_threads. Then you could clean up the seq file
    #    after each seq finishes, thus avoiding potentially having thousands of file in a dir
    # 2) Split the seq into even slices and process these in parallel (which the same cleanup as in 1). For sequences smaller than the
    #    target slice size, bundle them up together into a single file. Vastly more complex, partially implemented in the splice_genome method
    #    Allows for more consistency with parallelisation (since there should be no large outliers). But require a mapping strategy for the
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

def multiprocess_generic(cmd):
    print(" ".join(cmd))
    subprocess.run(cmd)
def slice_genome(genome_file, target_dir, target_slice_size):
    # The below is sort of tested
    # Without the
    target_seq_length = 50000000
    min_seq_length = 1000
    current_header = ""
    current_seq = ""
    seq_dict = {}
    for line in seq:
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
            file_out.write(
                ">" + header + "_sli" + str(seq_iterator) + "\n" + seq + "\n"
            )
            file_out.close()
            file_number += 1
            file_name = "genome_file_" + str(file_number)



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

    
