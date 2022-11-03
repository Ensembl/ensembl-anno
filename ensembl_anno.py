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
import argparse
import gc
import glob
import io
import math
import multiprocessing
import os
import pathlib
import random
import re
import shutil
import signal
import subprocess
import tempfile

from typing import Dict, List, Union

# project imports
from utils import (
    add_log_file_handler,
    check_exe,
    check_file,
    check_gtf_content,
    create_dir,
    create_paired_paths,
    get_seq_region_lengths,
    list_to_string,
    logger,
    prlimit_command,
)


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

    logger.info("loading_cmd: %s" % list_to_string(loading_cmd))
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


def run_repeatmasker_regions(
    genome_file: Union[pathlib.Path, str],
    repeatmasker_path,
    library,
    species,
    main_output_dir,
    num_threads: int,
):
    if not repeatmasker_path:
        repeatmasker_path = "RepeatMasker"

    check_exe(repeatmasker_path)
    repeatmasker_output_dir = create_dir(main_output_dir, "repeatmasker_output")
    os.chdir(repeatmasker_output_dir)

    output_file = repeatmasker_output_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "repeat")
        if transcript_count > 0:
            logger.info("Repeatmasker gtf file already exists, skipping analysis")
            return

    logger.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, min_seq_length=5000)
    slice_ids = create_slice_ids(
        seq_region_lengths, slice_size=1_000_000, overlap=0, min_length=5000
    )

    if not library:
        if not species:
            species = "homo"

        generic_repeatmasker_cmd = [
            repeatmasker_path,
            "-nolow",
            "-species",
            species,
            "-engine",
            "crossmatch",
            "-dir",
            repeatmasker_output_dir,
        ]
    else:
        generic_repeatmasker_cmd = [
            repeatmasker_path,
            "-nolow",
            "-lib",
            library,
            "-engine",
            "crossmatch",
            "-dir",
            repeatmasker_output_dir,
        ]

    logger.info("Running RepeatMasker processes")
    pool = multiprocessing.Pool(num_threads)
    tasks = []
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_repeatmasker,
            args=(
                generic_repeatmasker_cmd,
                slice_id,
                genome_file,
                repeatmasker_output_dir,
            ),
        )

    pool.close()
    pool.join()
    slice_output_to_gtf(
        repeatmasker_output_dir, ".rm.gtf", 1, "repeat_id", "repeatmask"
    )


def multiprocess_repeatmasker(
    generic_repeatmasker_cmd,
    slice_id,
    genome_file: Union[pathlib.Path, str],
    repeatmasker_output_dir: Union[pathlib.Path, str],
):
    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing slice to find repeats with RepeatMasker: %s:%s:%s"
        % (region_name, start, end)
    )
    seq = get_sequence(
        seq_region=region_name,
        start=start,
        end=end,
        strand=1,
        fasta_file=genome_file,
        output_dir=repeatmasker_output_dir,
    )

    slice_file_name = f"{region_name}.rs{start}.re{end}"
    region_fasta_file_path = repeatmasker_output_dir / f"{slice_file_name}.fa"
    with open(region_fasta_file_path, "w+") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")

    region_results_file_path = f"{region_fasta_file_path}.rm.gtf"
    repeatmasker_output_file_path = f"{region_fasta_file_path}.out"
    repeatmasker_masked_file_path = f"{region_fasta_file_path}.masked"
    repeatmasker_tbl_file_path = f"{region_fasta_file_path}.tbl"
    repeatmasker_log_file_path = f"{region_fasta_file_path}.log"
    repeatmasker_cat_file_path = f"{region_fasta_file_path}.cat"

    repeatmasker_cmd = generic_repeatmasker_cmd.copy()
    repeatmasker_cmd.append(region_fasta_file_path)
    logger.info("repeatmasker_cmd: %s" % list_to_string(repeatmasker_cmd))
    subprocess.run(repeatmasker_cmd)

    create_repeatmasker_gtf(
        repeatmasker_output_file_path, region_results_file_path, region_name
    )

    os.remove(region_fasta_file_path)
    if os.path.exists(region_results_file_path):
        os.remove(region_results_file_path)
    if os.path.exists(repeatmasker_masked_file_path):
        os.remove(repeatmasker_masked_file_path)
    if os.path.exists(repeatmasker_tbl_file_path):
        os.remove(repeatmasker_tbl_file_path)
    if os.path.exists(repeatmasker_log_file_path):
        os.remove(repeatmasker_log_file_path)
    if os.path.exists(repeatmasker_cat_file_path):
        os.remove(repeatmasker_cat_file_path)


def create_repeatmasker_gtf(
    repeatmasker_output_file_path, region_results_file_path, region_name
):
    with open(repeatmasker_output_file_path, "r") as repeatmasker_in, open(
        region_results_file_path, "w+"
    ) as repeatmasker_out:
        repeat_count = 1
        for line in repeatmasker_in:
            result_match = re.search(r"^\s*\d+\s+", line)
            if result_match:
                results = line.split()
                if results[-1] == "*":
                    results.pop()
                if not len(results) == 15:
                    continue

                score = results[0]
                start = results[5]
                end = results[6]
                strand = results[8]
                repeat_name = results[9]
                repeat_class = results[10]
                if strand == "+":
                    repeat_start = results[11]
                    repeat_end = results[12]
                else:
                    repeat_start = results[13]
                    repeat_end = results[12]
                    strand = "-"

                gtf_line = f'{region_name}\tRepeatMasker\trepeat\t{start}\t{end}\t.\t{strand}\t.\trepeat_id "{repeat_count}"; repeat_name "{repeat_name}"; repeat_class "{repeat_class}"; repeat_start "{repeat_start}"; repeat_end "{repeat_end}"; score "{score}";\n'

                repeatmasker_out.write(gtf_line)
                repeat_count += 1


def run_eponine_regions(
    genome_file: Union[pathlib.Path, str],
    java_path,
    eponine_path,
    main_output_dir,
    num_threads: int,
):
    if not java_path:
        java_path = "java"

    if not eponine_path:
        eponine_path = "/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/opt/eponine/libexec/eponine-scan.jar"

    check_file(eponine_path)
    check_exe(java_path)

    eponine_output_dir = create_dir(main_output_dir, "eponine_output")

    output_file = eponine_output_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "simple_feature")
        if transcript_count > 0:
            logger.info("Eponine gtf file already exists, skipping analysis")
            return

    logger.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, min_seq_length=5000)
    slice_ids = create_slice_ids(
        seq_region_lengths, slice_size=1_000_000, overlap=0, min_length=5000
    )

    threshold = "0.999"
    generic_eponine_cmd = [
        java_path,
        "-jar",
        eponine_path,
        "-threshold",
        threshold,
        "-seq",
    ]
    logger.info("Running Eponine processes")
    pool = multiprocessing.Pool(num_threads)
    tasks = []
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_eponine,
            args=(
                generic_eponine_cmd,
                slice_id,
                genome_file,
                eponine_output_dir,
            ),
        )

    pool.close()
    pool.join()
    slice_output_to_gtf(eponine_output_dir, ".epo.gtf", 1, "feature_id", "eponine")


def multiprocess_eponine(
    generic_eponine_cmd,
    slice_id,
    genome_file: Union[pathlib.Path, str],
    eponine_output_dir: Union[pathlib.Path, str],
):
    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing slice to find repeats with Eponine: %s:%s:%s"
        % (region_name, start, end)
    )
    seq = get_sequence(
        seq_region=region_name,
        start=start,
        end=end,
        strand=1,
        fasta_file=genome_file,
        output_dir=eponine_output_dir,
    )

    slice_file_name = f"{region_name}.rs{start}.re{end}"
    region_fasta_file_path = eponine_output_dir / f"{slice_file_name}.fa"

    with open(region_fasta_file_path, "w+") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")

    region_results_file_path = eponine_output_dir / f"{slice_file_name}.epo.gtf"

    eponine_output_file_path = f"{region_fasta_file_path}.epo"
    eponine_out = open(eponine_output_file_path, "w+")

    eponine_cmd = generic_eponine_cmd.copy()
    eponine_cmd.append(region_fasta_file_path)

    logger.info("eponine_cmd: %s" % list_to_string(eponine_cmd))
    subprocess.run(eponine_cmd, stdout=eponine_out)
    eponine_out.close()

    create_eponine_gtf(eponine_output_file_path, region_results_file_path, region_name)
    os.remove(eponine_output_file_path)
    os.remove(region_fasta_file_path)


def create_eponine_gtf(eponine_output_file_path, region_results_file_path, region_name):
    with open(eponine_output_file_path, "r") as eponine_in, open(
        region_results_file_path, "w+"
    ) as eponine_out:
        feature_count = 1
        for line in eponine_in:
            result_match = re.search(r"^" + region_name, line)

            if result_match:
                results = line.split()
                start = int(results[3])
                end = int(results[4])
                score = float(results[5])
                strand = results[6]

                # There's a one base offset on the reverse strand
                if strand == "-":
                    start -= 1
                    end -= 1

                gtf_line = f'{region_name}\tEponine\tsimple_feature\t{start}\t{end}\t.\t{strand}\t.\tfeature_id "{feature_count}"; score "{score}";\n'
                eponine_out.write(gtf_line)
                feature_count += 1


def run_cpg_regions(
    genome_file: Union[pathlib.Path, str], cpg_path, main_output_dir, num_threads: int
):
    if not cpg_path:
        cpg_path = "cpg_lh"

    check_exe(cpg_path)
    cpg_output_dir = create_dir(main_output_dir, "cpg_output")

    output_file = cpg_output_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "simple_feature")
        if transcript_count > 0:
            logger.info("Cpg gtf file already exists, skipping analysis")
            return

    logger.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, min_seq_length=5000)
    slice_ids = create_slice_ids(
        seq_region_lengths, slice_size=1_000_000, overlap=0, min_length=5000
    )

    logger.info("Running CpG processes")
    pool = multiprocessing.Pool(num_threads)
    tasks = []
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_cpg,
            args=(
                cpg_path,
                slice_id,
                genome_file,
                cpg_output_dir,
            ),
        )

    pool.close()
    pool.join()
    slice_output_to_gtf(cpg_output_dir, ".cpg.gtf", 1, "feature_id", "cpg")


def multiprocess_cpg(
    cpg_path, slice_id, genome_file, cpg_output_dir: Union[pathlib.Path, str]
):
    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing slice to find CpG islands with cpg_lh: %s:%s:%s"
        % (region_name, start, end)
    )
    seq = get_sequence(
        seq_region=region_name,
        start=start,
        end=end,
        strand=1,
        fasta_file=genome_file,
        output_dir=cpg_output_dir,
    )

    slice_file_name = f"{region_name}.rs{start}.re{end}"
    region_fasta_file_path = cpg_output_dir / f"{slice_file_name}.fa"

    with open(region_fasta_file_path, "w+") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")

    region_results_file_path = cpg_output_dir / f"{slice_file_name}.cpg.gtf"

    cpg_output_file_path = f"{region_fasta_file_path}.cpg"
    cpg_out = open(cpg_output_file_path, "w+")

    cpg_cmd = [cpg_path, region_fasta_file_path]
    logger.info("cpg_cmd: %s" % list_to_string(cpg_cmd))
    subprocess.run(cpg_cmd, stdout=cpg_out)
    cpg_out.close()

    create_cpg_gtf(cpg_output_file_path, region_results_file_path, region_name)
    os.remove(cpg_output_file_path)
    os.remove(region_fasta_file_path)


def create_cpg_gtf(cpg_output_file_path, region_results_file_path, region_name):
    cpg_min_length = 400
    cpg_min_gc_content = 50
    cpg_min_oe = 0.6

    with open(cpg_output_file_path, "r") as cpg_in, open(
        region_results_file_path, "w+"
    ) as cpg_out:
        feature_count = 1
        for line in cpg_in:
            result_match = re.search(r"^" + region_name, line)
            if result_match:
                results = line.split()
                start = int(results[1])
                end = int(results[2])
                length = end - start + 1
                score = float(results[3])
                gc_content = float(results[6])
                oe = results[7]

                if oe == "-" or oe == "inf":
                    oe = 0
                else:
                    oe = float(oe)

                if (
                    length >= cpg_min_length
                    and gc_content >= cpg_min_gc_content
                    and oe >= cpg_min_oe
                ):
                    gtf_line = f'{region_name}\tCpG\tsimple_feature\t{start}\t{end}\t.\t+\t.\tfeature_id "{feature_count}"; score "{score}";\n'
                    cpg_out.write(gtf_line)
                    feature_count += 1


def run_trnascan_regions(
    genome_file: Union[pathlib.Path, str],
    trnascan_path,
    trnascan_filter_path,
    main_output_dir,
    num_threads: int,
):
    if not trnascan_path:
        trnascan_path = "/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/bin/tRNAscan-SE"
    check_exe(trnascan_path)
    logger.info("trnascan_path: %s" % trnascan_path)

    if not trnascan_filter_path:
        trnascan_filter_path = "/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/bin/EukHighConfidenceFilter"
    # check_exe(trnascan_filter_path)
    check_file(trnascan_filter_path)
    logger.info("trnascan_filter_path: %s" % trnascan_filter_path)

    trnascan_output_dir = create_dir(main_output_dir, "trnascan_output")

    output_file = trnascan_output_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Trnascan gtf file already exists, skipping analysis")
            return

    logger.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, min_seq_length=5000)
    slice_ids = create_slice_ids(
        seq_region_lengths, slice_size=1_000_000, overlap=0, min_length=5000
    )

    generic_trnascan_cmd = [
        trnascan_path,
        None,
        "-o",
        None,
        "-f",
        None,
        "-H",
        "-q",
        "--detail",
        "-Q",
    ]
    logger.info("Running tRNAscan-SE processes")
    pool = multiprocessing.Pool(num_threads)
    tasks = []
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_trnascan,
            args=(
                generic_trnascan_cmd,
                slice_id,
                genome_file,
                trnascan_filter_path,
                trnascan_output_dir,
            ),
        )

    pool.close()
    pool.join()
    slice_output_to_gtf(trnascan_output_dir, ".trna.gtf", 1, None, None)


def multiprocess_trnascan(
    generic_trnascan_cmd,
    slice_id,
    genome_file: Union[pathlib.Path, str],
    trnascan_filter_path,
    trnascan_output_dir: Union[pathlib.Path, str],
):
    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing slice to find tRNAs using tRNAscan-SE: %s:%s:%s"
        % (region_name, start, end)
    )
    seq = get_sequence(
        seq_region=region_name,
        start=start,
        end=end,
        strand=1,
        fasta_file=genome_file,
        output_dir=trnascan_output_dir,
    )

    slice_file_name = f"{region_name}.rs{start}.re{end}"
    region_fasta_file_name = f"{slice_file_name}.fa"
    region_fasta_file_path = trnascan_output_dir / region_fasta_file_name
    with open(region_fasta_file_path, "w+") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")

    region_results_file_path = trnascan_output_dir / f"{slice_file_name}.trna.gtf"

    trnascan_output_file_path = f"{region_fasta_file_path}.trna"
    trnascan_ss_output_file_path = f"{trnascan_output_file_path}.ss"

    # The filter takes an output dir and a prefix and then uses those to make a path to a .out file
    trnascan_filter_file_prefix = f"{region_fasta_file_name}.filt"
    trnascan_filter_file_name = f"{trnascan_filter_file_prefix}.out"
    trnascan_filter_log_file_name = f"{trnascan_filter_file_prefix}.log"
    trnascan_filter_ss_file_name = f"{trnascan_filter_file_prefix}.ss"
    trnascan_filter_file_path = trnascan_output_dir / trnascan_filter_file_name
    trnascan_filter_log_file_path = trnascan_output_dir / trnascan_filter_log_file_name
    trnascan_filter_ss_file_path = trnascan_output_dir / trnascan_filter_ss_file_name

    trnascan_cmd = generic_trnascan_cmd.copy()
    trnascan_cmd[1] = region_fasta_file_path
    trnascan_cmd[3] = trnascan_output_file_path
    trnascan_cmd[5] = trnascan_ss_output_file_path

    logger.info("tRNAscan-SE command:\n%s" % list_to_string(trnascan_cmd))
    subprocess.run(trnascan_cmd)

    # If we have a blank output file at this point we want to stop and remove whatever files
    # have been created instead of moving onto the filter
    if os.stat(trnascan_output_file_path).st_size == 0:
        os.remove(trnascan_output_file_path)
        os.remove(region_fasta_file_path)
        if os.path.exists(trnascan_ss_output_file_path):
            os.remove(trnascan_ss_output_file_path)
        return

    filter_cmd = [
        trnascan_filter_path,
        "--result",
        trnascan_output_file_path,
        "--ss",
        trnascan_ss_output_file_path,
        "--output",
        trnascan_output_dir,
        "--prefix",
        trnascan_filter_file_prefix,
    ]
    logger.info("tRNAscan-SE filter command:\n%s" % list_to_string(filter_cmd))
    subprocess.run(filter_cmd)

    create_trnascan_gtf(
        region_results_file_path, trnascan_filter_file_path, region_name
    )
    if os.path.exists(trnascan_output_file_path):
        os.remove(trnascan_output_file_path)
    if os.path.exists(trnascan_ss_output_file_path):
        os.remove(trnascan_ss_output_file_path)
    if os.path.exists(trnascan_filter_file_path):
        os.remove(trnascan_filter_file_path)
    if os.path.exists(trnascan_filter_log_file_path):
        os.remove(trnascan_filter_log_file_path)
    if os.path.exists(trnascan_filter_ss_file_path):
        os.remove(trnascan_filter_ss_file_path)
    if os.path.exists(region_fasta_file_path):
        os.remove(region_fasta_file_path)


def create_trnascan_gtf(
    region_results_file_path, trnascan_filter_file_path, region_name
):
    with open(trnascan_filter_file_path, "r") as trna_in, open(
        region_results_file_path, "w+"
    ) as trna_out:
        gene_counter = 1
        for line in trna_in:
            result_match = re.search(r"^" + region_name, line)
            if result_match:
                results = line.split()
                start = int(results[2])
                end = int(results[3])
                trna_type = results[4]
                score = results[8]

                strand = "+"
                if start > end:
                    strand = "-"
                    temp_end = start
                    start = end
                    end = temp_end

                biotype = "tRNA_pseudogene"
                high_confidence_match = re.search(r"high confidence set", line)
                if high_confidence_match:
                    biotype = "tRNA"

                transcript_string = f'{region_name}\ttRNAscan\ttranscript\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gene_counter}"; transcript_id "{gene_counter}"; biotype "{biotype}";\n'
                exon_string = f'{region_name}\ttRNAscan\texon\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gene_counter}"; transcript_id "{gene_counter}"; exon_number "1"; biotype "{biotype}";\n'

                trna_out.write(transcript_string)
                trna_out.write(exon_string)
                gene_counter += 1


def run_dust_regions(
    genome_file: Union[pathlib.Path, str], dust_path, main_output_dir, num_threads: int
):
    if not dust_path:
        dust_path = "dustmasker"

    check_exe(dust_path)
    dust_output_dir = create_dir(main_output_dir, "dust_output")

    output_file = dust_output_dir / "annotation.gtf"
    logger.info("output_file: %s" % output_file)
    if output_file.is_file():
        transcript_count = check_gtf_content(output_file, "repeat")
        if transcript_count > 0:
            logger.info("Dust gtf file already exists, skipping analysis")
            return 0

    logger.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, min_seq_length=5000)
    slice_ids = create_slice_ids(
        seq_region_lengths, slice_size=1_000_000, overlap=0, min_length=5000
    )

    generic_dust_cmd = [dust_path, "-in"]
    logger.info("Running Dust processes")
    pool = multiprocessing.Pool(num_threads)
    tasks = []
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_dust,
            args=(
                generic_dust_cmd,
                slice_id,
                genome_file,
                dust_output_dir,
            ),
        )

    pool.close()
    pool.join()
    slice_output_to_gtf(dust_output_dir, ".dust.gtf", 1, "repeat_id", "dust")


def multiprocess_dust(
    generic_dust_cmd,
    slice_id,
    genome_file: Union[pathlib.Path, str],
    dust_output_dir: pathlib.Path,
):
    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing slice to find low complexity regions with Dust: %s:%s:%s"
        % (region_name, start, end)
    )
    seq = get_sequence(
        seq_region=region_name,
        start=start,
        end=end,
        strand=1,
        fasta_file=genome_file,
        output_dir=dust_output_dir,
    )

    slice_file_name = f"{region_name}.rs{start}.re{end}"
    region_fasta_file_path = dust_output_dir / f"{slice_file_name}.fa"
    with open(region_fasta_file_path, "w+") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")

    region_results_file_path = dust_output_dir / f"{slice_file_name}.dust.gtf"

    dust_output_file_path = f"{region_fasta_file_path}.dust"
    dust_out = open(dust_output_file_path, "w+")
    dust_cmd = generic_dust_cmd.copy()
    dust_cmd.append(region_fasta_file_path)
    logger.info("dust_cmd" % list_to_string(dust_cmd))
    subprocess.run(dust_cmd, stdout=dust_out)
    dust_out.close()

    create_dust_gtf(dust_output_file_path, region_results_file_path, region_name)
    os.remove(dust_output_file_path)
    os.remove(region_fasta_file_path)


def create_dust_gtf(dust_output_file_path, region_results_file_path, region_name):
    with open(dust_output_file_path, "r") as dust_in, open(
        region_results_file_path, "w+"
    ) as dust_out:
        repeat_count = 1
        for line in dust_in:
            result_match = re.search(r"(\d+)\ - (\d+)", line)
            if result_match:
                start = int(result_match.group(1)) + 1
                end = int(result_match.group(2)) + 1
                gtf_line = f'{region_name}\tDust\trepeat\t{start}\t{end}\t.\t+\t.\trepeat_id "{repeat_count}";\n'
                dust_out.write(gtf_line)
                repeat_count += 1


def run_trf_repeats(
    genome_file: Union[pathlib.Path, str], trf_path, main_output_dir, num_threads: int
):
    if not trf_path:
        trf_path = "trf"

    check_exe(trf_path)
    trf_output_dir = create_dir(main_output_dir, "trf_output")

    output_file = trf_output_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "repeat")
        if transcript_count > 0:
            logger.info("Trf gtf file already exists, skipping analysis")
            return

    logger.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, min_seq_length=5000)
    slice_ids = create_slice_ids(
        seq_region_lengths, slice_size=1_000_000, overlap=0, min_length=5000
    )

    match_score = 2
    mismatch_score = 5
    delta = 7
    pm = 80
    pi = 10
    minscore = 40
    maxperiod = 500

    trf_output_extension = (
        f".{match_score}.{mismatch_score}.{delta}.{pm}.{pi}.{minscore}.{maxperiod}.dat"
    )

    generic_trf_cmd = [
        trf_path,
        None,
        str(match_score),
        str(mismatch_score),
        str(delta),
        str(pm),
        str(pi),
        str(minscore),
        str(maxperiod),
        "-d",
        "-h",
    ]
    logger.info("Running TRF processes")
    pool = multiprocessing.Pool(num_threads)
    tasks = []
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_trf,
            args=(
                generic_trf_cmd,
                slice_id,
                genome_file,
                trf_output_dir,
                trf_output_extension,
            ),
        )

    pool.close()
    pool.join()
    slice_output_to_gtf(trf_output_dir, ".trf.gtf", 1, "repeat_id", "trf")


def multiprocess_trf(
    generic_trf_cmd,
    slice_id,
    genome_file: Union[pathlib.Path, str],
    trf_output_dir: pathlib.Path,
    trf_output_extension,
):
    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing slice to find tandem repeats with TRF: %s:%s:%s"
        % (region_name, start, end)
    )
    seq = get_sequence(
        seq_region=region_name,
        start=start,
        end=end,
        strand=1,
        fasta_file=genome_file,
        output_dir=trf_output_dir,
    )

    slice_file_name = f"{region_name}.rs{start}.re{end}"
    region_fasta_file_path = trf_output_dir / f"{slice_file_name}.fa"
    with open(region_fasta_file_path, "w+") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")

    region_results_file_path = trf_output_dir / f"{slice_file_name}.trf.gtf"

    # TRF writes to the current dir, so switch to the output dir for it
    os.chdir(trf_output_dir)
    trf_output_file_path = f"{region_fasta_file_path}{trf_output_extension}"
    trf_cmd = generic_trf_cmd.copy()
    trf_cmd[1] = region_fasta_file_path
    logger.info("trf_cmd: %s" % list_to_string(trf_cmd))
    subprocess.run(trf_cmd)
    create_trf_gtf(trf_output_file_path, region_results_file_path, region_name)
    os.remove(trf_output_file_path)
    os.remove(region_fasta_file_path)


def create_trf_gtf(trf_output_file_path, region_results_file_path, region_name):
    with open(trf_output_file_path, "r") as trf_in, open(
        region_results_file_path, "w+"
    ) as trf_out:
        repeat_count = 1
        for line in trf_in:
            result_match = re.search(r"^\d+", line)
            if result_match:
                results = line.split()
                if not len(results) == 15:
                    continue
                start = results[0]
                end = results[1]
                period = float(results[2])
                copy_number = float(results[3])
                percent_matches = float(results[5])
                score = float(results[7])
                repeat_consensus = results[13]
                if (
                    score < 50
                    and percent_matches >= 80
                    and copy_number > 2
                    and period < 10
                ) or (copy_number >= 2 and percent_matches >= 70 and score >= 50):
                    gtf_line = f'{region_name}\tTRF\trepeat\t{start}\t{end}\t.\t+\t.\trepeat_id "{repeat_count}"; score "{score}"; repeat_consensus "{repeat_consensus}";\n'
                    trf_out.write(gtf_line)
                    repeat_count += 1


def run_cmsearch_regions(
    genome_file: Union[pathlib.Path, str],
    cmsearch_path,
    rfam_cm_db_path,
    rfam_seeds_file_path,
    rfam_accession_file,
    main_output_dir: pathlib.Path,
    num_threads: int,
):
    if not cmsearch_path:
        cmsearch_path = "cmsearch"

    check_exe(cmsearch_path)
    rfam_output_dir = create_dir(main_output_dir, "rfam_output")

    os.chdir(rfam_output_dir)

    output_file = rfam_output_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("CMsearch GTF file already exists, skipping analysis")
            return

    rfam_dbname = "Rfam"
    rfam_user = "rfamro"
    rfam_host = "mysql-rfam-public.ebi.ac.uk"
    rfam_port = "4497"
    # rfam_accession_query_cmd = ["mysql -h", rfam_host,"-u",rfam_user,"-P",rfam_port,"-NB -e",rfam_dbname,"'select rfam_acc FROM (SELECT DISTINCT f.rfam_acc, f.rfam_id, f.type, f.description, f.gathering_cutoff, f.trusted_cutoff FROM full_region fr, rfamseq rf, taxonomy tx, family f WHERE rf.ncbi_id = tx.ncbi_id AND f.rfam_acc = fr.rfam_acc AND fr.rfamseq_acc = rf.rfamseq_acc AND LOWER(tx.tax_string) LIKE \'%" + clade + "%\' AND (f.type LIKE \'%snRNA%\' OR f.type LIKE \'%rRNA%\' OR LOWER(f.rfam_id) LIKE \'%rnase%\' OR LOWER(f.rfam_id) LIKE \'%vault%\' OR LOWER(f.rfam_id) LIKE \'%y_rna%\' OR f.rfam_id LIKE \'%Metazoa_SRP%\') AND is_significant = 1) AS TEMP WHERE rfam_id NOT LIKE \'%bacteria%\' AND rfam_id NOT LIKE \'%archaea%\' AND rfam_id NOT LIKE \'%microsporidia%\';'"]

    # mysql -hmysql-rfam-public.ebi.ac.uk -urfamro -P4497 Rfam -NB -e "select rfam_acc FROM (SELECT DISTINCT f.rfam_acc, f.rfam_id, f.type, f.description, f.gathering_cutoff, f.trusted_cutoff FROM full_region fr, rfamseq rf, taxonomy tx, family f WHERE rf.ncbi_id = tx.ncbi_id AND f.rfam_acc = fr.rfam_acc AND fr.rfamseq_acc = rf.rfamseq_acc AND LOWER(tx.tax_string) LIKE '%insect%' AND (f.type LIKE '%snRNA%' OR f.type LIKE '%rRNA%' OR LOWER(f.rfam_id) LIKE '%rnase%' OR LOWER(f.rfam_id) LIKE '%vault%' OR LOWER(f.rfam_id) LIKE '%y_rna%' OR f.rfam_id LIKE '%Metazoa_SRP%') AND is_significant = 1) AS TEMP WHERE rfam_id NOT LIKE '%bacteria%' AND rfam_id NOT LIKE '%archaea%' AND rfam_id NOT LIKE '%microsporidia%';"

    # rfam_accession_file = '/hps/nobackup2/production/ensembl/fergal/production/test_runs/non_verts/butterfly/rfam_insect_ids.txt'
    # rfam_accession_file = main_output_dir / rfam_accessions.txt'
    rfam_cm_db_path = (
        "/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.0/Rfam.cm"
    )
    rfam_seeds_file_path = (
        "/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.0/Rfam.seed"
    )
    rfam_selected_models_file = rfam_output_dir / "rfam_models.cm"
    with open(rfam_accession_file) as rfam_accessions_in:
        rfam_accessions = rfam_accessions_in.read().splitlines()

    with open(rfam_cm_db_path, "r") as rfam_cm_in:
        rfam_data = rfam_cm_in.read()

    rfam_models = rfam_data.split("//\n")
    with open(rfam_selected_models_file, "w+") as rfam_cm_out:
        for model in rfam_models:
            # The Rfam.cm file has INFERNAL and HMMR models, both are needed at this point
            # Later we just want the INFERNAL ones for looking at thresholds
            match = re.search(r"(RF\d+)", model)
            if match:
                model_accession = match.group(1)
                if model_accession in rfam_accessions:
                    rfam_cm_out.write(f"{model}//\n")

    seed_descriptions = get_rfam_seed_descriptions(rfam_seeds_file_path)
    cv_models = extract_rfam_metrics(rfam_selected_models_file)

    logger.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, min_seq_length=5000)
    slice_ids = create_slice_ids(
        seq_region_lengths, slice_size=1_000_000, overlap=0, min_length=5000
    )

    generic_cmsearch_cmd = [
        cmsearch_path,
        "--rfam",
        "--cpu",
        "1",
        "--nohmmonly",
        "--cut_ga",
        "--tblout",
    ]
    logger.info("Running Rfam")
    pool = multiprocessing.Pool(num_threads)
    results = []
    failed_slice_ids = []
    memory_limit = 3 * 1024**3
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_cmsearch,
            args=(
                generic_cmsearch_cmd,
                slice_id,
                genome_file,
                rfam_output_dir,
                rfam_selected_models_file,
                cv_models,
                seed_descriptions,
                memory_limit,
            ),
        )
    pool.close()
    pool.join()

    # Need to figure something more automated out here. At the moment it's just limiting to 5 cores and 5GB vram
    # Ideally we could look at the amount of mem requested and put something like 10GB per core and then figure
    # out how many cores to use (obviously not using more than the amount specified)
    memory_limit = 5 * 1024**3
    if num_threads > 5:
        num_threads = 5
    pool = multiprocessing.Pool(num_threads)
    for exception_file_path in rfam_output_dir.glob("*.rfam.except"):
        logger.info("Running himem job for failed region:\n%s" % exception_file_path)
        exception_file_name = exception_file_path.name
        match = re.search(r"(.+)\.rs(\d+)\.re(\d+)\.", exception_file_name)
        if match:
            except_region = match.group(1)
            except_start = match.group(2)
            except_end = match.group(3)
            except_slice_id = [except_region, except_start, except_end]
            pool.apply_async(
                multiprocess_cmsearch,
                args=(
                    generic_cmsearch_cmd,
                    except_slice_id,
                    genome_file,
                    rfam_output_dir,
                    rfam_selected_models_file,
                    cv_models,
                    seed_descriptions,
                    memory_limit,
                ),
            )
    pool.close()
    pool.join()

    slice_output_to_gtf(rfam_output_dir, ".rfam.gtf", 1, None, None)


def multiprocess_cmsearch(
    generic_cmsearch_cmd,
    slice_id,
    genome_file: Union[pathlib.Path, str],
    rfam_output_dir: pathlib.Path,
    rfam_selected_models_file,
    cv_models: Dict[str, Dict],
    seed_descriptions: Dict,
    memory_limit: int,
):
    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing Rfam data using cmsearch against slice: %s:%s:%s"
        % (region_name, start, end)
    )
    seq = get_sequence(
        seq_region=region_name,
        start=start,
        end=end,
        strand=1,
        fasta_file=genome_file,
        output_dir=rfam_output_dir,
    )

    slice_file_name = f"{region_name}.rs{start}.re{end}"

    region_fasta_file_path = rfam_output_dir / f"{slice_file_name}.fa"
    with open(region_fasta_file_path, "w+") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")

    region_tblout_file_path = rfam_output_dir / f"{slice_file_name}.tblout"
    region_results_file_path = rfam_output_dir / f"{slice_file_name}.rfam.gtf"

    exception_results_file_path = rfam_output_dir / f"{slice_file_name}.rfam.except"

    cmsearch_cmd = generic_cmsearch_cmd.copy()
    cmsearch_cmd.extend(
        [region_tblout_file_path, rfam_selected_models_file, region_fasta_file_path]
    )

    if memory_limit is not None:
        cmsearch_cmd = prlimit_command(cmsearch_cmd, memory_limit)

    logger.info("cmsearch_cmd: %s" % list_to_string(cmsearch_cmd))

    return_value = None
    try:
        return_value = subprocess.check_output(cmsearch_cmd)
    except subprocess.CalledProcessError as ex:
        # Note that writing to file was the only option here. If return_value was passed back, eventually it would clog the
        # tiny pipe that is used by the workers to send info back. That would mean that it would eventually just be one
        # worked running at a time
        logger.error(
            "Issue processing the following region with cmsearch: %s %s-%s. Return value: %s"
            % (region_name, start, end, return_value)
        )
        with open(exception_results_file_path, "w+") as exception_out:
            exception_out.write(f"{region_name} {start} {end}\n")
        os.remove(region_fasta_file_path)
        os.remove(region_tblout_file_path)
        return

    initial_table_results = parse_rfam_tblout(region_tblout_file_path, region_name)
    unique_table_results = remove_rfam_overlap(initial_table_results)
    filtered_table_results = filter_rfam_results(unique_table_results, cv_models)
    create_rfam_gtf(
        filtered_results=filtered_table_results,
        cm_models=cv_models,
        descriptions=seed_descriptions,
        region_name=region_name,
        region_results_file_path=region_results_file_path,
        genome_file=genome_file,
        rfam_output_dir=rfam_output_dir,
    )
    os.remove(region_fasta_file_path)
    os.remove(region_tblout_file_path)
    gc.collect()


def get_rfam_seed_descriptions(rfam_seeds_path) -> Dict:
    descriptions = {}

    # NOTE: for some reason the decoder breaks on the seeds file, so I have made this ignore errors
    with open(rfam_seeds_path, encoding="utf-8", errors="ignore") as rfam_seeds_in:
        for line in rfam_seeds_in:
            seed = line.rstrip()

            domain_match = re.search("^\#=GF AC   (.+)", seed)
            if domain_match:
                domain = domain_match.group(1)
                descriptions[domain] = {}
                continue

            description_match = re.search("^\#=GF DE   (.+)", seed)
            if description_match:
                description = description_match.group(1)
                descriptions[domain]["description"] = description
                continue

            name_match = re.search("^\#=GF ID   (.+)", seed)
            if name_match:
                name = name_match.group(1)
                descriptions[domain]["name"] = name
                continue

            type_match = re.search("^\#=GF TP   Gene; (.+)", seed)
            if type_match:
                rfam_type = type_match.group(1)
                descriptions[domain]["type"] = rfam_type
                continue

    return descriptions


def extract_rfam_metrics(rfam_selected_models_file: pathlib.Path) -> Dict[str, Dict]:
    with open(rfam_selected_models_file, "r") as rfam_cm_in:
        rfam_data = rfam_cm_in.read()

    rfam_models = rfam_data.split("//\n")
    parsed_cm_data = {}
    for model in rfam_models:
        temp = model.split("\n")
        model_name_match = re.search(r"NAME\s+(\S+)", model)
        match_infernal = re.search(r"INFERNAL", model)
        if model_name_match and match_infernal:
            model_name = model_name_match.group(1)
            parsed_cm_data[model_name] = {}
            for line in temp:
                name_match = re.search(r"^NAME\s+(\S+)", line)
                if name_match:
                    parsed_cm_data[model_name]["-name"] = name_match.group(1)
                    continue

                description_match = re.search(r"^DESC\s+(\S+)", line)
                if description_match:
                    parsed_cm_data[model_name][
                        "-description"
                    ] = description_match.group(1)
                    continue

                length_match = re.search(r"^CLEN\s+(\d+)", line)
                if length_match:
                    parsed_cm_data[model_name]["-length"] = length_match.group(1)
                    continue

                max_length_match = re.search(r"^W\s+(\d+)", line)
                if max_length_match:
                    parsed_cm_data[model_name]["-maxlength"] = max_length_match.group(1)
                    continue

                accession_match = re.search(r"^ACC\s+(\S+)", line)
                if accession_match:
                    parsed_cm_data[model_name]["-accession"] = accession_match.group(1)
                    continue

                threshold_match = re.search(r"^GA\s+(\d+)", line)
                if threshold_match:
                    parsed_cm_data[model_name]["-threshold"] = threshold_match.group(1)
                    continue

    return parsed_cm_data


def parse_rfam_tblout(
    region_tblout_file_path: pathlib.Path, region_name: str
) -> List[Dict]:
    """
    # original comment bellow this function:
    # NOTE some of the code above and the code commented out here is to do with creating
    # a DAF. As we don't have a python concept of this I'm leaving it out for the moment
    # but the code below is a reference
    #    my $daf = Bio::EnsEMBL::DnaDnaAlignFeature->new(
    #      -slice          => $self->queries,
    #      -start          => $strand == 1 ? $start : $end,
    #      -end            => $strand == 1 ? $end : $start,
    #      -strand         => $strand,
    #      -hstart         => $hstart,
    #      -hend           => $hend,
    #      -hstrand        => $strand,
    #      -score          => $score,
    #      -hseqname       => length($target_name) > 39 ? substr($target_name, 0, 39) : $target_name,,
    #      -p_value  => $evalue,
    #      -align_type => 'ensembl',
    #      -cigar_string  => abs($hend - $hstart) . "M",
    #   );
    """
    parsed_results = []
    with open(region_tblout_file_path, "r") as rfam_tbl_in:
        for result in rfam_tbl_in:
            parsed_tbl_data = {}
            if not re.match(r"^" + region_name, result):
                continue

            hit = result.split()
            accession = hit[3]
            # target_name = hit[0]
            query_name = hit[2]
            # hstart = int(hit[5])
            # hend = int(hit[6])
            start = int(hit[7])
            end = int(hit[8])
            if hit[9] == "+":
                strand = 1
            else:
                strand = -1
            # evalue = float(hit[15])
            score = float(hit[14])

            parsed_tbl_data["accession"] = accession
            parsed_tbl_data["start"] = start
            parsed_tbl_data["end"] = end
            parsed_tbl_data["strand"] = strand
            parsed_tbl_data["query_name"] = query_name
            parsed_tbl_data["score"] = score
            parsed_results.append(parsed_tbl_data)

    return parsed_results


def remove_rfam_overlap(parsed_tbl_data: List[Dict]) -> List[Dict]:
    excluded_structures = {}
    chosen_structures = []
    for structure_x in parsed_tbl_data:
        chosen_structure = structure_x
        structure_x_start = structure_x["start"]
        structure_x_end = structure_x["end"]
        structure_x_score = structure_x["score"]
        structure_x_accession = structure_x["accession"]
        structure_x_string = f"{structure_x_start}:{structure_x_end}:{structure_x_score}:{structure_x_accession}"
        for structure_y in parsed_tbl_data:
            structure_y_start = structure_y["start"]
            structure_y_end = structure_y["end"]
            structure_y_score = structure_y["score"]
            structure_y_accession = structure_y["accession"]
            structure_y_string = f"{structure_y_start}:{structure_y_end}:{structure_y_score}:{structure_y_accession}"
            if structure_y_string in excluded_structures:
                continue

            if (
                structure_x_start <= structure_y_end
                and structure_x_end >= structure_y_start
            ):
                if structure_x_score < structure_y_score:
                    chosen_structure = structure_y
                    excluded_structures[structure_x_string] = 1
                else:
                    excluded_structures[structure_y_string] = 1

        chosen_structures.append(chosen_structure)

    return chosen_structures


def filter_rfam_results(
    unfiltered_tbl_data: List[Dict], cv_models: Dict[str, Dict]
) -> List[Dict]:
    filtered_results = []
    for structure in unfiltered_tbl_data:
        query = structure["query_name"]
        if query in cv_models:
            threshold = float(cv_models[query]["-threshold"])

            if query == "LSU_rRNA_eukarya":
                threshold = 1700
            elif query == "LSU_rRNA_archaea":
                continue
            elif query == "LSU_rRNA_bacteria":
                continue
            elif query == "SSU_rRNA_eukarya":
                threshold = 1600
            elif query == "5_8S_rRNA":
                threshold = 85
            elif query == "5S_rRNA":
                threshold = 75

            if threshold and structure["score"] >= threshold:
                filtered_results.append(structure)

    return filtered_results


# NOTE: The below are notes from the perl code (which has extra code) about possible improvements that are not implemented there
# Although not included in RefSeq filters, additional filters that consider sizes and score_to_size ratios can be applied
# in future work to further exclude FPs
#
# my $is_valid_size = $mapping_length > $min_length && $mapping_length < $max_length ? 1 : 0;
# my $score_size_ratio = $result->{'score'} / $mapping_length;


def create_rfam_gtf(
    filtered_results: List[Dict],
    cm_models: Dict[str, Dict],
    descriptions,
    region_name,
    region_results_file_path,
    genome_file: Union[pathlib.Path, str],
    rfam_output_dir: pathlib.Path,
):
    if not filtered_results:
        return

    with open(region_results_file_path, "w+") as rfam_gtf_out:
        gene_counter = 1
        for structure in filtered_results:
            query = structure["query_name"]
            accession = structure["accession"]
            if query in cm_models:
                model = cm_models[query]
                if accession in descriptions:
                    description = descriptions[accession]
                    if "type" in description:
                        rfam_type = description["type"]
                    else:
                        description = None
                        rfam_type = "misc_RNA"
                domain = structure["query_name"]
                padding = model["-length"]
                # rfam_type = description['type']
                gtf_strand = structure["strand"]
                rnafold_strand = structure["strand"]
                if gtf_strand == 1:
                    start = structure["start"]
                    end = structure["end"]
                    gtf_strand = "+"
                else:
                    start = structure["end"]
                    end = structure["start"]
                    score = structure["score"]
                    gtf_strand = "-"
                    rnafold_strand = -1

                biotype = "misc_RNA"
                if re.match(r"^snRNA;", rfam_type):
                    biotype = "snRNA"
                if re.match(r"^snRNA; snoRNA", rfam_type):
                    biotype = "snoRNA"
                if re.match(r"^snRNA; snoRNA; scaRNA;", rfam_type):
                    biotype = "scaRNA"
                if re.match(r"rRNA;", rfam_type):
                    biotype = "rRNA"
                if re.match(r"antisense;", rfam_type):
                    biotype = "antisense"
                if re.match(r"antitoxin;", rfam_type):
                    biotype = "antitoxin"
                if re.match(r"ribozyme;", rfam_type):
                    biotype = "ribozyme"
                if re.match(r"" + domain, rfam_type):
                    biotype = domain
                if re.match(r"" + domain, rfam_type):
                    biotype = "Vault_RNA"
                if re.match(r"" + domain, rfam_type):
                    biotype = "Y_RNA"

                rna_seq = get_sequence(
                    seq_region=region_name,
                    start=start,
                    end=end,
                    strand=rnafold_strand,
                    fasta_file=genome_file,
                    output_dir=rfam_output_dir,
                )
                valid_structure = check_rnafold_structure(rna_seq, rfam_output_dir)

                if not valid_structure:
                    continue

                transcript_string = f'{region_name}\tRfam\ttranscript\t{start}\t{end}\t.\t{gtf_strand}\t.\tgene_id "{gene_counter}"; transcript_id "{gene_counter}"; biotype "{biotype}";\n'
                exon_string = f'{region_name}\tRfam\texon\t{start}\t{end}\t.\t{gtf_strand}\t.\tgene_id "{gene_counter}"; transcript_id "{gene_counter}"; exon_number "1"; biotype "{biotype}";\n'

                rfam_gtf_out.write(transcript_string)
                rfam_gtf_out.write(exon_string)
                gene_counter += 1


def check_rnafold_structure(seq: str, rfam_output_dir: pathlib.Path):
    # Note there's some extra code in the RNAfold Perl module for encoding the structure into an attrib
    # Could consider implementing this when running for loading into an Ensembl db
    structure = 0
    with tempfile.NamedTemporaryFile(
        mode="w+t", delete=False, dir=rfam_output_dir
    ) as rna_temp_in:
        rna_temp_in.write(f">seq1\n{seq}\n")
        rna_in_file_path = rna_temp_in.name

    rnafold_cmd = ["RNAfold", "--infile", rna_in_file_path]
    rnafold_output = subprocess.Popen(rnafold_cmd, stdout=subprocess.PIPE)
    for line in io.TextIOWrapper(rnafold_output.stdout, encoding="utf-8"):
        match = re.search(r"([().]+)\s\(\s*(-*\d+.\d+)\)\n$", line)
        if match:
            structure = match.group(1)
            # score = match.group(2)
            break
    rna_temp_in.close()
    os.remove(rna_in_file_path)

    return structure


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


def run_red(
    red_path,
    genome_file: pathlib.Path,
    main_output_dir: Union[pathlib.Path, str],
) -> pathlib.Path:
    if not red_path:
        red_path = "Red"

    check_exe(red_path)
    red_dir = create_dir(main_output_dir, "red_output")
    red_mask_dir = create_dir(red_dir, "mask_output")
    red_repeat_dir = create_dir(red_dir, "repeat_output")
    red_genome_dir = create_dir(red_dir, "genome_dir")

    genome_file_name = genome_file.name
    red_genome_file = red_genome_dir / genome_file_name
    genome_file_stem = genome_file.stem
    masked_genome_file = red_mask_dir / f"{genome_file_stem}.msk"
    repeat_coords_file = red_repeat_dir / f"{genome_file_stem}.rpt"
    gtf_output_file_path = red_dir / "annotation.gtf"

    if masked_genome_file.exists():
        logger.warning(
            "Masked Genome file already found on the path to the Red mask output dir. Will not create a new file"
        )
        create_red_gtf(repeat_coords_file, gtf_output_file_path)
        return masked_genome_file

    if red_genome_file.exists():
        logger.warning(
            "Unmasked genome file already found on the path to the Red genome dir, will not create a sym link"
        )
    else:
        logger.info(
            "Creating sym link of the genome file to the Red genome dir:\n%s\nto\n%s"
            % (genome_file, red_genome_file)
        )
        red_genome_file.symlink_to(genome_file)

    if not red_genome_file.exists():
        logger.error(
            "Could not find the genome file or sym link to the original file in the Red genome dir at:\n%s"
            % red_genome_file
        )

    logger.info("Running Red, this may take some time depending on the genome size")
    subprocess.run(
        [red_path, "-gnm", red_genome_dir, "-msk", red_mask_dir, "-rpt", red_repeat_dir]
    )

    logger.info("Completed running Red")

    create_red_gtf(repeat_coords_file, gtf_output_file_path)

    return masked_genome_file


def create_red_gtf(repeat_coords_file, gtf_output_file_path):
    with open(repeat_coords_file, "r") as red_in, open(
        gtf_output_file_path, "w+"
    ) as red_out:
        for repeat_id, line in enumerate(red_in, start=1):
            result_match = re.search(r"^\>(.+)\:(\d+)\-(\d+)", line)
            if result_match:
                region_name = result_match.group(1)
                # Note that Red is 0-based, so add 1
                start = int(result_match.group(2)) + 1
                end = int(result_match.group(3)) + 1
                gtf_line = f'{region_name}\tRed\trepeat\t{start}\t{end}\t.\t+\t.\trepeat_id "{repeat_id}";\n'
                red_out.write(gtf_line)


def run_genblast_align(
    genblast_path,
    convert2blastmask_path,
    makeblastdb_path,
    genblast_dir: pathlib.Path,
    protein_file,
    masked_genome_file: pathlib.Path,
    max_intron_length,
    num_threads: int,
    genblast_timeout_secs,
    main_script_dir: pathlib.Path,
):
    if not genblast_path:
        genblast_path = "genblast"

    check_exe(genblast_path)

    if not convert2blastmask_path:
        convert2blastmask_path = "convert2blastmask"

    check_exe(convert2blastmask_path)

    if not makeblastdb_path:
        makeblastdb_path = "makeblastdb"

    check_exe(makeblastdb_path)

    genblast_dir = create_dir(genblast_dir)

    output_file = genblast_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Genblast GTF file exists, skipping analysis")
            return

    # genblast_output_file = genblast_dir / "genblast"

    asnb_file = f"{masked_genome_file}.asnb"
    logger.info("ASNB file: %s" % asnb_file)

    if not os.path.exists("alignscore.txt"):
        original_alignscore_path = main_script_dir / "support_files" / "alignscore.txt"
        subprocess.run(["cp", original_alignscore_path, "./"])

    if not masked_genome_file.exists():
        raise FileNotFoundError(
            "Masked genome file does not exist: %s" % masked_genome_file
        )

    if not os.path.exists(protein_file):
        raise FileNotFoundError("Protein file does not exist: %s" % protein_file)

    if not os.path.exists(asnb_file):
        run_convert2blastmask(convert2blastmask_path, masked_genome_file, asnb_file)
    else:
        logger.info("Found an existing asnb, so will skip convert2blastmask")

    if not os.path.exists(asnb_file):
        raise FileNotFoundError("asnb file does not exist: %s" % asnb_file)

    run_makeblastdb(makeblastdb_path, masked_genome_file, asnb_file)

    batched_protein_files = split_protein_file(
        protein_file=protein_file, protein_output_dir=genblast_dir, batch_size=20
    )

    pool = multiprocessing.Pool(num_threads)
    for batched_protein_file in batched_protein_files:
        pool.apply_async(
            multiprocess_genblast,
            args=(
                batched_protein_file,
                masked_genome_file,
                genblast_path,
                genblast_timeout_secs,
                max_intron_length,
            ),
        )
    pool.close()
    pool.join()

    logger.info("Completed running GenBlast")
    logger.info("Combining output into single GTF")
    generate_genblast_gtf(genblast_dir)


def multiprocess_genblast(
    batched_protein_file: pathlib.Path,
    masked_genome_file,
    genblast_path,
    genblast_timeout_secs,
    max_intron_length,
):
    # batch_num = os.path.splitext(batched_protein_file)[0]
    # batch_dir = os.path.dirname(batched_protein_file)
    logger.info("Running GenBlast on %s:" % batched_protein_file)

    genblast_cmd = [
        genblast_path,
        "-p",
        "genblastg",
        "-q",
        batched_protein_file,
        "-t",
        masked_genome_file,
        "-g",
        "T",
        "-pid",
        "-r",
        "1",
        "-P",
        "blast",
        "-gff",
        "-e",
        "1e-1",
        "-c",
        "0.8",
        "-W",
        "3",
        "-softmask",
        "-scodon",
        "50",
        "-i",
        "30",
        "-x",
        "10",
        "-n",
        "30",
        "-d",
        str(max_intron_length),
        "-o",
        batched_protein_file,
    ]

    logger.info("genblast_cmd: %s" % list_to_string(genblast_cmd))
    # Using the child process termination as described here:
    # https://alexandra-zaharia.github.io/posts/kill-subprocess-and-its-children-on-timeout-python/
    try:
        p = subprocess.Popen(genblast_cmd, start_new_session=True)
        p.wait(timeout=genblast_timeout_secs)
    except subprocess.TimeoutExpired:
        logger.error("Timeout reached for file:\n%s" % batched_protein_file)
        subprocess.run(["touch", f"{batched_protein_file}.except"])
        os.killpg(os.getpgid(p.pid), signal.SIGTERM)

    files_to_delete = glob.glob(f"{batched_protein_file}*msk.blast*")
    files_to_delete.append(batched_protein_file)
    for file_to_delete in files_to_delete:
        subprocess.run(["rm", file_to_delete])


def generate_genblast_gtf(genblast_dir: pathlib.Path):
    genblast_extension = "_1.1c_2.3_s1_0_16_1"
    file_out_name = genblast_dir / "annotation.gtf"
    with open(file_out_name, "w+") as file_out:
        for path in genblast_dir.glob("**/*"):
            if path.suffix == ".gff":
                gtf_string = convert_gff_to_gtf(path)
                file_out.write(gtf_string)
            elif (
                str(path).endswith(".fa.blast")
                or str(path).endswith(".fa.blast.report")
                or str(path).endswith(genblast_extension)
            ):
                subprocess.run(["rm", path])


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


def split_protein_file(
    protein_file, protein_output_dir: pathlib.Path, batch_size: int = 20
) -> List[pathlib.Path]:
    for i in range(0, 10):
        create_dir(protein_output_dir, f"bin_{i}")

    batched_protein_files = []
    with open(protein_file) as file_in:
        seq_count = 0
        batch_count = 0
        current_record = ""
        initial_seq = 1
        for line in file_in:
            match = re.search(r">(.+)$", line)
            if match and not initial_seq and seq_count % batch_size == 0:
                num_dir = random.randint(0, 9)
                file_out_name = (
                    protein_output_dir / f"bin_{num_dir}" / f"{batch_count}.fa"
                )
                with open(file_out_name, "w+") as file_out:
                    file_out.write(current_record)

                batch_count += 1
                seq_count += 1
                current_record = line
                batched_protein_files.append(file_out_name)
            elif match:
                current_record += line
                initial_seq = 0
                seq_count += 1
            else:
                current_record += line

    if current_record:
        num_dir = random.randint(0, 9)
        file_out_name = protein_output_dir / f"bin_{num_dir}" / f"{batch_count}.fa"
        with open(file_out_name, "w+") as file_out:
            file_out.write(current_record)
        batched_protein_files.append(file_out_name)

    return batched_protein_files


def run_convert2blastmask(convert2blastmask_path, masked_genome_file, asnb_file):
    asnb_file = f"{masked_genome_file}.asnb"
    cmd = [
        convert2blastmask_path,
        "-in",
        masked_genome_file,
        "-parse_seqids",
        "-masking_algorithm",
        "other",
        "-masking_options",
        '"REpeatDetector, default"',
        "-outfmt",
        "maskinfo_asn1_bin",
        "-out",
        asnb_file,
    ]
    logger.info(
        "Running convert2blastmask prior to GenBlast:\n%s" % list_to_string(cmd)
    )
    subprocess.run(cmd)
    logger.info("Completed running convert2blastmask")


def run_makeblastdb(makeblastdb_path, masked_genome_file, asnb_file):
    cmd = [
        makeblastdb_path,
        "-in",
        masked_genome_file,
        "-dbtype",
        "nucl",
        "-parse_seqids",
        "-mask_data",
        asnb_file,
        "-max_file_sz",
        "10000000000",
    ]
    logger.info("Running makeblastdb prior to GenBlast:\n%s" % list_to_string(cmd))
    subprocess.run(cmd)
    logger.info("Completed running makeblastdb")


def run_trimming(
    main_output_dir: Union[pathlib.Path, str],
    short_read_fastq_dir,
    delete_pre_trim_fastq,
    num_threads: int,
):
    trim_galore_path = "trim_galore"
    check_exe(trim_galore_path)

    # TODO
    # update type upstream
    short_read_fastq_dir = pathlib.Path(short_read_fastq_dir)

    trim_dir = create_dir(main_output_dir, "trim_galore_output")

    file_types = ["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"]
    fastq_file_list = []
    for file_type in file_types:
        fastq_file_list.extend(short_read_fastq_dir.glob(file_type))

    fastq_paired_paths = create_paired_paths(fastq_file_list)

    for fastq_paired_files in fastq_paired_paths:
        logger.info(
            "FASTQ file path(s): %s" % list_to_string(fastq_paired_files, separator=",")
        )

    generic_trim_galore_cmd = [
        trim_galore_path,
        "--illumina",
        "--quality",
        "20",
        "--length",
        "50",
        "--output_dir",
        trim_dir,
    ]

    pool = multiprocessing.Pool(num_threads)
    for fastq_paired_files in fastq_paired_paths:
        pool.apply_async(
            multiprocess_trim_galore,
            args=(
                generic_trim_galore_cmd,
                fastq_paired_files,
                trim_dir,
            ),
        )
        if delete_pre_trim_fastq:
            for fastq_path in fastq_paired_files:
                logger.info(
                    "Removing original fastq file post trimming:\n%s" % fastq_path
                )
                subprocess.run(["rm", fastq_path])

    pool.close()
    pool.join()

    for trimmed_fastq_path in trim_dir.glob("*.fq.gz"):
        logger.info("Trimmed file path:\n%s" % trimmed_fastq_path)
        sub_patterns = r"|".join(("_val_1.fq", "_val_2.fq", "_trimmed.fq"))
        updated_file_path = re.sub(sub_patterns, ".fq", str(trimmed_fastq_path))
        logger.info("Updated file path:\n%s" % updated_file_path)
        subprocess.run(["mv", trimmed_fastq_path, updated_file_path])

        files_to_delete_list = []
        for file_type in file_types:
            files_to_delete_list.extend(short_read_fastq_dir.glob(file_type))


def multiprocess_trim_galore(
    generic_trim_galore_cmd,
    fastq_paired_files: List[pathlib.Path],
    trim_dir: pathlib.Path,
):
    fastq_file = fastq_files[0]

    if len(fastq_files) == 2:
        fastq_file_pair = fastq_files[1]
    # len(fastq_files) == 1
    else:
        fastq_file_pair = None

    trim_galore_cmd = generic_trim_galore_cmd
    if fastq_file_pair:
        trim_galore_cmd.append("--paired")

    trim_galore_cmd.append(fastq_file)

    if fastq_file_pair:
        trim_galore_cmd.append(fastq_file_pair)

    logger.info(
        "Running Trim Galore with the following command:\n%s"
        % list_to_string(trim_galore_cmd)
    )
    subprocess.run(trim_galore_cmd)


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
    main_script_dir: pathlib.Path,
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

    if (
        subsample_script_path is None
        or not pathlib.Path(subsample_script_path).exists()
    ):
        subsample_script_path = (
            main_script_dir / "support_scripts" / "subsample_fastq.py"
        )

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
    fastq_paired_paths = create_paired_paths(fastq_file_list)

    # Subsamples in parallel if there's a value set
    if max_reads_per_sample:
        pool = multiprocessing.Pool(num_threads)
        for fastq_files in fastq_paired_paths:
            fastq_file = fastq_files[0]

            if len(fastq_files) == 2:
                fastq_file_pair = fastq_files[1]
            # len(fastq_files) == 1
            else:
                fastq_file_pair = None

            fastq_file_subsampled = fastq_file.parent / f"{fastq_file.name}.sub"
            fastq_file_pair_subsampled = (
                fastq_file_pair.parent / f"{fastq_file_pair.name}.sub"
            )
            if (
                fastq_file_pair
                and fastq_file_subsampled.exists()
                and fastq_file_pair_subsampled.exists()
            ):
                logger.info(
                    "Found existing .sub files on the fastq path for both members of the pair, will use those instead of subsampling again:\n%s\n%s"
                    % (fastq_file_subsampled, fastq_file_pair_subsampled)
                )
            elif fastq_file_pair:
                pool.apply_async(
                    run_subsample_script,
                    args=(
                        subsample_script_path,
                        fastq_file,
                        fastq_file_pair,
                    ),
                )
            elif fastq_file_subsampled.exists():
                logger.info(
                    "Found an existing .sub file on the fastq path, will use that instead:\n%s"
                    % f"{fastq_file}.sub"
                )
            else:
                pool.apply_async(
                    run_subsample_script,
                    args=(
                        subsample_script_path,
                        fastq_file,
                        fastq_file_pair,
                    ),
                )

        pool.close()
        pool.join()

    fastq_paired_paths = check_for_fastq_subsamples(fastq_paired_paths)

    if not fastq_paired_paths:
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
    for fastq_files in fastq_paired_paths:
        # NOTE
        # There might be a bug here. The paired FASTQ files elements of the output list
        # of check_for_fastq_subsamples were strings containing paths separated by comma,
        # but the following code originally was handled these elements as if they were
        # a single path. The next command was added to select the first file from the pair,
        # but maybe there is something more to be done here.
        fastq_file = fastq_files[0]
        logger.info("fastq_file: %s" % fastq_file)
        fastq_file_name = fastq_file.name
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

        logger.info("Processing %s" % fastq_file)
        # star_command = [star_path,'--outFilterIntronMotifs','RemoveNoncanonicalUnannotated','--outSAMstrandField','intronMotif','--runThreadN',str(num_threads),'--twopassMode','Basic','--runMode','alignReads','--genomeDir',star_dir,'--readFilesIn',fastq_file,'--outFileNamePrefix',(star_dir + '/'),'--outTmpDir',star_tmp_dir,'--outSAMtype','SAM','--alignIntronMax',str(max_intron_length),'--outSJfilterIntronMaxVsReadN','5000','10000','25000','40000','50000','50000','50000','50000','50000','100000']

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
            fastq_file,
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


def run_subsample_script(subsample_script_path, fastq_file, fastq_file_pair):
    subsample_script_cmd = [
        "python3",
        subsample_script_path,
        "--fastq_file",
        fastq_file,
    ]
    if fastq_file_pair:
        subsample_script_cmd.extend(["--fastq_file_pair", fastq_file_pair])
    subprocess.run(subsample_script_cmd)


def check_for_fastq_subsamples(
    fastq_paired_paths: List[List[pathlib.Path]],
) -> List[List[pathlib.Path]]:
    """
    Replace FASTQ files with their corresponding subsampled versions if they exist.

    Original comment:
    This should probably removed at some point as it is needlessly complicated
    Would be better to just build into the previous step
    Mainly just about making sure that if you have subsamples they're used and if you have pairs they're paired

    Args:
        fastq_paired_paths: list of lists containing the FASTQ paired file paths
    Returns:
        the input list with FASTQ file paths replaced with their corresponding
        subsampled versions
    """
    for index, fastq_files in enumerate(fastq_paired_paths):
        fastq_file = fastq_files[0]
        fastq_file_subsampled = fastq_file.parent / f"{fastq_file.name}.sub"

        if len(fastq_files) == 2:
            fastq_file_pair = fastq_files[1]
            fastq_file_pair_subsampled = (
                fastq_file_pair.parent / f"{fastq_file_pair.name}.sub"
            )
        # len(fastq_files) == 1
        else:
            fastq_file_pair = None
            fastq_file_pair_subsampled = None

        if fastq_file_subsampled.exists():
            logger.info(
                "Found a subsampled file extension, will use that instead of the original file:\n%s"
                % fastq_file_subsampled
            )
            fastq_paired_paths[index] = [fastq_file_subsampled]
        else:
            fastq_paired_paths[index] = [fastq_file]

        if fastq_file_pair_subsampled.exists():
            logger.info(
                "Found a subsampled paired file extension, will use that instead of the original file:\n%s"
                % fastq_file_pair_subsampled
            )
            fastq_paired_paths[index].append(fastq_file_pair_subsampled)
        elif fastq_file_pair:
            fastq_paired_paths[index].append(fastq_file_pair)

        logger.info("Entry at current index:\n%s" % fastq_paired_paths[index])

    return fastq_paired_paths


def run_minimap2_align(
    minimap2_path,
    paftools_path,
    main_output_dir,
    long_read_fastq_dir,
    genome_file: pathlib.Path,
    max_intron_length,
    num_threads: int,
):
    if not minimap2_path:
        minimap2_path = "minimap2"

    check_exe(minimap2_path)

    if not paftools_path:
        paftools_path = "paftools.js"

    check_exe(paftools_path)

    minimap2_output_dir = create_dir(main_output_dir, "minimap2_output")

    output_file = minimap2_output_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Minimap2 GTF file exists skipping analysis")
            return

    genome_file_name = genome_file.name
    genome_file_index = f"{genome_file_name}.mmi"
    minimap2_index_file = minimap2_output_dir / genome_file_index
    minimap2_hints_file = minimap2_output_dir / "minimap2_hints.gff"

    fastq_file_list = []
    for fastq_file in glob.glob(f"{long_read_fastq_dir}/*.fastq"):
        fastq_file_list.append(fastq_file)

    for fastq_file in glob.glob(f"{long_read_fastq_dir}/*.fq"):
        fastq_file_list.append(fastq_file)

    if not fastq_file_list:
        # NOTE: should update this to have a param that says it's okay to be empty if there's an override param
        # This is because it's important that an external user would have an exception if they provided a long read dir
        # but there was nothing in it. Whereas in the Ensembl pipeline there might not be long read data, but we don't
        # know this when the command line is being constructed. This is also true for the short read data. An alternative
        # would be to put in another analysis into the Ensembl pipeline to construct this commandline after the transcriptomic
        # data has been searched for
        return
        # raise IndexError('The list of fastq files is empty. Fastq dir:\n%s' % long_read_fastq_dir)

    if not minimap2_index_file.exists():
        logger.info("Did not find an index file for minimap2. Will create now")
        subprocess.run(
            [
                minimap2_path,
                "-t",
                str(num_threads),
                "-d",
                minimap2_index_file,
                genome_file,
            ]
        )

    if not minimap2_index_file.exists():
        raise FileNotFoundError(
            "Failed to create the minimap2 index file at:\n%s" % minimap2_index_file
        )

    logger.info("Running minimap2 on the files in the long read fastq dir")
    for fastq_file_path in fastq_file_list:
        fastq_file_name = fastq_file_path.name
        sam_file = minimap2_output_dir / f"{fastq_file_name}.sam"
        bed_file = minimap2_output_dir / f"{fastq_file_name}.bed"
        bed_file_out = open(bed_file, "w+")
        logger.info("Processing %s" % fastq_file)
        subprocess.run(
            [
                minimap2_path,
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
                fastq_file_path,
                "-o",
                sam_file,
            ]
        )
        logger.info("Creating bed file from SAM")
        subprocess.run([paftools_path, "splice2bed", sam_file], stdout=bed_file_out)
        bed_file_out.close()

    bed_to_gtf(minimap2_output_dir)

    logger.info("Completed running minimap2")


def bed_to_gtf(minimap2_output_dir: pathlib.Path):
    gtf_file_path = minimap2_output_dir / "annotation.gtf"
    with open(gtf_file_path, "w+") as gtf_out:
        gene_id = 1
        for bed_file in minimap2_output_dir.glob("*.bed"):
            logger.info("Converting bed to GTF:\n%s" % bed_file)
            with open(bed_file) as bed_in:
                for line in bed_in:
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
                    exons = bed_to_exons(block_sizes, block_starts, offset)
                    transcript_line = [
                        seq_region_name,
                        "minimap",
                        "transcript",
                        0,
                        0,
                        ".",
                        strand,
                        ".",
                        f'gene_id "minimap_{gene_id}"; transcript_id "minimap_{gene_id}"',
                    ]
                    transcript_start = None
                    transcript_end = None
                    exon_records = []
                    for index, exon_coords in enumerate(exons, start=1):
                        if (
                            transcript_start is None
                            or exon_coords[0] < transcript_start
                        ):
                            transcript_start = exon_coords[0]

                        if transcript_end is None or exon_coords[1] > transcript_end:
                            transcript_end = exon_coords[1]

                        exon_line = [
                            seq_region_name,
                            "minimap",
                            "exon",
                            str(exon_coords[0]),
                            str(exon_coords[1]),
                            ".",
                            strand,
                            ".",
                            f'gene_id "minimap_{gene_id}"; transcript_id "minimap_{gene_id}"; exon_number "{index}";',
                        ]

                        exon_records.append(exon_line)

                    transcript_line[3] = str(transcript_start)
                    transcript_line[4] = str(transcript_end)

                    gtf_out.write("\t".join(transcript_line) + "\n")
                    for exon_line in exon_records:
                        gtf_out.write("\t".join(exon_line) + "\n")

                    gene_id += 1


def bed_to_exons(block_sizes, block_starts, offset):
    exons = []
    for i, element in enumerate(block_sizes):
        block_start = offset + int(block_starts[i]) + 1
        block_end = block_start + int(block_sizes[i]) - 1

        if block_end < block_start:
            logger.warning("Warning: block end is less than block start, skipping exon")
            continue

        exon_coords = [str(block_start), str(block_end)]
        exons.append(exon_coords)

    return exons


def check_transcriptomic_output(main_output_dir: pathlib.Path):
    # This will check across the various transcriptomic dirs and check there's actually some data
    transcriptomic_dirs = ["scallop_output", "stringtie_output", "minimap2_output"]
    total_lines = 0
    min_lines = 100_000
    for transcriptomic_dir in transcriptomic_dirs:
        full_file_path = main_output_dir / transcriptomic_dir / "annotation.gtf"
        if not full_file_path.exists():
            logger.warning(
                'Warning, no annotation.gtf found for "%s". This might be fine, e.g. no long read data were provided'
                % transcriptomic_dir
            )
            continue
        num_lines = sum(1 for line in open(full_file_path))
        total_lines += num_lines
        logger.info(
            'For "%s" found a total of %s in the annotation.gtf file'
            % (transcriptomic_dir, num_lines)
        )
    if total_lines == 0:
        raise IOError(
            "Anno was run with transcriptomic mode enabled, but the transcriptomic annotation files are empty"
        )
    elif total_lines <= min_lines:
        raise IOError(
            "Anno was run with transcriptomic mode enabled, but the total number of lines in the output files were less than the min expected value\nFound: %s\nMin allowed: %s"
            % (total_lines, min_lines)
        )
    else:
        logger.info(
            "Found %s total lines across the transcriptomic files. Checks passed"
            % total_lines
        )


# start gene g1
# 1       AUGUSTUS        gene    1       33908   1       +       .       g1
# 1       AUGUSTUS        transcript      1       33908   .       +       .       g1.t1
# 1       AUGUSTUS        CDS     3291    3585    .       +       2       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        exon    3291    3585    .       +       .       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        CDS     11377   11510   .       +       1       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        exon    11377   11510   .       +       .       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        CDS     30726   30871   .       +       2       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        exon    30726   30871   .       +       .       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        CDS     32975   33502   .       +       0       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        exon    32975   33908   .       +       .       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        stop_codon      33500   33502   .       +       0       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        tts     33908   33908   .       +       .       transcript_id "g1.t1"; gene_id "g1";
# protein sequence = [GGRGEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEKKEKKEEEEEEEEEEEEEEEK
# EEGEGMEGRDIMGIRGKQKQPRKQPELSSSFPTFLLTQKTPYCPESIKLLLILRNNIQTIVFYKGPKGVINDWRKFKLESEDGDSIPPSKKEILRQMS
# SPQSRDDSKERMSRKMSIQEYELIHQDKEDESCLRKYRRQCMQDMHQKLSFGPRYGFVYELETGEQFLETIEKEQKVTTIVVNIYEDGVRGCDALNSS
# LACLAVEYPMVKFCKIKASNTGAGDRFSTDVLPTLLVYKGGELISNFISVAEQFAEEFFAVDVESFLNEYGLLPEREIHDLEQTNMEDEDIE]
# Evidence for and against this transcript:
# % of transcript supported by hints (any source): 0
# CDS exons: 0/4
# CDS introns: 0/4
# 5'UTR exons and introns: 0/0
# 3'UTR exons and introns: 0/1
# hint groups fully obeyed: 0
# incompatible hint groups: 0
# end gene g1
###


def augustus_output_to_gtf(
    augustus_output_dir: pathlib.Path,
    augustus_genome_dir: pathlib.Path,
):
    gtf_file_path = augustus_output_dir / "annotation.gtf"
    with open(gtf_file_path, "w+") as gtf_out:
        record_count = 1
        for gff_file_path in augustus_genome_dir.glob("*.aug"):
            gff_file_name = gff_file_path.name
            match = re.search(r"\.rs(\d+)\.re(\d+)\.", gff_file_name)
            start_offset = int(match.group(1))

            exon_number = 1
            current_exon_hints_total = 0
            current_exon_hints_match = 0
            current_intron_hints_total = 0
            current_intron_hints_match = 0
            current_record = []
            with open(gff_file_path, "r") as gff_in:
                for line in gff_in:
                    match = re.search(r"# CDS exons\: (\d+)\/(\d+)", line)
                    if match:
                        current_exon_hints_match = match.group(1)
                        current_exon_hints_total = match.group(2)

                    match = re.search(r"# CDS introns\: (\d+)\/(\d+)", line)
                    if match:
                        current_introns_hints_match = match.group(1)
                        current_introns_hints_total = match.group(2)

                    if re.search(r"# end gene", line):
                        for output_line in current_record:
                            gtf_out.write(output_line)
                        current_record = []
                        current_exon_hints_total = 0
                        current_exon_hints_match = 0
                        current_intron_hints_total = 0
                        current_intron_hints_match = 0
                        record_count += 1
                        exon_number = 1

                    values = line.split("\t")
                    if (
                        len(values) == 9
                        and values[1] == "AUGUSTUS"
                        and (values[2] == "transcript" or values[2] == "exon")
                    ):
                        values[3] = str(int(values[3]) + (start_offset - 1))
                        values[4] = str(int(values[4]) + (start_offset - 1))
                        # fmt: off
                        values[8] = f'gene_id "aug{record_count}"; transcript_id "aug{record_count}";'
                        # fmt: on
                        if values[2] == "exon":
                            values[8] += f' exon_number "{exon_number}";'
                            exon_number += 1
                        values[8] += "\n"
                        current_record.append("\t".join(values))


def run_augustus_predict(
    augustus_path,
    main_output_dir: pathlib.Path,
    masked_genome_file,
    num_threads: int,
):
    min_seq_length = 1000

    if not augustus_path:
        augustus_path = (
            "/hps/nobackup2/production/ensembl/jma/src/Augustus/bin/augustus"
        )
    check_exe(augustus_path)

    bam2hints_path = "/homes/fergal/bin/bam2hints"
    bam2wig_path = "/homes/fergal/bin/bam2wig"
    wig2hints_path = "/homes/fergal/bin/wig2hints"

    # Run bam2hints, bam2wig, wig2hints, then combine the hints into a single file
    # Multiprocess with all three steps in the MP as that would be fastest

    augustus_dir = create_dir(main_output_dir, "augustus_output")
    augustus_hints_dir = create_dir(augustus_dir, "hints")
    augustus_genome_dir = create_dir(augustus_dir, "genome_dir")
    augustus_evidence_dir = create_dir(augustus_dir, "evidence")
    augustus_hints_file = augustus_evidence_dir / "augustus_hints.gff"
    star_dir = main_output_dir / "star_output"
    # minimap2_output_dir = main_output_dir / "minimap2_output"

    if star_dir.exists():
        logger.info("Found a STAR output dir, generating hints from any .sj.tab files")
        generate_hints(
            bam2hints_path,
            bam2wig_path,
            wig2hints_path,
            augustus_hints_dir,
            star_dir,
            num_threads,
        )
        with open(augustus_hints_file, "w+") as hints_out:
            for gff_file in augustus_hints_dir.glob("*.bam.hints.gff"):
                with open(gff_file, "r") as gff_in:
                    for line in gff_in:
                        hints_out.write(line)

    seq_region_lengths = get_seq_region_lengths(genome_file, min_seq_length=5000)
    slice_ids = create_slice_ids(
        seq_region_lengths, slice_size=1_000_000, overlap=100_000, min_length=5000
    )

    generic_augustus_cmd = [
        augustus_path,
        "--species=human",
        "--UTR=on",
        (
            "--extrinsicCfgFile="
            + "/hps/nobackup2/production/ensembl/jma/src/Augustus/config/extrinsic/extrinsic.M.RM.E.W.P.cfg"
        ),
    ]

    pool = multiprocessing.Pool(num_threads)
    tasks = []

    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_augustus_id,
            args=(
                generic_augustus_cmd,
                slice_id,
                masked_genome_file,
                augustus_hints_file,
                augustus_genome_dir,
            ),
        )

    pool.close()
    pool.join()

    augustus_output_to_gtf(augustus_dir, augustus_genome_dir)


def create_slice_ids(
    seq_region_lengths,
    slice_size: int = 1_000_000,
    overlap: int = 0,
    min_length: int = 0,
) -> List[List]:
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


def generate_hints(
    bam2hints_path,
    bam2wig_path,
    wig2hints_path,
    augustus_hints_dir: Union[pathlib.Path, str],
    star_dir: pathlib.Path,
    num_threads: int,
):
    pool = multiprocessing.Pool(num_threads)
    for bam_file in star_dir.glob("*.bam"):
        pool.apply_async(
            multiprocess_augustus_hints,
            args=(
                bam2hints_path,
                bam2wig_path,
                wig2hints_path,
                bam_file,
                augustus_hints_dir,
            ),
        )
    pool.close()
    pool.join()


def multiprocess_augustus_hints(
    bam2hints_path,
    bam2wig_path,
    wig2hints_path,
    bam_file: pathlib.Path,
    augustus_hints_dir: pathlib.Path,
):
    bam_file_name = bam_file.name
    logger.info("Processing %s for Augustus hints" % bam_file_name)

    bam2hints_file_name = f"{bam_file_name}.hints.gff"
    bam2hints_file_path = augustus_hints_dir / bam2hints_file_name
    bam2hints_cmd = [
        bam2hints_path,
        f"--in={bam_file}",
        f"--out={bam2hints_file_path}",
        "--maxintronlen=100000",
    ]
    logger.info("bam2hints command:\n%s" % list_to_string(bam2hints_cmd))
    subprocess.run(bam2hints_cmd)

    # bam2wig_cmd = [bam2wig_path, "-D", augustus_hints_dir, bam_file]
    # logger.info("bam2wig command:\n%s" % list_to_string(bam2wig_cmd))
    # subprocess.run(bam2wig_cmd)

    # wig2hints is odd in that it runs directly off STDIN and then just prints to STDOUT,
    # so the code below is implemented in steps as it's not good practice to use pipes and
    # redirects in a subprocess command
    # wig_file_name = re.sub(".bam", ".wig", bam_file_name)
    # wig_file_path = augustus_hints_dir / wig_file_name
    # wig_hints_file_name = f"{wig_file_name}.hints.gff"
    # wig_hints_file_path = augustus_hints_dir / wig_hints_file_name
    # logger.info("Writing wig file info to hints file:\n%s" % wig_hints_file_name)
    # wig2hints_out = open(wig_hints_file_path, "w+")
    # wigcat = subprocess.Popen(("cat", wig_file_path), stdout=subprocess.PIPE)
    # subprocess.run(wig2hints_path, stdin=wigcat.stdout, stdout=wig2hints_out)
    # wig2hints_out.close()


def multiprocess_augustus_id(
    cmd,
    slice_id,
    genome_file: Union[pathlib.Path, str],
    hints_file,
    output_dir: pathlib.Path,
):
    region = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]
    seq = get_sequence(
        seq_region=region,
        start=start,
        end=end,
        strand=1,
        fasta_file=genome_file,
        output_dir=output_dir,
    )

    region_fasta_file_name = f"{region}.rs{start}.re{end}.fa"
    region_fasta_file_path = output_dir / region_fasta_file_name
    region_augustus_file_path = output_dir / f"{region_fasta_file_name}.aug"

    with open(region_fasta_file_path, "w+") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")

    region_hints_file = create_slice_hints_file(
        region, start, end, hints_file, region_fasta_file_path
    )

    aug_out = open(region_augustus_file_path, "w+")

    augustus_forward = cmd.copy()
    augustus_forward.extend(
        [f"--hintsfile={region_hints_file}", "--strand=forward", region_fasta_file_path]
    )
    subprocess.run(augustus_forward, stdout=aug_out)

    augustus_backward = cmd.copy()
    augustus_forward.extend(
        [
            f"--hintsfile={region_hints_file}",
            "--strand=backward",
            region_fasta_file_path,
        ]
    )
    subprocess.run(augustus_backward, stdout=aug_out)

    aug_out.close()


def create_slice_hints_file(region, start, end, hints_file, region_fasta_file_path):
    """
    Note this is trying to be memory and file efficient at the cost of speed
    So files are only created as needed and the hints are being read line by line as written as needed
    This comes with the downside of being slow, but it's only a very small amount of time relative
    to how slow the step is in total. Given that this step in general eats up a low of memory, saving as much
    as possible here is not a bad thing even if it's adding in an overhead by continuously reading the hints file
    """
    region_hints_file_path = f"{region_fasta_file_path}.gff"
    with open(hints_file) as hints_in, open(region_hints_file_path, "w+") as hints_out:
        for hint_line in hints_in:
            hint_line_values = hint_line.split("\t")
            if not len(hint_line_values) == 9:
                continue

            hint_region = hint_line_values[0]
            hint_region_start = int(hint_line_values[3])
            hint_region_end = int(hint_line_values[4])

            if (
                hint_region == region
                and hint_region_start >= start
                and hint_region_end <= end
            ):
                hint_line_values[3] = str(int(hint_line_values[3]) - (start - 1))
                hint_line_values[4] = str(int(hint_line_values[4]) - (start - 1))
                hints_out.write("\t".join(hint_line_values))

    return region_hints_file_path


def run_stringtie_assemble(
    stringtie_path, samtools_path, main_output_dir, num_threads: int
):
    if not stringtie_path:
        stringtie_path = shutil.which("stringtie")
    check_exe(stringtie_path)

    if not samtools_path:
        samtools_path = shutil.which("samtools")
    check_exe(samtools_path)

    stringtie_dir = create_dir(main_output_dir, "stringtie_output")

    output_file = stringtie_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("StringTie GTF file exists skipping analysis")
            return

    stringtie_merge_input_file = stringtie_dir / "stringtie_assemblies.txt"
    stringtie_merge_output_file = stringtie_dir / "annotation.gtf"
    star_dir = main_output_dir / "star_output"

    if star_dir.exists():
        logger.info("Found a STAR output dir, will load sam file")

    sorted_bam_files = list(star_dir.glob("*.bam"))

    if not sorted_bam_files:
        raise IndexError(
            "The list of sorted bam files is empty, expected them in STAR output dir:\n%s"
            % star_dir
        )

    # Don't know why this isn't multiprocessed, probably cos it was fast enough in serial. But consider multiprocessing if
    # the mem usage is low
    for sorted_bam_file in sorted_bam_files:
        sorted_bam_file_name = sorted_bam_file.name
        transcript_file_name = re.sub(".bam", ".stringtie.gtf", sorted_bam_file_name)
        transcript_file_path = stringtie_dir / transcript_file_name

        if transcript_file_path.exists():
            logger.info(
                "Found an existing stringtie gtf file, will not overwrite:\n%s"
                % transcript_file_path
            )
        else:
            logger.info(
                "Running Stringtie on: %s, writing output to:\n%s"
                % (sorted_bam_file_name, transcript_file_path)
            )
            subprocess.run(
                [
                    stringtie_path,
                    sorted_bam_file,
                    "-o",
                    transcript_file_path,
                    "-p",
                    str(num_threads),
                    "-t",
                    "-a",
                    "15",
                ]
            )

    # Now need to merge
    logger.info("Creating Stringtie merge input file: %s" % stringtie_merge_input_file)
    with open(stringtie_merge_input_file, "w+") as gtf_list_out:
        for gtf_file in stringtie_dir.glob("*.stringtie.gtf"):
            transcript_count = check_gtf_content(gtf_file, "transcript")
            if transcript_count > 0:
                gtf_list_out.write(f"{gtf_file}\n")
            else:
                logger.warning(
                    "Warning, skipping file with no transcripts:\n%s" % gtf_file
                )

    logger.info("Merging Stringtie results to:\n%s" % stringtie_merge_output_file)
    subprocess.run(
        [
            stringtie_path,
            "--merge",
            "-o",
            stringtie_merge_output_file,
            stringtie_merge_input_file,
        ]
    )


def run_scallop_assemble(scallop_path, stringtie_path, main_output_dir: pathlib.Path):
    if not scallop_path:
        scallop_path = shutil.which("scallop")
    check_exe(scallop_path)

    if not stringtie_path:
        stringtie_path = shutil.which("stringtie")
    check_exe(stringtie_path)

    scallop_dir = create_dir(main_output_dir, "scallop_output")

    output_file = scallop_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Scallop gtf file exists, skipping analysis")
            return

    stringtie_merge_input_file = scallop_dir / "scallop_assemblies.txt"
    stringtie_merge_output_file = scallop_dir / "annotation.gtf"
    star_dir = main_output_dir / "star_output"
    memory_limit = 40 * 1024**3

    if star_dir.exists():
        logger.info("Found a STAR output dir, will load bam files")

    sorted_bam_files = list(star_dir.glob("*.bam"))

    if not sorted_bam_files:
        raise IndexError(
            "Empty list of sorted bam files, expected them in STAR output dir:\n%s"
            % star_dir
        )

    # Don't know why this isn't multiprocessed, probably cos it was fast enough in serial.
    # But consider multiprocessing if the mem usage is low.
    for sorted_bam_file in sorted_bam_files:
        sorted_bam_file_name = sorted_bam_file.name
        transcript_file_name = re.sub(".bam", ".scallop.gtf", sorted_bam_file_name)
        transcript_file_path = scallop_dir / transcript_file_name

        if transcript_file_path.exists():
            logger.info(
                "Found an existing scallop gtf file, will not overwrite:\n%s"
                % transcript_file_path
            )
        else:
            logger.info(
                "Running Scallop on: %s, writing output to:\n%s"
                % (sorted_bam_file_name, transcript_file_path)
            )
            scallop_cmd = [
                scallop_path,
                "-i",
                sorted_bam_file,
                "-o",
                transcript_file_path,
                "--min_flank_length",
                "10",
            ]
            if memory_limit is not None:
                scallop_cmd = prlimit_command(scallop_cmd, memory_limit)

            return_value = None
            try:
                return_value = subprocess.check_output(scallop_cmd)
            except subprocess.CalledProcessError as ex:
                logger.error(
                    "Issue processing the following region with scallop\nReturn value: %s"
                    % return_value
                )

    #      subprocess.run([scallop_path,'-i',sorted_bam_file,'-o',transcript_file_path,'--min_flank_length','10'])

    # Now need to merge
    logger.info("Creating Stringtie merge input file: %s" % stringtie_merge_input_file)

    with open(stringtie_merge_input_file, "w+") as gtf_list_out:
        for gtf_file in scallop_dir.glob("*.scallop.gtf"):
            transcript_count = check_gtf_content(gtf_file, "transcript")
            if transcript_count > 0:
                gtf_list_out.write(f"{gtf_file}\n")
            else:
                logger.warning(
                    "Warning, skipping file with no transcripts:\n%s" % gtf_file
                )

    logger.info("Merging Scallop results to:\n%s" % stringtie_merge_output_file)
    subprocess.run(
        [
            stringtie_path,
            "--merge",
            "-o",
            stringtie_merge_output_file,
            stringtie_merge_input_file,
        ]
    )


def run_finalise_geneset(
    main_script_dir: pathlib.Path,
    main_output_dir: pathlib.Path,
    genome_file: Union[pathlib.Path, str],
    seq_region_names,
    validation_type,
    diamond_validation_db,
    num_threads: int,
):
    if validation_type is None:
        validation_type = "relaxed"
    logger.info("Setting validation type to %s" % validation_type)

    final_annotation_dir = create_dir(main_output_dir, "annotation_output")
    region_annotation_dir = create_dir(final_annotation_dir, "initial_region_gtfs")
    final_region_annotation_dir = create_dir(final_annotation_dir, "final_region_gtfs")
    utr_region_annotation_dir = create_dir(final_annotation_dir, "utr_region_gtfs")
    validation_dir = create_dir(final_annotation_dir, "cds_validation")
    seq_region_lengths = get_seq_region_lengths(genome_file, min_seq_length=0)

    # This used to be a list of output dirs and a loop which was neat, I'm converting to a list of conditions as
    # it's more straightforward with the renaming and having to merge scallop and stringtie
    protein_annotation_raw = main_output_dir / "genblast_output" / "annotation.gtf"
    minimap2_annotation_raw = main_output_dir / "minimap2_output" / "annotation.gtf"
    stringtie_annotation_raw = main_output_dir / "stringtie_output" / "annotation.gtf"
    scallop_annotation_raw = main_output_dir / "scallop_output" / "annotation.gtf"
    busco_annotation_raw = main_output_dir / "busco_output" / "annotation.gtf"

    transcript_selector_script = (
        main_script_dir / "support_scripts_perl" / "select_best_transcripts.pl"
    )
    finalise_geneset_script = (
        main_script_dir / "support_scripts_perl" / "finalise_geneset.pl"
    )
    clean_geneset_script = main_script_dir / "support_scripts_perl" / "clean_geneset.pl"
    clean_utrs_script = (
        main_script_dir / "support_scripts_perl" / "clean_utrs_and_lncRNAs.pl"
    )
    gtf_to_seq_script = main_script_dir / "support_scripts_perl" / "gtf_to_seq.pl"

    transcriptomic_annotation_raw = final_annotation_dir / "transcriptomic_raw.gtf"
    with open(transcriptomic_annotation_raw, "w+") as file_out:
        for transcriptomic_file in [
            minimap2_annotation_raw,
            scallop_annotation_raw,
            stringtie_annotation_raw,
        ]:
            if not transcriptomic_file.exists():
                logger.info(
                    'No annotation.gtf file found in "%s", skipping'
                    % transcriptomic_file
                )
                continue

            with open(transcriptomic_file) as file_in:
                for line in file_in:
                    file_out.write(line.rstrip())

    # Copy the raw files into the annotation dir, this is not needed as such,
    # but collecting them in one place and relabelling is helpful for a user
    if busco_annotation_raw.exists():
        subprocess.run(
            [
                "cp",
                busco_annotation_raw,
                final_annotation_dir / "busco_raw.gtf",
            ]
        )

    if protein_annotation_raw.exists():
        subprocess.run(
            [
                "cp",
                protein_annotation_raw,
                final_annotation_dir / "protein_raw.gtf",
            ]
        )

    gtf_files = ["transcriptomic_raw.gtf", "protein_raw.gtf", "busco_raw.gtf"]
    generic_select_cmd = [
        "perl",
        transcript_selector_script,
        "-genome_file",
        genome_file,
    ]
    pool = multiprocessing.Pool(num_threads)
    for seq_region_name in seq_region_names:
        # The selection script needs different params depending on whether the seqs are from transcriptomic data or not
        seq_region_length = seq_region_lengths[seq_region_name]
        region_details = f"{seq_region_name}.rs1.re{seq_region_length}"
        transcriptomic_region_gtf_path = (
            region_annotation_dir / f"{region_details}.trans.gtf"
        )
        busco_region_gtf_path = region_annotation_dir / f"{region_details}.busco.gtf"
        protein_region_gtf_path = (
            region_annotation_dir / f"{region_details}.protein.gtf"
        )

        if transcriptomic_annotation_raw.exists():
            logger.info("Finalising transcriptomic data for: %s" % seq_region_name)
            transcriptomic_annotation_select = re.sub(
                "_raw.gtf", "_sel.gtf", str(transcriptomic_annotation_raw)
            )
            cmd = generic_select_cmd.copy()
            cmd.extend(
                [
                    "-region_details",
                    region_details,
                    "-input_gtf_file",
                    transcriptomic_annotation_raw,
                    "-output_gtf_file",
                    transcriptomic_region_gtf_path,
                    "-cds_search",
                    "-final_biotype",
                    "transcriptomic",
                ]
            )
            pool.apply_async(subprocess_run_and_log, args=(cmd,))

        if busco_annotation_raw.exists():
            logger.info("Finalising BUSCO data for: %s" % seq_region_name)
            busco_annotation_select = re.sub(
                "_raw.gtf", "_sel.gtf", str(busco_annotation_raw)
            )
            cmd = generic_select_cmd.copy()
            cmd.extend(
                [
                    "-region_details",
                    region_details,
                    "-input_gtf_file",
                    busco_annotation_raw,
                    "-output_gtf_file",
                    busco_region_gtf_path,
                    "-all_cds_exons",
                    "-final_biotype",
                    "busco",
                ]
            )
            pool.apply_async(subprocess_run_and_log, args=(cmd,))

        if protein_annotation_raw.exists():
            logger.info("Finalising protein data for: %s" % seq_region_name)
            protein_annotation_select = re.sub(
                "_raw.gtf", "_sel.gtf", str(protein_annotation_raw)
            )
            cmd = generic_select_cmd.copy()
            cmd.extend(
                [
                    "-region_details",
                    region_details,
                    "-input_gtf_file",
                    protein_annotation_raw,
                    "-output_gtf_file",
                    protein_region_gtf_path,
                    "-clean_transcripts",
                    "-all_cds_exons",
                    "-final_biotype",
                    "protein",
                ]
            )
            pool.apply_async(subprocess_run_and_log, args=(cmd,))

    pool.close()
    pool.join()

    # At this point we will have the region files for all the,

    merge_finalise_output_files(
        final_annotation_dir=final_annotation_dir,
        region_annotation_dir=region_annotation_dir,
        extension=".trans.gtf",
        id_label="transcriptomic",
    )
    merge_finalise_output_files(
        final_annotation_dir=final_annotation_dir,
        region_annotation_dir=region_annotation_dir,
        extension=".busco.gtf",
        id_label="busco",
    )
    merge_finalise_output_files(
        final_annotation_dir=final_annotation_dir,
        region_annotation_dir=region_annotation_dir,
        extension=".protein.gtf",
        id_label="protein",
    )

    # Create a single GTF file with all the selected transcripts now that they have proper ids
    fully_merged_gtf_path = final_annotation_dir / "all_selected_transcripts.gtf"
    fully_merged_gtf_out = open(fully_merged_gtf_path, "w+")

    merge_gtf_cmd = ["cat"]
    merge_gtf_cmd.extend(list(final_annotation_dir.glob("*_sel.gtf")))
    subprocess.run(merge_gtf_cmd, stdout=fully_merged_gtf_out)
    fully_merged_gtf_out.close()

    # Now collapse the gene set
    generic_finalise_cmd = [
        "perl",
        finalise_geneset_script,
        "-genome_file",
        genome_file,
    ]

    pool = multiprocessing.Pool(num_threads)
    for seq_region_name in seq_region_names:
        seq_region_length = seq_region_lengths[seq_region_name]
        region_details = f"{seq_region_name}.rs1.re{seq_region_length}"
        final_region_gtf_path = (
            final_region_annotation_dir / f"{region_details}.final.gtf"
        )

        cmd = generic_finalise_cmd.copy()
        cmd.extend(
            [
                "-region_details",
                region_details,
                "-input_gtf_file",
                fully_merged_gtf_path,
                "-output_gtf_file",
                final_region_gtf_path,
            ]
        )
        pool.apply_async(subprocess_run_and_log, args=(cmd,))

    pool.close()
    pool.join()

    merge_finalise_output_files(
        final_annotation_dir=final_annotation_dir,
        region_annotation_dir=final_region_annotation_dir,
        extension=".final.gtf",
        id_label="prevalidation",
    )
    merged_gtf_file = final_annotation_dir / "prevalidation_sel.gtf"
    merged_cdna_file = final_annotation_dir / "prevalidation_sel.cdna.fa"
    merged_amino_acid_file = final_annotation_dir / "prevalidation_sel.prot.fa"
    updated_gtf_lines = validate_coding_transcripts(
        merged_cdna_file,
        merged_amino_acid_file,
        validation_dir,
        validation_type,
        diamond_validation_db,
        merged_gtf_file,
        num_threads,
    )
    postvalidation_gtf_file = final_annotation_dir / "postvalidation.gtf"
    with open(postvalidation_gtf_file, "w+") as file_out:
        for line in updated_gtf_lines:
            file_out.write(line)

    cleaned_initial_gtf_file = final_annotation_dir / "cleaned_pre_utr.gtf"
    cleaned_utr_gtf_file = final_annotation_dir / "annotation.gtf"

    cleaning_cmd = [
        "perl",
        clean_geneset_script,
        "-genome_file",
        genome_file,
        "-gtf_file",
        postvalidation_gtf_file,
        "-output_gtf_file",
        cleaned_initial_gtf_file,
    ]
    logger.info("Cleaning initial set:\n%s" % list_to_string(cleaning_cmd))
    subprocess.run(cleaning_cmd)

    # Clean UTRs
    generic_clean_utrs_cmd = [
        "perl",
        clean_utrs_script,
        "-genome_file",
        genome_file,
        "-input_gtf_file",
        cleaned_initial_gtf_file,
    ]
    pool = multiprocessing.Pool(num_threads)
    for seq_region_name in seq_region_names:
        seq_region_length = seq_region_lengths[seq_region_name]
        region_details = f"{seq_region_name}.rs1.re{seq_region_length}"
        utr_region_gtf_path = utr_region_annotation_dir / f"{region_details}.utr.gtf"

        cmd = generic_clean_utrs_cmd.copy()
        cmd.extend(
            [
                "-region_details",
                region_details,
                "-input_gtf_file",
                cleaned_initial_gtf_file,
                "-output_gtf_file",
                utr_region_gtf_path,
            ]
        )
        pool.apply_async(subprocess_run_and_log, args=(cmd,))
    pool.close()
    pool.join()

    merge_finalise_output_files(
        final_annotation_dir=final_annotation_dir,
        region_annotation_dir=utr_region_annotation_dir,
        extension=".utr.gtf",
        id_label="annotation",
    )
    subprocess.run(
        [
            "mv",
            final_annotation_dir / "annotation_sel.gtf",
            cleaned_utr_gtf_file,
        ]
    )

    dumping_cmd = [
        "perl",
        gtf_to_seq_script,
        "-genome_file",
        genome_file,
        "-gtf_file",
        cleaned_utr_gtf_file,
    ]
    logger.info(
        "Dumping transcript and translation sequences:\n%s"
        % list_to_string(dumping_cmd)
    )
    subprocess.run(dumping_cmd)

    logger.info("Finished creating gene set")


def validate_coding_transcripts(
    cdna_file,
    amino_acid_file,
    validation_dir: pathlib.Path,
    validation_type,
    diamond_validation_db,
    gtf_file,
    num_threads: int,
):
    logger.info("Running CDS validation with RNAsamba and CPC2")
    rnasamba_weights = "/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/rnasamba_data/full_length_weights.hdf5"
    rnasamba_output_path = validation_dir / "rnasamba.tsv.txt"
    cpc2_output_path = validation_dir / "cpc2.tsv"
    rnasamba_volume = f"{validation_dir}/:/app:rw"
    rnasamba_cmd = [
        "singularity",
        "exec",
        "--bind",
        rnasamba_volume,
        "/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/rnasamba_latest.sif",
        "rnasamba",
        "classify",
        rnasamba_output_path,
        cdna_file,
        rnasamba_weights,
    ]
    logger.info("rnasamba_cmd: %s" % list_to_string(rnasamba_cmd))
    subprocess.run(rnasamba_cmd)
    cpc2_volume = f"{validation_dir}/:/app:rw"
    cpc2_cmd = [
        "singularity",
        "exec",
        "--bind",
        cpc2_volume,
        "/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/test_cpc2.sif",
        "python3",
        "/CPC2_standalone-1.0.1/bin/CPC2.py",
        "-i",
        cdna_file,
        "--ORF",
        "-o",
        cpc2_output_path,
    ]
    logger.info("cpc2_cmd: %s" % list_to_string(cpc2_cmd))
    subprocess.run(cpc2_cmd)
    cpc2_output_path = f"{cpc2_output_path}.txt"

    check_file(rnasamba_output_path)
    check_file(cpc2_output_path)

    logger.info("diamond validation")
    diamond_results = None
    if diamond_validation_db is not None:
        diamond_output_dir = create_dir(validation_dir, "diamond_output")
        diamond_validation(
            diamond_validation_db, amino_acid_file, diamond_output_dir, num_threads
        )
        diamond_results = read_diamond_results(diamond_output_dir)

    logger.info("read results")
    rnasamba_results = read_rnasamba_results(rnasamba_output_path)
    cpc2_results = read_cpc2_results(cpc2_output_path)
    combined_results = combine_results(rnasamba_results, cpc2_results, diamond_results)
    logger.info("read gtf genes")
    parsed_gtf_genes = read_gtf_genes(gtf_file)
    updated_gtf_lines = update_gtf_genes(
        parsed_gtf_genes, combined_results, validation_type
    )

    return updated_gtf_lines


def diamond_validation(
    diamond_validation_db,
    amino_acid_file,
    diamond_output_dir: pathlib.Path,
    num_threads: int,
):
    batched_protein_files = split_protein_file(
        protein_file=amino_acid_file,
        protein_output_dir=diamond_output_dir,
        batch_size=100,
    )

    pool = multiprocessing.Pool(num_threads)
    for batched_protein_file in batched_protein_files:
        pool.apply_async(
            multiprocess_diamond,
            args=(
                batched_protein_file,
                diamond_output_dir,
                diamond_validation_db,
            ),
        )
    pool.close()
    pool.join()


def multiprocess_diamond(
    batched_protein_file: pathlib.Path,
    diamond_output_dir,
    diamond_validation_db,
):
    # batch_num = os.path.splitext(batched_protein_file)[0]
    # batch_dir = os.path.dirname(batched_protein_file)
    diamond_output_file = f"{batched_protein_file}.dmdout"
    diamond_cmd = [
        "diamond",
        "blastp",
        "--query",
        batched_protein_file,
        "--db",
        diamond_validation_db,
        "--out",
        diamond_output_file,
    ]

    logger.info(
        "Running diamond on %s:\n%s"
        % (batched_protein_file, list_to_string(diamond_cmd))
    )
    subprocess.run(diamond_cmd)
    subprocess.run(["mv", diamond_output_file, diamond_output_dir])


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
                    f'; biotype "{biotype}";',
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
                        f'; biotype "{biotype}";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    rnasamba_coding_potential == "coding"
                    and cpc2_coding_potential == "coding"
                    and peptide_length >= min_single_exon_pep_length
                ):
                    transcript_line = re.sub(
                        f'; biotype "{biotype}";',
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
                        f'; biotype "{biotype}";',
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
                        f'; biotype "{biotype}";',
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
                        f'; biotype "{biotype}";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                else:
                    continue
            else:
                if diamond_e_value is not None:
                    transcript_line = re.sub(
                        f'; biotype "{biotype}";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    rnasamba_coding_potential == "coding"
                    and cpc2_coding_potential == "coding"
                    and peptide_length >= min_multi_exon_pep_length
                ):
                    transcript_line = re.sub(
                        f'; biotype "{biotype}";',
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
                        f'; biotype "{biotype}";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif transcript_length >= 200:
                    transcript_line = re.sub(
                        f'; biotype "{biotype}";',
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


def read_rnasamba_results(file_path: Union[pathlib.Path, str]):
    results = []
    with open(file_path) as file_in:
        for line in file_in:
            line = line.rstrip()
            match = re.search(r"^sequence_name", line)
            if match:
                continue

            eles = line.split("\t")
            if not len(eles) == 3:
                continue

            transcript_id = eles[0]
            coding_probability = eles[1]
            coding_potential = eles[2]
            results.append([transcript_id, coding_probability, coding_potential])

    return results


def read_cpc2_results(file_path):
    results = []
    with open(file_path) as file_in:
        for line in file_in:
            line = line.rstrip()
            match = re.search(r"^#ID", line)
            if match:
                continue

            eles = line.split("\t")
            if not len(eles) == 9:
                continue

            transcript_id = eles[0]
            transcript_length = eles[1]
            peptide_length = eles[2]
            coding_probability = eles[7]
            coding_potential = eles[8]
            results.append(
                [
                    transcript_id,
                    coding_probability,
                    coding_potential,
                    transcript_length,
                    peptide_length,
                ]
            )

    return results


def read_diamond_results(diamond_output_dir: pathlib.Path):
    results = []
    for diamond_file in diamond_output_dir.glob("*.dmdout"):
        with open(diamond_file) as file_in:
            for line in file_in:
                line = line.rstrip()

                eles = line.split("\t")
                if not len(eles) == 12:
                    continue

                transcript_id = eles[0]
                e_value = eles[10]
                results.append([transcript_id, e_value])

    return results


def combine_results(rnasamba_results, cpc2_results, diamond_results):
    transcript_ids = {}

    for result in rnasamba_results:
        transcript_id = result[0]
        coding_probability = result[1]
        coding_potential = result[2]

        if transcript_id not in transcript_ids:
            transcript_ids[transcript_id] = [coding_probability, coding_potential]

    for result in cpc2_results:
        transcript_id = result[0]
        coding_probability = result[1]
        coding_potential = result[2]
        transcript_length = result[3]
        peptide_length = result[4]
        transcript_ids[transcript_id].extend(
            [coding_probability, coding_potential, transcript_length, peptide_length]
        )

    if diamond_results is not None:
        for result in diamond_results:
            transcript_id = result[0]
            e_value = result[1]
            # There seems to be an issue where there are a small number of sequences that don't make it into the cpc2/rnasamba output
            # Should code in a system for this, but it would be good to understand why it happens to begin with. Seems to be the same
            # number of missing seqs in both, so maybe a shared cut-off
            if transcript_id in transcript_ids:
                transcript_ids[transcript_id].extend([e_value])

    return transcript_ids


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


def merge_finalise_output_files(
    final_annotation_dir: pathlib.Path,
    region_annotation_dir: pathlib.Path,
    extension: str,
    id_label: str,
):
    """
    The below is not great, it's a bit messy because there might be some cases where there aren't
    translations. So it's not as straightforward as reading the records across all three files
    in parallel. The solution is to just load the seqs into memory and index them on the current
    header, which should correspond to a transcript/gene id in the GTF. When writing the results
    into the single merged files the ids will be updated to be unique and consistent across the
    three file types
    """
    merged_gtf_file = final_annotation_dir / f"{id_label}_sel.gtf"
    merged_cdna_file = final_annotation_dir / f"{id_label}_sel.cdna.fa"
    merged_amino_acid_file = final_annotation_dir / f"{id_label}_sel.prot.fa"

    gene_id_counter = 0
    transcript_id_counter = 0
    with open(merged_gtf_file, "w+") as gtf_out, open(
        merged_cdna_file, "w+"
    ) as cdna_out, open(merged_amino_acid_file, "w+") as amino_acid_out:
        gtf_files = region_annotation_dir.glob(f"*{extension}")
        for gtf_file in gtf_files:
            logger.info("GTF file: %s" % gtf_file)
            cdna_file = f"{gtf_file}.cdna"
            with open(cdna_file) as cdna_in:
                cdna_seq_index = fasta_to_dict(cdna_in.readlines())
            amino_acid_file = f"{gtf_file}.prot"
            with open(amino_acid_file) as amino_acid_in:
                amino_acid_seq_index = fasta_to_dict(amino_acid_in.readlines())

            current_gene_id = ""
            with open(gtf_file) as gtf_in:
                for line in gtf_in:
                    if re.search(r"^#", line):
                        continue

                    eles = line.split("\t")
                    if not len(eles) == 9:
                        continue

                    match = re.search(
                        r'gene_id "([^"]+)".+transcript_id "([^"]+)"', line
                    )
                    if match and eles[2] == "transcript":
                        transcript_id_counter += 1

                    gene_id = match.group(1)
                    transcript_id = match.group(2)

                    if not current_gene_id:
                        gene_id_counter += 1
                        current_gene_id = gene_id

                    if not gene_id == current_gene_id:
                        gene_id_counter += 1
                        current_gene_id = gene_id

                    new_gene_id = f"{id_label}_{gene_id_counter}"
                    new_transcript_id = f"{id_label}_{transcript_id_counter}"
                    line = re.sub(
                        f'gene_id "{gene_id}"', f'gene_id "{new_gene_id}"', line
                    )
                    line = re.sub(
                        f'transcript_id "{transcript_id}"',
                        f'transcript_id "{new_transcript_id}"',
                        line,
                    )
                    gtf_out.write(line)

                    if eles[2] == "transcript":
                        new_header = f">{new_transcript_id}\n"
                        cdna_out.write(new_header + cdna_seq_index[transcript_id])

                        if transcript_id in amino_acid_seq_index:
                            amino_acid_out.write(
                                new_header + amino_acid_seq_index[transcript_id]
                            )


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
    logger.info("subprocess_run_and_log command: %s" % list_to_string(command))
    subprocess.run(command)


def get_sequence(
    seq_region,
    start: int,
    end: int,
    strand: int,
    fasta_file,
    output_dir: Union[pathlib.Path, str],
) -> str:
    """
    TODO:
    The functionality of this function could be better implemented by using pybedtools
    or even better pyfaidx.
    """
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


def main():
    """
    main function
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output_dir",
        type=str,
        help="Path where the output and temp files will write to. Uses current dir by default",
    )
    parser.add_argument(
        "--genome_file", type=str, required=True, help="Path to the fasta genome file"
    )
    parser.add_argument(
        "--num_threads", type=int, default=1, help="Number of threads to use"
    )
    parser.add_argument(
        "--run_masking",
        action="store_true",
        help="Run Red to find repeats and softmask the genome. Otherwise provide a softmasked genome",
    )
    parser.add_argument(
        "--red_path",
        type=str,
        help="Path to Red executable. See http://toolsmith.ens.utulsa.edu",
    )
    parser.add_argument(
        "--genblast_path",
        type=str,
        help="Path to GenBlast executable. See http://genome.sfu.ca/genblast/download.html",
    )
    parser.add_argument(
        "--convert2blastmask_path",
        type=str,
        help="Path to convert2blastmask executable",
    )
    parser.add_argument(
        "--makeblastdb_path", type=str, help="Path to makeblastdb executable"
    )
    parser.add_argument(
        "--run_genblast",
        action="store_true",
        help="Run GenBlast to align protein sequences",
    )
    parser.add_argument(
        "--genblast_timeout",
        type=int,
        default=10800,
        help="GenBlast timeout in seconds",
    )
    parser.add_argument(
        "--run_busco",
        action="store_true",
        help="Run GenBlast to align BUSCO protein sequences",
    )
    parser.add_argument(
        "--protein_file",
        type=str,
        help="Path to a fasta file with protein sequences",
    )
    parser.add_argument(
        "--busco_protein_file",
        type=str,
        help="Path to a fasta file with BUSCO protein sequences",
    )
    parser.add_argument(
        "--rfam_accessions_file",
        type=str,
        help="Path to a file with Rfam CM accessions, one accession per line, to use with cmsearch",
    )
    parser.add_argument(
        "--run_star", action="store_true", help="Run STAR for short read alignment"
    )
    parser.add_argument(
        "--star_path", type=str, help="Path to STAR for short read alignment"
    )
    parser.add_argument(
        "--max_reads_per_sample",
        nargs="?",
        type=int,
        default=0,
        help="The maximum number of reads to use per sample. Default=0 (unlimited)",
    )
    parser.add_argument(
        "--max_total_reads",
        nargs="?",
        type=int,
        default=0,
        help="The maximum total number of reads. Default=0 (unlimited)",
    )
    parser.add_argument(
        "--short_read_fastq_dir",
        help="Path to short read fastq dir for running with STAR",
    )
    parser.add_argument(
        "--max_intron_length",
        nargs="?",
        type=int,
        default=100_000,
        help="The maximum intron size for alignments. Default=100,000",
    )
    parser.add_argument(
        "--run_minimap2",
        action="store_true",
        help="Run minimap2 for long read alignment",
    )
    parser.add_argument(
        "--minimap2_path",
        type=str,
        help="Path to minimap2 for long read alignment",
    )
    parser.add_argument(
        "--paftools_path",
        type=str,
        help="Path to paftools for SAM to BED conversion",
    )
    parser.add_argument(
        "--long_read_fastq_dir",
        type=str,
        help="Path to long read fastq dir for running with minimap2",
    )
    parser.add_argument(
        "--run_augustus",
        action="store_true",
        help="Run Augustus with hints for gene/transcript prediction",
    )
    parser.add_argument("--augustus_path", type=str, help="Path to Augustus")
    parser.add_argument(
        "--run_stringtie",
        action="store_true",
        help="Run Stringtie on the results from the STAR alignments",
    )
    parser.add_argument(
        "--run_scallop",
        action="store_true",
        help="Run Scallop on the results from the STAR alignments",
    )
    parser.add_argument("--stringtie_path", type=str, help="Path to Stringtie")
    parser.add_argument("--scallop_path", type=str, help="Path to Scallop")
    parser.add_argument(
        "--subsample_script_path",
        type=str,
        help="Path to ensembl-anno subsampling script",
    )
    parser.add_argument("--samtools_path", type=str, help="Path to subsampling script")
    parser.add_argument(
        "--finalise_geneset",
        action="store_true",
        help="Used to finalise the gene set from the various GTF files generated",
    )
    parser.add_argument(
        "--db_details",
        type=str,
        help="A comma separated string of dbname,host,port,user,pass",
    )
    parser.add_argument(
        "--run_cmsearch",
        action="store_true",
        help="Search for sncRNA structures using Rfam and cmsearch",
    )
    parser.add_argument(
        "--run_trf", action="store_true", help="Run TRF to find tandem repeats"
    )
    parser.add_argument("--trf_path", type=str, help="Path to TRF")
    parser.add_argument(
        "--run_dust",
        action="store_true",
        help="Run Dust to find low complexity regions",
    )
    parser.add_argument("--dust_path", type=str, help="Path to Dust")
    parser.add_argument(
        "--run_repeatmasker",
        action="store_true",
        help="Run RepeatMasker to find repeat regions",
    )
    parser.add_argument(
        "--repeatmasker_path", action="store_true", help="Path to RepeatMasker"
    )
    parser.add_argument(
        "--run_trnascan", action="store_true", help="Run tRNAscan-SE to find tRNAs"
    )
    parser.add_argument("--trnascan_path", type=str, help="Path to tRNAscan-SE")
    parser.add_argument(
        "--trnascan_filter_path",
        type=str,
        help="Path to tRNAscan-SE high confidence filter",
    )
    parser.add_argument(
        "--run_cpg", action="store_true", help="Run cpg_lh to find CpG islands"
    )
    parser.add_argument("--cpg_path", type=str, help="Path to cpg_lh")
    parser.add_argument(
        "--run_eponine",
        action="store_true",
        help="Run Eponine to find transcription start sites",
    )
    parser.add_argument("--eponine_path", type=str, help="Path to Eponine jar file")
    parser.add_argument(
        "--java_path", type=str, help="Path to Java for use with Eponine"
    )
    parser.add_argument(
        "--run_full_annotation",
        action="store_true",
        help="Run a full annotation, will automatically check for input data and run tools based on that",
    )
    parser.add_argument("--run_repeats", action="store_true", help="Run Red, Dust, TRF")
    parser.add_argument(
        "--run_simple_features", action="store_true", help="Run CpG, Eponine"
    )
    parser.add_argument(
        "--run_sncrnas", action="store_true", help="Run Rfam, tRNAscan-SE"
    )
    parser.add_argument(
        "--run_transcriptomic",
        action="store_true",
        help="Run STAR, Stringtie2, Scallop, minimap2 (if short_read_fastq_dir and/or long_read_fastq_dir are provided)",
    )
    parser.add_argument(
        "--run_proteins",
        action="store_true",
        help="Run GenBlast if protein_file and/or busco_protein_file have also been set",
    )
    parser.add_argument(
        "--diamond_validation_db",
        type=str,
        help="Use a Diamond db with blastp mode to help validate cds sequences",
    )
    parser.add_argument(
        "--validation_type",
        type=str,
        help='The strength of evidence needed to validate an ORF as protein coding, can be "relaxed" or "moderate"',
    )
    parser.add_argument(
        "--load_to_ensembl_db",
        # TODO
        # resolve boolean vs string argument type
        action="store_true",
        help="Load results to an Ensembl db, must also provide the db_details flag",
    )
    parser.add_argument(
        "--trim_fastq",
        action="store_true",
        help="Trim the short read files using Trim Galore",
    )
    parser.add_argument(
        "--delete_pre_trim_fastq",
        action="store_true",
        help="Delete the original fastq files after trimming",
    )
    parser.add_argument(
        "--repeatmasker_library", type=str, help="Specify library for repeatmasker"
    )
    parser.add_argument(
        "--repeatmasker_species",
        type=str,
        # TODO
        # does a default value have to been defined here?
        # run_repeatmasker_regions function doesn't assume an existing value for species
        # default="homo",
        help="Specify species for repeatmasker (default homo)",
    )
    args = parser.parse_args()

    work_dir = args.output_dir
    genome_file = args.genome_file
    num_threads = args.num_threads
    masked_genome_file = genome_file  # This will be updated later if Red is run
    run_masking = args.run_masking
    red_path = args.red_path
    genblast_path = args.genblast_path
    convert2blastmask_path = args.convert2blastmask_path
    makeblastdb_path = args.makeblastdb_path
    run_genblast = args.run_genblast
    genblast_timeout = args.genblast_timeout
    run_busco = args.run_busco
    protein_file = args.protein_file
    busco_protein_file = args.busco_protein_file
    rfam_accessions_file = args.rfam_accessions_file
    run_star = args.run_star
    star_path = args.star_path
    short_read_fastq_dir = args.short_read_fastq_dir
    max_intron_length = args.max_intron_length
    max_reads_per_sample = args.max_reads_per_sample
    max_total_reads = args.max_total_reads
    run_minimap2 = args.run_minimap2
    minimap2_path = args.minimap2_path
    paftools_path = args.paftools_path
    long_read_fastq_dir = args.long_read_fastq_dir
    run_augustus = args.run_augustus
    augustus_path = args.augustus_path
    run_stringtie = args.run_stringtie
    run_scallop = args.run_scallop
    stringtie_path = args.stringtie_path
    scallop_path = args.scallop_path
    subsample_script_path = args.subsample_script_path
    samtools_path = args.samtools_path
    finalise_geneset = args.finalise_geneset
    db_details = args.db_details
    run_cmsearch = args.run_cmsearch
    run_trf = args.run_trf
    trf_path = args.trf_path
    run_dust = args.run_dust
    dust_path = args.dust_path
    run_trnascan = args.run_trnascan
    trnascan_path = args.trnascan_path
    trnascan_filter_path = args.trnascan_filter_path
    run_cpg = args.run_cpg
    cpg_path = args.cpg_path
    run_eponine = args.run_eponine
    eponine_path = args.eponine_path
    java_path = args.java_path
    run_repeatmasker = args.run_repeatmasker
    repeatmasker_path = args.repeatmasker_path
    run_full_annotation = args.run_full_annotation
    run_repeats = args.run_repeats
    run_simple_features = args.run_simple_features
    run_sncrnas = args.run_sncrnas
    run_transcriptomic = args.run_transcriptomic
    run_proteins = args.run_proteins
    diamond_validation_db = args.diamond_validation_db
    validation_type = args.validation_type
    load_to_ensembl_db = args.load_to_ensembl_db
    trim_fastq = args.trim_fastq
    delete_pre_trim_fastq = args.delete_pre_trim_fastq
    library = args.repeatmasker_library
    species = args.repeatmasker_species

    main_script_dir = pathlib.Path(os.path.realpath(__file__)).parent

    genome_file = pathlib.Path(genome_file)
    if not genome_file.exists():
        raise FileNotFoundError("File does not exist: %s" % genome_file)

    if not work_dir:
        work_dir = pathlib.Path.cwd()
    work_dir = pathlib.Path(work_dir)

    # save log to a file
    log_file_path = work_dir / "ensembl_anno.log"
    add_log_file_handler(logger, log_file_path)

    logger.info("work directory: %s" % work_dir)

    if not work_dir.exists():
        logger.info("Work dir does not exist, creating:")
        create_dir(work_dir)

    if num_threads == 1:
        logger.warning(
            "Thread count is set to the default value 1; this might be slow."
        )

    # If the run_full_annotation flag is set then we want to set a standardised set of analyses
    if run_full_annotation:
        run_repeats = True
        run_simple_features = True
        run_sncrnas = True
        run_transcriptomic = True
        run_proteins = True
        finalise_geneset = True

    # These are subsets of the analyses that can be run, group by type
    if run_repeats:
        run_masking = True
        run_dust = True
        run_trf = True

    if run_simple_features:
        run_cpg = True
        run_eponine = True

    if run_sncrnas:
        if rfam_accessions_file:
            run_cmsearch = True
        run_trnascan = True

    if run_transcriptomic:
        if short_read_fastq_dir:
            run_star = True
            run_scallop = True
            run_stringtie = True
        if long_read_fastq_dir:
            run_minimap2 = True

    if run_proteins:
        if protein_file:
            run_genblast = True
        if busco_protein_file:
            run_busco = True

    # Collect a list of seq region names, most useful for multiprocessing regions
    seq_region_names = get_seq_region_names(genome_file)
    logger.debug("seq region names:\n%s" % seq_region_names)

    # Repeat analyses
    ############################################################################
    if run_masking:
        logger.info("Running masking via Red")
        masked_genome_file = run_red(red_path, genome_file, work_dir)
        logger.info("Masked genome file: %s" % masked_genome_file)
    else:
        logger.info("Not running masking, presuming the genome file is softmasked")

    if run_dust:
        logger.info("Annotating low complexity regions")
        run_dust_regions(genome_file, dust_path, work_dir, num_threads)

    if run_trf:
        logger.info("Annotating tandem repeats")
        run_trf_repeats(genome_file, trf_path, work_dir, num_threads)

    if run_repeatmasker:
        logger.info("Annotating repeats with RepeatMasker")
        run_repeatmasker_regions(
            genome_file, repeatmasker_path, library, species, work_dir, num_threads
        )
    ############################################################################

    # Simple feature analyses
    ############################################################################
    if run_cpg:
        logger.info("Annotating CpG islands")
        run_cpg_regions(genome_file, cpg_path, work_dir, num_threads)

    if run_eponine:
        logger.info("Running Eponine to find transcription start sites")
        run_eponine_regions(genome_file, java_path, eponine_path, work_dir, num_threads)
    ############################################################################

    # sncRNA analyses
    ############################################################################
    # Search Rfam with cmsearch
    if run_cmsearch:
        logger.info("Annotating sncRNAs")
        run_cmsearch_regions(
            genome_file, None, None, None, rfam_accessions_file, work_dir, num_threads
        )

    if run_trnascan:
        logger.info("Annotating tRNAs")
        run_trnascan_regions(
            genome_file, trnascan_path, trnascan_filter_path, work_dir, num_threads
        )
    ############################################################################

    # Transcriptomic analyses
    ############################################################################
    if trim_fastq:
        run_trimming(work_dir, short_read_fastq_dir, delete_pre_trim_fastq, num_threads)

    # Run STAR
    if run_star:
        logger.info("Running STAR")
        run_star_align(
            star_path=star_path,
            trim_fastq=trim_fastq,
            subsample_script_path=subsample_script_path,
            main_output_dir=work_dir,
            short_read_fastq_dir=short_read_fastq_dir,
            genome_file=genome_file,
            max_reads_per_sample=max_reads_per_sample,
            max_total_reads=max_total_reads,
            max_intron_length=max_intron_length,
            num_threads=num_threads,
            main_script_dir=main_script_dir,
        )

    # Run Scallop
    if run_scallop:
        logger.info("Running Scallop")
        run_scallop_assemble(scallop_path, stringtie_path, work_dir)

    # Run Stringtie
    if run_stringtie:
        logger.info("Running StringTie")
        run_stringtie_assemble(stringtie_path, samtools_path, work_dir, num_threads)

    # Run minimap2
    if run_minimap2:
        logger.info("Running minimap2")
        run_minimap2_align(
            minimap2_path,
            paftools_path,
            work_dir,
            long_read_fastq_dir,
            genome_file,
            max_intron_length,
            num_threads,
        )

    if run_transcriptomic:
        check_transcriptomic_output(work_dir)
    ############################################################################

    # Protein analyses
    ############################################################################
    # Run GenBlast
    if run_genblast:
        logger.info("Running GenBlast")
        run_genblast_align(
            genblast_path,
            convert2blastmask_path,
            makeblastdb_path,
            work_dir / "genblast_output",
            protein_file,
            masked_genome_file,
            max_intron_length,
            num_threads,
            genblast_timeout,
            main_script_dir,
        )

    # Run GenBlast on BUSCO set, gives higher priority when creating the final genes in cases where transcriptomic data are missing or fragmented
    if run_busco:
        logger.info("Running GenBlast of BUSCO proteins")
        run_genblast_align(
            genblast_path,
            convert2blastmask_path,
            makeblastdb_path,
            work_dir / "busco_output",
            busco_protein_file,
            masked_genome_file,
            max_intron_length,
            num_threads,
            genblast_timeout,
            main_script_dir,
        )
    ############################################################################

    # Finalisation analyses
    ############################################################################
    # Do some magic
    if finalise_geneset:
        logger.info("Finalise geneset")
        run_finalise_geneset(
            main_script_dir,
            work_dir,
            genome_file,
            seq_region_names,
            validation_type,
            diamond_validation_db,
            num_threads,
        )
    ############################################################################

    # Other analyses
    ############################################################################
    # Run Augustus
    if run_augustus:
        logger.info("Running Augustus")
        run_augustus_predict(augustus_path, work_dir, masked_genome_file, num_threads)
    ############################################################################

    if load_to_ensembl_db:
        load_results_to_ensembl_db(
            main_script_dir,
            load_to_ensembl_db,
            genome_file,
            work_dir,
            db_details,
            num_threads,
        )


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.info("Interrupted with CTRL-C, exiting...")
