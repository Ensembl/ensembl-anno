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
import sys
import errno
import logging

from typing import List, Union
from utils import (
    add_log_file_handler,
    check_exe,
    check_file,
    check_gtf_content,
    create_dir,
    get_seq_region_lengths,
    logger,
    prlimit_command,
    load_results_to_ensembl_db,
    generic_load_records_to_ensembl_db,
    multiprocess_load_records_to_ensembl_db,
    batch_gtf_records,
    run_find_orfs,
    find_orf_phased_region,
    slice_output_to_gtf,
    convert_gff_to_gtf,
    set_attributes,
    create_slice_ids,
    update_gtf_genes,
    read_gtf_genes,
    fasta_to_dict,
    splice_junction_to_gff,
    split_genome,
    multiprocess_generic,
    reverse_complement,
    get_seq_region_names,
    slice_genome,
    subprocess_run_and_log,
    get_sequence,
    seq_region_names,
    list_to_string,
#    coallate_results,
)

def run_genblast_align(
    genblast_path,
    convert2blastmask_path,
    makeblastdb_path,
    genblast_dir: pathlib.Path,
    protein_file,
    masked_genome_file,
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

    if not os.path.exists(masked_genome_file):
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

