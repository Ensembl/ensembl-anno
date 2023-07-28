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
GenBlast identifies homologous gene sequences in genomic databases.
One of the key features of GenBlast is its flexibility to handle
comparative genomics tasks and accurately identify homologs even when
the sequences have undergone significant evolutionary changes.
This capability makes it a valuable resource for researchers studying gene
evolution, gene families, and gene function across diverse species.

GenBlast has been widely used in various genomic analyses and is available as a standalone command-line tool or as part of different bioinformatics pipelines. Researchers in the field of comparative genomics and gene function analysis often rely on GenBlast to perform sensitive homology searches and obtain valuable insights into the evolutionary relationships and functional conservation of genes in different organisms.


She, R., Chu, J.S., Uyar, B., Wang, J., Wang, K., and Chen, N. (2011).
GenBlastA: enabling BLAST to identify homologous gene sequences.
Genome Res., 21(5): 936-949.
"""
__all__ = ["run_trnascan"]

import logging
import logging.config
import multiprocessing
from os import PathLike
from pathlib import Path
import re
import subprocess
from typing import List
import argschema

from ensembl.tools.anno.utils._utils import (
    check_exe,
    check_file,
    create_dir,
    check_gtf_content,
    get_seq_region_length,
    get_slice_id,
    slice_output_to_gtf,
    get_sequence,
)

logger = logging.getLogger(__name__)

def run_genblast_align(
        masked_genome:Path,
        output_dir:Path,
        protein_dataset:Path,
        genblast_timeout_secs:int,
        max_intron_length:int,
        genblast_bin : Path=Path("genblast"),
        convert2blastmask_bin: Path=Path("convert2blastmask"),
        makeblastdb_bin:Path=Path("makeblastdb"),
        genblast_timeout_secs:int,
        num_threads:int=1,
)->None:
    """
    Executes GenBlast on genomic slices
    Args:
            masked_genome : Masked genome file path.
            output_dir: Working directory path.
            protein_dataset: Protein dataset (Uniprot/OrthoDb) path.
            genblast_timeout_secs: Time for timeout (sec).
            max_intron_length: Maximum intron length,
            genblast_bin : Path=Path("genblast"),
            convert2blastmask_bin: Path=Path("convert2blastmask"),
            makeblastdb_bin:Path=Path("makeblastdb"),
            genblast_timeout_secs:int,
            num_threads: int, number of threads.
    """

    check_exe(genblast_bin)
    check_exe(convert2blastmask_bin)
    check_exe(makeblastdb_bin)
    genblast_dir = create_dir(output_dir, "genblast_output")
    output_file = genblast_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Genblast gtf file exists, skipping analysis")
            return

    genblast_output_file = genblast_dir / "genblast"

    asnb_file = f"{masked_genome}.asnb"
    if not os.path.exists("alignscore.txt"):
        # Get the repo directory
        repo_root_dir = Path(__file__).parent.parent.parent.parent.parent
        shutil.copy(f"{repo_root_dir}/data/alignscore.txt", "./"
        )
    #        subprocess.run(
    #            [
    #                "cp",
    #                os.environ["ENSCODE"] + "/ensembl-anno/support_files/alignscore.txt",
    #                "./",
    #            ]
    #        )

    if not os.path.exists(masked_genome_file):
        raise IOError("Masked genome file does not exist: %s" % masked_genome_file)

    if not os.path.exists(protein_file):
        raise IOError("Protein file does not exist: %s" % protein_file)

    if not os.path.exists(asnb_file):
        run_convert2blastmask(convert2blastmask_path, masked_genome_file, asnb_file)
    else:
        logger.info("Found an existing asnb, so will skip convert2blastmask")

    if not os.path.exists(asnb_file):
        raise IOError("asnb file does not exist: %s" % asnb_file)

    run_makeblastdb(makeblastdb_path, masked_genome_file, asnb_file)

    batched_protein_files = split_protein_file(protein_file, genblast_dir, 20)

    pool = multiprocessing.Pool(int(num_threads))
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
    batched_protein_file,
    masked_genome_file,
    genblast_path,
    genblast_timeout_secs,
    max_intron_length,
):

    batch_num = os.path.splitext(batched_protein_file)[0]
    batch_dir = os.path.dirname(batched_protein_file)
    logger.info("Running GenBlast on " + batched_protein_file + ":")

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

    logger.info(" ".join(genblast_cmd))
    # Using the child process termination as described here:
    # https://alexandra-zaharia.github.io/posts/kill-subprocess
    # -and-its-children-on-timeout-python/
    try:
        p = subprocess.Popen(genblast_cmd, start_new_session=True)
        p.wait(timeout=genblast_timeout_secs)
    except subprocess.TimeoutExpired:
        logger.error("Timeout reached for file:\n" + batched_protein_file)
        subprocess.run(["touch", (batched_protein_file + ".except")])
        os.killpg(os.getpgid(p.pid), signal.SIGTERM)

    files_to_delete = glob.glob(batched_protein_file + "*msk.blast*")
    files_to_delete.append(batched_protein_file)
    for file_to_delete in files_to_delete:
        subprocess.run(["rm", file_to_delete])


def generate_genblast_gtf(genblast_dir):
    logger.info("generate_genblast_gtf")
    file_out_name = os.path.join(genblast_dir, "annotation.gtf")
    file_out = open(file_out_name, "w+")
    genblast_extension = "_1.1c_2.3_s1_0_16_1"
    for root, dirs, files in os.walk(genblast_dir):
        for genblast_file in files:
            genblast_file = os.path.join(root, genblast_file)
            if genblast_file.endswith(".gff"):
                gtf_string = convert_gff_to_gtf(genblast_file)
                file_out.write(gtf_string)
            elif (
                genblast_file.endswith(".fa.blast")
                or genblast_file.endswith(".fa.blast.report")
                or genblast_file.endswith(genblast_extension)
            ):
                # subprocess.run(["rm", genblast_file])
                pathlib.Path(genblast_file).unlink()
    file_out.close()

def split_protein_file(protein_file, protein_output_dir, batch_size):
    if batch_size is None:
        batch_size = 20

    batched_protein_files = []

    for i in range(0, 10):
        utils.create_dir(protein_output_dir, ("bin_" + str(i)))

    file_in = open(protein_file)
    line = file_in.readline()
    seq_count = 0
    batch_count = 0
    current_record = ""
    initial_seq = 1
    while line:
        num_dir = random.randint(0, 9)
        match = re.search(r">(.+)$", line)
        if match and not initial_seq and seq_count % batch_size == 0:
            file_out_name = os.path.join(
                protein_output_dir,
                ("bin_" + str(random.randint(0, 9))),
                (str(batch_count) + ".fa"),
            )
            file_out = open(file_out_name, "w+")
            file_out.write(current_record)
            file_out.close()
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
        line = file_in.readline()
    file_in.close()

    if current_record:
        file_out_name = os.path.join(
            protein_output_dir,
            ("bin_" + str(random.randint(0, 9))),
            (str(batch_count) + ".fa"),
        )
        file_out = open(file_out_name, "w+")
        file_out.write(current_record)
        file_out.close()
        batched_protein_files.append(file_out_name)

    return batched_protein_files


def run_convert2blastmask(convert2blastmask_path, masked_genome_file, asnb_file):

    asnb_file = masked_genome_file + ".asnb"
    logger.info("Running convert2blastmask prior to GenBlast:")
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
    logger.info(" ".join(cmd))
    subprocess.run(cmd)
    logger.info("Completed running convert2blastmask")


def run_makeblastdb(makeblastdb_path, masked_genome_file, asnb_file):

    logger.info("Running makeblastdb prior to GenBlast")
    subprocess.run(
        [
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
    )
    logger.info("Completed running makeblastdb")

