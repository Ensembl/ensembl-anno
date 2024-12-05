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
GenBlast has been widely used in various genomic analyses and is available as
a standalone command-line tool or as part of different bioinformatics pipelines.
Researchers in the field of comparative genomics and gene function analysis
often rely on GenBlast to perform sensitive homology searches and obtain
valuable insights into the evolutionary relationships and functional conservation
of genes in different organisms.


She, R., Chu, J.S., Uyar, B., Wang, J., Wang, K., and Chen, N. (2011).
GenBlastA: enabling BLAST to identify homologous gene sequences.
Genome Res., 21(5): 936-949.
"""
__all__ = ["run_genblast"]

import logging
import logging.config
import multiprocessing
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import argparse

from ensembl.tools.anno.utils._utils import (
    check_exe,
    create_dir,
    check_gtf_content,
    split_protein_file
)

logger = logging.getLogger(__name__)


def run_genblast(#pylint:disable=dangerous-default-value
    masked_genome: Path,
    output_dir: Path,
    protein_dataset: Path,
    max_intron_length: int,
    genblast_timeout_secs: int = 10800,
    genblast_bin: Path = Path("genblast"),
    convert2blastmask_bin: Path = Path("convert2blastmask"),
    makeblastdb_bin: Path = Path("makeblastdb"),
    num_threads: int = 1,
    protein_set: str = "uniprot",
) -> None:
    """
    
    Executes GenBlast on genomic slices
    
            :param masked_genome: Masked genome file path.
            :type masked_genome: Path
            :param output_dir: Working directory path.
            :type output_dir: Path
            :param protein_dataset: Protein dataset (Uniprot/OrthoDb) path.
            :type protein_dataset: Path
            :param genblast_timeout_secs: Time for timeout (sec).
            :type genblast_timeout_secs: int, default 10800
            :param max_intron_length: Maximum intron length.
            :type max_intron_length: int 
            :param genblast_bin: Software path.
            :type genblast_bin: Path, default genblast
            :param convert2blastmask_bin: Software path.
            :type convert2blastmask_bin: Path, default convert2blastmask
            :param makeblastdb_bin: Software path.
            :type makeblastdb_bin: Path, default makeblastdb
            :param genblast_timeout: seconds
            :type genblast_timeout: int, default 1
            :param num_threads: int, number of threads.
            :type num_threads: int, default 1 
            :param protein_set: Source 
            :type str: ["uniprot", "orthodb"]
            
            :return: None
            :rtype: None
            
    """

    check_exe(genblast_bin)
    check_exe(convert2blastmask_bin)
    check_exe(makeblastdb_bin)
    if protein_set not in {"uniprot", "orthodb"}:
        raise ValueError("protein_set must be either 'uniprot' or 'orthodb'")
    if protein_set == "uniprot":
        genblast_dir = create_dir(output_dir, "uniprot_output")
    elif protein_set == "orthodb":
        genblast_dir = create_dir(output_dir, "orthodb_output")
    output_file = genblast_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Genblast gtf file exists, skipping analysis")
            return
    logging.info(Path(f"{output_dir}/alignscore.txt"))
    if not Path(f"{genblast_dir}/alignscore.txt").exists():
        # Get the repo directory
        repo_root_dir = Path(__file__).parents[6]
        shutil.copy(Path(f"{repo_root_dir}/data/alignscore.txt"), genblast_dir)

    if not masked_genome.exists():
        raise IOError(f"Masked genome file does not exist: {masked_genome}")
    if not protein_dataset.exists():
        raise IOError(f"Protein file does not exist: {protein_dataset}")
    asnb_file = Path(f"{masked_genome}.asnb")
    if asnb_file.exists():
        logger.info("Found an existing asnb, so will skip convert2blastmask")
    else:
        _run_convert2blastmask(convert2blastmask_bin, masked_genome, asnb_file)
    _run_makeblastdb(makeblastdb_bin, masked_genome, asnb_file)
    batched_protein_files = split_protein_file(
        protein_dataset, genblast_dir, num_threads
    )
    pool = multiprocessing.Pool(num_threads)  # pylint:disable=consider-using-with
    for batched_protein_file in batched_protein_files:
        pool.apply_async(
            _multiprocess_genblast,
            args=(
                batched_protein_file,
                masked_genome,
                genblast_bin,
                genblast_timeout_secs,
                max_intron_length,
            ),
        )
    pool.close()
    pool.join()
    _generate_genblast_gtf(genblast_dir)
    for i in range(0, 10):
        shutil.rmtree(genblast_dir / f"bin_{i}")
    logger.info("Completed running GenBlast")


def _multiprocess_genblast(
    protein_file: Path,
    masked_genome: Path,
    genblast_bin: Path,
    genblast_timeout: int,
    max_intron_length: int,
):
    """
    Executes GenBlast on genomic slice
    Args:
            protein_file: Path of a single batched file.
            masked_genome : Masked genome file path.
            genblast_bin : Software path.
            genblast_timeout_secs: Time for timeout (sec).
            max_intron_length: Maximum intron length.
            Command line options:
            -P	Search program used to produce HSPs,
                can be either "blast" or "wublast", default is "blast",
                optional
            -p	specifies the program option of genBlast: genblasta or genblastg
            -q	List of query sequences to blast, must be in fasta format,
                required
            -t	The target database of genomic sequences in fasta format,
                required
            -g	parameter for blast: Perform gapped alignment (T/F)
                [default: F], optional
            -d	parameter for genBlast: maximum allowed distance between HSPs
                within the same gene, a non-negative integer [default: 100000],
                optional
            -r	parameter for genBlast: number of ranks in the output,
                a positive integer, optional
            -e	parameter for blast: The e-value, [default: 1e-2],
                optional
            -c	parameter for genBlast: minimum percentage of query gene
                coverage in the output, between 0 and 1 (e.g. for 50%
                gene coverage, use "0.5"), optional
            -W	parameter for blast: Set word size, 0 means using blast default [default: 0],
                optional
            -scodon The number of base pairs to search for start codon within the region of HSP
                        group (inside the first HSP). If not specified, default is 15.
            -i	parameter for genBlastG: minimum intron length, optional.
                If not specified, the default value is 15.
            -x	parameter for genBlastG: minimum internal exon length, optional.
                If not specified, default is 20.
            -n	parameter for genBlastG: maximum number of splice sites per region, optional.
                If not specified, default is 20.
            -gff	output options: turn on GFF output
            -o	output filename, optional. If not specified, the output
                will be the same as the query filename with ".gblast"
                extension.
            -pid turn on final alignment PID computation (global alignment between predicted
                gene and query) in output.
            -softmask	With this option NCBI blast will create a masking library,
                you need to use it when blasting against a whole genome
    """
    logger.info("Running GenBlast on : %s", protein_file)

    genblast_cmd = [
        str(genblast_bin),
        "-p",
        "genblastg",
        "-q",
        str(protein_file),
        "-t",
        str(masked_genome),
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
        str(protein_file),
    ]

    logger.info(" ".join(genblast_cmd))
    # Using the child process termination as described here:
    # https://alexandra-zaharia.github.io/posts/kill-subprocess
    # -and-its-children-on-timeout-python/
    try:
        p = subprocess.Popen(# pylint:disable=consider-using-with
            genblast_cmd, start_new_session=True
        )
        p.wait(timeout=genblast_timeout)
    except subprocess.TimeoutExpired:
        logger.error("Timeout reached for file: %s \n", protein_file)
        subprocess.run(# pylint:disable=subprocess-run-check
            ["touch", (Path(f"{protein_file}.except"))]
        )
        os.killpg(os.getpgid(p.pid), signal.SIGTERM)


def _generate_genblast_gtf(genblast_dir: Path) -> None:
    """
    Collect output from geneblast and create the final gtf file
    genblast_dir: Working directory path.
    """

    output_file = genblast_dir / "annotation.gtf"
    with open(output_file, "w+", encoding="utf8") as file_out:
        genblast_extension = "_1.1c_2.3_s1_0_16_1"
        for path in genblast_dir.rglob("*"):
            # for root, dirs, files in os.walk(genblast_dir):
            # for genblast_file in files:
            # genblast_file = os.path.join(root, genblast_file)
            if path.is_file() and path.suffix == ".gff":
                gtf_string = _convert_genblast_gff_to_gtf(path)
                file_out.write(gtf_string)
            elif path.is_file() and path.suffix in (
                ".fa.blast",
                ".fa.blast.report",
                genblast_extension,
            ):
                path.unlink()



def _run_convert2blastmask(
    convert2blastmask_bin: Path, masked_genome: Path, asnb_file: Path
) -> None:
    """
    Convert masking information in lower-case masked FASTA input to file
    formats suitable for makeblastdb.
    convert2blastmask_bin : Software path.
    masked_genome: Path of masked genome file.
    asnb_file: Path of assembly file.
    """
    logger.info("Running convert2blastmask prior to GenBlast:")
    cmd = [
        str(convert2blastmask_bin),
        "-in",
        str(masked_genome),
        "-parse_seqids",
        "-masking_algorithm",  # mask_program_name
        "other",
        "-masking_options",  # mask_program_options
        '"REpeatDetector, default"',
        "-outfmt",  # output_format
        "maskinfo_asn1_bin",
        "-out",
        str(asnb_file),
    ]
    logger.info(" ".join(cmd))
    subprocess.run(cmd, check=True)
    logger.info("Completed running convert2blastmask")


def _run_makeblastdb(makeblastdb_bin: Path, masked_genome: Path, asnb_file: Path) -> None:
    """
    Application to create BLAST databases.
    makeblastdb_bin : Software path.
    masked_genome: Path of masked genome file.
    asnb_file: Path of assembly file.
    """
    logger.info("Running makeblastdb prior to GenBlast")
    subprocess.run(  # pylint:disable=subprocess-run-check
        [
            str(makeblastdb_bin),
            "-in",
            str(masked_genome),
            "-dbtype",  # molecule_type
            "nucl",
            "-parse_seqids",
            "-mask_data",
            str(asnb_file),
            "-max_file_sz",  # number_of_bytes
            "10000000000",
        ]
    )
    logger.info("Completed running makeblastdb")


def _convert_genblast_gff_to_gtf(gff_file: Path) -> str:
    """
    Convert the content of gtf file in gff format
    gff_file: Path for the gff file
    """
    gtf_string = ""
    with open(gff_file, "r", encoding="utf8") as file_in:
        for line in file_in:
            results = line.split()
            if len(results) == 9:
                results[2] = "exon" if results[2] == "coding_exon" else results[2]
                attributes = _set_genblast_attributes(str(results[8]), str(results[2]))
                results[8] = attributes
                converted_line = "\t".join(results)
                gtf_string += converted_line + "\n"
    return gtf_string


def _set_genblast_attributes(attributes: str, feature_type: str) -> str:
    """
    Given the list of attributes in the genblast output,
    define the new attributes for the gtf file.
    attributes: GenBlast attribute list
    feature_type: transcript or exon
    Example genBlast output #pylint: disable=line-too-long, trailing-whitespace
    1       genBlastG       transcript      131128674       131137049       252.729 -       .       ID=259447-R1-1-A1;Name=259447;PID=84.65;Coverage=94.22;Note=PID:84.65-Cover:94.22
    1       genBlastG       coding_exon     131137031       131137049       .       -       .       ID=259447-R1-1-A1-E1;Parent=259447-R1-1-A1
    1       genBlastG       coding_exon     131136260       131136333       .       -       .       ID=259447-R1-1-A1-E2;Parent=259447-R1-1-A1
    1       genBlastG       coding_exon     131128674       131130245       .       -       .       ID=259447-R1-1-A1-E3;Parent=259447-R1-1-A1
    """
    converted_attributes = ""
    split_attributes = attributes.split(";")
    if feature_type == "transcript":
        match = re.search(r"Name\=(.+)$", split_attributes[1])
        assert match
        name = match.group(1)
        converted_attributes = f'gene_id "{name}"; transcript_id "{name}";'
    elif feature_type == "exon":
        match = re.search(r"\-E(\d+);Parent\=(.+)\-R\d+\-\d+\-", attributes)
        assert match
        exon_rank = match.group(1)
        name = match.group(2)
        converted_attributes = (
            f'gene_id "{name}"; transcript_id "{name}"; exon_number "{exon_rank}";'
        )

    return converted_attributes


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Genblast arguments")
    parser.add_argument(
        "--masked_genome_file", required=True, help="Masked genome file path"
    )
    parser.add_argument("--output_dir", required=True, help="Output directory path")
    parser.add_argument("--protein_file", required=True, help="Path for the protein dataset")
    parser.add_argument(
        "--genblast_timeout_secs", type=int, default=10800, help="Genblast timeout period"
    )
    parser.add_argument(
        "--max_intron_length", type=int, required=True, help="Maximum intron length"
    )
    parser.add_argument(
        "--genblast_bin",
        default="genblast",
        help="Genblast executable path",
    )
    parser.add_argument(
        "--convert2blastmask_bin",
        default="convert2blastmask",
        help="convert2blastmask executable path",
    )
    parser.add_argument(
        "--makeblastdb_bin",
        default="makeblastdb",
        help="makeblastdb executable path",
    )
    parser.add_argument("--num_threads", type=int, default=1, help="Number of threads")
    parser.add_argument(
        "--protein_set",
        required=True,
        choices=["uniprot", "orthodb"],
        help="Protein set [uniprot, orthodb]",
    )
    return parser.parse_args()

def main():
    """Genblast's entry-point."""
    args = parse_args()

    log_file_path = create_dir(args.output_dir, "log") / "genblast.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"

    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )

    run_genblast(
        Path(args.masked_genome_file),
        Path(args.output_dir),
        Path(args.protein_file),
        args.max_intron_length,
        args.genblast_timeout_secs,
        Path(args.genblast_bin),
        Path(args.convert2blastmask_bin),
        Path(args.makeblastdb_bin),
        args.num_threads,
        args.protein_set,
    )

if __name__ == "__main__":
    main()
