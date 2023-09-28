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
tRNAscan-SE identifies 99-100% of transfer RNA genes in DNA sequence while
giving less than one false positive per 15 gigabases.
Lowe TM, Eddy SR: tRNAscan-SE: a program for improved detection of transfer
RNA genes in genomic sequence.
Nucleic Acids Res. 1997, 25(5):955-64. [PMID: 9023104]
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


def run_trnascan(
    genome_file: PathLike,
    output_dir: Path,
    trnascan_bin: Path = Path("tRNAscan-SE"),
    trnascan_filter: Path = Path("EukHighConfidenceFilter"),
    num_threads: int = 1,
) -> None:
    """
    Executes tRNAscan-SE on genomic slices
        :param genome_file: Genome file path.
        :type genome_file: PathLike 
        :param output_dir:  working directory path.
        :type output_dir: Path  
        :param trnascan_bin: tRNAscan-SE software path.
        :type trnascan_bin: Path, default tRNAscan-SE
        :param trnascan_filter: tRNAscan-SE filter set path.
        :type trnascan_filter: Path, default EukHighConfidenceFilter
        :param num_threads: int, number of threads.
        :type num_threads: int, default 1 
                            
        :return: None
        :rtype: None
    """
    check_exe(trnascan_bin)
    check_file(trnascan_filter)
    trnascan_dir = create_dir(output_dir, "trnascan_output")
    output_file = trnascan_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Trnascan gtf file exists, skipping analysis")
            return
    logger.info("Creating list of genomic slices")
    seq_region_to_length = get_seq_region_length(genome_file, 5000)
    slice_ids_per_region = get_slice_id(seq_region_to_length, 1000000, 0, 5000)
    trnascan_cmd = [
        str(trnascan_bin),
        None,
        "-o",
        None,
        "-f",
        None,
        "-H",  # show both primary and secondary structure components to covariance model bit scores
        "-q",  # quiet mode
        "--detail",
        "-Q",
    ]
    logger.info("Running tRNAscan-SE")
    pool = multiprocessing.Pool(num_threads)  # pylint: disable=consider-using-with
    for slice_id in slice_ids_per_region:
        pool.apply_async(
            _multiprocess_trnascan,
            args=(
                trnascan_cmd,
                slice_id,
                genome_file,
                trnascan_filter,
                trnascan_dir,
            ),
        )

    pool.close()
    pool.join()
    slice_output_to_gtf(
        output_dir=trnascan_dir, unique_ids=True, file_extension=".trna.gtf"
    )
    for gtf_file in trnascan_dir.glob("*.trna.gtf"):
        gtf_file.unlink()


def _multiprocess_trnascan(
    trnascan_cmd: List[str],
    slice_id: List[str],
    genome_file: Path,
    trnascan_filter: Path,
    trnascan_dir: Path,
) -> None:
    """
    Run tRNAscan-SE on multiprocess on genomic slices
    Args:
        trnascan_cmd: tRNAscan-SE command to execute.
        slice_id: Slice Id to run tRNAscan-SE on.
        genome_file : Genome file.
        trnascan_dir : tRNAscan-SE output dir.
        trnascan_filter: tRNAscan-SE filter set.
    """
    region_name, start, end = slice_id
    logger.info(
        "Processing slice to find tRNAs using tRNAscan-SE:%s:%s:%s",
        region_name,
        start,
        end,
    )
    seq = get_sequence(region_name, int(start), int(end), 1, genome_file, trnascan_dir)
    slice_name = f"{region_name}.rs{start}.re{end}"
    slice_file = trnascan_dir / f"{slice_name}.fa"
    with open(slice_file, "w+", encoding="utf8") as region_out:
        region_out.write(f">{region_name}\n{seq}\n")
    # trnscan output
    region_results = trnascan_dir / f"{slice_name}.trna.gtf"
    output_file = Path(f"{slice_file}.trna")
    ss_output_file = Path(f"{output_file}.ss")
    # filtering
    filter_prefix_file = f"{slice_name}.filt"
    filter_output_file = trnascan_dir / f"{filter_prefix_file}.out"
    filter_log_file = trnascan_dir / f"{filter_prefix_file}.log"
    filter_ss_file = trnascan_dir / f"{filter_prefix_file}.ss"
    # trnascan_cmd = generic_trnascan_cmd.copy()
    trnascan_cmd[1], trnascan_cmd[3], trnascan_cmd[5] = (
        str(slice_file),
        str(output_file),
        str(ss_output_file),
    )
    logger.info("tRNAscan-SE command: %s", " ".join(trnascan_cmd))
    subprocess.run(trnascan_cmd, check=True)
    # If the trnascan output is empty there is no need to go on with filtering
    if output_file.stat().st_size == 0:
        output_file.unlink()
        slice_file.unlink()
        ss_output_file.unlink(missing_ok=True)
        return

    filter_cmd = [
        str(trnascan_filter),
        "--result",  # tRNAscan-SE output file used as input
        str(output_file),
        "--ss",  # tRNAscan-SE secondary structure file used as input
        str(ss_output_file),
        "--output",
        str(trnascan_dir),
        "--prefix",
        str(filter_prefix_file),
    ]
    logger.info(
        "tRNAscan-SE filter command: %s", " ".join(str(item) for item in filter_cmd)
    )
    subprocess.run(filter_cmd)#pylint:disable=subprocess-run-check
    _create_trnascan_gtf(region_results, filter_output_file, region_name)
    output_file.unlink(missing_ok=True)
    slice_file.unlink(missing_ok=True)
    ss_output_file.unlink(missing_ok=True)
    Path(filter_prefix_file).unlink(missing_ok=True)
    filter_log_file.unlink(missing_ok=True)
    filter_ss_file.unlink(missing_ok=True)
    filter_output_file.unlink(missing_ok=True)


def _create_trnascan_gtf(
    region_results: Path, filter_output_file: Path, region_name: str
) -> None:
    """
    Read the fasta file and save the content in gtf format
    All the genomic slices are collected in a single gtf output
    Args:
        region_results : GTF file with the results per region.
        filter_file : GTF file with the filtered results per region.
        region_name :Coordinates of genomic slice.

    tRNAscan-SE output format:
    col0: GtRNAdb Gene Symbol - gene ID in corresponding genome
    col1: tRNAscan-SE ID - tRNA ID in tRNAscan-SE prediction results
    col2-3: Locus - Genomic coordinates of predicted gene
    col4: Isotype (from Anticodon) - tRNA isotype determined by anticodon
    col5: Anticodon - anticodon of predicted tRNA gene
    col6-7: Intron boundaries
    col8: General tRNA Model Score - covariance model bit score from tRNAscan-SE results
    col9: Best Isotype Model - best matching (highest scoring) isotype determined
    by isotype-specific covariance model classification
    col10-11-12: Anticodon and Isotype Model Agreement - consistency between anticodon
    from predicted gene sequence and best isotype model
    col13: Features - special gene features that may include gene set categorization,
    number of introns, possible pseudogenes, possible truncation, or base-pair mismatches
    """
    with open(filter_output_file, "r", encoding="utf8") as trna_in, open(
        region_results, "w+", encoding="utf8"
    ) as trna_out:
        gene_counter = 1
        for line in trna_in:
            result_match = re.search(r"^" + region_name, line)
            if result_match:
                results = line.split()
                start = int(results[2])
                end = int(results[3])
                strand = "+"
                if start > end:
                    strand = "-"
                    start, end = end, start
                biotype = (
                    "tRNA"
                    if re.search(r"high confidence set", line)
                    else "tRNA_pseudogene"
                )
                transcript_string = (
                    f"{region_name}\ttRNAscan\ttranscript\t{start}\t{end}\t.\t"
                    f'{strand}\t.\tgene_id "{gene_counter}"; transcript_id '
                    f'"{gene_counter}"; biotype "{biotype}";\n'
                )
                exon_string = (
                    f"{region_name}\ttRNAscan\texon\t{start}\t{end}\t.\t"
                    f'{strand}\t.\tgene_id "{gene_counter}"; transcript_id '
                    f'"{gene_counter}"; exon_number "1"; biotype "{biotype}";\n'
                )
                trna_out.write(transcript_string)
                trna_out.write(exon_string)
                trna_out.flush()
                gene_counter += 1


class InputSchema(argschema.ArgSchema):
    """Input arguments expected to run tRNAscan-SE."""

    genome_file = argschema.fields.InputFile(
        required=True, description="Genome file path"
    )
    trnascan_bin = argschema.fields.String(
        required=False,
        default="tRNAscan-SE",
        description="tRNAscan-SE executable path",
    )
    trnascan_filter = argschema.fields.String(
        required=False,
        default="/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/bin/EukHighConfidenceFilter",
        description="tRNAscan-SE filter path",
    )
    output_dir = argschema.fields.OutputDir(
        required=True, description="Output directory path"
    )
    num_threads = argschema.fields.Integer(
        required=False, default=1, description="Number of threads"
    )


def main() -> None:
    """tRNAscan-SE's entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    log_file_path = create_dir(mod.args["output_dir"], "log") / "trnascan.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"
    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )
    run_trnascan(
        mod.args["genome_file"],
        mod.args["output_dir"],
        mod.args["trnascan_bin"],
        Path(mod.args["trnascan_filter"]),
        mod.args["num_threads"],
    )


if __name__ == "__main__":
    main()
