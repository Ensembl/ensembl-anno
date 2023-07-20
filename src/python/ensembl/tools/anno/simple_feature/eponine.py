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
Eponine is a probabilistic method for detecting transcription start sites (TSS)
in mammalian genomic sequence, with good specificity and excellent positional accuracy.
Down TA, Hubbard TJ. Computational detection and location of transcription start sites
in mammalian genomic DNA. Genome Res. 2002 Mar;12(3):458-61. doi: 10.1101/gr.216102.
PMID: 11875034; PMCID: PMC155284.
"""
__all__ = ["run_eponine"]

import logging
import logging.config
import multiprocessing
import os
from os import PathLike
from pathlib import Path
import re
import subprocess
from tempfile import TemporaryDirectory
from typing import List
import argschema

from ensembl.tools.anno.utils._utils import (
    check_exe,
    check_file,
    create_dir,
    check_gtf_content,
    get_sequence,
    get_seq_region_length,
    get_slice_id,
    slice_output_to_gtf,
)

logger = logging.getLogger("__name__")


def run_eponine(
    genome_file: PathLike,
    output_dir: Path,
    num_threads: int = 1,
    java_bin: Path = Path("java"),
    eponine_bin: Path = Path(
        "/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/opt/eponine/libexec/eponine-scan.jar"
    ),
    eponine_threshold: float = 0.999,
) -> None:
    """
    Run Eponine on genomic slices
    Args:
        genome_file : Genome file path.
        output_dir : Working directory path.
        java_bin : Java path.
        eponine_bin : Eponine software path
        num_threads: Number of threads.
    """
    check_file(eponine_bin)
    check_exe(java_bin)
    eponine_dir = create_dir(output_dir, "eponine_output")
    # os.chdir(str(eponine_dir))
    output_file = eponine_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "simple_feature")
        if transcript_count > 0:
            logger.info("Eponine gtf file exists, skipping analysis")
            return
    logger.info("Creating list of genomic slices")
    seq_region_to_length = get_seq_region_length(genome_file, 5000)
    slice_ids_per_region = get_slice_id(
        seq_region_to_length, slice_size=1000000, overlap=0, min_length=5000
    )

    eponine_cmd = [
        str(java_bin),
        "-jar",
        str(eponine_bin),
        "-threshold",
        str(eponine_threshold),
        "-seq",
    ]
    logger.info("Running Eponine")
    pool = multiprocessing.Pool(int(num_threads))  # pylint:disable=consider-using-with
    for slice_id in slice_ids_per_region:
        pool.apply_async(
            _multiprocess_eponine,
            args=(
                eponine_cmd,
                slice_id,
                eponine_dir,
                Path(genome_file),
            ),
        )
    pool.close()
    pool.join()
    slice_output_to_gtf(eponine_dir, "feature_id", "eponine", True, ".epo.gtf")
    for gtf_file in eponine_dir.glob("*.epo.gtf"):
        gtf_file.unlink()


def _multiprocess_eponine(
    eponine_cmd: List[str],
    slice_id: List[str],
    eponine_dir: Path,
    genome_file: Path,
) -> None:
    """
    Run Eponine on multiprocess on genomic slices
    Args:
        eponine_cmd: Eponine command to execute.
        slice_id: List of slice IDs.
        eponine_dir : Eponine output directory path.
        genome_file : Genome file.
    """
    region_name, start, end = slice_id
    logger.info(
        "Processing slice to find transcription start sites with Eponine: %s:%s:%s",
        region_name,
        start,
        end,
    )
    seq = get_sequence(region_name, int(start), int(end), 1, genome_file, eponine_dir)
    slice_name = f"{region_name}.rs{start}.re{end}"
    #with tempfile.TemporaryDirectory(dir=eponine_dir) as tmpdirname:
    slice_file = eponine_dir / f"{slice_name}.fa"
    with open(slice_file, "w+", encoding="utf8") as region_out:
        region_out.write(f">{region_name}\n{seq}\n")
    region_results = eponine_dir / f"{slice_name}.epo.gtf"
    output_file = Path(f"{slice_file}.epo")
    eponine_cmd = eponine_cmd.copy()
    eponine_cmd.append(str(slice_file))
    logging.info(eponine_cmd)
    with open(output_file, "w+") as eponine_out:
        subprocess.run(eponine_cmd, stdout=eponine_out, check=True)
    _create_eponine_gtf(output_file, region_results, region_name)
    slice_file.unlink()
    output_file.unlink()


def _create_eponine_gtf(
    output_file: Path,
    region_results: Path,
    region_name: str,
) -> None:
    """
    Read the fasta file and save the content in gtf format
    All the genomic slices are collected in a single gtf output
    Args:
        output_file: GTF file with final results.
        region_results: GTF file with the results per region.
        region_name: Coordinates of genomic slice.
    """
    with open(output_file, "r", encoding="utf8") as eponine_in, open(
        region_results, "w+", encoding="utf8"
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
                logging.info(results)
                # There's a one base offset on the reverse strand
                if strand == "-":
                    start -= 1
                    end -= 1

                gtf_line = (
                    f"{region_name}\tEponine\tsimple_feature\t"
                    f"{start}\t{end}\t.\t{strand}\t.\t"
                    f'feature_id "{feature_count}"; score "{score}";\n'
                )
                eponine_out.write(gtf_line)
                feature_count += 1


class InputSchema(argschema.ArgSchema):
    """Input arguments expected to run Eponine."""

    genome_file = argschema.fields.InputFile(
        required=True, description="Genome file path"
    )
    output_dir = argschema.fields.OutputDir(
        required=True, description="Output directory path"
    )
    num_threads = argschema.fields.Integer(
        required=False, default=1, description="Number of threads"
    )
    java_bin = argschema.fields.String(
        required=False,
        default="java",
        description="Java executable path",
    )
    eponine_bin = argschema.fields.String(
        required=False,
        default="/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/opt/eponine/libexec/eponine-scan.jar",  # pylint:disable=line-too-long
        description="Java executable path",
    )
    eponine_threashold = argschema.fields.Float(
        required=False, default=0.999, description="Eponine threashold"
    )


def main() -> None:
    """Eponine's entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    log_file_path = create_dir(mod.args["output_dir"], "log") / "eponine.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"
    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )
    run_eponine(
        mod.args["genome_file"],
        mod.args["output_dir"],
        mod.args["num_threads"],
        Path(mod.args["java_bin"]),
        Path(mod.args["eponine_bin"]),
        mod.args["eponine_threashold"],
    )


if __name__ == "__main__":
    main()
