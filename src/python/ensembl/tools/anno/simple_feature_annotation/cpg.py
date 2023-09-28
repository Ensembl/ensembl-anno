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
Set of discriminant functions that can recognize structural and compositional features
such as CpG islands, promoter regions and first splice-donor sites.
Davuluri RV, Grosse I, Zhang MQ: Computational identification of promoters and
first exons in the human genome. Nat Genet. 2001, 29(4):412-417. [PMID: 11726928]
"""
__all__ = ["run_cpg"]
import logging
import logging.config
import multiprocessing
from os import PathLike
from pathlib import Path
import re
import subprocess
from tempfile import TemporaryDirectory
from typing import List,Union
import argschema

from ensembl.tools.anno.utils._utils import (
    check_exe,
    create_dir,
    check_gtf_content,
    get_seq_region_length,
    get_slice_id,
    slice_output_to_gtf,
    get_sequence,
)

logger = logging.getLogger(__name__)


def run_cpg(
    genome_file: PathLike,
    output_dir: Path,
    cpg_bin: Path = Path("cpg_lh"),
    cpg_min_length: int = 400,
    cpg_min_gc_content: int = 50,
    cpg_min_oe: float = 0.6,
    num_threads: int = 1,
) -> None:
    """
    Run CpG islands on genomic slices
    Args:
        :param genome_file : Genome file path.
        :type genome_file: PathLike
        :param output_dir : Working directory path
        :param cpg_bin : CpG software path.
        :param cpg_min_length : Min length of CpG islands
        :param cpg_min_gc_content : Min GC frequency percentage
        :param cpg_min_oe :  Min ratio of the observed to expected number of CpG (CpGo/e)
        :param num_threads: int, number of threads.
        
        :return: None
        :rtype: None
    """

    check_exe(cpg_bin)
    cpg_dir = create_dir(output_dir, "cpg_output")
    output_file = cpg_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "simple_feature")
        if transcript_count > 0:
            logger.info("Cpg gtf file exists")
            return
    logger.info("Creating list of genomic slices")
    seq_region_to_length = get_seq_region_length(genome_file, 5000)
    slice_ids_per_region = get_slice_id(
        seq_region_to_length, slice_size=1000000, overlap=0, min_length=5000
    )
    logger.info("Running CpG")
    pool = multiprocessing.Pool(int(num_threads))  # pylint:disable=consider-using-with
    for slice_id in slice_ids_per_region:
        pool.apply_async(
            _multiprocess_cpg,
            args=(
                cpg_bin,
                slice_id,
                genome_file,
                cpg_dir,
                cpg_min_length,
                cpg_min_gc_content,
                cpg_min_oe,
            ),
        )

    pool.close()
    pool.join()
    slice_output_to_gtf(cpg_dir, "feature_id", "cpg", True, ".cpg.gtf")
    for gtf_file in cpg_dir.glob("*.cpg.gtf"):
        gtf_file.unlink()


def _multiprocess_cpg(
    cpg_bin: Path,
    slice_id: List[str],
    genome_file: Path,
    cpg_dir: Path,
    cpg_min_length: int = 400,
    cpg_min_gc_content: int = 50,
    cpg_min_oe: float = 0.6,
) -> None:
    """
    Annotation of CpG islands on multiprocess on genomic slices
    Args:
        cpg_bin: CpG software path.
        slice_id: Slice id to run CpG on.
        genome_file : Genome file.
        cpg_dir : Output dir.
        cpg_min_length : Min length of CpG islands
        cpg_min_gc_content : Min GC frequency percentage
        cpg_min_oe :  Min ratio of the observed to expected number of CpG (CpGo/e)
    """
    region_name, start, end = slice_id
    logger.info(
        "Processing slice to find CpG islands with cpg_lh: %s:%s:%s",
        region_name,
        start,
        end,
    )
    seq = get_sequence(region_name, int(start), int(end), 1, genome_file, cpg_dir)
    slice_name = f"{region_name}.rs{start}.re{end}"
    #with TemporaryDirectory(dir=cpg_dir) as tmpdirname:
    slice_file = cpg_dir / f"{slice_name}.fa"
    with open(slice_file, "w+", encoding="utf8") as region_out:
        region_out.write(f">{region_name}\n{seq}\n")
    region_results = cpg_dir / f"{slice_file}.cpg.gtf"
    output_file = Path(f"{slice_file}.cpg")
    cpg_cmd = [str(cpg_bin), str(slice_file)]
    with open(output_file, "w+", encoding="utf8") as cpg_out:
        subprocess.run(cpg_cmd, stdout=cpg_out, check=True)
        _create_cpg_gtf(
            output_file,
            region_results,
            region_name,
            cpg_min_length,
            cpg_min_gc_content,
            cpg_min_oe,
        )
    slice_file.unlink()
    output_file.unlink()


def _create_cpg_gtf(
    output_file: Path,
    region_results: Path,
    region_name: str,
    cpg_min_length: int = 400,
    cpg_min_gc_content: int = 50,
    cpg_min_oe: float = 0.6,
) -> None:
    """
    Read the fasta file and save the content in gtf format
    All the genomic slices are collected in a single gtf output
    Args:
        output_file : GTF file with final results.
        region_results : GTF file with the results per region.
        region_name :Coordinates of genomic slice.
        cpg_dir : Output dir.
        cpg_min_length : Min length of CpG islands
        cpg_min_gc_content : Min GC frequency percentage
        cpg_min_oe :  Min ratio of the observed to expected number of CpG (CpGo/e)
    """
    with open(output_file, "r", encoding="utf8") as cpg_in, open(region_results, "w+", encoding="utf8") as cpg_out:
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
                oe_score_str = results[7]
                oe_score: Union[float, int]
                if oe_score_str in ("-", "inf"):
                    oe_score=0
                else:
                    oe_score=float(oe_score_str)
                if (
                    int(length) >= int(cpg_min_length)
                    and gc_content >= int(cpg_min_gc_content)
                    and oe_score >= float(cpg_min_oe)
                ):
                    gtf_line = (
                        f"{region_name}\tCpG\tsimple_feature\t{start}\t"
                        f'{end}\t.\t+\t.\tfeature_id "{feature_count}"; score "{score}";\n'
                    )
                    cpg_out.write(gtf_line)


class InputSchema(argschema.ArgSchema):
    """Input arguments expected to run CpG software."""

    genome_file = argschema.fields.InputFile(
        required=True, description="Genome file path"
    )
    output_dir = argschema.fields.OutputDir(
        required=True, description="Output directory path"
    )
    cpg_bin = argschema.fields.String(
        required=False,
        default="cpg_lh",
        description="CpG executable path",
    )
    cpg_min_length = argschema.fields.Integer(
        required=False,
        default="400",
        description="Min length of CpG islands",
    )
    cpg_min_gc_content = argschema.fields.Integer(
        required=False,
        default="50",
        description="Min GC frequency percentage",
    )
    cpg_min_oe = argschema.fields.Float(
        required=False,
        default="0.6",
        description="Min ratio of the observed to expected number of CpG (CpGo/e)",
    )
    num_threads = argschema.fields.Integer(
        required=False, default=1, description="Number of threads"
    )


def main() -> None:
    """CpG's entry-point."""
    mod = argschema.ArgSchemaParser(schema_type=InputSchema)
    log_file_path = create_dir(mod.args["output_dir"], "log") / "cpg.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"
    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )
    run_cpg(
        mod.args["genome_file"],
        mod.args["output_dir"],
        mod.args["cpg_bin"],
        mod.args["cpg_min_length"],
        mod.args["cpg_min_gc_content"],
        mod.args["cpg_min_oe"],
        mod.args["num_threads"],
    )
