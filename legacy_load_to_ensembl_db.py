# See the NOTICE file distributed with this work for additional information #pylint: disable=missing-module-docstring
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
Script to load annotation results into an Ensembl database.
"""

import logging
import re
import subprocess
import tempfile
from multiprocessing import Pool
from pathlib import Path
from typing import Any

from src.python.ensembl.tools.anno.utils import _utils

logger = logging.getLogger(__name__)


def _load_gtf_to_db(  # pylint: disable=too-many-arguments, too-many-locals
    gtf_file: Path,
    description: str,
    load_type: str,
    analysis_name: str,
    batch_size: int,
    *,
    make_single_transcript_genes: bool,
    db_loading_script: Path,
    genome_file: Path,
    db_details: dict[str, Any],
    db_loading_dir: Path,
    run_repeatmasker_analysis: int,
    num_threads: int,
) -> None:
    """Load a GTF file into the Ensembl database."""

    if not gtf_file.exists():
        logger.error(
            "Did not find %s annotation file, so not loading. Path checked: %s",
            description,
            gtf_file,
        )
        return

    logger.info("Loading %s to db", description)

    gtf_records = batch_gtf_records(  # pylint: disable=too-many-function-args
        gtf_file,
        batch_size,
        load_type,
    )

    generic_load_records_to_ensembl_db(
        make_single_transcript_genes,
        run_repeatmasker_analysis,
        db_loading_script,
        genome_file,
        db_details,
        db_loading_dir,
        load_type,
        analysis_name,
        gtf_records,
        num_threads,
    )


def load_results_to_ensembl_db(  # pylint: disable=too-many-arguments, too-many-positional-arguments, too-many-locals
    main_script_dir: Path,
    make_single_transcript_genes: bool,
    genome_file: Path,
    main_output_dir: Path,
    db_details: dict[str, Any],
    num_threads: int,
    repeatmasker_analysis: str,
) -> None:
    """
    Load annotation results into an Ensembl database.
    Args:
        main_script_dir: Directory containing the main script.
        make_single_transcript_genes: Whether to load results to Ensembl database.
        genome_file: Path to the genome file.
        main_output_dir: Directory containing the output of annotation analyses.
        db_details: Database connection details.
        num_threads: Number of threads to use for loading.
        repeatmasker_analysis: Name of the RepeatMasker analysis to use in the database.
    """
    db_loading_script = main_script_dir / "support_scripts_perl" / "load_gtf_ensembl.pl"

    db_loading_dir = _utils.create_dir(
        main_output_dir,
        "db_loading",
    )
    analyses = (
        ("annotation_output", "main geneset", "gene", "ensembl", 200, 0),
        ("rfam_output", "Rfam-based sncRNA genes", "gene", "ncrna", 500, 0),
        ("trnascan_output", "tRNAScan-SE tRNA genes", "gene", "ncrna", 500, 0),
        ("dust_output", "Dust repeats", "single_line_feature", "dust", 500, 0),
        ("red_output", "Red repeats", "single_line_feature", "repeatdetector", 500, 0),
        ("trf_output", "TRF repeats", "single_line_feature", "trf", 500, 0),
        (
            "repeatmasker_output",
            "RepeatMasker repeats",
            "single_line_feature",
            repeatmasker_analysis,
            500,
            1,
        ),  # pylint: disable=line-too-long
        ("cpg_output", "CpG islands", "single_line_feature", "cpg", 500, 0),
        ("eponine_output", "Eponine features", "single_line_feature", "eponine", 500, 0),
    )

    for (
        output_dir,
        description,
        load_type,
        analysis_name,
        batch_size,
        run_repeatmasker_analysis,
    ) in analyses:  # pylint: disable=line-too-long
        _load_gtf_to_db(
            gtf_file=(main_output_dir / output_dir / "annotation.gtf"),
            description=description,
            load_type=load_type,
            analysis_name=analysis_name,
            batch_size=batch_size,
            make_single_transcript_genes=make_single_transcript_genes,
            db_loading_script=db_loading_script,
            genome_file=genome_file,
            db_details=db_details,
            db_loading_dir=db_loading_dir,
            num_threads=num_threads,
            run_repeatmasker_analysis=run_repeatmasker_analysis,
        )

    logger.info("Finished loading records to db")


def generic_load_records_to_ensembl_db(  # pylint: disable=too-many-arguments, too-many-positional-arguments, too-many-locals
    make_single_transcript_genes: bool,
    run_repeatmasker_analysis: int,
    db_loading_script: Path,
    genome_file: Path,
    db_details: dict[str, Any],
    db_loading_dir: Path,
    load_type: str,
    analysis_name: str,
    gtf_records: list[list[str]],
    num_threads: int,
) -> None:
    """Load batches of GTF records in parallel.
    Args:
        make_single_transcript_genes: Whether to load results as single transcript genes.
        db_loading_script: Path to the database loading script.
        genome_file: Path to the genome file.
        db_details: Database connection details.
        db_loading_dir: Directory to use for temporary files during loading.
        load_type: Type of features being loaded (e.g., "gene", "single_line_feature").
        analysis_name: Name of the analysis to use in the database.
        gtf_records: List of batches of GTF records to load.
        num_threads: Number of threads to use for loading.
    """

    with Pool(processes=num_threads) as pool:
        for record_batch in gtf_records:
            pool.apply_async(
                multiprocess_load_records_to_ensembl_db,
                args=(
                    make_single_transcript_genes,
                    db_loading_script,
                    genome_file,
                    db_details,
                    db_loading_dir,
                    load_type,
                    analysis_name,
                    run_repeatmasker_analysis,
                    record_batch,
                ),
            )

        pool.close()
        pool.join()

    logger.info("Finished loading records to db")


def multiprocess_load_records_to_ensembl_db(  # pylint: disable=too-many-arguments, too-many-branches, too-many-locals, too-many-positional-arguments
    make_single_transcript_genes: bool,
    db_loading_script: Path,
    genome_file: Path,
    db_details: str,
    output_dir: Path,
    load_type: str,
    analysis_name: str,
    run_repeatmasker_analysis: int,
    record_batch: list[str],
) -> None:
    """Load a single GTF batch into the Ensembl database."""

    with tempfile.NamedTemporaryFile(
        mode="w",
        delete=False,
        dir=output_dir,
        encoding="utf-8",
    ) as gtf_temp_out:
        gtf_temp_out.writelines(record_batch)
        gtf_temp_file_path = Path(gtf_temp_out.name)

    db_name, db_host, db_port, db_user, db_pass = db_details.split(",")

    loading_cmd: list[str] = [
        "perl",
        str(db_loading_script),
        "-genome_file",
        str(genome_file),
        "-dbname",
        db_name,
        "-host",
        db_host,
        "-port",
        db_port,
        "-user",
        db_user,
        "-pass",
        db_pass,
        "-gtf_file",
        str(gtf_temp_file_path),
        "-analysis_name",
        analysis_name,
        "-load_type",
        load_type,
        "-run_repeatmasker_analysis",
        str(run_repeatmasker_analysis),
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

        if make_single_transcript_genes is True:
            loading_cmd.append(
                "-make_single_transcript_genes",
            )

    logger.info(" ".join(loading_cmd))

    try:
        subprocess.run(
            loading_cmd,
            check=True,
        )
    finally:
        gtf_temp_file_path.unlink(missing_ok=True)

    logger.info(
        "Finished: %s",
        gtf_temp_file_path,
    )


def batch_gtf_records(  # pylint: disable=too-many-branches, too-many-arguments
    input_gtf_file: Path,
    batch_size: int,
    record_type: str,
) -> list[list[str]]:
    """
    Batch GTF records for database loading.

    Args:
        input_gtf_file: Input GTF file.
        batch_size: Number of records per batch.
        record_type: Either "gene" or "single_line_feature".

    Returns:
        List of batches containing GTF lines.
    """
    records: list[list[str]] = []
    current_batch: list[str] = []

    with input_gtf_file.open(encoding="utf-8") as gtf_in:
        if record_type == "gene":
            current_gene_id = ""
            record_counter = 0
            current_record_batch: list[str] = []

            for line in gtf_in:
                if line.startswith("#"):
                    continue
                if len(line.split("\t")) != 9:
                    continue
                match = re.search(r'gene_id "([^"]+)"', line)
                if not match:
                    logger.warning("Could not find gene_id in line: %s", line)
                    continue
                gene_id = match.group(1)

                if not current_gene_id:
                    record_counter += 1
                    current_gene_id = gene_id

                elif gene_id != current_gene_id:
                    record_counter += 1
                    if record_counter % batch_size == 0 and current_record_batch:
                        records.append(current_record_batch)
                        current_record_batch = []
                    current_gene_id = gene_id
                current_record_batch.append(line)
        elif record_type == "single_line_feature":
            record_counter = 0
            current_record_batch = []
            current_gene_id = ""
            for line in gtf_in:
                if line.startswith("#"):
                    continue
                if len(line.split("\t")) != 9:
                    continue
                record_counter += 1
                if record_counter % batch_size == 0 and current_record_batch:
                    records.append(current_record_batch)
                    current_record_batch = []
                current_record_batch.append(line)
            records.append(current_record_batch)
        else:
            raise ValueError(f"Unsupported record_type: {record_type}")

    if current_batch:
        records.append(current_batch)

    return records
