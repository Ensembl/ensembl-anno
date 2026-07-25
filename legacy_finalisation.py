# pylint: disable=too-many-arguments, too-many-locals, too-many-branches, too-many-statements, line-too-long, too-many-lines, too-many-public-methods, too-many-instance-attributes, subprocess-run-check,consider-using-with
"""Legacy finalisation module"""

import logging
import json
import multiprocessing
from multiprocessing.pool import AsyncResult
import os
from pathlib import Path
import re
import subprocess
import shutil
from typing import Any, Dict, Union, cast

from src.python.ensembl.tools.anno.utils._utils import (
    check_file,
    create_dir,
    check_gtf_content,
    get_seq_region_lengths,
    split_protein_file,
)

# Use same config path layout as monolithic script


def load_json(path: Union[str, Path]) -> Dict[str, Any]:
    """Load a JSON file and return its contents as a dictionary.

    Args:
        path (Union[str, Path]): _description_

    Returns:
        Dict[str, Any]: _description_
    """
    path = Path(path)
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


config = load_json(Path(os.environ["ENSCODE"]) / "ensembl-anno" / "conf" / "config.json")


def file_ok(path: Union[str, Path], min_size: int = 1) -> bool:
    """Check if a file exists, is a file, and is above a minimum size.
    Args:
        path (Union[str, Path]): The path to the file to check.
        min_size (int, optional): The minimum size in bytes for the file. Defaults to 1.
    Returns:
        bool: True if the file exists, is a file, and is above the minimum size; False otherwise.
    """
    path = Path(path)
    return path.exists() and path.is_file() and path.stat().st_size >= min_size


def gtf_ok(path: Union[str, Path]) -> bool:
    """Check if a GTF file exists, is a file, and contains at least one transcript.
    Args:
        path (Union[str, Path]): The path to the GTF file to check.
    Returns:
        bool: True if the GTF file exists, is a file, and contains at
        least one transcript; False otherwise.
    """
    if not file_ok(path):
        return False
    try:
        return check_gtf_content(Path(path), "transcript") > 0
    except Exception:  # pylint: disable=broad-except
        return False


def skip_if_exists(desc: str, path: Union[str, Path], check_fn=file_ok) -> bool:
    """Skip a step if a file already exists and passes a given check function.
    Args:
        desc (str): A description of the file being checked, used for logging purposes.
        path (Union[str, Path]): The path to the file to check.
        check_fn (_type_, optional): The check function to use. Defaults to file_ok.

    Returns:
        bool: True if the file exists and passes the check function, False otherwise.
    """
    if check_fn(path):
        logging.info("%s exists (%s), skipping", desc, path)
        return True
    return False


def run_finalise_geneset(  # pylint: disable=too-many-arguments, too-many-positional-arguments
    main_script_dir: Path,
    main_output_dir: Path,
    genome_file: Path,
    seq_region_names: list[str],
    validation_type: str | None,
    diamond_validation_db: Path | None,
    num_threads: int,
) -> None:
    """
    Create the final annotation geneset.
    Asociated steps include:
    - Selecting the best transcripts for each evidence type (protein, transcriptomic,
    busco)
    - Merging the selected transcripts
    - Validating coding transcripts using RNAsamba, CPC2 and optionally Diamond
    - Cleaning the annotation and UTRs using custom scripts
    - Dumping transcript and translation sequences from the final annotation
    - The validation step can be run under different levels of strictness, controlled
    by the validation_type parameter. This mainly affects how single exon coding
    transcripts are validated, with more relaxed validation allowing more single exon
    coding transcripts to be retained.
    - The final annotation is written to a file called annotation.gtf in the final
    annotation directory.
    Args:
    main_script_dir: Path to the main script directory, used to locate support scripts.
        main_output_dir: Path to the main output directory, where intermediate results
        are stored.
        genome_file: Path to the genome file, used for various steps including validation
        and cleaning.
        seq_region_names: List of sequence region names to process.
        validation_type: The level of strictness for validating coding transcripts.
        Can be "relaxed", "moderate" or "strict". Affects mainly the validation of
        single exon coding transcripts.
        diamond_validation_db: Optional path to a Diamond database for validating
        coding transcripts. If not provided, Diamond validation will be skipped.
        num_threads: Number of threads to use for parallel processing.
    """

    validation_type = validation_type or "relaxed"
    logging.info(
        "Setting validation type to %s",
        validation_type,
    )

    final_annotation_dir = create_dir(
        main_output_dir,
        "annotation_output",
    )
    region_annotation_dir = create_dir(
        final_annotation_dir,
        "initial_region_gtfs",
    )
    final_region_annotation_dir = create_dir(
        final_annotation_dir,
        "final_region_gtfs",
    )
    utr_region_annotation_dir = create_dir(
        final_annotation_dir,
        "utr_region_gtfs",
    )
    validation_dir = create_dir(
        final_annotation_dir,
        "cds_validation",
    )

    seq_region_lengths = get_seq_region_lengths(
        genome_file,
        0,
    )

    output_file = final_annotation_dir / "annotation.gtf"

    logging.info(
        "Checking whether final annotation already exists",
    )

    if gtf_ok(output_file):
        logging.info(
            "Final annotation already exists, skipping",
        )
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logging.info("Final gtf file exists")
            return

    logging.info(
        "No final annotation found, continuing",
    )

    genblast_uniprot_annotation_raw = main_output_dir / "genblast_uniprot_output" / "annotation.gtf"

    genblast_orthodb_annotation_raw = main_output_dir / "genblast_orthodb_output" / "annotation.gtf"

    miniprot_uniprot_annotation_raw = main_output_dir / "miniprot_uniprot_output" / "annotation.gtf"

    miniprot_orthodb_annotation_raw = main_output_dir / "miniprot_orthodb_output" / "annotation.gtf"

    minimap2_annotation_raw = main_output_dir / "minimap2_output" / "annotation.gtf"

    stringtie_annotation_raw = main_output_dir / "stringtie_output" / "annotation.gtf"

    scallop_annotation_raw = main_output_dir / "scallop_output" / "annotation.gtf"

    transcript_selector_script = main_script_dir / "support_scripts_perl" / "select_best_transcripts.pl"

    finalise_geneset_script = main_script_dir / "support_scripts_perl" / "finalise_geneset.pl"

    clean_geneset_script = main_script_dir / "support_scripts_perl" / "clean_geneset.pl"

    clean_utrs_script = main_script_dir / "support_scripts_perl" / "clean_utrs_and_lncRNAs.pl"

    gtf_to_seq_script = main_script_dir / "support_scripts_perl" / "gtf_to_seq.pl"

    transcriptomic_annotation_raw = final_annotation_dir / "transcriptomic_raw.gtf"

    busco_annotation_raw = final_annotation_dir / "busco_raw.gtf"

    protein_annotation_raw = final_annotation_dir / "protein_raw.gtf"

    if not skip_if_exists(
        "Transcriptomic raw GTF",
        transcriptomic_annotation_raw,
    ):
        with transcriptomic_annotation_raw.open(
            "w",
            encoding="utf-8",
        ) as out_handle:
            for transcriptomic_file in (
                minimap2_annotation_raw,
                scallop_annotation_raw,
                stringtie_annotation_raw,
            ):
                if not transcriptomic_file.exists():
                    logging.info(
                        "Missing %s, skipping",
                        transcriptomic_file,
                    )
                    continue

                with transcriptomic_file.open(
                    encoding="utf-8",
                ) as in_handle:
                    shutil.copyfileobj(
                        in_handle,
                        out_handle,
                    )

    if genblast_uniprot_annotation_raw.exists() or miniprot_uniprot_annotation_raw.exists():
        # shutil.copy2(
        #    busco_annotation_raw,
        #    final_annotation_dir / "busco_raw.gtf",
        # )
        with busco_annotation_raw.open(
            "w",
            encoding="utf-8",
        ) as out_handle:
            for protein_file in (
                genblast_uniprot_annotation_raw,
                miniprot_uniprot_annotation_raw,
            ):
                if not protein_file.exists():
                    logging.info(
                        "Missing %s, skipping",
                        protein_file,
                    )
                    continue

                with protein_file.open(
                    encoding="utf-8",
                ) as in_handle:
                    shutil.copyfileobj(
                        in_handle,
                        out_handle,
                    )
    if genblast_orthodb_annotation_raw.exists() or miniprot_orthodb_annotation_raw.exists():
        # shutil.copy2(
        #    protein_annotation_raw,
        #    final_annotation_dir / "busco_raw.gtf",
        # )
        with protein_annotation_raw.open(
            "w",
            encoding="utf-8",
        ) as out_handle:
            for protein_file in (
                genblast_orthodb_annotation_raw,
                miniprot_orthodb_annotation_raw,
            ):
                if not protein_file.exists():
                    logging.info(
                        "Missing %s, skipping",
                        protein_file,
                    )
                    continue

                with protein_file.open(
                    encoding="utf-8",
                ) as in_handle:
                    shutil.copyfileobj(
                        in_handle,
                        out_handle,
                    )

    # if protein_annotation_raw.exists():
    #    shutil.copy2(
    #        protein_annotation_raw,
    #        final_annotation_dir / "protein_raw.gtf",
    #    )

    generic_select_cmd = [
        "perl",
        transcript_selector_script,
        "-genome_file",
        genome_file,
    ]

    def run_biotype(
        biotype_name: str,
        annotation_raw: Path,
        output_suffix: str,
        extra_args: list[str],
        threads: int | None = None,
    ) -> None:
        """
        Select the best transcripts for a given evidence type.
        """

        if not annotation_raw.exists():
            logging.info(
                "Missing annotation file: %s",
                annotation_raw,
            )
            return

        # pool_size = threads or num_threads

        async_results: list[tuple[str, str, AsyncResult]] = []

        # with multiprocessing.Pool(
        #    processes=pool_size,
        # ) as pool:
        pool_size = int(threads) if threads is not None else int(num_threads)
        pool = multiprocessing.Pool(pool_size)
        for seq_region_name in seq_region_names:
            region_details = f"{seq_region_name}.rs1" f".re{seq_region_lengths[seq_region_name]}"

            region_gtf_path = region_annotation_dir / f"{region_details}.{output_suffix}.gtf"

            if gtf_ok(region_gtf_path):
                logging.info(
                    "%s region GTF exists, skipping: %s",
                    biotype_name,
                    region_gtf_path,
                )
                continue

            logging.info(
                "Finalising %s data for: %s",
                biotype_name,
                seq_region_name,
            )

            cmd = [
                *generic_select_cmd,
                "-region_details",
                region_details,
                "-input_gtf_file",
                str(annotation_raw),
                "-output_gtf_file",
                str(region_gtf_path),
                "-final_biotype",
                biotype_name,
                *extra_args,
            ]

            result = pool.apply_async(
                multiprocess_finalise_geneset,
                args=(cmd,),
            )

            async_results.append(
                (
                    seq_region_name,
                    biotype_name,
                    result,
                )
            )
        pool.close()
        pool.join()
        for (
            region_name,
            biotype,
            result,
        ) in async_results:
            try:
                logging.info(
                    "Waiting for %s (%s)...",
                    region_name,
                    biotype,
                )
                result.get()
                logging.info(
                    "Finished %s (%s)",
                    region_name,
                    biotype,
                )
            except Exception as exc:
                logging.error(
                    "Job exception for %s (%s): %s",
                    region_name,
                    biotype,
                    exc,
                )
                raise

    # Transcriptomic
    biotype_configs = (
        (
            "transcriptomic",
            transcriptomic_annotation_raw,
            "trans",
            ["-cds_search"],
            num_threads,
        ),
        # BUSCO / OrthoDB
        (
            "busco",
            busco_annotation_raw,
            "busco",
            ["-all_cds_exons"],
            num_threads,
        ),
        # Protein / UniProt
        (
            "protein",
            protein_annotation_raw,
            "protein",
            [
                "-clean_transcripts",
                "-all_cds_exons",
            ],
            num_threads,
        ),
    )

    for (
        biotype_name,
        annotation_raw,
        output_suffix,
        extra_args,
        num_threads,  # pylint: disable=redefined-argument-from-local
    ) in biotype_configs:
        run_biotype(
            biotype_name=biotype_name,
            annotation_raw=annotation_raw,
            output_suffix=output_suffix,
            extra_args=extra_args,
            threads=num_threads,
        )

        logging.info(
            "%s select best transcript finished",
            biotype_name.capitalize(),
        )
    for suffix, label in (
        (".trans.gtf", "transcriptomic"),
        (".busco.gtf", "busco"),
        (".protein.gtf", "protein"),
    ):
        merge_finalise_output_files(
            final_annotation_dir,
            region_annotation_dir,
            suffix,
            label,
        )
    fully_merged_gtf_path = final_annotation_dir / "all_selected_transcripts.gtf"

    with fully_merged_gtf_path.open(
        "w",
        encoding="utf-8",
    ) as out_handle:
        for gtf_file in sorted(final_annotation_dir.glob("*_sel.gtf")):
            with gtf_file.open(
                encoding="utf-8",
            ) as in_handle:
                shutil.copyfileobj(
                    in_handle,
                    out_handle,
                )
    generic_finalise_cmd: list[str] = [
        "perl",
        str(finalise_geneset_script),
        "-genome_file",
        str(genome_file),
    ]
    # with multiprocessing.Pool(
    #    processes=num_threads,
    # ) as pool:
    pool = multiprocessing.Pool(int(num_threads))
    for seq_region_name in seq_region_names:
        region_details = f"{seq_region_name}.rs1" f".re{seq_region_lengths[seq_region_name]}"

        final_region_gtf_path = final_region_annotation_dir / f"{region_details}.final.gtf"

        if gtf_ok(final_region_gtf_path):
            logging.info(
                "Final region GTF exists, skipping: %s",
                final_region_gtf_path,
            )
            continue

        cmd = [
            *generic_finalise_cmd,
            "-region_details",
            region_details,
            "-input_gtf_file",
            str(fully_merged_gtf_path),
            "-output_gtf_file",
            str(final_region_gtf_path),
        ]

        pool.apply_async(
            multiprocess_finalise_geneset,
            args=(cmd,),
        )
    pool.close()
    pool.join()

    merge_finalise_output_files(
        final_annotation_dir,
        final_region_annotation_dir,
        ".final.gtf",
        "prevalidation",
    )
    merged_gtf_file = final_annotation_dir / "prevalidation_sel.gtf"

    merged_cdna_file = final_annotation_dir / "prevalidation_sel.cdna.fa"

    merged_amino_acid_file = final_annotation_dir / "prevalidation_sel.prot.fa"

    postvalidation_gtf_file = final_annotation_dir / "postvalidation.gtf"
    if gtf_ok(postvalidation_gtf_file):
        logging.info(
            "Post-validation GTF exists, skipping validation",
        )
    else:
        updated_gtf_lines = validate_coding_transcripts(
            merged_cdna_file,
            merged_amino_acid_file,
            validation_dir,
            validation_type,
            diamond_validation_db,
            merged_gtf_file,
            num_threads,
        )

        with postvalidation_gtf_file.open(
            "w",
            encoding="utf-8",
        ) as file_out:
            file_out.writelines(
                updated_gtf_lines,
            )
    cleaned_initial_gtf_file = final_annotation_dir / "cleaned_pre_utr.gtf"

    cleaned_utr_gtf_file = final_annotation_dir / "annotation.gtf"
    cleaning_cmd: list[str] = [
        "perl",
        str(clean_geneset_script),
        "-genome_file",
        str(genome_file),
        "-gtf_file",
        str(postvalidation_gtf_file),
        "-output_gtf_file",
        str(cleaned_initial_gtf_file),
    ]
    if gtf_ok(cleaned_utr_gtf_file):
        logging.info(
            "Cleaned annotation exists, skipping",
        )
    else:
        subprocess.run(
            cleaning_cmd,
            check=True,
        )

    generic_clean_utrs_cmd: list[str] = [
        "perl",
        str(clean_utrs_script),
        "-genome_file",
        str(genome_file),
        "-input_gtf_file",
        str(cleaned_initial_gtf_file),
    ]

    # with multiprocessing.Pool(
    #    processes=num_threads,
    # ) as pool:
    pool = multiprocessing.Pool(int(num_threads))
    for seq_region_name in seq_region_names:
        region_details = f"{seq_region_name}.rs1" f".re{seq_region_lengths[seq_region_name]}"

        utr_region_gtf_path = utr_region_annotation_dir / f"{region_details}.utr.gtf"

        cmd = [
            *generic_clean_utrs_cmd,
            "-region_details",
            region_details,
            "-input_gtf_file",
            str(cleaned_initial_gtf_file),
            "-output_gtf_file",
            str(utr_region_gtf_path),
        ]

        pool.apply_async(
            multiprocess_generic,
            args=(cmd,),
        )

    pool.close()
    pool.join()

    merge_finalise_output_files(
        final_annotation_dir,
        utr_region_annotation_dir,
        ".utr.gtf",
        "annotation",
    )

    annotation_sel_gtf = final_annotation_dir / "annotation_sel.gtf"

    shutil.move(
        annotation_sel_gtf,
        cleaned_utr_gtf_file,
    )

    logging.info(
        "Dumping transcript and translation sequences",
    )

    dumping_cmd: list[str] = [
        "perl",
        str(gtf_to_seq_script),
        "-genome_file",
        str(genome_file),
        "-gtf_file",
        str(cleaned_utr_gtf_file),
    ]

    logging.info(
        " ".join(dumping_cmd),
    )

    subprocess.run(
        dumping_cmd
        # check=True,
    )

    logging.info(
        "Finished creating gene set",
    )


def validate_coding_transcripts(  # pylint: disable=too-many-arguments, too-many-positional-arguments
    cdna_file: Path,
    amino_acid_file: Path,
    validation_dir: Path,
    validation_type: str,
    diamond_validation_db: Path | None,
    gtf_file: Path,
    num_threads: int,
) -> list[str]:
    """
    Validate coding transcripts using RNAsamba, CPC2 and optionally Diamond.
    Validation results are combined and used to update the biotypes of transcripts in
    the GTF file. The updated GTF lines are returned as a list of strings.
    Asociated steps include:
    - Running RNAsamba and CPC2 on the merged transcriptome to get coding potential
    predictions for each transcript
    - Optionally running Diamond to validate transcripts against a protein database,
    which can provide strong evidence for coding potential
    - Combining the results from RNAsamba, CPC2 and Diamond to get an overall validation
    result for each transcript
    - Reading the original GTF file and updating the biotypes of transcripts based on
    the combined validation results, with different levels of strictness for single exon
    coding transcripts depending on the validation_type parameter
    Args:
    cdna_file: Path to the FASTA file containing the cDNA sequences of the transcripts.
    amino_acid_file: Path to the FASTA file containing the amino acid sequences of the
    transcripts.
    validation_dir: Path to the directory where validation results will be stored.
    validation_type: The level of strictness for validating coding transcripts.
    Can be "relaxed", "moderate" or "strict". Affects mainly the validation of single
    exon coding transcripts.
    diamond_validation_db: Optional path to a Diamond database for validating coding
    transcripts. If not provided, Diamond validation will be skipped.
    gtf_file: Path to the GTF file containing the original transcript annotations.
    num_threads: Number of threads to use for parallel processing.
    """

    validation_dir = Path(validation_dir)

    rnasamba_output = validation_dir / "rnasamba.tsv.txt"
    cpc2_output = validation_dir / "cpc2.tsv"
    cpc2_output_file = validation_dir / "cpc2.tsv.txt"

    _run_rnasamba(
        cdna_file=cdna_file,
        validation_dir=validation_dir,
        output_file=rnasamba_output,
    )

    _run_cpc2(
        cdna_file=cdna_file,
        validation_dir=validation_dir,
        output_file=cpc2_output,
        expected_output=cpc2_output_file,
    )

    diamond_results: list[list[str]] | None = None

    if diamond_validation_db is not None:
        logging.info("Running Diamond validation")

        diamond_output_dir = create_dir(
            validation_dir,
            "diamond_output",
        )

        diamond_validation(
            diamond_validation_db,
            amino_acid_file,
            diamond_output_dir,
            num_threads,
        )

        diamond_results = read_diamond_results(
            diamond_output_dir,
        )

    logging.info("Reading validation results")

    rnasamba_results = read_rnasamba_results(
        rnasamba_output,
    )

    cpc2_results = read_cpc2_results(
        cpc2_output_file,
    )

    combined_results = combine_results(
        rnasamba_results,
        cpc2_results,
        diamond_results,
    )

    parsed_gtf_genes = read_gtf_genes(
        gtf_file,
    )

    return update_gtf_genes(
        parsed_gtf_genes,
        combined_results,
        validation_type,
    )


def _run_rnasamba(
    cdna_file: Path,
    validation_dir: Path,
    output_file: Path,
) -> None:
    """
    Run RNAsamba if results are not present.
    Args:        cdna_file: Path to the FASTA file containing the cDNA sequences of
    the transcripts.
        validation_dir: Path to the directory where validation results will be stored.
        output_file: Path to the file where RNAsamba results will be written.
    """

    if output_file.exists():
        logging.info(
            "Found RNAsamba output. Skipping.",
        )
        return

    logging.info(
        "Running CDS validation with RNAsamba",
    )

    cmd: list[str] = [
        "singularity",
        "exec",
        "--bind",
        f"{validation_dir}:/app:rw",
        config["rnasamba"]["software"],
        "rnasamba",
        "classify",
        str(output_file),
        str(cdna_file),
        config["rnasamba"]["weights"],
    ]

    logging.info(" ".join(cmd))

    subprocess.run(
        cmd
        # check=True,
    )

    check_file(output_file)

    logging.info(
        "RNAsamba completed successfully",
    )


def _run_cpc2(
    cdna_file: Path,
    validation_dir: Path,
    output_file: Path,
    expected_output: Path,
) -> None:
    """
    Run CPC2 if results are not present.
    Args:
    cdna_file: Path to the FASTA file containing the cDNA sequences of the transcripts.
    validation_dir: Path to the directory where validation results will be stored.
    output_file: Path to the file where CPC2 results will be written.
    expected_output: Path to the file that should be produced by CPC2. This is used to
    check whether CPC2 has already been run, as the actual output file may have a
    different name or extension.
    """

    if expected_output.exists():
        logging.info(
            "Found CPC2 output. Skipping.",
        )
        return

    logging.info(
        "Running CDS validation with CPC2",
    )

    cmd: list[str] = [
        "singularity",
        "exec",
        "--bind",
        f"{validation_dir}:/app:rw",
        config["cpc2"]["software"],
        "python3",
        "/CPC2_standalone-1.0.1/bin/CPC2.py",
        "-i",
        str(cdna_file),
        "--ORF",
        "-o",
        str(output_file),
    ]

    logging.info(" ".join(cmd))

    subprocess.run(
        cmd
        # check=True,
    )

    check_file(expected_output)

    logging.info(
        "CPC2 completed successfully",
    )


def diamond_validation(
    diamond_validation_db: Path,
    amino_acid_file: Path,
    diamond_output_dir: Path,
    num_threads: int,
) -> None:
    """
    Validate proteins against a Diamond database.
    Asociated steps include:
    - Splitting the amino acid FASTA file into batches for parallel processing
    - Running Diamond blastp for each batch against the provided Diamond database
    Args:
    diamond_validation_db: Path to the Diamond database to use for validation.
    amino_acid_file: Path to the FASTA file containing the amino acid sequences of the
    transcripts.
    diamond_output_dir: Path to the directory where Diamond output files will be stored.
    num_threads: Number of threads to use for parallel processing.
    """

    logging.info(
        "Starting Diamond validation",
    )
    logging.info(
        "Diamond DB: %s",
        diamond_validation_db,
    )
    logging.info(
        "Amino acid file: %s",
        amino_acid_file,
    )
    logging.info(
        "Diamond output dir: %s",
        diamond_output_dir,
    )
    logging.info(
        "Threads: %s",
        num_threads,
    )

    batched_protein_files = split_protein_file(
        amino_acid_file,
        diamond_output_dir,
        100,
    )

    logging.info(
        "Split protein file into %d batches",
        len(batched_protein_files),
    )

    if not batched_protein_files:
        raise RuntimeError(
            "No batched protein files were created",
        )

    async_results: list[tuple[Path, AsyncResult]] = []

    # with multiprocessing.Pool(
    #    processes=num_threads,
    # ) as pool:
    pool = multiprocessing.Pool(int(num_threads))
    for batch_file in batched_protein_files:
        logging.info(
            "Submitting Diamond job for batch: %s",
            batch_file,
        )

        result = pool.apply_async(
            multiprocess_diamond,
            args=(
                Path(batch_file),
                diamond_output_dir,
                diamond_validation_db,
            ),
        )

        async_results.append(
            (
                Path(batch_file),
                result,
            )
        )

    for batch_file, result in async_results:
        try:
            result.get()
        except Exception as exc:
            logging.error(
                "Diamond failed for batch %s: %s",
                batch_file,
                exc,
            )
            raise
    pool.close()
    pool.join()
    logging.info(
        "Diamond validation finished",
    )


def multiprocess_diamond(
    batched_protein_file: Path,
    diamond_output_dir: Path,
    diamond_validation_db: Path,
) -> None:
    """
    Run Diamond on a protein batch.
    Args:
    batched_protein_file: Path to the FASTA file containing the amino acid sequences
    for     this batch.
    diamond_output_dir: Path to the directory where Diamond output files will be stored.
    diamond_validation_db: Path to the Diamond database to use for validation.
    """

    if not batched_protein_file.exists():
        raise FileNotFoundError(f"Missing batch file: " f"{batched_protein_file}")

    diamond_output_file = diamond_output_dir / f"{batched_protein_file.name}.dmdout"

    diamond_cmd: list[str] = [
        "diamond",
        "blastp",
        "--query",
        str(batched_protein_file),
        "--db",
        str(diamond_validation_db),
        "--out",
        str(diamond_output_file),
        # "--outfmt",
        # "6",
        # "--threads",
        # "1",
    ]

    logging.info(
        "Running Diamond for %s",
        batched_protein_file.name,
    )

    subprocess.run(
        diamond_cmd
        # check=True,
        # capture_output=True,
        # text=True,
    )

    if not diamond_output_file.exists():
        raise FileNotFoundError(f"Diamond produced no output file: " f"{diamond_output_file}")

    if diamond_output_file.stat().st_size == 0:
        raise RuntimeError(f"Diamond output file is empty: " f"{diamond_output_file}")

    logging.info(
        "Diamond completed successfully: %s",
        batched_protein_file.name,
    )


def update_gtf_genes(
    parsed_gtf_genes: dict[str, dict[str, dict[str, str]]],
    combined_results: dict[str, list[str]],
    validation_type: str,
) -> list[str]:
    """
    Update GTF gene annotations based on validation results.

    Args:
    parsed_gtf_genes: A nested dictionary containing parsed GTF gene annotations.
    The structure is
    {gene_id: {transcript_id: {"transcript": transcript_line, "exons": [exon_lines]}}}.
    combined_results: A dictionary containing combined validation results for each
    transcript. The
    structure is {transcript_id: [rnasamba_coding_probability, rnasamba_coding_potential,
    cpc2_coding_probability, cpc2_coding_potential, transcript_length, peptide_length, (optional) diamond_e_value]}.
    validation_type: The level of strictness for validating coding transcripts. Can be
    "relaxed" or "strict".

    Returns:
    A list of strings representing the updated GTF lines for the genes, with biotypes
    updated based on the validation results.
    """

    output_lines = []

    for gene_id in parsed_gtf_genes.keys():
        transcript_ids = parsed_gtf_genes[gene_id].keys()
        for transcript_id in transcript_ids:
            # if transcript_id not in combined_results:
            # ANNA skip/delete this transcript TEMPORARY FIX FOR 3BP CDS
            #    continue
            transcript_line = parsed_gtf_genes[gene_id][transcript_id]["transcript"]
            single_cds_exon_transcript = 0
            translation_match = re.search(r'; translation_coords "([^"]+)";', transcript_line)
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

            avg_coding_probability = (rnasamba_coding_probability + cpc2_coding_probability) / 2
            max_coding_probability = max(rnasamba_coding_probability, cpc2_coding_probability)

            match = re.search(r'; biotype "([^"]+)";', transcript_line)
            if match:
                biotype = match.group(1)  # pylint: disable=possibly-used-before-assignment
                if biotype in ("busco", "protein"):
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
                    (rnasamba_coding_potential == "coding" or cpc2_coding_potential == "coding")
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
                if diamond_e_value is not None and peptide_length >= min_single_exon_pep_length:
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    (rnasamba_coding_potential == "coding" and cpc2_coding_potential == "coding")
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
                    (rnasamba_coding_potential == "coding" or cpc2_coding_potential == "coding")
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
                    transcript_line = re.sub(' translation_coords "[^"]+";', "", transcript_line)
                else:
                    continue

            output_lines.append(transcript_line)
            output_lines.extend(exon_lines)

    return output_lines


def read_rnasamba_results(
    file_path: Path,
) -> list[list[str]]:
    """Read RNAsamba results.
    The output file has three columns: sequence_name, coding_probability, coding_potential.
    The coding_probability is a value between 0 and 1, where higher values indicate
    a higher likelihood of being a protein-coding transcript.
    The coding_potential is a categorical prediction of "coding" or "noncoding" based
    on the coding_probability and the model's internal thresholds.
    A typical threshold for classifying a transcript as coding might be around 0.5, but
    this can vary depending on the specific model and dataset used for training.
    Args:    file_path: Path to the RNAsamba output file.
    Returns: A list of lists, where each inner list contains the transcript ID, coding
    probability and coding potential for
    """

    results: list[list[str]] = []

    with file_path.open(
        encoding="utf-8",
    ) as file_handle:
        for line in file_handle:
            line = line.rstrip()

            if line.startswith("sequence_name"):
                continue

            fields = line.split("\t")

            if len(fields) != 3:
                continue

            transcript_id = fields[0]
            coding_probability = fields[1]
            coding_potential = fields[2]
            results.append([transcript_id, coding_probability, coding_potential])

    return results


def read_cpc2_results(
    file_path: Path,
) -> list[list[str]]:
    """Read CPC2 results.
    The output file has nine columns: ID, Length, ORF length, Fickett score, Hexamer
    score, ORF integrity, ORF coverage, coding probability and coding potential.
    The coding probability is a value between 0 and 1, where higher values indicate a
    higher likelihood of being a protein-coding transcript.
    The coding potential is a categorical prediction of "coding" or "noncoding" based
    on the coding probability and the model's internal thresholds.
    A typical threshold for classifying a transcript as coding might be around 0.5,
    but this can vary depending on the specific model and dataset used for training.
    Args:    file_path: Path to the CPC2 output file.
    Returns: A list of lists, where each inner list contains the transcript ID, coding
    probability, coding potential, transcript length and peptide length for each transcript.
    """

    results: list[list[str]] = []

    with file_path.open(
        encoding="utf-8",
    ) as file_handle:
        for line in file_handle:
            line = line.rstrip()

            if line.startswith("#ID"):
                continue

            fields = line.split("\t")

            if len(fields) != 9:
                continue

            transcript_id = fields[0]
            transcript_length = fields[1]
            peptide_length = fields[2]
            coding_probability = fields[7]
            coding_potential = fields[8]
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


def read_diamond_results(
    diamond_output_dir: Path,
) -> list[list[str]]:
    """Read Diamond results.
    The output file has 12 columns: query id, subject id, percentage identity, alignment
    length, mismatches, gap opens, q. start, q. end, s. start, s. end, e-value and bit
    score.
    The e-value is a measure of the statistical significance of the match, with lower
    values indicating more significant matches. A common threshold for considering a
    match to be significant is an e-value of 1e-5 or lower, but this can vary depending
    on the specific analysis and dataset.
    Args:    diamond_output_dir: Path to the directory containing the Diamond output files.
    Returns: A list of lists, where each inner list contains the transcript ID and e-value
    for each significant match found by Diamond.
    """

    results: list[list[str]] = []

    for file_path in diamond_output_dir.glob("*.dmdout"):
        with file_path.open(
            encoding="utf-8",
        ) as file_handle:
            for line in file_handle:
                fields = line.rstrip().split("\t")

                if len(fields) != 12:
                    continue

                transcript_id = fields[0]
                e_value = fields[10]
                results.append([transcript_id, e_value])

    return results


def combine_results(
    rnasamba_results: list[list[str]],
    cpc2_results: list[list[str]],
    diamond_results: list[list[str]] | None,
) -> dict[str, list[str]]:
    """Combine validation results.
    Asociated steps include:
    - Creating a dictionary to store the combined results for each transcript,
    indexed by transcript ID
    - Populating the dictionary with RNAsamba results, using the transcript ID
    as the key
    Adding CPC2 results to the dictionary, matching transcripts by ID and extending
    the existing entries with the new results
    - If Diamond results are available, adding them to the dictionary in the same
    way, ensuring that only transcripts that have results from both RNAsamba and
    CPC2 are included in the final combined results, as these are the transcripts
    for which we have the most comprehensive validation information. Diamond results
    can be added to these transcripts if available, but transcripts without RNAsamba
    or CPC2 results will be excluded from the final combined results regardless of
    whether they have Diamond matches or not.
    Args:
    rnasamba_results: A list of lists containing the RNAsamba results for each transcript
    cpc2_results: A list of lists containing the CPC2 results for each transcript
    diamond_results: An optional list of lists containing the Diamond results for each
    transcript. If not provided, Diamond validation will be skipped and the combined
    results will only include RNAsamba and CPC2 results.
    Returns:
    A dictionary containing the combined validation results for each transcript,
    indexed by transcript ID. Each entry in the dictionary is a list of strings containing
    the RNAsamba coding probability, RNAsamba coding potential, CPC2 coding probability,
    CPC2 coding potential, transcript length, peptide length and optionally the Diamond
    e-value for each transcript.
    """

    transcript_ids: dict[
        str,
        list[str],
    ] = {}

    for (
        transcript_id,
        coding_probability,
        coding_potential,
    ) in rnasamba_results:
        if transcript_id not in transcript_ids:
            transcript_ids.setdefault(
                transcript_id,
                [
                    coding_probability,
                    coding_potential,
                ],
            )

    for (
        transcript_id,
        coding_probability,
        coding_potential,
        transcript_length,
        peptide_length,
    ) in cpc2_results:
        if transcript_id not in transcript_ids:
            continue

        transcript_ids[transcript_id].extend(
            [
                coding_probability,
                coding_potential,
                transcript_length,
                peptide_length,
            ]
        )

    if diamond_results:
        for transcript_id, e_value in diamond_results:
            # There seems to be an issue where there are a small number
            # of sequences that don't make it into the cpc2/rnasamba output
            # Should code in a system for this, but it would be good to
            # understand why it happens to begin with. Seems to be the same
            # number of missing seqs in both, so maybe a shared cut-off
            if transcript_id in transcript_ids:
                transcript_ids[transcript_id].append(
                    e_value,
                )

    return transcript_ids


def read_gtf_genes(
    gtf_file: Path,
):
    """Read transcript and exon records from a GTF.
    The GTF file is parsed to extract gene and transcript information, which is
    stored in a nested dictionary structure. The outer dictionary is keyed by
    gene ID, and each value is another dictionary keyed by transcript ID. Each
    transcript entry contains the original GTF line for the transcript feature
    and a list of GTF lines for the exon features associated with that transcript.
    Args:
        gtf_file: Path to the GTF file to be parsed.
    Returns:
        A nested dictionary containing the parsed gene and transcript information
        from the GTF file. The structure is
        {gene_id: {transcript_id: {"transcript": transcript_line, "exons": [exon_lines]}}}.
    """

    gtf_genes: dict[str, dict[str, dict[str, str | list[str]]]] = {}

    with gtf_file.open(
        encoding="utf-8",
    ) as file_handle:
        for line in file_handle:
            fields = line.split("\t")

            if len(fields) != 9:
                continue

            match = re.search(r'gene_id "([^"]+)".+transcript_id "([^"]+)"', line)

            if match is None:
                continue

            gene_id = match.group(1)
            transcript_id = match.group(2)
            feature_type = fields[2]

            if gene_id not in gtf_genes:
                gtf_genes[gene_id] = {}
            if feature_type == "transcript":
                gtf_genes[gene_id][transcript_id] = {}
                gtf_genes[gene_id][transcript_id]["transcript"] = line
                gtf_genes[gene_id][transcript_id]["exons"] = []
            elif feature_type == "exon":
                cast(list[str], gtf_genes[gene_id][transcript_id]["exons"]).append(line)

    return gtf_genes


def fasta_to_dict(
    fasta_lines: list[str],
) -> dict[str, str]:
    """
    Convert FASTA records into a dictionary.
    The FASTA file is read line by line, and sequences are stored in a dictionary
    where the keys are the sequence headers (without the leading '>') and the values
    are the corresponding sequences. The sequences are concatenated if they span
    multiple lines in the FASTA file. Each sequence is also appended with a newline
    character at the end.
    Args:
    fasta_lines: A list of strings representing the lines of a FASTA file.
    Returns:
    A dictionary where the keys are sequence headers (without the leading '>') and
    the values are the corresponding sequences, concatenated if they span multiple
    lines and with a newline character appended at the end
    """

    sequences: dict[str, str] = {}

    current_header: str | None = None
    current_sequence: list[str] = []

    for line in fasta_lines:
        line = line.rstrip()

        if line.startswith(">"):
            if current_header is not None:
                sequences[current_header] = "".join(current_sequence) + "\n"

            current_header = line[1:]
            current_sequence = []

        else:
            current_sequence.append(line)

    if current_header is not None:
        sequences[current_header] = "".join(current_sequence) + "\n"

    return sequences


def merge_finalise_output_files(
    final_annotation_dir: Path,
    region_annotation_dir: Path,
    extension: str,
    id_label: str,
) -> None:
    """
        Merge regional GTF, cDNA and protein files.
        This function takes the regional GTF files and their corresponding cDNA and
        protein FASTA files, and merges them into single GTF, cDNA and protein files
        in the final annotation directory. During the merging process, gene and transcript
        IDs are updated to be unique and consistent across all three file types.
        The function reads each regional GTF file, extracts the gene and transcript
        information, updates the IDs, and writes the merged records to the output files.
        The cDNA and protein sequences are indexed in memory to allow for efficient
        retrieval when writing the merged FASTA files.
        Args:
        final_annotation_dir: Path to the directory where the merged output files will
        be written.
        region_annotation_dir: Path to the directory containing the regional GTF files
        and their corresponding c
    DNA and protein FASTA files.
        extension: The file extension of the regional GTF files to be merged (e.g., ".utr.gtf").
        id_label: A label to prefix the gene and transcript IDs with in the merged
        output files, to ensure uniqueness and consistency across the merged GTF, cDNA
        and protein files.
    """

    gtf_files = sorted(region_annotation_dir.glob(f"*{extension}"))

    merged_gtf_file = final_annotation_dir / f"{id_label}_sel.gtf"

    merged_cdna_file = final_annotation_dir / f"{id_label}_sel.cdna.fa"

    merged_protein_file = final_annotation_dir / f"{id_label}_sel.prot.fa"
    # The below is not great, it's a bit messy because there might be
    # some cases where there aren't translations. So it's not as
    # straightforward as reading the records across all three files
    # in parallel. The solution is to just load the seqs into
    # memory and index them on the current header, which should
    # correspond to a transcript/gene id in the GTF. When writing the
    # results into the single merged files the ids will be updated to
    # be unique and consistent across the header,  three file types

    gene_id_counter = 0
    transcript_id_counter = 0

    with (
        merged_gtf_file.open(
            "w",
            encoding="utf-8",
        ) as gtf_out,
        merged_cdna_file.open(
            "w",
            encoding="utf-8",
        ) as cdna_out,
        merged_protein_file.open(
            "w",
            encoding="utf-8",
        ) as protein_out,
    ):
        for gtf_file in gtf_files:
            logging.info(
                "Processing %s",
                gtf_file,
            )

            cdna_file = Path(f"{gtf_file}.cdna")

            protein_file = Path(f"{gtf_file}.prot")

            with cdna_file.open(
                encoding="utf-8",
            ) as handle:
                cdna_index = fasta_to_dict(handle.readlines())

            with protein_file.open(
                encoding="utf-8",
            ) as handle:
                protein_index = fasta_to_dict(handle.readlines())

            current_gene_id = ""

            with gtf_file.open(
                encoding="utf-8",
            ) as gtf_in:
                for line in gtf_in:
                    if line.startswith("#"):
                        continue

                    fields = line.split("\t")

                    if len(fields) != 9:
                        continue

                    match = re.search(r'gene_id "([^"]+)".+transcript_id "([^"]+)"', line)

                    if match is None:
                        continue

                    gene_id = match.group(1)
                    transcript_id = match.group(2)

                    feature_type = fields[2]

                    if feature_type == "transcript":
                        transcript_id_counter += 1

                    if not current_gene_id:
                        gene_id_counter += 1
                        current_gene_id = gene_id

                    elif gene_id != current_gene_id:
                        gene_id_counter += 1
                        current_gene_id = gene_id

                    new_gene_id = f"{id_label}_" f"{gene_id_counter}"

                    new_transcript_id = f"{id_label}_" f"{transcript_id_counter}"

                    line = line.replace(
                        f'gene_id "{gene_id}"',
                        (f"gene_id " f'"{new_gene_id}"'),
                    )

                    line = line.replace(
                        (f"transcript_id " f'"{transcript_id}"'),
                        (f"transcript_id " f'"{new_transcript_id}"'),
                    )

                    gtf_out.write(line)

                    if feature_type != "transcript":
                        continue

                    new_header = f">{new_transcript_id}\n"

                    cdna_out.write(new_header)

                    cdna_out.write(cdna_index[transcript_id])

                    protein_seq = protein_index.get(transcript_id)

                    if protein_seq:
                        protein_out.write(new_header + protein_seq)

                        # protein_out.write(protein_seq)


def multiprocess_generic(cmd):
    """Run command and fail hard on error"""
    logging.info("Worker executing: %s", " ".join(str(x) for x in cmd))
    subprocess.run(cmd)


def multiprocess_finalise_geneset(cmd):
    """Run command and fail hard on error, with logging."""
    logging.info("Worker executing: %s", " ".join(str(x) for x in cmd))

    result = subprocess.run(
        cmd
        # capture_output=True,
        # text=True,
        # check=False,
    )

    if result.returncode != 0:
        logging.error(
            "Command failed with return code %d: %s",
            result.returncode,
            " ".join(str(x) for x in cmd),
        )
        logging.error("STDERR:\n%s", result.stderr)
        if result.stdout:
            logging.error("STDOUT:\n%s", result.stdout)

        # propagate failure to parent
        raise RuntimeError(f"Selector failed (return code {result.returncode})")

    logging.info("Command succeeded: %s", " ".join(str(x) for x in cmd))
    return result
