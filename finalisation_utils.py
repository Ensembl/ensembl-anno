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
    protein_annotation_raw = main_output_dir / "genblast_output", "annotation.gtf"
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
                "_raw.gtf", "_sel.gtf", transcriptomic_annotation_raw
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
                "_raw.gtf", "_sel.gtf", busco_annotation_raw
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
                "_raw.gtf", "_sel.gtf", protein_annotation_raw
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
    logger.info("Cleaning initial set:\n%s" % " ".join(cleaning_cmd))
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
        "Dumping transcript and translation sequences:\n%s" % " ".join(dumping_cmd)
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
    logger.info("rnasamba_cmd: %s" % " ".join(rnasamba_cmd))
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
    logger.info("cpc2_cmd: %s" % " ".join(cpc2_cmd))
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
        "Running diamond on %s:\n%s" % (batched_protein_file, " ".join(diamond_cmd))
    )
    subprocess.run(diamond_cmd)
    subprocess.run(["mv", diamond_output_file, diamond_output_dir])


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

