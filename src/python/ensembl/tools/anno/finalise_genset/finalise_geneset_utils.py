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
    main_script_dir,
    main_output_dir,
    genome_file,
    seq_region_names,
    validation_type,
    diamond_validation_db,
    num_threads,
):

    if validation_type is None:
        logger.info("Setting validation type to relaxed")
    else:
        logger.info("Setting validation type to " + validation_type)

    final_annotation_dir = utils.create_dir(main_output_dir, "annotation_output")
    region_annotation_dir = utils.create_dir(final_annotation_dir, "initial_region_gtfs")
    final_region_annotation_dir = utils.create_dir(
        final_annotation_dir, "final_region_gtfs"
    )
    utr_region_annotation_dir = utils.create_dir(final_annotation_dir, "utr_region_gtfs")
    validation_dir = utils.create_dir(final_annotation_dir, "cds_validation")
    seq_region_lengths = utils.get_seq_region_lengths(genome_file, 0)

    logger.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(final_annotation_dir, "annotation.gtf")
    if os.path.exists(output_file):
        logger.info("final_annotation_dir exists")
        transcript_count = utils.check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Final gtf file exists")
            return
    else:
        logger.info("No gtf file, go on with the analysis")
    # This used to be a list of output dirs and a loop which was neat,
    # I'm coverting to a list of conditions as
    # it's more straightforward with the renaming
    # and having to merge scallop and stringtie
    protein_annotation_raw = os.path.join(
        main_output_dir, "genblast_output", "annotation.gtf"
    )
    minimap2_annotation_raw = os.path.join(
        main_output_dir, "minimap2_output", "annotation.gtf"
    )
    stringtie_annotation_raw = os.path.join(
        main_output_dir, "stringtie_output", "annotation.gtf"
    )
    scallop_annotation_raw = os.path.join(
        main_output_dir, "scallop_output", "annotation.gtf"
    )
    busco_annotation_raw = os.path.join(main_output_dir, "busco_output", "annotation.gtf")

    transcript_selector_script = os.path.join(
        main_script_dir, "support_scripts_perl", "select_best_transcripts.pl"
    )
    finalise_geneset_script = os.path.join(
        main_script_dir, "support_scripts_perl", "finalise_geneset.pl"
    )
    clean_geneset_script = os.path.join(
        main_script_dir, "support_scripts_perl", "clean_geneset.pl"
    )
    clean_utrs_script = os.path.join(
        main_script_dir, "support_scripts_perl", "clean_utrs_and_lncRNAs.pl"
    )
    gtf_to_seq_script = os.path.join(
        main_script_dir, "support_scripts_perl", "gtf_to_seq.pl"
    )

    transcriptomic_annotation_raw = os.path.join(
        final_annotation_dir, "transcriptomic_raw.gtf"
    )
    file_out = open(transcriptomic_annotation_raw, "w+")
    for transcriptomic_file in [
        minimap2_annotation_raw,
        scallop_annotation_raw,
        stringtie_annotation_raw,
    ]:

        if not os.path.exists(transcriptomic_file):
            logger.info(
                "No annotation.gtf file found in " + transcriptomic_file + ", skipping"
            )
            continue

        file_in = open(transcriptomic_file)
        line = file_in.readline()
        while line:
            print(line.rstrip(), file=file_out)
            line = file_in.readline()
        file_in.close()
    file_out.close()

    # Copy the raw files into the annotation dir, this is not needed
    # as such, but collecting them in one place and relabelling is
    # helpful for a user
    if os.path.exists(busco_annotation_raw):
        subprocess.run(
            [
                "cp",
                busco_annotation_raw,
                os.path.join(final_annotation_dir, "busco_raw.gtf"),
            ]
        )

    if os.path.exists(protein_annotation_raw):
        subprocess.run(
            [
                "cp",
                protein_annotation_raw,
                os.path.join(final_annotation_dir, "protein_raw.gtf"),
            ]
        )

    gtf_files = ["transcriptomic_raw.gtf", "protein_raw.gtf", "busco_raw.gtf"]
    generic_select_cmd = [
        "perl",
        transcript_selector_script,
        "-genome_file",
        genome_file,
    ]
    pool = multiprocessing.Pool(int(num_threads))
    for seq_region_name in seq_region_names:
        # The selection script needs different params depending on
        # whether the seqs are from transcriptomic data or not
        region_details = (
            seq_region_name + ".rs1" + ".re" + str(seq_region_lengths[seq_region_name])
        )
        transcriptomic_region_gtf_path = os.path.join(
            region_annotation_dir, (region_details + ".trans.gtf")
        )
        busco_region_gtf_path = os.path.join(
            region_annotation_dir, (region_details + ".busco.gtf")
        )
        protein_region_gtf_path = os.path.join(
            region_annotation_dir, (region_details + ".protein.gtf")
        )

        if os.path.exists(transcriptomic_annotation_raw):
            logger.info("Finalising transcriptomic data for: " + seq_region_name)
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
            pool.apply_async(multiprocess_finalise_geneset, args=(cmd,))

        if os.path.exists(busco_annotation_raw):
            logger.info("Finalising BUSCO data for: " + seq_region_name)
            busco_annotation_select = re.sub("_raw.gtf", "_sel.gtf", busco_annotation_raw)
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
            pool.apply_async(multiprocess_finalise_geneset, args=(cmd,))

        if os.path.exists(protein_annotation_raw):
            logger.info("Finalising protein data for: " + seq_region_name)
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
            pool.apply_async(multiprocess_finalise_geneset, args=(cmd,))

    pool.close()
    pool.join()

    # At this point we will have the region files for all the,
    merge_finalise_output_files(
        final_annotation_dir,
        region_annotation_dir,
        ".trans.gtf",
        "transcriptomic",
    )
    merge_finalise_output_files(
        final_annotation_dir, region_annotation_dir, ".busco.gtf", "busco"
    )
    merge_finalise_output_files(
        final_annotation_dir, region_annotation_dir, ".protein.gtf", "protein"
    )

    # Create a single GTF file with all the selected transcripts
    # now that they have proper ids
    fully_merged_gtf_path = os.path.join(
        final_annotation_dir, "all_selected_transcripts.gtf"
    )
    fully_merged_gtf_out = open(fully_merged_gtf_path, "w+")

    merge_gtf_cmd = ["cat"]
    merge_gtf_cmd.extend(glob.glob(final_annotation_dir + "/*_sel.gtf"))
    subprocess.run(merge_gtf_cmd, stdout=fully_merged_gtf_out)
    fully_merged_gtf_out.close()

    # Now collapse the gene set
    generic_finalise_cmd = [
        "perl",
        finalise_geneset_script,
        "-genome_file",
        genome_file,
    ]

    pool = multiprocessing.Pool(int(num_threads))
    for seq_region_name in seq_region_names:
        region_details = (
            seq_region_name + ".rs1" + ".re" + str(seq_region_lengths[seq_region_name])
        )
        final_region_gtf_path = os.path.join(
            final_region_annotation_dir, (region_details + ".final.gtf")
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
        pool.apply_async(multiprocess_finalise_geneset, args=(cmd,))

    pool.close()
    pool.join()

    merge_finalise_output_files(
        final_annotation_dir,
        final_region_annotation_dir,
        ".final.gtf",
        "prevalidation",
    )
    merged_gtf_file = os.path.join(final_annotation_dir, ("prevalidation_sel.gtf"))
    merged_cdna_file = os.path.join(final_annotation_dir, ("prevalidation_sel.cdna.fa"))
    merged_amino_acid_file = os.path.join(
        final_annotation_dir, ("prevalidation_sel.prot.fa")
    )
    updated_gtf_lines = validate_coding_transcripts(
        merged_cdna_file,
        merged_amino_acid_file,
        validation_dir,
        validation_type,
        diamond_validation_db,
        merged_gtf_file,
        num_threads,
    )
    postvalidation_gtf_file = os.path.join(final_annotation_dir, ("postvalidation.gtf"))
    file_out = open(postvalidation_gtf_file, "w+")
    for line in updated_gtf_lines:
        file_out.write(line)
    file_out.close()

    cleaned_initial_gtf_file = os.path.join(final_annotation_dir, ("cleaned_pre_utr.gtf"))
    cleaned_utr_gtf_file = os.path.join(final_annotation_dir, ("annotation.gtf"))

    logger.info("Cleaning initial set")
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
    logger.info(" ".join(cleaning_cmd))
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
    pool = multiprocessing.Pool(int(num_threads))
    for seq_region_name in seq_region_names:
        region_details = (
            seq_region_name + ".rs1" + ".re" + str(seq_region_lengths[seq_region_name])
        )
        utr_region_gtf_path = os.path.join(
            utr_region_annotation_dir, (region_details + ".utr.gtf")
        )

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
        pool.apply_async(multiprocess_generic, args=(cmd,))
    pool.close()
    pool.join()

    merge_finalise_output_files(
        final_annotation_dir,
        utr_region_annotation_dir,
        ".utr.gtf",
        "annotation",
    )
    subprocess.run(
        [
            "mv",
            os.path.join(final_annotation_dir, "annotation_sel.gtf"),
            cleaned_utr_gtf_file,
        ]
    )

    logger.info("Dumping transcript and translation sequences")
    dumping_cmd = [
        "perl",
        gtf_to_seq_script,
        "-genome_file",
        genome_file,
        "-gtf_file",
        cleaned_utr_gtf_file,
    ]
    logger.info(" ".join(dumping_cmd))
    subprocess.run(dumping_cmd)

    logger.info("Finished creating gene set")


def validate_coding_transcripts(
    cdna_file,
    amino_acid_file,
    validation_dir,
    validation_type,
    diamond_validation_db,
    gtf_file,
    num_threads,
):

    logger.info("Running CDS validation with RNAsamba and CPC2")
    rnasamba_weights = config["rnasamba"]["weights"]
    rnasamba_output_path = os.path.join(validation_dir, "rnasamba.tsv.txt")
    cpc2_output_path = os.path.join(validation_dir, "cpc2.tsv")
    rnasamba_volume = validation_dir + "/:/app:rw"
    rnasamba_cmd = [
        "singularity",
        "exec",
        "--bind",
        rnasamba_volume,
        config["rnasamba"]["software"],
        "rnasamba",
        "classify",
        rnasamba_output_path,
        cdna_file,
        rnasamba_weights,
    ]
    logger.info(" ".join(rnasamba_cmd))
    subprocess.run(rnasamba_cmd)
    cpc2_volume = validation_dir + "/:/app:rw"
    cpc2_cmd = [
        "singularity",
        "exec",
        "--bind",
        cpc2_volume,
        config["cpc2"]["software"],
        "python3",
        "/CPC2_standalone-1.0.1/bin/CPC2.py",
        "-i",
        cdna_file,
        "--ORF",
        "-o",
        cpc2_output_path,
    ]
    logger.info(" ".join(cpc2_cmd))
    subprocess.run(cpc2_cmd)
    cpc2_output_path = cpc2_output_path + ".txt"

    check_file(rnasamba_output_path)
    check_file(cpc2_output_path)

    logger.info("diamond validation")
    diamond_results = None
    if diamond_validation_db is not None:
        diamond_output_dir = utils.create_dir(validation_dir, "diamond_output")
        diamond_validation(
            diamond_validation_db,
            amino_acid_file,
            diamond_output_dir,
            num_threads,
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
    diamond_validation_db, amino_acid_file, diamond_output_dir, num_threads
):

    batched_protein_files = split_protein_file(amino_acid_file, diamond_output_dir, 100)

    pool = multiprocessing.Pool(int(num_threads))
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
    batched_protein_file,
    diamond_output_dir,
    diamond_validation_db,
):

    batch_num = os.path.splitext(batched_protein_file)[0]
    batch_dir = os.path.dirname(batched_protein_file)
    diamond_output_file = batched_protein_file + ".dmdout"
    logger.info("Running diamond on " + batched_protein_file + ":")

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

    logger.info(" ".join(diamond_cmd))
    subprocess.run(diamond_cmd)
    subprocess.run(["mv", diamond_output_file, diamond_output_dir])
def read_rnasamba_results(file_path):

    results = []

    file_in = open(file_path)
    line = file_in.readline()
    while line:
        line = line.rstrip()
        match = re.search(r"^sequence_name", line)
        if match:
            line = file_in.readline()
            continue

        eles = line.split("\t")
        if not len(eles) == 3:
            line = file_in.readline()
            continue

        transcript_id = eles[0]
        coding_probability = eles[1]
        coding_potential = eles[2]
        results.append([transcript_id, coding_probability, coding_potential])
        line = file_in.readline()
    file_in.close()

    return results


def read_cpc2_results(file_path):

    results = []

    file_in = open(file_path)
    line = file_in.readline()
    while line:
        line = line.rstrip()
        match = re.search(r"^#ID", line)
        if match:
            line = file_in.readline()
            continue

        eles = line.split("\t")
        if not len(eles) == 9:
            line = file_in.readline()
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
        line = file_in.readline()
    file_in.close()

    return results


def read_diamond_results(diamond_output_dir):

    results = []
    diamond_files = glob.glob(diamond_output_dir + "/*.dmdout")
    for file_path in diamond_files:
        file_in = open(file_path)
        line = file_in.readline()
        while line:
            line = line.rstrip()

            eles = line.split("\t")
            if not len(eles) == 12:
                line = file_in.readline()
                continue

            transcript_id = eles[0]
            e_value = eles[10]
            results.append([transcript_id, e_value])
            line = file_in.readline()
    file_in.close()

    return results


def combine_results(rnasamba_results, cpc2_results, diamond_results):

    transcript_ids = {}

    for result in rnasamba_results:
        transcript_id = result[0]
        coding_probability = result[1]
        coding_potential = result[2]

        if transcript_id not in transcript_ids:
            transcript_ids[transcript_id] = [
                coding_probability,
                coding_potential,
            ]

    for result in cpc2_results:
        transcript_id = result[0]
        coding_probability = result[1]
        coding_potential = result[2]
        transcript_length = result[3]
        peptide_length = result[4]
        transcript_ids[transcript_id].extend(
            [
                coding_probability,
                coding_potential,
                transcript_length,
                peptide_length,
            ]
        )

    if diamond_results is not None:
        for result in diamond_results:
            transcript_id = result[0]
            e_value = result[1]
            # There seems to be an issue where there are a small number
            # of sequences that don't make it into the cpc2/rnasamba output
            # Should code in a system for this, but it would be good to
            # understand why it happens to begin with. Seems to be the same
            # number of missing seqs in both, so maybe a shared cut-off
            if transcript_id in transcript_ids:
                transcript_ids[transcript_id].extend([e_value])

    return transcript_ids

def merge_finalise_output_files(
    final_annotation_dir, region_annotation_dir, extension, id_label
):

    gtf_files = glob.glob(region_annotation_dir + "/*" + extension)

    merged_gtf_file = os.path.join(final_annotation_dir, (id_label + "_sel.gtf"))
    merged_cdna_file = os.path.join(final_annotation_dir, (id_label + "_sel.cdna.fa"))
    merged_amino_acid_file = os.path.join(
        final_annotation_dir, (id_label + "_sel.prot.fa")
    )

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
    gtf_out = open(merged_gtf_file, "w+")
    cdna_out = open(merged_cdna_file, "w+")
    amino_acid_out = open(merged_amino_acid_file, "w+")
    for gtf_file in gtf_files:
        logger.info("GTF file: " + gtf_file)
        cdna_seq_index = {}
        amino_acid_seq_index = {}
        cdna_file = gtf_file + ".cdna"
        amino_acid_file = gtf_file + ".prot"
        cdna_in = open(cdna_file)
        amino_acid_in = open(amino_acid_file)
        cdna_seq_index = fasta_to_dict(cdna_in.readlines())
        amino_acid_seq_index = fasta_to_dict(amino_acid_in.readlines())
        cdna_in.close()
        amino_acid_in.close()

        current_gene_id = ""
        gtf_in = open(gtf_file)
        line = gtf_in.readline()
        while line:
            if re.search(r"^#", line):
                line = gtf_in.readline()
                continue

            eles = line.split("\t")
            if not len(eles) == 9:
                line = gtf_in.readline()
                continue

            match = re.search(r'gene_id "([^"]+)".+transcript_id "([^"]+)"', line)
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

            new_gene_id = id_label + "_" + str(gene_id_counter)
            new_transcript_id = id_label + "_" + str(transcript_id_counter)
            line = re.sub(
                'gene_id "' + gene_id + '"',
                ('gene_id "' + new_gene_id + '"'),
                line,
            )
            line = re.sub(
                'transcript_id "' + transcript_id + '"',
                ('transcript_id "' + new_transcript_id + '"'),
                line,
            )
            gtf_out.write(line)
            line = gtf_in.readline()

            if eles[2] == "transcript":
                new_header = ">" + new_transcript_id + "\n"
                cdna_out.write(new_header + cdna_seq_index[transcript_id])

                if transcript_id in amino_acid_seq_index:
                    amino_acid_out.write(new_header + amino_acid_seq_index[transcript_id])

    gtf_out.close()
    cdna_out.close()
    amino_acid_out.close()


def multiprocess_finalise_geneset(cmd):

    print(" ".join(cmd))
    subprocess.run(cmd)

