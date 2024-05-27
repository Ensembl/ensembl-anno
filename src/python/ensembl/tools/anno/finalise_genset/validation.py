def validate_coding_transcripts(
    cdna_file,
    amino_acid_file,
    validation_dir,
    validation_type,
    diamond_validation_db,
    gtf_file,
    num_threads,
):

    logging.info("Running CDS validation with RNAsamba and CPC2")
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
    logging.info(" ".join(rnasamba_cmd))
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
    logging.info(" ".join(cpc2_cmd))
    subprocess.run(cpc2_cmd)
    cpc2_output_path = cpc2_output_path + ".txt"

    check_file(rnasamba_output_path)
    check_file(cpc2_output_path)

    logging.info("diamond validation")
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

    logging.info("read results")
    rnasamba_results = read_rnasamba_results(rnasamba_output_path)
    cpc2_results = read_cpc2_results(cpc2_output_path)
    combined_results = combine_results(rnasamba_results, cpc2_results, diamond_results)
    logging.info("read gtf genes")
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
    logging.info("Running diamond on " + batched_protein_file + ":")

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

    logging.info(" ".join(diamond_cmd))
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
def update_gtf_genes(parsed_gtf_genes, combined_results, validation_type):
    output_lines = []

    for gene_id in parsed_gtf_genes.keys():
        transcript_ids = parsed_gtf_genes[gene_id].keys()
        for transcript_id in transcript_ids:
            transcript_line = parsed_gtf_genes[gene_id][transcript_id]["transcript"]
            single_cds_exon_transcript = 0
            translation_match = re.search(
                r'; translation_coords "([^"]+)";', transcript_line
            )
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

            avg_coding_probability = (
                rnasamba_coding_probability + cpc2_coding_probability
            ) / 2
            max_coding_probability = max(
                rnasamba_coding_probability, cpc2_coding_probability
            )

            match = re.search(r'; biotype "([^"]+)";', transcript_line)
            biotype = match.group(1)
            if biotype == "busco" or biotype == "protein":
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
                    (
                        rnasamba_coding_potential == "coding"
                        or cpc2_coding_potential == "coding"
                    )
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
                if (
                    diamond_e_value is not None
                    and peptide_length >= min_single_exon_pep_length
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    (
                        rnasamba_coding_potential == "coding"
                        and cpc2_coding_potential == "coding"
                    )
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
                    (
                        rnasamba_coding_potential == "coding"
                        or cpc2_coding_potential == "coding"
                    )
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
                    transcript_line = re.sub(
                        ' translation_coords "[^"]+";', "", transcript_line
                    )
                else:
                    continue

            output_lines.append(transcript_line)
            output_lines.extend(exon_lines)

    return output_lines
def read_gtf_genes(gtf_file):
    gtf_genes = {}
    gtf_in = open(gtf_file)
    line = gtf_in.readline()
    while line:
        eles = line.split("\t")
        if not len(eles) == 9:
            line = gtf_in.readline()
            continue

        match = re.search(r'gene_id "([^"]+)".+transcript_id "([^"]+)"', line)

        if not match:
            line = gtf_in.readline()
            continue

        gene_id = match.group(1)
        transcript_id = match.group(2)
        feature_type = eles[2]
        if gene_id not in gtf_genes:
            gtf_genes[gene_id] = {}
        if feature_type == "transcript":
            gtf_genes[gene_id][transcript_id] = {}
            gtf_genes[gene_id][transcript_id]["transcript"] = line
            gtf_genes[gene_id][transcript_id]["exons"] = []
        elif feature_type == "exon":
            gtf_genes[gene_id][transcript_id]["exons"].append(line)
        line = gtf_in.readline()
    gtf_in.close()

    return gtf_genes
