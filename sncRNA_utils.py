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

def run_trnascan_regions(
    genome_file: Union[pathlib.Path, str],
    trnascan_path,
    trnascan_filter_path,
    main_output_dir,
    num_threads: int,
):
    if not trnascan_path:
        trnascan_path = "/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/bin/tRNAscan-SE"
    check_exe(trnascan_path)
    logger.info("trnascan_path: %s" % trnascan_path)

    if not trnascan_filter_path:
        trnascan_filter_path = "/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/bin/EukHighConfidenceFilter"
    # check_exe(trnascan_filter_path)
    check_file(trnascan_filter_path)
    logger.info("trnascan_filter_path: %s" % trnascan_filter_path)

    trnascan_output_dir = create_dir(main_output_dir, "trnascan_output")

    output_file = trnascan_output_dir / "annotation.gtf"
    if os.path.exists(output_file):
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Trnascan gtf file already exists, skipping analysis")
            return

    logger.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, min_seq_length=5000)
    slice_ids = create_slice_ids(
        seq_region_lengths, slice_size=1_000_000, overlap=0, min_length=5000
    )

    generic_trnascan_cmd = [
        trnascan_path,
        None,
        "-o",
        None,
        "-f",
        None,
        "-H",
        "-q",
        "--detail",
        "-Q",
    ]
    logger.info("Running tRNAscan-SE processes")
    pool = multiprocessing.Pool(num_threads)
    tasks = []
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_trnascan,
            args=(
                generic_trnascan_cmd,
                slice_id,
                genome_file,
                trnascan_filter_path,
                trnascan_output_dir,
            ),
        )

    pool.close()
    pool.join()
    slice_output_to_gtf(trnascan_output_dir, ".trna.gtf", 1, None, None)


def multiprocess_trnascan(
    generic_trnascan_cmd,
    slice_id,
    genome_file: Union[pathlib.Path, str],
    trnascan_filter_path,
    trnascan_output_dir: Union[pathlib.Path, str],
):
    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing slice to find tRNAs using tRNAscan-SE: %s:%s:%s"
        % (region_name, start, end)
    )
    seq = get_sequence(
        seq_region=region_name,
        start=start,
        end=end,
        strand=1,
        fasta_file=genome_file,
        output_dir=trnascan_output_dir,
    )

    slice_file_name = f"{region_name}.rs{start}.re{end}"
    region_fasta_file_name = f"{slice_file_name}.fa"
    region_fasta_file_path = trnascan_output_dir / region_fasta_file_name
    with open(region_fasta_file_path, "w+") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")

    region_results_file_name = f"{slice_file_name}.trna.gtf"
    region_results_file_path = trnascan_output_dir / region_results_file_name

    trnascan_output_file_path = f"{region_fasta_file_path}.trna"
    trnascan_ss_output_file_path = f"{trnascan_output_file_path}.ss"

    # The filter takes an output dir and a prefix and then uses those to make a path to a .out file
    trnascan_filter_file_prefix = f"{region_fasta_file_name}.filt"
    trnascan_filter_file_name = f"{trnascan_filter_file_prefix}.out"
    trnascan_filter_log_file_name = f"{trnascan_filter_file_prefix}.log"
    trnascan_filter_ss_file_name = f"{trnascan_filter_file_prefix}.ss"
    trnascan_filter_file_path = trnascan_output_dir / trnascan_filter_file_name
    trnascan_filter_log_file_path = trnascan_output_dir / trnascan_filter_log_file_name
    trnascan_filter_ss_file_path = trnascan_output_dir / trnascan_filter_ss_file_name

    trnascan_cmd = generic_trnascan_cmd.copy()
    trnascan_cmd[1] = region_fasta_file_path
    trnascan_cmd[3] = trnascan_output_file_path
    trnascan_cmd[5] = trnascan_ss_output_file_path

    logger.info("tRNAscan-SE command:\n%s" % " ".join(trnascan_cmd))
    subprocess.run(trnascan_cmd)

    # If we have a blank output file at this point we want to stop and remove whatever files
    # have been created instead of moving onto the filter
    if os.stat(trnascan_output_file_path).st_size == 0:
        os.remove(trnascan_output_file_path)
        os.remove(region_fasta_file_path)
        if os.path.exists(trnascan_ss_output_file_path):
            os.remove(trnascan_ss_output_file_path)
        return

    filter_cmd = [
        trnascan_filter_path,
        "--result",
        trnascan_output_file_path,
        "--ss",
        trnascan_ss_output_file_path,
        "--output",
        trnascan_output_dir,
        "--prefix",
        trnascan_filter_file_prefix,
    ]
    logger.info("tRNAscan-SE filter command:\n%s" % " ".join(filter_cmd))
    subprocess.run(filter_cmd)

    create_trnascan_gtf(
        region_results_file_path, trnascan_filter_file_path, region_name
    )
    if os.path.exists(trnascan_output_file_path):
        os.remove(trnascan_output_file_path)
    if os.path.exists(trnascan_ss_output_file_path):
        os.remove(trnascan_ss_output_file_path)
    if os.path.exists(trnascan_filter_file_path):
        os.remove(trnascan_filter_file_path)
    if os.path.exists(trnascan_filter_log_file_path):
        os.remove(trnascan_filter_log_file_path)
    if os.path.exists(trnascan_filter_ss_file_path):
        os.remove(trnascan_filter_ss_file_path)
    if os.path.exists(region_fasta_file_path):
        os.remove(region_fasta_file_path)


def create_trnascan_gtf(
    region_results_file_path, trnascan_filter_file_path, region_name
):
    with open(trnascan_filter_file_path, "r") as trna_in, open(
        region_results_file_path, "w+"
    ) as trna_out:
        gene_counter = 1
        for line in trna_in:
            result_match = re.search(r"^" + region_name, line)
            if result_match:
                results = line.split()
                start = int(results[2])
                end = int(results[3])
                trna_type = results[4]
                score = results[8]

                strand = "+"
                if start > end:
                    strand = "-"
                    temp_end = start
                    start = end
                    end = temp_end

                biotype = "tRNA_pseudogene"
                high_confidence_match = re.search(r"high confidence set", line)
                if high_confidence_match:
                    biotype = "tRNA"

                transcript_string = f'{region_name}\ttRNAscan\ttranscript\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gene_counter}"; transcript_id "{gene_counter}"; biotype "{biotype}";\n'
                exon_string = f'{region_name}\ttRNAscan\texon\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gene_counter}"; transcript_id "{gene_counter}"; exon_number "1"; biotype "{biotype}";\n'

                trna_out.write(transcript_string)
                trna_out.write(exon_string)
                gene_counter += 1

def run_cmsearch_regions(
    genome_file: Union[pathlib.Path, str],
    cmsearch_path,
    rfam_cm_db_path,
    rfam_seeds_file_path,
    rfam_accession_file,
    main_output_dir,
    num_threads: int,
):
    if not cmsearch_path:
        cmsearch_path = "cmsearch"

    check_exe(cmsearch_path)
    rfam_output_dir = create_dir(main_output_dir, "rfam_output")

    os.chdir(rfam_output_dir)

    output_file = rfam_output_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("CMsearch GTF file already exists, skipping analysis")
            return

    rfam_dbname = "Rfam"
    rfam_user = "rfamro"
    rfam_host = "mysql-rfam-public.ebi.ac.uk"
    rfam_port = "4497"
    #  rfam_accession_query_cmd = ["mysql -h", rfam_host,"-u",rfam_user,"-P",rfam_port,"-NB -e",rfam_dbname,"'select rfam_acc FROM (SELECT DISTINCT f.rfam_acc, f.rfam_id, f.type, f.description, f.gathering_cutoff, f.trusted_cutoff FROM full_region fr, rfamseq rf, taxonomy tx, family f WHERE rf.ncbi_id = tx.ncbi_id AND f.rfam_acc = fr.rfam_acc AND fr.rfamseq_acc = rf.rfamseq_acc AND LOWER(tx.tax_string) LIKE \'%" + clade + "%\' AND (f.type LIKE \'%snRNA%\' OR f.type LIKE \'%rRNA%\' OR LOWER(f.rfam_id) LIKE \'%rnase%\' OR LOWER(f.rfam_id) LIKE \'%vault%\' OR LOWER(f.rfam_id) LIKE \'%y_rna%\' OR f.rfam_id LIKE \'%Metazoa_SRP%\') AND is_significant = 1) AS TEMP WHERE rfam_id NOT LIKE \'%bacteria%\' AND rfam_id NOT LIKE \'%archaea%\' AND rfam_id NOT LIKE \'%microsporidia%\';'"]

    # mysql -hmysql-rfam-public.ebi.ac.uk -urfamro -P4497 Rfam -NB -e "select rfam_acc FROM (SELECT DISTINCT f.rfam_acc, f.rfam_id, f.type, f.description, f.gathering_cutoff, f.trusted_cutoff FROM full_region fr, rfamseq rf, taxonomy tx, family f WHERE rf.ncbi_id = tx.ncbi_id AND f.rfam_acc = fr.rfam_acc AND fr.rfamseq_acc = rf.rfamseq_acc AND LOWER(tx.tax_string) LIKE '%insect%' AND (f.type LIKE '%snRNA%' OR f.type LIKE '%rRNA%' OR LOWER(f.rfam_id) LIKE '%rnase%' OR LOWER(f.rfam_id) LIKE '%vault%' OR LOWER(f.rfam_id) LIKE '%y_rna%' OR f.rfam_id LIKE '%Metazoa_SRP%') AND is_significant = 1) AS TEMP WHERE rfam_id NOT LIKE '%bacteria%' AND rfam_id NOT LIKE '%archaea%' AND rfam_id NOT LIKE '%microsporidia%';"

    #  rfam_accession_file = '/hps/nobackup2/production/ensembl/fergal/production/test_runs/non_verts/butterfly/rfam_insect_ids.txt'
    # rfam_accession_file = os.path.join(main_output_dir,'rfam_accessions.txt')
    rfam_cm_db_path = (
        "/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.0/Rfam.cm"
    )
    rfam_seeds_file_path = (
        "/hps/nobackup/flicek/ensembl/genebuild/blastdb/ncrna/Rfam_14.0/Rfam.seed"
    )
    rfam_selected_models_file = rfam_output_dir / "rfam_models.cm"
    with open(rfam_accession_file) as rfam_accessions_in:
        rfam_accessions = rfam_accessions_in.read().splitlines()

    with open(rfam_cm_db_path, "r") as rfam_cm_in:
        rfam_data = rfam_cm_in.read()

    rfam_models = rfam_data.split("//\n")
    with open(rfam_selected_models_file, "w+") as rfam_cm_out:
        for model in rfam_models:
            # The Rfam.cm file has INFERNAL and HMMR models, both are needed at this point
            # Later we just want the INFERNAL ones for looking at thresholds
            match = re.search(r"(RF\d+)", model)
            if match:
                model_accession = match.group(1)
                if model_accession in rfam_accessions:
                    rfam_cm_out.write(f"{model}//\n")

    seed_descriptions = get_rfam_seed_descriptions(rfam_seeds_file_path)
    cv_models = extract_rfam_metrics(rfam_selected_models_file)

    logger.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, min_seq_length=5000)
    slice_ids = create_slice_ids(
        seq_region_lengths, slice_size=1_000_000, overlap=0, min_length=5000
    )

    generic_cmsearch_cmd = [
        cmsearch_path,
        "--rfam",
        "--cpu",
        "1",
        "--nohmmonly",
        "--cut_ga",
        "--tblout",
    ]
    logger.info("Running Rfam")
    pool = multiprocessing.Pool(num_threads)
    results = []
    failed_slice_ids = []
    memory_limit = 3 * 1024**3
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_cmsearch,
            args=(
                generic_cmsearch_cmd,
                slice_id,
                genome_file,
                rfam_output_dir,
                rfam_selected_models_file,
                cv_models,
                seed_descriptions,
                memory_limit,
            ),
        )
    pool.close()
    pool.join()

    # Need to figure something more automated out here. At the moment it's just limiting to 5 cores and 5GB vram
    # Ideally we could look at the amount of mem requested and put something like 10GB per core and then figure
    # out how many cores to use (obviously not using more than the amount specified)
    memory_limit = 5 * 1024**3
    if num_threads > 5:
        num_threads = 5
    pool = multiprocessing.Pool(num_threads)
    for exception_file_path in rfam_output_dir.glob("*.rfam.except"):
        logger.info("Running himem job for failed region:\n%s" % exception_file_path)
        exception_file_name = exception_file_path.name
        match = re.search(r"(.+)\.rs(\d+)\.re(\d+)\.", exception_file_name)
        if match:
            except_region = match.group(1)
            except_start = match.group(2)
            except_end = match.group(3)
            except_slice_id = [except_region, except_start, except_end]
            pool.apply_async(
                multiprocess_cmsearch,
                args=(
                    generic_cmsearch_cmd,
                    except_slice_id,
                    genome_file,
                    rfam_output_dir,
                    rfam_selected_models_file,
                    cv_models,
                    seed_descriptions,
                    memory_limit,
                ),
            )
    pool.close()
    pool.join()

    slice_output_to_gtf(rfam_output_dir, ".rfam.gtf", 1, None, None)


def multiprocess_cmsearch(
    generic_cmsearch_cmd,
    slice_id,
    genome_file: Union[pathlib.Path, str],
    rfam_output_dir: pathlib.Path,
    rfam_selected_models_file,
    cv_models,
    seed_descriptions,
    memory_limit,
):
    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing Rfam data using cmsearch against slice: %s:%s:%s"
        % (region_name, start, end)
    )
    seq = get_sequence(
        seq_region=region_name,
        start=start,
        end=end,
        strand=1,
        fasta_file=genome_file,
        output_dir=rfam_output_dir,
    )

    slice_file_name = f"{region_name}.rs{start}.re{end}"
    region_fasta_file_name = f"{slice_file_name}.fa"
    region_fasta_file_path = rfam_output_dir / region_fasta_file_name
    with open(region_fasta_file_path, "w+") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")

    region_tblout_file_name = f"{slice_file_name}.tblout"
    region_tblout_file_path = rfam_output_dir / region_tblout_file_name
    region_results_file_name = f"{slice_file_name}.rfam.gtf"
    region_results_file_path = rfam_output_dir / region_results_file_name

    exception_results_file_name = f"{slice_file_name}.rfam.except"
    exception_results_file_path = rfam_output_dir / exception_results_file_name

    cmsearch_cmd = generic_cmsearch_cmd.copy()
    cmsearch_cmd.append(region_tblout_file_path)
    cmsearch_cmd.append(rfam_selected_models_file)
    cmsearch_cmd.append(region_fasta_file_path)

    if memory_limit is not None:
        cmsearch_cmd = prlimit_command(cmsearch_cmd, memory_limit)

    logger.info("cmsearch_cmd: %s" % " ".join(cmsearch_cmd))

    return_value = None
    try:
        return_value = subprocess.check_output(cmsearch_cmd)
    except subprocess.CalledProcessError as ex:
        # Note that writing to file was the only option here. If return_value was passed back, eventually it would clog the
        # tiny pipe that is used by the workers to send info back. That would mean that it would eventually just be one
        # worked running at a time
        logger.error(
            "Issue processing the following region with cmsearch: %s %s-%s. Return value: %s"
            % (region_name, start, end, return_value)
        )
        with open(exception_results_file_path, "w+") as exception_out:
            exception_out.write(f"{region_name} {start} {end}\n")
        os.remove(region_fasta_file_path)
        os.remove(region_tblout_file_path)
        return

    initial_table_results = parse_rfam_tblout(region_tblout_file_path, region_name)
    unique_table_results = remove_rfam_overlap(initial_table_results)
    filtered_table_results = filter_rfam_results(unique_table_results, cv_models)
    if filtered_table_results:
        create_rfam_gtf(
            filtered_results=filtered_table_results,
            cm_models=cv_models,
            descriptions=seed_descriptions,
            region_name=region_name,
            region_results_file_path=region_results_file_path,
            genome_file=genome_file,
            rfam_output_dir=rfam_output_dir,
        )
    os.remove(region_fasta_file_path)
    os.remove(region_tblout_file_path)
    gc.collect()


def get_rfam_seed_descriptions(rfam_seeds_path):
    descriptions = {}

    # NOTE: for some reason the decoder breaks on the seeds file, so I have made this ignore errors
    with open(rfam_seeds_path, encoding="utf-8", errors="ignore") as rfam_seeds_in:
        for line in rfam_seeds_in:
            seed = line.rstrip()

            domain_match = re.search("^\#=GF AC   (.+)", seed)
            if domain_match:
                domain = domain_match.group(1)
                descriptions[domain] = {}
                continue

            description_match = re.search("^\#=GF DE   (.+)", seed)
            if description_match:
                description = description_match.group(1)
                descriptions[domain]["description"] = description
                continue

            name_match = re.search("^\#=GF ID   (.+)", seed)
            if name_match:
                name = name_match.group(1)
                descriptions[domain]["name"] = name
                continue

            type_match = re.search("^\#=GF TP   Gene; (.+)", seed)
            if type_match:
                rfam_type = type_match.group(1)
                descriptions[domain]["type"] = rfam_type
                continue

    return descriptions


def extract_rfam_metrics(rfam_selected_models_file):
    with open(rfam_selected_models_file, "r") as rfam_cm_in:
        rfam_data = rfam_cm_in.read()

    rfam_models = rfam_data.split("//\n")
    parsed_cm_data = {}
    for model in rfam_models:
        temp = model.split("\n")
        model_name_match = re.search(r"NAME\s+(\S+)", model)
        match_infernal = re.search(r"INFERNAL", model)
        if model_name_match and match_infernal:
            model_name = model_name_match.group(1)
            parsed_cm_data[model_name] = {}
            for line in temp:
                name_match = re.search(r"^NAME\s+(\S+)", line)
                if name_match:
                    parsed_cm_data[model_name]["-name"] = name_match.group(1)
                    continue

                description_match = re.search(r"^DESC\s+(\S+)", line)
                if description_match:
                    parsed_cm_data[model_name][
                        "-description"
                    ] = description_match.group(1)
                    continue

                length_match = re.search(r"^CLEN\s+(\d+)", line)
                if length_match:
                    parsed_cm_data[model_name]["-length"] = length_match.group(1)
                    continue

                max_length_match = re.search(r"^W\s+(\d+)", line)
                if max_length_match:
                    parsed_cm_data[model_name]["-maxlength"] = max_length_match.group(1)
                    continue

                accession_match = re.search(r"^ACC\s+(\S+)", line)
                if accession_match:
                    parsed_cm_data[model_name]["-accession"] = accession_match.group(1)
                    continue

                threshold_match = re.search(r"^GA\s+(\d+)", line)
                if threshold_match:
                    parsed_cm_data[model_name]["-threshold"] = threshold_match.group(1)
                    continue

    return parsed_cm_data


def parse_rfam_tblout(region_tblout_file_path, region_name):
    parsed_results = []
    with open(region_tblout_file_path, "r") as rfam_tbl_in:
        for result in rfam_tbl_in:
            parsed_tbl_data = {}
            if not re.match(r"^" + region_name, result):
                continue

            hit = result.split()
            accession = hit[3]
            target_name = hit[0]
            query_name = hit[2]
            hstart = hit[5]
            hend = hit[6]
            start = hit[7]
            end = hit[8]
            if hit[9] == "+":
                strand = 1
            else:
                strand = -1
            evalue = hit[15]
            score = hit[14]

            parsed_tbl_data["accession"] = accession
            parsed_tbl_data["start"] = start
            parsed_tbl_data["end"] = end
            parsed_tbl_data["strand"] = strand
            parsed_tbl_data["query_name"] = query_name
            parsed_tbl_data["score"] = score
            parsed_results.append(parsed_tbl_data)

    return parsed_results


# NOTE some of the code above and the code commented out here is to do with creating
# a DAF. As we don't have a python concept of this I'm leaving it out for the moment
# but the code below is a reference
#    my $daf = Bio::EnsEMBL::DnaDnaAlignFeature->new(
#      -slice          => $self->queries,
#      -start          => $strand == 1 ? $start : $end,
#      -end            => $strand == 1 ? $end : $start,
#      -strand         => $strand,
#      -hstart         => $hstart,
#      -hend           => $hend,
#      -hstrand        => $strand,
#      -score          => $score,
#      -hseqname       => length($target_name) > 39 ? substr($target_name, 0, 39) : $target_name,,
#      -p_value  => $evalue,
#      -align_type => 'ensembl',
#      -cigar_string  => abs($hend - $hstart) . "M",
#   );


def remove_rfam_overlap(parsed_tbl_data):
    excluded_structures = {}
    chosen_structures = []
    for structure_x in parsed_tbl_data:
        chosen_structure = structure_x
        structure_x_start = int(structure_x["start"])
        structure_x_end = int(structure_x["end"])
        structure_x_score = float(structure_x["score"])
        structure_x_accession = structure_x["accession"]
        structure_x_string = "{structure_x_start}:{structure_x_end}:{structure_x_score}:{structure_x_accession}"
        for structure_y in parsed_tbl_data:
            structure_y_start = int(structure_y["start"])
            structure_y_end = int(structure_y["end"])
            structure_y_score = float(structure_y["score"])
            structure_y_accession = structure_y["accession"]
            structure_y_string = "{structure_y_start}:{structure_y_end}:{structure_y_score}:{structure_y_accession}"
            if structure_y_string in excluded_structures:
                continue

            if (
                structure_x_start <= structure_y_end
                and structure_x_end >= structure_y_start
            ):
                if structure_x_score < structure_y_score:
                    chosen_structure = structure_y
                    excluded_structures[structure_x_string] = 1
                else:
                    excluded_structures[structure_y_string] = 1

        chosen_structures.append(chosen_structure)

    return chosen_structures


def filter_rfam_results(unfiltered_tbl_data, cv_models):
    filtered_results = []
    for structure in unfiltered_tbl_data:
        query = structure["query_name"]
        if query in cv_models:
            threshold = cv_models[query]["-threshold"]
            if query == "LSU_rRNA_eukarya":
                threshold = 1700

            elif query == "LSU_rRNA_archaea":
                continue

            elif query == "LSU_rRNA_bacteria":
                continue

            elif query == "SSU_rRNA_eukarya":
                threshold = 1600

            elif query == "5_8S_rRNA":
                threshold = 85

            elif query == "5S_rRNA":
                threshold = 75

            if threshold and float(structure["score"]) >= float(threshold):
                filtered_results.append(structure)

    return filtered_results


# NOTE: The below are notes from the perl code (which has extra code) about possible improvements that are not implemented there
# Although not included in RefSeq filters, additional filters that consider sizes and score_to_size ratios can be applied
# in future work to further exclude FPs
#
# my $is_valid_size = $mapping_length > $min_length && $mapping_length < $max_length ? 1 : 0;
# my $score_size_ratio = $result->{'score'} / $mapping_length;


def create_rfam_gtf(
    filtered_results,
    cm_models,
    descriptions,
    region_name,
    region_results_file_path,
    genome_file: Union[pathlib.Path, str],
    rfam_output_dir,
):
    if not filtered_results:
        return

    with open(region_results_file_path, "w+") as rfam_gtf_out:
        gene_counter = 1
        for structure in filtered_results:
            query = structure["query_name"]
            accession = structure["accession"]
            if query in cm_models:
                model = cm_models[query]
                if accession in descriptions:
                    description = descriptions[accession]
                    if "type" in description:
                        rfam_type = description["type"]
                    else:
                        description = None
                        rfam_type = "misc_RNA"
                domain = structure["query_name"]
                padding = model["-length"]
                # rfam_type = description['type']
                gtf_strand = structure["strand"]
                rnafold_strand = structure["strand"]
                if gtf_strand == 1:
                    start = structure["start"]
                    end = structure["end"]
                    gtf_strand = "+"
                else:
                    start = structure["end"]
                    end = structure["start"]
                    score = structure["score"]
                    gtf_strand = "-"
                    rnafold_strand = -1

                biotype = "misc_RNA"
                if re.match(r"^snRNA;", rfam_type):
                    biotype = "snRNA"
                if re.match(r"^snRNA; snoRNA", rfam_type):
                    biotype = "snoRNA"
                if re.match(r"^snRNA; snoRNA; scaRNA;", rfam_type):
                    biotype = "scaRNA"
                if re.match(r"rRNA;", rfam_type):
                    biotype = "rRNA"
                if re.match(r"antisense;", rfam_type):
                    biotype = "antisense"
                if re.match(r"antitoxin;", rfam_type):
                    biotype = "antitoxin"
                if re.match(r"ribozyme;", rfam_type):
                    biotype = "ribozyme"
                if re.match(r"" + domain, rfam_type):
                    biotype = domain
                if re.match(r"" + domain, rfam_type):
                    biotype = "Vault_RNA"
                if re.match(r"" + domain, rfam_type):
                    biotype = "Y_RNA"

                rna_seq = get_sequence(
                    seq_region=region_name,
                    start=start,
                    end=end,
                    strand=rnafold_strand,
                    fasta_file=genome_file,
                    output_dir=rfam_output_dir,
                )
                valid_structure = check_rnafold_structure(rna_seq, rfam_output_dir)

                if not valid_structure:
                    continue

                transcript_string = f'{region_name}\tRfam\ttranscript\t{start}\t{end}\t.\t{gtf_strand}\t.\tgene_id "{gene_counter}"; transcript_id "{gene_counter}"; biotype "{biotype}";\n'
                exon_string = f'{region_name}\tRfam\texon\t{start}\t{end}\t.\t{gtf_strand}\t.\tgene_id "{gene_counter}"; transcript_id "{gene_counter}"; exon_number "1"; biotype "{biotype}";\n'

                rfam_gtf_out.write(transcript_string)
                rfam_gtf_out.write(exon_string)
                gene_counter += 1


def check_rnafold_structure(seq, rfam_output_dir):
    # Note there's some extra code in the RNAfold Perl module for encoding the structure into an attrib
    # Could consider implementing this when running for loading into an Ensembl db
    structure = 0
    with tempfile.NamedTemporaryFile(
        mode="w+t", delete=False, dir=rfam_output_dir
    ) as rna_temp_in:
        rna_temp_in.write(f">seq1\n{seq}\n")
        rna_in_file_path = rna_temp_in.name

    rnafold_cmd = ["RNAfold", "--infile", rna_in_file_path]
    rnafold_output = subprocess.Popen(rnafold_cmd, stdout=subprocess.PIPE)
    for line in io.TextIOWrapper(rnafold_output.stdout, encoding="utf-8"):
        match = re.search(r"([().]+)\s\(\s*(-*\d+.\d+)\)\n$", line)
        if match:
            structure = match.group(1)
            score = match.group(2)
            break
    rna_temp_in.close()
    os.remove(rna_in_file_path)

    return structure
