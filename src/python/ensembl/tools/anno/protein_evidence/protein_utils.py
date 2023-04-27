def run_genblast_align(
    genblast_path,
    convert2blastmask_path,
    makeblastdb_path,
    genblast_dir,
    protein_file,
    masked_genome_file,
    max_intron_length,
    num_threads,
    genblast_timeout_secs,
):

    if not genblast_path:
        genblast_path = config["genblast"]["software"]

    utils.check_exe(genblast_path)

    if not convert2blastmask_path:
        convert2blastmask_path = config["convert2blastmask"]["software"]

    utils.check_exe(convert2blastmask_path)

    if not makeblastdb_path:
        makeblastdb_path = config["makeblastdb"]["software"]

    utils.check_exe(makeblastdb_path)

    utils.create_dir(genblast_dir, None)

    logger.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(genblast_dir, "annotation.gtf")
    if os.path.exists(output_file):
        transcript_count = utils.check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Genblast gtf file exists")
            return
    else:
        logger.info("No gtf file, go on with the analysis")

    genblast_output_file = os.path.join(genblast_dir, "genblast")

    asnb_file = masked_genome_file + ".asnb"
    logger.info("ASNB file: %s" % asnb_file)

    if not os.path.exists("alignscore.txt"):
        shutil.copy(
            os.environ["ENSCODE"] + "/ensembl-anno/support_files/alignscore.txt", "./"
        )
    #        subprocess.run(
    #            [
    #                "cp",
    #                os.environ["ENSCODE"] + "/ensembl-anno/support_files/alignscore.txt",
    #                "./",
    #            ]
    #        )

    if not os.path.exists(masked_genome_file):
        raise IOError("Masked genome file does not exist: %s" % masked_genome_file)

    if not os.path.exists(protein_file):
        raise IOError("Protein file does not exist: %s" % protein_file)

    if not os.path.exists(asnb_file):
        run_convert2blastmask(convert2blastmask_path, masked_genome_file, asnb_file)
    else:
        logger.info("Found an existing asnb, so will skip convert2blastmask")

    if not os.path.exists(asnb_file):
        raise IOError("asnb file does not exist: %s" % asnb_file)

    run_makeblastdb(makeblastdb_path, masked_genome_file, asnb_file)

    batched_protein_files = split_protein_file(protein_file, genblast_dir, 20)

    pool = multiprocessing.Pool(int(num_threads))
    for batched_protein_file in batched_protein_files:
        pool.apply_async(
            multiprocess_genblast,
            args=(
                batched_protein_file,
                masked_genome_file,
                genblast_path,
                genblast_timeout_secs,
                max_intron_length,
            ),
        )
    pool.close()
    pool.join()

    logger.info("Completed running GenBlast")
    logger.info("Combining output into single GTF")
    generate_genblast_gtf(genblast_dir)


def multiprocess_genblast(
    batched_protein_file,
    masked_genome_file,
    genblast_path,
    genblast_timeout_secs,
    max_intron_length,
):

    batch_num = os.path.splitext(batched_protein_file)[0]
    batch_dir = os.path.dirname(batched_protein_file)
    logger.info("Running GenBlast on " + batched_protein_file + ":")

    genblast_cmd = [
        genblast_path,
        "-p",
        "genblastg",
        "-q",
        batched_protein_file,
        "-t",
        masked_genome_file,
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
        batched_protein_file,
    ]

    logger.info(" ".join(genblast_cmd))
    # Using the child process termination as described here:
    # https://alexandra-zaharia.github.io/posts/kill-subprocess
    # -and-its-children-on-timeout-python/
    try:
        p = subprocess.Popen(genblast_cmd, start_new_session=True)
        p.wait(timeout=genblast_timeout_secs)
    except subprocess.TimeoutExpired:
        logger.error("Timeout reached for file:\n" + batched_protein_file)
        subprocess.run(["touch", (batched_protein_file + ".except")])
        os.killpg(os.getpgid(p.pid), signal.SIGTERM)

    files_to_delete = glob.glob(batched_protein_file + "*msk.blast*")
    files_to_delete.append(batched_protein_file)
    for file_to_delete in files_to_delete:
        subprocess.run(["rm", file_to_delete])


def generate_genblast_gtf(genblast_dir):
    logger.info("generate_genblast_gtf")
    file_out_name = os.path.join(genblast_dir, "annotation.gtf")
    file_out = open(file_out_name, "w+")
    genblast_extension = "_1.1c_2.3_s1_0_16_1"
    for root, dirs, files in os.walk(genblast_dir):
        for genblast_file in files:
            genblast_file = os.path.join(root, genblast_file)
            if genblast_file.endswith(".gff"):
                gtf_string = convert_gff_to_gtf(genblast_file)
                file_out.write(gtf_string)
            elif (
                genblast_file.endswith(".fa.blast")
                or genblast_file.endswith(".fa.blast.report")
                or genblast_file.endswith(genblast_extension)
            ):
                # subprocess.run(["rm", genblast_file])
                pathlib.Path(genblast_file).unlink()
    file_out.close()

def split_protein_file(protein_file, protein_output_dir, batch_size):
    if batch_size is None:
        batch_size = 20

    batched_protein_files = []

    for i in range(0, 10):
        utils.create_dir(protein_output_dir, ("bin_" + str(i)))

    file_in = open(protein_file)
    line = file_in.readline()
    seq_count = 0
    batch_count = 0
    current_record = ""
    initial_seq = 1
    while line:
        num_dir = random.randint(0, 9)
        match = re.search(r">(.+)$", line)
        if match and not initial_seq and seq_count % batch_size == 0:
            file_out_name = os.path.join(
                protein_output_dir,
                ("bin_" + str(random.randint(0, 9))),
                (str(batch_count) + ".fa"),
            )
            file_out = open(file_out_name, "w+")
            file_out.write(current_record)
            file_out.close()
            batch_count += 1
            seq_count += 1
            current_record = line
            batched_protein_files.append(file_out_name)
        elif match:
            current_record += line
            initial_seq = 0
            seq_count += 1
        else:
            current_record += line
        line = file_in.readline()
    file_in.close()

    if current_record:
        file_out_name = os.path.join(
            protein_output_dir,
            ("bin_" + str(random.randint(0, 9))),
            (str(batch_count) + ".fa"),
        )
        file_out = open(file_out_name, "w+")
        file_out.write(current_record)
        file_out.close()
        batched_protein_files.append(file_out_name)

    return batched_protein_files


def run_convert2blastmask(convert2blastmask_path, masked_genome_file, asnb_file):

    asnb_file = masked_genome_file + ".asnb"
    logger.info("Running convert2blastmask prior to GenBlast:")
    cmd = [
        convert2blastmask_path,
        "-in",
        masked_genome_file,
        "-parse_seqids",
        "-masking_algorithm",
        "other",
        "-masking_options",
        '"REpeatDetector, default"',
        "-outfmt",
        "maskinfo_asn1_bin",
        "-out",
        asnb_file,
    ]
    logger.info(" ".join(cmd))
    subprocess.run(cmd)
    logger.info("Completed running convert2blastmask")


def run_makeblastdb(makeblastdb_path, masked_genome_file, asnb_file):

    logger.info("Running makeblastdb prior to GenBlast")
    subprocess.run(
        [
            makeblastdb_path,
            "-in",
            masked_genome_file,
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-mask_data",
            asnb_file,
            "-max_file_sz",
            "10000000000",
        ]
    )
    logger.info("Completed running makeblastdb")

