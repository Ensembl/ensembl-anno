def run_trf_repeats(  # pylint: disable=too-many-locals
    genome_file: typing.Union[pathlib.Path, str],
    trf_path: str,
    main_output_dir: str,
    num_threads: int,
):
    """
    Run trf on genomic slices

    Args:
        genome_file : pathlib.Path
        trf_path : str
        main_output_dir : str
        num_threads: int

    Return:
        gtfs with the masked sequence for each genomic slice
    """

    if not trf_path:
        trf_path = config["trf"]["software"]

    check_exe(trf_path)
    trf_output_dir = pathlib.Path(create_dir(main_output_dir, "trf_output"))
    os.chdir(str(trf_output_dir))
    output_file = trf_output_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "repeat")
        if transcript_count > 0:
            logger.info("Trf gtf file exists, skipping analysis")
            return

    logger.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, 5000)
    slice_ids = create_slice_ids(
        seq_region_lengths, slice_size=1000000, overlap=0, min_length=5000
    )
    match_score = config["trf"]["match_score"]
    mismatch_score = config["trf"]["mismatch_score"]
    delta = config["trf"]["delta"]
    pm = config["trf"]["pm"]  # pylint: disable=invalid-name
    pi = config["trf"]["pi"]  # pylint: disable=invalid-name
    minscore = config["trf"]["minscore"]
    maxperiod = config["trf"]["maxperiod"]

    trf_output_extension = (
        f".{match_score}.{mismatch_score}.{delta}.{pm}.{pi}.{minscore}.{maxperiod}.dat"
    )

    generic_trf_cmd = [
        trf_path,
        None,
        str(match_score),
        str(mismatch_score),
        str(delta),
        str(pm),
        str(pi),
        str(minscore),
        str(maxperiod),
        "-d",
        "-h",
    ]
    logger.info("Running TRF")
    pool = multiprocessing.Pool(int(num_threads))
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_trf,
            args=(
                generic_trf_cmd,
                slice_id,
                genome_file,
                trf_output_dir,
                trf_output_extension,
            ),
        )
    pool.close()
    pool.join()
    slice_output_to_gtf(str(trf_output_dir), ".trf.gtf", 1, "repeat_id", "trf")
    for gtf_file in pathlib.Path(trf_output_dir).glob("*.trf.gtf"):
        gtf_file.unlink()

def multiprocess_trf(  # pylint: disable=too-many-locals
    generic_trf_cmd: list,
    slice_id: str,
    genome_file: pathlib.Path,
    trf_output_dir: pathlib.Path,
    trf_output_extension: str,
):
    """
    Run TRF on multiprocess on genomic slices
    Args:
        generic_trf_cmd:list
        slice_id:str
        genome_file : pathlib.Path
        trf_output_dir : pathlib.Path
        trf_output_extension: str
    """
    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing slice to find tandem repeats with TRF:%s:%s:%s",
        region_name,
        start,
        end,
    )
    seq = get_sequence(region_name, start, end, 1, genome_file, str(trf_output_dir))

    slice_file_name = f"{region_name}.rs{start}.re{end}"
    with tempfile.TemporaryDirectory(dir=trf_output_dir) as tmpdirname:
        region_fasta_file_path = trf_output_dir / tmpdirname / f"{slice_file_name}.fa"
        with open(region_fasta_file_path, "w+") as region_fasta_out:
            region_fasta_out.write(f">{region_name}\n{seq}\n")

        region_results_file_path = trf_output_dir / f"{slice_file_name}.trf.gtf"
        # TRF writes to the current dir, so swtich to the output dir for it
        # os.chdir(str(trf_output_dir))
        trf_output_file_path = pathlib.Path(
            f"{region_fasta_file_path}{trf_output_extension}"
        )
        trf_cmd = generic_trf_cmd.copy()
        trf_cmd[1] = str(region_fasta_file_path)
        logger.info("trf_cmd: %s", trf_cmd)
        # with open(trf_output_file_path, "w+") as trf_out:
        subprocess.run(  # pylint: disable=subprocess-run-check
            trf_cmd, cwd=trf_output_dir / tmpdirname
        )  # pylint: disable=subprocess-run-check
        create_trf_gtf(trf_output_file_path, region_results_file_path, region_name)
        # trf_output_file_path.unlink()
        # region_fasta_file_path.unlink()

def create_trf_gtf(
    trf_output_file_path: pathlib.Path,
    region_results_file_path: pathlib.Path,
    region_name: str,
):  # pylint: disable=too-many-locals
    """
        Read the fasta file and save the content in gtf format

        All the genomic slices are collected in a single gtf output
        Args:
            trf_output_file_path : pathlib.Path
            region_results_file_path : pathlib.Path
            region_name :str

        trf_output_file_path is txt file space delimited where the colummns are
        cols 1+2:  Indices of the repeat relative to the start of the sequence
    col 3:     Period size of the repeat
    col 4:     Number of copies aligned with the consensus pattern
    col 5:     Size of consensus pattern (may differ slightly from the period size)
    col 6:     Percent of matches between adjacent copies overall
    col 7:     Percent of indels between adjacent copies overall
    col 8:     Alignment score
    cols 9-12: Percent composition for each of the four nucleotides
    col 13:    Entropy measure based on percent composition
    col 14:    Consensus sequence
    col 15:    Repeat sequence
    """
    with open(trf_output_file_path, "r") as trf_in, open(
        region_results_file_path, "w+"
    ) as trf_out:
        repeat_count = 1
        for line in trf_in:
            result_match = re.search(r"^\d+", line)
            if result_match:
                results = line.split()
                if not len(results) == 15:
                    continue
                start = results[0]
                end = results[1]
                period = float(results[2])
                copy_number = float(results[3])
                percent_matches = float(results[5])
                score = float(results[7])
                repeat_consensus = results[13]
                if (  # pylint: disable=too-many-boolean-expressions
                    score < 50
                    and percent_matches >= 80
                    and copy_number > 2
                    and period < 10
                ) or (copy_number >= 2 and percent_matches >= 70 and score >= 50):
                    gtf_line = (
                        f"{region_name}\tTRF\trepeat\t{start}\t{end}\t.\t+\t.\t"
                        f'repeat_id "{repeat_count}"; score "{score}"; '
                        f'repeat_consensus "{repeat_consensus}";\n'
                    )

                    trf_out.write(gtf_line)
                    repeat_count += 1

