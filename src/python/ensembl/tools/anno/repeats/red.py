def run_red(
    red_path: str, main_output_dir: str, genome_file: typing.Union[pathlib.Path, str]
) -> str:
    """
    Run Red on genome file

    Args:
        red_path : str
        main_output_dir : pathlib.Path
        genome_file : pathlib.Path

    Return:
        masked genome file
    """
    if not red_path:
        red_path = config["red"]["software"]
    genome_file = pathlib.Path(genome_file)
    check_exe(red_path)
    red_dir = pathlib.Path(create_dir(main_output_dir, "red_output"))
    red_mask_dir = pathlib.Path(create_dir(red_dir, "mask_output"))
    red_repeat_dir = pathlib.Path(create_dir(red_dir, "repeat_output"))
    red_genome_dir = pathlib.Path(create_dir(red_dir, "genome_dir"))

    sym_link_genome_cmd = "ln -s " + str(genome_file) + " " + str(red_genome_dir)
    logger.info(genome_file)
    genome_file_name = genome_file.name
    red_genome_file = red_genome_dir / genome_file_name
    genome_file_stem = genome_file.stem
    masked_genome_file = red_mask_dir / f"{genome_file_stem}.msk"
    repeat_coords_file = red_repeat_dir / f"{genome_file_stem}.rpt"
    gtf_output_file_path = red_dir / "annotation.gtf"

    if masked_genome_file.exists():
        logger.warning(
            "Masked Genome file already found on the path to the Red mask output dir. \
            Will not create a new file"
        )
        create_red_gtf(repeat_coords_file, gtf_output_file_path)
        return str(masked_genome_file)
    if red_genome_file.exists():
        logger.warning(
            "Unmasked genome file already found on the path to the Red genome dir, \
            will not create a sym link"
        )

    else:
        logger.info(
            "Preparing to sym link the genome file to the Red genome dir. Cmd\n %s",
            sym_link_genome_cmd,
        )
        # subprocess.run(["ln", "-s", genome_file, red_genome_dir])
        red_genome_file.symlink_to(genome_file)
    if not red_genome_file.exists():
        logger.error(
            "Could not find the genome file in the Red genome dir or sym link \
            to the original file. Path expected:\n%s",
            red_genome_file,
        )

    logger.info("Running Red, this may take some time depending on the genome size")
    subprocess.run(
        [
            red_path,
            "-gnm",
            red_genome_dir,
            "-msk",
            red_mask_dir,
            "-rpt",
            red_repeat_dir,
        ],
        check=True,
    )

    logger.info("Completed running Red")

    create_red_gtf(repeat_coords_file, gtf_output_file_path)
    reveal_locals()
    return str(masked_genome_file)
def create_red_gtf(repeat_coords_file: pathlib.Path, gtf_output_file_path: pathlib.Path):
    """
    Create Red gtf file  from masked  genome file

    Args:
        repeat_coords_file: pathlib.Path
        gtf_output_file_path : pathlib.Path
    """
    with open(repeat_coords_file, "r") as red_in, open(
        gtf_output_file_path, "w+"
    ) as red_out:
        for repeat_id, line in enumerate(red_in, start=1):
            result_match = re.search(r"^\>(.+)\:(\d+)\-(\d+)", line)
            if result_match:
                region_name = result_match.group(1)
                # Note that Red is 0-based, so add 1
                start = int(result_match.group(2)) + 1
                end = int(result_match.group(3)) + 1
                gtf_line = (
                    f"{region_name}\tRed\trepeat\t{start}\t"
                    f'{end}\t.\t+\t.\trepeat_id "{repeat_id}";\n'
                )
                red_out.write(gtf_line)

