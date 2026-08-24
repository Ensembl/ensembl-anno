import argparse
import logging

def check_transcriptomic_output(gtf_dict, logfile_path) -> None:
    """
    Validate transcriptomic annotation outputs.

    Args:
        main_output_dir: Main output directory.

    Raises:
        OSError: If no transcriptomic annotations are found or the total
            number of lines is below the expected threshold.
    """
    logger = logging.getLogger(__name__)
    logging.basicConfig(filename=logfile_path, encoding='utf-8', level=logging.DEBUG)
    total_lines = 0
    min_lines = 100000
    for tool, annotation_file in gtf_dict.items():
        if annotation_file is None:
            logger.warning(
                "No annotation gtf was generated for %s. This may be expected "
                "(e.g. no long-read data were provided).",
                tool,
            )
            continue
        with annotation_file.open(encoding="utf-8") as file_handle:
            num_lines = sum(1 for _ in file_handle)

        total_lines += num_lines

        logger.info(
            "Found %s lines for %s gtf",
            num_lines,
            tool,
        )
    if total_lines == 0:
        raise OSError(
            "Transcriptomic mode was enabled, but all transcriptomic annotation files are empty."
        )  

    if total_lines <= min_lines:
        raise OSError(
            f"Transcriptomic mode was enabled, but the total number of "
            f"lines across annotation files is below the expected "
            f"minimum. Found: {total_lines}. Minimum allowed: {min_lines}."
        )
    logger.info(
        "Found %s total lines across transcriptomic annotation files. Checks passed.",
        total_lines,
    )



def parse_args(tools):
    parser = argparse.ArgumentParser(description="Arguments for script to check contents of transcriptomic gtfs")
    parser.add_argument("--logfile", help="Path to output logfile recording status of each gtf")
    for tool in tools:
        parser.add_argument(f"--{tool}", help=f"Path to {tool} gtf", required=False, default=None)
    args = parser.parse_args()
    return args

def make_gtf_dict(tools, args):
    gtf_dict = {}
    for tool in tools:
            gtf_dict[tool] = args.tool
    return gtf_dict
    

if __name__ == "__main__":
    tools = ['scallop', 'stringtie', 'minimap2']
    args = parse_args(tools)
    gtf_dict = make_gtf_dict(tools, args)
    check_transcriptomic_output(gtf_dict, args.logfile)
