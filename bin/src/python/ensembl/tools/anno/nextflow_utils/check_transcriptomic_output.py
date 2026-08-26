import argparse
import logging

def check_transcriptomic_output(args):
    """
    Check transcriptomic annotation outputs.

    Args:
        main_output_dir: Main output directory.

    Raises:
        OSError: If no transcriptomic annotations are found or the total
            number of lines is below the expected threshold.
    """
    logger = logging.getLogger(__name__)
    logging.basicConfig(filename=args.logfile, encoding='utf-8', level=logging.DEBUG)
    total_lines = 0
    min_lines = args.min_lines
    gtfs = args.gtfs
    for annotation_file in gtfs: 
        with open(annotation_file, 'r', encoding="utf-8") as file_handle:
            num_lines = sum(1 for _ in file_handle)

        total_lines += num_lines

        logger.info(
            "Found %s lines for %s",
            num_lines,
            annotation_file,
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

    with open(args.checks_passed_file, 'w') as success_handle:
        success_handle.write('all checks passed, success')


def parse_args():
    parser = argparse.ArgumentParser(description="Arguments for script to check contents of transcriptomic gtfs")
    parser.add_argument("--gtfs", nargs='+', help="Path to gtfs")
    parser.add_argument("--logfile", help="Path to output logfile recording status of each gtf")
    parser.add_argument("--min_lines", type=int, default=100000, help="Minimum total number of gtf lines expected across all tools")
    parser.add_argument("--checks_passed_file", help="Path to write file to if all checks pass. Useful to ensure nextflow behaves appropriately.")
    args = parser.parse_args()
    return args
    

if __name__ == "__main__":
    args = parse_args()
    check_transcriptomic_output(args)
