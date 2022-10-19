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


# standard library
import errno
import logging
import os
import pathlib
import re
import shutil
import sys

from typing import Union


# logging formats
logging_formatter_time_message = logging.Formatter(
    fmt="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

# set up base logger
logger = logging.getLogger("main_logger")
logger.setLevel(logging.DEBUG)
logger.propagate = False
# create console handler and add to logger
console_handler = logging.StreamHandler(sys.stderr)
console_handler.setLevel(logging.DEBUG)
console_handler.setFormatter(logging_formatter_time_message)
logger.addHandler(console_handler)


def add_log_file_handler(
    logger: logging.Logger,
    log_file_path: Union[pathlib.Path, str],
    logging_formatter: logging.Formatter = logging_formatter_time_message,
):
    """
    Create file handler and add to logger.
    """
    file_handler = logging.FileHandler(log_file_path, mode="a+")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(logging_formatter)
    logger.addHandler(file_handler)


def check_exe(exe_path):
    if not shutil.which(exe_path):
        raise OSError('Executable file not found at "%s"' % exe_path)


def check_file(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file_path)


def check_gtf_content(gtf_file: Union[pathlib.Path, str], content_obj):
    """
    This just checks how many transcript lines are in a GTF
    """
    transcript_count = 0
    with open(gtf_file) as gtf_in:
        for line in gtf_in:
            eles = line.split("\t")
            if not len(eles) == 9:
                continue
            if eles[2] == content_obj:
                transcript_count += 1
    logger.info("%s GTF transcript count: %s" % (gtf_file, transcript_count))
    return transcript_count


def create_dir(main_output_dir: Union[pathlib.Path, str], dir_name: str = None):
    """
    Create directory or subdirectory and log operations.

    Args:
        main_output_dir: main output directory path
        dir_name: optional subdirectory to be created
    Returns:
        created directory Path object
    """
    main_output_dir = pathlib.Path(main_output_dir)

    if dir_name:
        target_dir = main_output_dir / dir_name
    else:
        target_dir = main_output_dir

    try:
        target_dir.mkdir()
    except FileExistsError:
        logger.warning('Directory "%s" already exists' % target_dir)
    except OSError:
        logger.error('Failed to create directory "%s"' % target_dir)
        sys.exit()
    else:
        logger.info('Successfully created directory "%s"' % target_dir)

    return target_dir


def create_paired_paths(fastq_file_paths):
    """
    TODO
    standardize to str or pathlib.Path paths
    """
    path_dict = {}
    final_list = []

    for path in fastq_file_paths:
        match = re.search(r"(.+)_\d+\.(fastq|fq)", str(path))
        if not match:
            logger.error(
                "Could not find _1 or _2 at the end of the prefix for file. Assuming file is not paired: %s"
                % path
            )
            final_list.append([path])
            continue

        prefix = match.group(1)
        if prefix in path_dict:
            # path_dict[prefix] = path_dict[prefix] + ',' + path
            path_dict[prefix].append(path)
        else:
            path_dict[prefix] = [path]

    for pair in path_dict:
        final_list.append(path_dict[pair])

    return final_list


def get_seq_region_lengths(genome_file: Union[pathlib.Path, str], min_seq_length: int):
    current_header = ""
    current_seq = ""

    seq_regions = {}
    with open(genome_file) as file_in:
        for line in file_in:
            match = re.search(r">(.+)$", line)
            if match and current_header:
                if len(current_seq) > min_seq_length:
                    seq_regions[current_header] = len(current_seq)

                current_seq = ""
                current_header = match.group(1)
            elif match:
                current_header = match.group(1)
            else:
                current_seq += line.rstrip()

        if len(current_seq) > min_seq_length:
            seq_regions[current_header] = len(current_seq)

    return seq_regions


def prlimit_command(command_list: list, virtual_memory_limit: int):
    """
    Uses the `prlimit` program to set a memory limit for a command list to be run with subprocess.

    prlimit - get and set process resource limits
    -v, --as[=limits]
        Address space limit.

    Args:
        command_list: original subprocess command list
        virtual_memory_limit: virtual memory limit in bytes
    Returns:
        memory limited subprocess command list
    """
    return ["prlimit", f"-v{virtual_memory_limit}"] + command_list


def list_to_string(original_list: List) -> str:
    """
    Create a string with the original list elements string representations concatenated
    with spaces between them.

    Args:
        original_list: original list
    Returns:
        generated list element string
    """
    return str.join(" ", [str(element) for element in original_list])
