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
import argparse
import gzip
import multiprocessing
import pathlib
import random
import re


def subsample(
    fastq_files: list,
    output_files: list,
    subsample_read_limit: int,
    num_processes: int,
    compressed: bool,
):
    fastq_file = fastq_files[0]
    fastq_file_pair = fastq_files[1]

    output_file = output_files[0]
    output_file_pair = output_files[1]

    check_compression = re.search(r"\.gz$", fastq_file)
    if check_compression:
        print("Found a .gz extension, so assuming compression")
        compressed = True

    # Count the file to begin with
    if compressed:
        num_lines = sum(1 for line in gzip.open(fastq_file))
    else:
        num_lines = sum(1 for line in open(fastq_file))

    range_limit = int(num_lines / 4)

    if range_limit <= subsample_read_limit:
        print(
            f"Number of reads ({range_limit}) is less than the max allowed read count ({subsample_read_limit}), no need to subsample"
        )
        return

    random_indices = {}

    rand_list = random.sample(range(0, range_limit - 1), subsample_read_limit)
    for idx, item in enumerate(rand_list):
        random_indices[rand_list[idx] * 4] = 1

    # Note that because of the checks in main this will only be true if there is a paired file that exists
    if num_processes == 2:
        print("Processing paired files in parallel")
        pool = multiprocessing.Pool(num_processes)
        pool.apply_async(
            print_subsample,
            args=(
                fastq_file,
                output_file,
                random_indices,
                compressed,
            ),
        )
        pool.apply_async(
            print_subsample,
            args=(
                fastq_file_pair,
                output_file_pair,
                random_indices,
                compressed,
            ),
        )
        pool.close()
        pool.join()

    else:
        print_subsample(fastq_file, output_file, random_indices, compressed)
        if fastq_file_pair:
            print_subsample(
                fastq_file_pair, output_file_pair, random_indices, compressed
            )


def print_subsample(fastq_file, output_file, random_indices, compressed: bool):
    if compressed:
        file_in = gzip.open(fastq_file, "rt")
    else:
        file_in = open(fastq_file)

    line_index = 0
    with open(output_file, "w+") as file_out:
        l1 = file_in.readline()
        l2 = file_in.readline()
        l3 = file_in.readline()
        l4 = file_in.readline()

        if line_index in random_indices:
            file_out.write(l1)
            file_out.write(l2)
            file_out.write(l3)
            file_out.write(l4)
        line_index += 4

        while l4:
            if line_index in random_indices:
                file_out.write(l1)
                file_out.write(l2)
                file_out.write(l3)
                file_out.write(l4)

            l1 = file_in.readline()
            l2 = file_in.readline()
            l3 = file_in.readline()
            l4 = file_in.readline()
            line_index += 4

    file_in.close()


def main():
    """
    main function
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq_file", type=str, required=True, help="fastq file path")
    parser.add_argument("--fastq_file_pair", type=str, help="paired file path")
    parser.add_argument(
        "--output_file",
        type=str,
        help='Designated output file. Defaults to appending ".sub" to input file',
    )
    parser.add_argument(
        "--output_file_pair",
        type=str,
        help='Designated output file for the paired file. Defaults to appending ".sub" to input file',
    )
    parser.add_argument(
        "--subsample_read_limit",
        type=int,
        default=100_000_000,
        help="Maximum number of reads to subsample. Defaults to 100,000,000",
    )
    parser.add_argument(
        "--parallelize",
        action="store_true",
        help="Process files in parallel if a paired file is provided",
    )
    parser.add_argument(
        "--compressed",
        action="store_true",
        help="Set if the files are compressed with gzip",
    )

    args = parser.parse_args()

    fastq_file = args.fastq_file
    fastq_file_pair = args.fastq_file_pair
    output_file = args.output_file
    output_file_pair = args.output_file_pair
    subsample_read_limit = args.subsample_read_limit
    parallelize = args.parallelize
    compressed = args.compressed

    if not pathlib.Path(fastq_file).exists():
        raise FileNotFoundError("Fastq file not found: %s" % fastq_file)

    if fastq_file_pair and not pathlib.Path(fastq_file_pair).exists():
        raise FileNotFoundError("Paired fastq file not found: %s" % fastq_file_path)

    if not output_file:
        output_file = f"{fastq_file}.sub"
        print("Using output file:\n%s" % output_file)

    if fastq_file_pair and not output_file_pair:
        output_file_pair = f"{fastq_file_pair}.sub"
        print("Using output file for the paired file:\n%s" % output_file_pair)

    if fastq_file_pair and parallelize:
        num_processes = 2
    else:
        num_processes = 1
    print("num_processes: %s" % num_processes)

    fastq_files = [fastq_file, fastq_file_pair]
    output_files = [output_file, output_file_pair]

    subsample(
        fastq_files, output_files, subsample_read_limit, num_processes, compressed
    )


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.info("Interrupted with CTRL-C, exiting...")
