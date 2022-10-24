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


def run_repeatmasker_regions(
    genome_file: Union[pathlib.Path, str],
    repeatmasker_path,
    library,
    species,
    main_output_dir,
    num_threads: int,
):
    if not repeatmasker_path:
        repeatmasker_path = "RepeatMasker"

    check_exe(repeatmasker_path)
    repeatmasker_output_dir = create_dir(main_output_dir, "repeatmasker_output")
    os.chdir(repeatmasker_output_dir)

    output_file = repeatmasker_output_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "repeat")
        if transcript_count > 0:
            logger.info("Repeatmasker gtf file already exists, skipping analysis")
            return

    logger.info("Creating list of genomic slices")
    seq_region_lengths = get_seq_region_lengths(genome_file, min_seq_length=5000)
    slice_ids = create_slice_ids(
        seq_region_lengths, slice_size=1_000_000, overlap=0, min_length=5000
    )

    if not library:
        if not species:
            species = "homo"

        generic_repeatmasker_cmd = [
            repeatmasker_path,
            "-nolow",
            "-species",
            species,
            "-engine",
            "crossmatch",
            "-dir",
            repeatmasker_output_dir,
        ]
    else:
        generic_repeatmasker_cmd = [
            repeatmasker_path,
            "-nolow",
            "-lib",
            library,
            "-engine",
            "crossmatch",
            "-dir",
            repeatmasker_output_dir,
        ]

    logger.info("Running RepeatMasker processes")
    pool = multiprocessing.Pool(num_threads)
    tasks = []
    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_repeatmasker,
            args=(
                generic_repeatmasker_cmd,
                slice_id,
                genome_file,
                repeatmasker_output_dir,
            ),
        )

    pool.close()
    pool.join()
    slice_output_to_gtf(
        repeatmasker_output_dir, ".rm.gtf", 1, "repeat_id", "repeatmask"
    )


def multiprocess_repeatmasker(
    generic_repeatmasker_cmd,
    slice_id,
    genome_file: Union[pathlib.Path, str],
    repeatmasker_output_dir: Union[pathlib.Path, str],
):
    region_name = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]

    logger.info(
        "Processing slice to find repeats with RepeatMasker: %s:%s:%s"
        % (region_name, start, end)
    )
    seq = get_sequence(
        seq_region=region_name,
        start=start,
        end=end,
        strand=1,
        fasta_file=genome_file,
        output_dir=repeatmasker_output_dir,
    )

    slice_file_name = f"{region_name}.rs{start}.re{end}"
    region_fasta_file_name = f"{slice_file_name}.fa"
    region_fasta_file_path = repeatmasker_output_dir / region_fasta_file_name
    with open(region_fasta_file_path, "w+") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")

    region_results_file_path = f"{region_fasta_file_path}.rm.gtf"
    repeatmasker_output_file_path = f"{region_fasta_file_path}.out"
    repeatmasker_masked_file_path = f"{region_fasta_file_path}.masked"
    repeatmasker_tbl_file_path = f"{region_fasta_file_path}.tbl"
    repeatmasker_log_file_path = f"{region_fasta_file_path}.log"
    repeatmasker_cat_file_path = f"{region_fasta_file_path}.cat"

    repeatmasker_cmd = generic_repeatmasker_cmd.copy()
    repeatmasker_cmd.append(region_fasta_file_path)
    logger.info("repeatmasker_cmd: %s" % " ".join(repeatmasker_cmd))
    subprocess.run(repeatmasker_cmd)

    create_repeatmasker_gtf(
        repeatmasker_output_file_path, region_results_file_path, region_name
    )

    os.remove(region_fasta_file_path)
    if os.path.exists(region_results_file_path):
        os.remove(region_results_file_path)
    if os.path.exists(repeatmasker_masked_file_path):
        os.remove(repeatmasker_masked_file_path)
    if os.path.exists(repeatmasker_tbl_file_path):
        os.remove(repeatmasker_tbl_file_path)
    if os.path.exists(repeatmasker_log_file_path):
        os.remove(repeatmasker_log_file_path)
    if os.path.exists(repeatmasker_cat_file_path):
        os.remove(repeatmasker_cat_file_path)


def create_repeatmasker_gtf(
    repeatmasker_output_file_path, region_results_file_path, region_name
):
    with open(repeatmasker_output_file_path, "r") as repeatmasker_in, open(
        region_results_file_path, "w+"
    ) as repeatmasker_out:
        repeat_count = 1
        for line in repeatmasker_in:
            result_match = re.search(r"^\s*\d+\s+", line)
            if result_match:
                results = line.split()
                if results[-1] == "*":
                    results.pop()
                if not len(results) == 15:
                    continue

                score = results[0]
                start = results[5]
                end = results[6]
                strand = results[8]
                repeat_name = results[9]
                repeat_class = results[10]
                if strand == "+":
                    repeat_start = results[11]
                    repeat_end = results[12]
                else:
                    repeat_start = results[13]
                    repeat_end = results[12]
                    strand = "-"

                gtf_line = f'{region_name}\tRepeatMasker\trepeat\t{start}\t{end}\t.\t{strand}\t.\trepeat_id "{repeat_count}"; repeat_name "{repeat_name}"; repeat_class "{repeat_class}"; repeat_start "{repeat_start}"; repeat_end "{repeat_end}"; score "{score}";\n'

                repeatmasker_out.write(gtf_line)
                repeat_count += 1


