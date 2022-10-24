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


