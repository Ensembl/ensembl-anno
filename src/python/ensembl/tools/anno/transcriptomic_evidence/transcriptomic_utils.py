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


def run_trimming(
    main_output_dir, short_read_fastq_dir, delete_pre_trim_fastq, num_threads
):

    trim_galore_path = "trim_galore"
    utils.check_exe(trim_galore_path)

    trim_dir = utils.create_dir(main_output_dir, "trim_galore_output")

    fastq_file_list = []
    file_types = ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz")
    for file_type in file_types:
        fastq_file_list.extend(glob.glob(os.path.join(short_read_fastq_dir, file_type)))

    fastq_file_list = create_paired_paths(fastq_file_list)

    for fastq_file_path in fastq_file_list:
        logger.info("fastaq file path" + fastq_file_path)

    generic_trim_galore_cmd = [
        trim_galore_path,
        "--illumina",
        "--quality",
        "20",
        "--length",
        "50",
        "--output_dir",
        trim_dir,
    ]

    pool = multiprocessing.Pool(int(num_threads))
    for fastq_files in fastq_file_list:
        pool.apply_async(
            multiprocess_trim_galore,
            args=(
                generic_trim_galore_cmd,
                fastq_files,
                trim_dir,
            ),
        )
        if delete_pre_trim_fastq:
            for file_path in fastq_files:
                logger.info("Removing original fastq file post trimming:\n" + file_path)
                subprocess.run(["rm", file_path])

    pool.close()
    pool.join()

    trimmed_fastq_list = glob.glob(os.path.join(trim_dir, "*.fq.gz"))
    for trimmed_fastq_path in trimmed_fastq_list:
        logger.info("Trimmed file path:\n" + trimmed_fastq_path)
        sub_patterns = r"|".join(("_val_1.fq", "_val_2.fq", "_trimmed.fq"))
        updated_file_path = re.sub(sub_patterns, ".fq", trimmed_fastq_path)
        logger.info("Updated file path:\n" + updated_file_path)
        subprocess.run(["mv", trimmed_fastq_path, updated_file_path])

        files_to_delete_list = []
        for file_type in file_types:
            files_to_delete_list.extend(
                glob.glob(os.path.join(short_read_fastq_dir, file_type))
            )


def multiprocess_trim_galore(generic_trim_galore_cmd, fastq_files, trim_dir):

    fastq_file = fastq_files[0]
    fastq_file_pair = None
    if len(fastq_files) == 2:
        fastq_file_pair = fastq_files[1]

    trim_galore_cmd = generic_trim_galore_cmd
    if fastq_file_pair:
        trim_galore_cmd.append("--paired")

    trim_galore_cmd.append(fastq_file)

    if fastq_file_pair:
        trim_galore_cmd.append(fastq_file_pair)

    logger.info("Running Trim Galore with the following command:")
    logger.info(" ".join(trim_galore_cmd))
    subprocess.run(trim_galore_cmd)


def run_star_align(
    star_path,
    trim_fastq,
    subsample_script_path,
    main_output_dir,
    short_read_fastq_dir,
    genome_file,
    max_reads_per_sample,
    max_total_reads,
    max_short_read_intron_length,
    num_threads,
):
    # !!! Need to add in samtools path above instead of just using 'samtools' in command

    if not star_path:
        star_path = config["star"]["software"]

    utils.check_exe(star_path)

    # If trimming has been enabled then switch the path for
    # short_read_fastq_dir from the original location to the trimmed fastq dir
    if trim_fastq:
        short_read_fastq_dir = os.path.join(main_output_dir, "trim_galore_output")

    #  if not os.path.exists(subsample_script_path):
    subsample_script_path = "subsample_fastq.py"

    star_dir = utils.create_dir(main_output_dir, "star_output")

    logger.info("Skip analysis if the final log file already exists")
    log_final_out = pathlib.Path(os.path.join(star_dir, "Log.final.out"))
    log_out = pathlib.Path(os.path.join(star_dir, "Log.out"))
    log_progress_out = pathlib.Path(os.path.join(star_dir, "Log.progress.out"))
    if log_final_out.is_file() and log_out.is_file() and log_progress_out.is_file():
        logger.info("Star gtf file exists")
        return
    else:
        logger.info("No log files, go on with the analysis")

    star_tmp_dir = os.path.join(star_dir, "tmp")
    if os.path.exists(star_tmp_dir):
        subprocess.run(["rm", "-rf", star_tmp_dir])

    star_index_file = os.path.join(star_dir, "SAindex")

    fastq_file_list = []
    file_types = ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz")
    for file_type in file_types:
        fastq_file_list.extend(glob.glob(os.path.join(short_read_fastq_dir, file_type)))

    # This works out if the files are paired or not
    fastq_file_list = create_paired_paths(fastq_file_list)

    # Subsamples in parallel if there's a value set
    if max_reads_per_sample:
        pool = multiprocessing.Pool(int(num_threads))
        for fastq_files in fastq_file_list:
            fastq_file = fastq_files[0]
            fastq_file_pair = ""
            if len(fastq_files) == 2:
                fastq_file_pair = fastq_files[1]

            if (
                fastq_file_pair
                and os.path.exists(fastq_file + ".sub")
                and os.path.exists(fastq_file_pair + ".sub")
            ):
                logger.info(
                    "Found an existing .sub files on the fastq path for both members of the pair, \
                    will use those instead of subsampling again. Files:"
                )
                logger.info(fastq_file + ".sub")
                logger.info(fastq_file_pair + ".sub")
            elif fastq_file_pair:
                pool.apply_async(
                    run_subsample_script,
                    args=(
                        fastq_file,
                        fastq_file_pair,
                        subsample_script_path,
                    ),
                )
            elif os.path.exists(fastq_file + ".sub"):
                logger.info(
                    "Found an existing .sub file on the fastq path, will use that instead. File:"
                )
                logger.info(fastq_file + ".sub")
            else:
                pool.apply_async(
                    run_subsample_script,
                    args=(
                        fastq_file,
                        fastq_file_pair,
                        subsample_script_path,
                    ),
                )

        pool.close()
        pool.join()

    fastq_file_list = check_for_fastq_subsamples(fastq_file_list)

    if not fastq_file_list:
        raise IndexError(
            "The list of fastq files is empty. Fastq dir:\n%s" % short_read_fastq_dir
        )

    if not os.path.exists(star_index_file):
        logger.info("Did not find an index file for Star. Will create now")
        seq_region_lengths = utils.get_seq_region_lengths(genome_file, 0)
        genome_size = sum(seq_region_lengths.values())
        index_bases = min(14, math.floor((math.log(genome_size, 2) / 2) - 1))
        subprocess.run(
            [
                star_path,
                "--runThreadN",
                str(num_threads),
                "--runMode",
                "genomeGenerate",
                "--outFileNamePrefix",
                (star_dir + "/"),
                "--genomeDir",
                star_dir,
                "--genomeSAindexNbases",
                str(index_bases),
                "--genomeFastaFiles",
                genome_file,
            ]
        )

    if not star_index_file:
        raise IOError(
            "The index file does not exist. Expected path:\n%s" % star_index_file
        )

    logger.info("Running Star on the files in the fastq dir")
    for fastq_file_path in fastq_file_list:
        logger.info(fastq_file_path)
        fastq_file_name = os.path.basename(fastq_file_path)
        check_compression = re.search(r".gz$", fastq_file_name)

        # If there's a tmp dir already, the most likely cause is that \
        # STAR failed on the previous input file(s)
        # In this case STAR would effectively break for all the rest
        # of the files, as it won't run if the tmp
        # dir exists. So clean it up and put out a warning. Another approach
        # would be to just name the tmp dir
        # uniquely. Also there should be code checking the return on
        # STAR anyway. So later when this is being
        # cleaned up to have proper tests, need to decide on the best implementation
        if os.path.exists(star_tmp_dir):
            logger.error(
                "Found an existing tmp dir, implies potential failure on \
                previous file. Removing tmp dir"
            )
            try:
                shutil.rmtree(star_tmp_dir)
            except OSError as e:
                logger.error("Error: %s - %s." % (e.filename, e.strerror))

        sam_file_path = os.path.join(star_dir, (fastq_file_name + ".sam"))
        junctions_file_path = os.path.join(star_dir, (fastq_file_name + ".sj.tab"))
        sam_file_name = os.path.basename(sam_file_path)
        sam_temp_file_path = os.path.join(star_dir, (sam_file_name + ".tmp"))
        bam_sort_file_path = os.path.join(star_dir, re.sub(".sam", ".bam", sam_file_name))

        if (
            os.path.isfile(bam_sort_file_path)
            and not os.stat(bam_sort_file_path).st_size == 0
        ):
            logger.info(
                "Found an existing bam file for the fastq file, \
                presuming the file has been processed, will skip"
            )
            continue

        logger.info("Processing %s" % fastq_file_path)
        #    star_command = [star_path,'--outFilterIntronMotifs',
        #'RemoveNoncanonicalUnannotated','--outSAMstrandField','intronMotif',
        #'--runThreadN',str(num_threads),'--twopassMode','Basic','--runMode',
        #'alignReads','--genomeDir',star_dir,'--readFilesIn',fastq_file_path,
        #'--outFileNamePrefix',(star_dir + '/'),'--outTmpDir',star_tmp_dir,
        #'--outSAMtype','SAM','--alignIntronMax',str(max_intron_length),
        #'--outSJfilterIntronMaxVsReadN','5000','10000','25000','40000',
        #'50000','50000','50000','50000','50000','100000']

        star_command = [
            star_path,
            "--outFilterIntronMotifs",
            "RemoveNoncanonicalUnannotated",
            "--outSAMstrandField",
            "intronMotif",
            "--runThreadN",
            str(num_threads),
            "--twopassMode",
            "Basic",
            "--runMode",
            "alignReads",
            "--genomeDir",
            star_dir,
            "--readFilesIn",
            fastq_file_path,
            "--outFileNamePrefix",
            (star_dir + "/"),
            "--outTmpDir",
            star_tmp_dir,
            "--outSAMtype",
            "SAM",
            "--alignIntronMax",
            str(max_intron_length),
        ]

        if check_compression:
            star_command.append("--readFilesCommand")
            star_command.append("gunzip")
            star_command.append("-c")

        subprocess.run(star_command)
        subprocess.run(["mv", os.path.join(star_dir, "Aligned.out.sam"), sam_file_path])
        subprocess.run(["mv", os.path.join(star_dir, "SJ.out.tab"), junctions_file_path])

        logger.info("Converting samfile into sorted bam file. Bam file:")
        logger.info(bam_sort_file_path)
        subprocess.run(
            [
                "samtools",
                "sort",
                "-@",
                str(num_threads),
                "-T",
                sam_temp_file_path,
                "-o",
                bam_sort_file_path,
                sam_file_path,
            ]
        )

        logger.info("Removing sam file")
        subprocess.run(["rm", sam_file_path])

    logger.info("Completed running STAR")


def run_subsample_script(fastq_file, fastq_file_pair, subsample_script_path):

    if fastq_file_pair:
        subprocess.run(
            [
                "python3",
                subsample_script_path,
                "--fastq_file",
                fastq_file,
                "--fastq_file_pair",
                fastq_file_pair,
            ]
        )
    else:
        subprocess.run(["python3", subsample_script_path, "--fastq_file", fastq_file])


def check_for_fastq_subsamples(fastq_file_list):
    # This should probably removed at some point as it is needlessly complicated
    # Would be better to just build into the previous step
    # Mainly just about making sure that if you have subsamples
    # they're used and if you have pairs they're paired
    for idx, fastq_files in enumerate(fastq_file_list):
        fastq_file = fastq_files[0]
        subsample_file = fastq_file + ".sub"

        fastq_file_pair = ""
        subsample_file_pair = ""
        if len(fastq_files) == 2:
            fastq_file_pair = fastq_files[1]
            subsample_file_pair = fastq_file_pair + ".sub"

        # This bit will replace the list entry with a string,
        # don't need a list after this function for each pair/file
        if os.path.exists(subsample_file):
            logger.info(
                "Found a subsampled file extension, \
                will use that instead of the original file. Path:"
            )
            logger.info(subsample_file)
            fastq_file_list[idx] = subsample_file
        else:
            fastq_file_list[idx] = fastq_file

        # This bit just concats the paired file (or subsampled paired file) if it exists
        if os.path.exists(subsample_file_pair):
            logger.info(
                "Found a subsampled paired file extension,\
                will use that instead of the original file. Path:"
            )
            logger.info(subsample_file_pair)
            fastq_file_list[idx] = subsample_file + "," + subsample_file_pair
        elif fastq_file_pair:
            fastq_file_list[idx] = fastq_file + "," + fastq_file_pair

        logger.info("Entry at current index:")
        logger.info(fastq_file_list[idx])

    return fastq_file_list


def run_minimap2_align(
    minimap2_path,
    paftools_path,
    main_output_dir,
    long_read_fastq_dir,
    genome_file,
    max_intron_length,
    num_threads,
):

    if not minimap2_path:
        minimap2_path = config["minimap2"]["software"]

    utils.check_exe(minimap2_path)

    if not paftools_path:
        paftools_path = config["minimap2"]["paftools_path"]

    utils.check_exe(paftools_path)

    minimap2_output_dir = utils.create_dir(main_output_dir, "minimap2_output")

    logger.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(minimap2_output_dir, "annotation.gtf")
    if os.path.exists(output_file):
        transcript_count = utils.check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Minimap2 gtf file exists")
            return
    else:
        logger.info("No gtf file, go on with the analysis")

    genome_file_name = os.path.basename(genome_file)
    genome_file_index = genome_file_name + ".mmi"
    minimap2_index_file = os.path.join(minimap2_output_dir, genome_file_index)
    minimap2_hints_file = os.path.join(minimap2_output_dir, "minimap2_hints.gff")

    fastq_file_list = []
    for fastq_file in glob.glob(long_read_fastq_dir + "/*.fastq"):
        fastq_file_list.append(fastq_file)

    for fastq_file in glob.glob(long_read_fastq_dir + "/*.fq"):
        fastq_file_list.append(fastq_file)

    if not fastq_file_list:
        # NOTE: should update this to have a param that says
        # it's okay to be empty if there's an override param
        # This is because it's important that an external user
        # would have an exception if they provided a long read dir
        # but there was nothing in it. Whereas in the Ensembl pipeline
        # there might not be long read data, but we don't
        # know this when the command line is being constructed.
        # This is also true for the short read data. An alternative
        # would be to put in another analysis into the Ensembl pipeline
        # to construct this commandline after the transcriptomic
        # data has been searched for
        return
        # raise IndexError('The list of fastq files is empty.
        # Fastq dir:\n%s' % long_read_fastq_dir)

    if not os.path.exists(minimap2_index_file):
        logger.info("Did not find an index file for minimap2. Will create now")
        subprocess.run(
            [
                minimap2_path,
                "-t",
                str(num_threads),
                "-d",
                os.path.join(minimap2_index_file),
                genome_file,
            ]
        )

    if not minimap2_index_file:
        raise IOError(
            "The minimap2 index file does not exist. Expected path:\n%s"
            % minimap2_index_file
        )

    logger.info("Running minimap2 on the files in the long read fastq dir")
    for fastq_file_path in fastq_file_list:
        fastq_file_name = os.path.basename(fastq_file_path)
        sam_file = os.path.join(minimap2_output_dir, (fastq_file_name + ".sam"))
        bed_file = os.path.join(minimap2_output_dir, (fastq_file_name + ".bed"))
        bed_file_out = open(bed_file, "w+")
        logger.info("Processing %s" % fastq_file)
        subprocess.run(
            [
                minimap2_path,
                "-G",
                str(max_intron_length),
                "-t",
                str(num_threads),
                "--cs",
                "--secondary=no",
                "-ax",
                "splice",
                "-u",
                "b",
                minimap2_index_file,
                fastq_file_path,
                "-o",
                sam_file,
            ]
        )
        logger.info("Creating bed file from SAM")
        subprocess.run([paftools_path, "splice2bed", sam_file], stdout=bed_file_out)
        bed_file_out.close()

    bed_to_gtf(minimap2_output_dir)

    logger.info("Completed running minimap2")

def create_paired_paths(fastq_file_paths):
    path_dict = {}
    final_list = []

    for path in fastq_file_paths:
        match = re.search(r"(.+)_\d+\.(fastq|fq)", path)
        if not match:
            logger.error(
                "Could not find _1 or _2 at the end of the prefix \
                for file. Assuming file is not paired:"
            )
            logger.error(path)
            final_list.append([path])
            continue

        prefix = match.group(1)
        if prefix in path_dict:
            #      path_dict[prefix] = path_dict[prefix] + ',' + path
            path_dict[prefix].append(path)
        else:
            path_dict[prefix] = [path]

    for pair in path_dict:
        final_list.append(path_dict[pair])

    return final_list

def check_transcriptomic_output(main_output_dir):

    # This will check across the various transcriptomic
    # dirs and check there's actually some data
    transcriptomic_dirs = [
        "scallop_output",
        "stringtie_output",
        "minimap2_output",
    ]
    total_lines = 0
    min_lines = 100000
    for transcriptomic_dir in transcriptomic_dirs:
        full_file_path = os.path.join(
            main_output_dir, transcriptomic_dir, "annotation.gtf"
        )
        if not os.path.exists(full_file_path):
            logger.warning(
                "Warning, no annotation.gtf found for "
                + transcriptomic_dir
                + ". This might be fine, e.g. no long read data were provided"
            )
            continue
        num_lines = sum(1 for line in open(full_file_path))
        total_lines = total_lines + num_lines
        logger.info(
            "For "
            + transcriptomic_dir
            + " found a total of "
            + str(num_lines)
            + " in the annotation.gtf file"
        )
    if total_lines == 0:
        raise IOError(
            "Anno was run with transcriptomic mode enabled,\
            but the transcriptomic annotation files are empty"
        )
    elif total_lines <= min_lines:
        raise IOError(
            "Anno was run with transcriptomic mode enabled, \
            but the total number of lines in the output \
            files were less than the min expected value"
            + "\n"
            "Found: " + str(total_lines) + "\n"
            "Min allowed: " + str(min_lines)
        )

    else:
        logger.info(
            "Found "
            + str(total_lines)
            + " total lines across the transcriptomic files. Checks passed"
        )


# start gene g1 #pylint: disable=line-too-long
# 1       AUGUSTUS        gene    1       33908   1       +       .       g1
# 1       AUGUSTUS        transcript      1       33908   .       +       .       g1.t1
# 1       AUGUSTUS        CDS     3291    3585    .       +       2       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        exon    3291    3585    .       +       .       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        CDS     11377   11510   .       +       1       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        exon    11377   11510   .       +       .       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        CDS     30726   30871   .       +       2       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        exon    30726   30871   .       +       .       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        CDS     32975   33502   .       +       0       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        exon    32975   33908   .       +       .       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        stop_codon      33500   33502   .       +       0       transcript_id "g1.t1"; gene_id "g1";
# 1       AUGUSTUS        tts     33908   33908   .       +       .       transcript_id "g1.t1"; gene_id "g1";
# protein sequence = [GGRGEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEKKEKKEEEEEEEEEEEEEEEK
# EEGEGMEGRDIMGIRGKQKQPRKQPELSSSFPTFLLTQKTPYCPESIKLLLILRNNIQTIVFYKGPKGVINDWRKFKLESEDGDSIPPSKKEILRQMS
# SPQSRDDSKERMSRKMSIQEYELIHQDKEDESCLRKYRRQCMQDMHQKLSFGPRYGFVYELETGEQFLETIEKEQKVTTIVVNIYEDGVRGCDALNSS
# LACLAVEYPMVKFCKIKASNTGAGDRFSTDVLPTLLVYKGGELISNFISVAEQFAEEFFAVDVESFLNEYGLLPEREIHDLEQTNMEDEDIE]
# Evidence for and against this transcript:
# % of transcript supported by hints (any source): 0
# CDS exons: 0/4
# CDS introns: 0/4
# 5'UTR exons and introns: 0/0
# 3'UTR exons and introns: 0/1
# hint groups fully obeyed: 0
# incompatible hint groups: 0
# end gene g1
###


def augustus_output_to_gtf(augustus_output_dir, augustus_genome_dir):

    gtf_file_path = os.path.join(augustus_output_dir, "annotation.gtf")
    gtf_out = open(gtf_file_path, "w+")
    record_count = 1
    for gff_file_path in glob.glob(augustus_genome_dir + "/*.aug"):
        gff_file_name = os.path.basename(gff_file_path)
        match = re.search(r"\.rs(\d+)\.re(\d+)\.", gff_file_name)
        start_offset = int(match.group(1))

        exon_number = 1
        current_exon_hints_total = 0
        current_exon_hints_match = 0
        current_intron_hints_total = 0
        current_intron_hints_match = 0
        current_record = []
        gff_in = open(gff_file_path, "r")
        line = gff_in.readline()
        while line:
            match = re.search(r"# CDS exons\: (\d+)\/(\d+)", line)
            if match:
                current_exon_hints_match = match.group(1)
                current_exon_hints_total = match.group(2)

            match = re.search(r"# CDS introns\: (\d+)\/(\d+)", line)
            if match:
                current_introns_hints_match = match.group(1)
                current_introns_hints_total = match.group(2)

            if re.search(r"# end gene", line):
                for output_line in current_record:
                    gtf_out.write(output_line)
                current_record = []
                current_exon_hints_total = 0
                current_exon_hints_match = 0
                current_intron_hints_total = 0
                current_intron_hints_match = 0
                record_count += 1
                exon_number = 1

            values = line.split("\t")
            if (
                len(values) == 9
                and values[1] == "AUGUSTUS"
                and (values[2] == "transcript" or values[2] == "exon")
            ):
                values[3] = str(int(values[3]) + (start_offset - 1))
                values[4] = str(int(values[4]) + (start_offset - 1))
                values[8] = (
                    'gene_id "aug'
                    + str(record_count)
                    + '"; transcript_id "aug'
                    + str(record_count)
                    + '";'
                )
                if values[2] == "exon":
                    values[8] = values[8] + ' exon_number "' + str(exon_number) + '";'
                    exon_number += 1
                values[8] = values[8] + "\n"
                current_record.append("\t".join(values))

            line = gff_in.readline()
        gff_in.close()
    gtf_out.close()


def run_augustus_predict(augustus_path, main_output_dir, masked_genome_file, num_threads):

    min_seq_length = 1000

    if not augustus_path:
        augustus_path = config["augustus"]["software"]

    bam2hints_path = config["augustus"]["bam2hints_path"]
    bam2wig_path = config["augustus"]["bam2wig_path"]
    wig2hints_path = config["augustus"]["wig2hints_path"]
    utils.check_exe(augustus_path)

    # Run bam2hints, bam2wig, wig2hints, then combine the hints into a single file
    # Multiprocess with all three steps in the MP as that would be fastest

    augustus_dir = utils.create_dir(main_output_dir, "augustus_output")
    augustus_hints_dir = utils.create_dir(augustus_dir, "hints")
    augustus_genome_dir = utils.create_dir(augustus_dir, "genome_dir")
    augustus_evidence_dir = utils.create_dir(augustus_dir, "evidence")
    augustus_hints_file = os.path.join(augustus_evidence_dir, "augustus_hints.gff")
    star_dir = os.path.join(main_output_dir, "star_output")
    minimap2_output_dir = os.path.join(main_output_dir, "minimap2_output")

    if os.path.exists(star_dir):
        logger.info("Found a Star output dir, generating hints from any .sj.tab files")
        generate_hints(
            bam2hints_path,
            bam2wig_path,
            wig2hints_path,
            augustus_hints_dir,
            star_dir,
            num_threads,
        )
        hints_out = open(augustus_hints_file, "w+")
        for gff_file in glob.glob(augustus_hints_dir + "/*.bam.hints.gff"):
            gff_in = open(gff_file, "r")
            line = gff_in.readline()
            while line:
                hints_out.write(line)
                line = gff_in.readline()
            gff_in.close()
        hints_out.close()

    seq_region_lengths = utils.get_seq_region_lengths(genome_file, 5000)
    slice_ids = utils.create_slice_ids(seq_region_lengths, 1000000, 100000, 5000)

    generic_augustus_cmd = [
        augustus_path,
        "--species=human",
        "--UTR=on",
        (
            "--extrinsicCfgFile=" + "/hps/nobackup2/production/ensembl/jma/src/Augustus/"
            "config/extrinsic/extrinsic.M.RM.E.W.P.cfg"
        ),
    ]

    pool = multiprocessing.Pool(int(num_threads))
    tasks = []

    for slice_id in slice_ids:
        pool.apply_async(
            multiprocess_augustus_id,
            args=(
                generic_augustus_cmd,
                slice_id,
                masked_genome_file,
                augustus_hints_file,
                augustus_genome_dir,
            ),
        )

    pool.close()
    pool.join()

    augustus_output_to_gtf(augustus_dir, augustus_genome_dir)


def generate_hints(
    bam2hints_path,
    bam2wig_path,
    wig2hints_path,
    augustus_hints_dir,
    star_dir,
    num_threads,
):

    pool = multiprocessing.Pool(int(num_threads))
    for bam_file in glob.glob(star_dir + "/*.bam"):
        pool.apply_async(
            multiprocess_augustus_hints,
            args=(
                bam2hints_path,
                bam2wig_path,
                wig2hints_path,
                bam_file,
                augustus_hints_dir,
            ),
        )
    pool.close()
    pool.join()


def multiprocess_augustus_hints(
    bam2hints_path, bam2wig_path, wig2hints_path, bam_file, augustus_hints_dir
):
    bam_file_name = os.path.basename(bam_file)
    logger.info("Processing " + bam_file_name + " for Augustus hints")

    bam2hints_file_name = bam_file_name + ".hints.gff"
    bam2hints_file_path = os.path.join(augustus_hints_dir, bam2hints_file_name)
    bam2hints_cmd = [
        bam2hints_path,
        ("--in=" + bam_file),
        ("--out=" + bam2hints_file_path),
        "--maxintronlen=100000",
    ]
    logger.info("bam2hints command:\n" + " ".join(bam2hints_cmd))
    subprocess.run(bam2hints_cmd)


#  bam2wig_cmd = [bam2wig_path,'-D',augustus_hints_dir,bam_file]
#  print("bam2wig command:\n" + ' '.join(bam2wig_cmd))
#  subprocess.run(bam2wig_cmd)


# wig2hints is odd in that it runs directly off STDIN and then just prints to STDOUT,
# so the code below is implemented in steps as it's not good practice to use pipes and
# redirects in a subprocess command
#  wig_file_name = re.sub('.bam','.wig',bam_file_name)
#  wig_file_path = os.path.join(augustus_hints_dir,wig_file_name)
#  wig_hints_file_name = (wig_file_name + '.hints.gff')
#  wig_hints_file_path =  os.path.join(augustus_hints_dir,wig_hints_file_name)
#  print("Writing wig file info to hints file:\n" + wig_hints_file_name)
#  wig2hints_out = open(wig_hints_file_path,'w+')
#  wigcat = subprocess.Popen(('cat',wig_file_path), stdout=subprocess.PIPE)
#  subprocess.run(wig2hints_path, stdin=wigcat.stdout, stdout=wig2hints_out)
#  wig2hints_out.close()


def multiprocess_augustus_id(cmd, slice_id, genome_file, hints_file, output_dir):

    region = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]
    seq = utils.get_sequence(region, start, end, 1, genome_file, output_dir)

    region_fasta_file_name = region + ".rs" + str(start) + ".re" + str(end) + ".fa"
    region_fasta_file_path = os.path.join(output_dir, region_fasta_file_name)
    region_augustus_file_path = os.path.join(
        output_dir, (region_fasta_file_name + ".aug")
    )

    region_fasta_out = open(region_fasta_file_path, "w+")
    region_fasta_out.write(">" + region + "\n" + seq + "\n")
    region_fasta_out.close()

    region_hints_file = create_slice_hints_file(
        region, start, end, hints_file, region_fasta_file_path
    )

    aug_out = open(region_augustus_file_path, "w+")

    augustus_forward = cmd.copy()
    augustus_forward.append(("--hintsfile=" + region_hints_file))
    augustus_forward.append("--strand=forward")
    augustus_forward.append(region_fasta_file_path)
    subprocess.run(augustus_forward, stdout=aug_out)

    augustus_backward = cmd.copy()
    augustus_backward.append(("--hintsfile=" + region_hints_file))
    augustus_backward.append("--strand=backward")
    augustus_backward.append(region_fasta_file_path)
    subprocess.run(augustus_backward, stdout=aug_out)

    aug_out.close()


def create_slice_hints_file(region, start, end, hints_file, region_fasta_file_path):

    # Note this is trying to be memory and file efficient at the cost of speed
    # So files are only created as needed and the hints are
    # being read line by line as written as needed
    # This comes with the downside of being slow, but it's
    # only a very small amount of time relative
    # to how slow the step is in total. Given that this step
    # in general eats up a low of memory, saving as much
    # as possible here is not a bad thing even if it's adding
    # in an overhead by continuously reading the hints file

    region_hints_file_path = region_fasta_file_path + ".gff"
    hints_in = open(hints_file)
    hints_out = open(region_hints_file_path, "w+")
    hint_line = hints_in.readline()
    while hint_line:
        hint_line_values = hint_line.split("\t")
        if not len(hint_line_values) == 9:
            hint_line = hints_in.readline()
            continue

        hint_region = hint_line_values[0]
        hint_region_start = int(hint_line_values[3])
        hint_region_end = int(hint_line_values[4])

        if (
            hint_region == region
            and hint_region_start >= start
            and hint_region_end <= end
        ):
            hint_line_values[3] = str(int(hint_line_values[3]) - (start - 1))
            hint_line_values[4] = str(int(hint_line_values[4]) - (start - 1))
            hints_out.write("\t".join(hint_line_values))

        hint_line = hints_in.readline()
    hints_in.close()
    hints_out.close()

    return region_hints_file_path


def run_stringtie_assemble(
    stringtie_path, samtools_path, main_output_dir, genome_file, num_threads
):

    if not stringtie_path:
        stringtie_path = shutil.which("stringtie")
    utils.check_exe(stringtie_path)

    if not samtools_path:
        samtools_path = shutil.which("samtools")
    utils.check_exe(samtools_path)

    stringtie_dir = utils.create_dir(main_output_dir, "stringtie_output")

    logger.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(stringtie_dir, "annotation.gtf")
    if os.path.exists(output_file):
        transcript_count = utils.check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Stringtie gtf file exists")
            return
    else:
        logger.info("No gtf file, go on with the analysis")

    stringtie_merge_input_file = os.path.join(stringtie_dir, "stringtie_assemblies.txt")
    stringtie_merge_output_file = os.path.join(stringtie_dir, "annotation.gtf")
    star_dir = os.path.join(main_output_dir, "star_output")

    if os.path.exists(star_dir):
        logger.info("Found a Star output dir, will load sam file")

    sorted_bam_files = []
    for bam_file in glob.glob(star_dir + "/*.bam"):
        sorted_bam_files.append(bam_file)

    if not sorted_bam_files:
        raise IndexError(
            "The list of sorted bam files is empty, expected \
            them in Star output dir. Star dir:\n%s"
            % star_dir
        )

    # Don't know why this isn't multiprocessed, probably cos it
    # was fast enough in serial. But consider multiprocessing if
    # the mem usage is low
    for sorted_bam_file in sorted_bam_files:
        sorted_bam_file_name = os.path.basename(sorted_bam_file)
        transcript_file_name = re.sub(".bam", ".stringtie.gtf", sorted_bam_file_name)
        transcript_file_path = os.path.join(stringtie_dir, transcript_file_name)

        if os.path.exists(transcript_file_path):
            logger.info(
                "Found an existing stringtie gtf file, will not overwrite. File found:"
            )
            logger.info(transcript_file_path)
        else:
            logger.info("Running Stringtie on: " + sorted_bam_file_name)
            logger.info("Writing output to: " + transcript_file_path)
            subprocess.run(
                [
                    stringtie_path,
                    sorted_bam_file,
                    "-o",
                    transcript_file_path,
                    "-p",
                    str(num_threads),
                    "-t",
                    "-a",
                    "15",
                ]
            )

    # Now need to merge
    logger.info("Creating Stringtie merge input file: " + stringtie_merge_input_file)

    # Now need to merge
    logger.info("Creating Stringtie merge input file: " + stringtie_merge_input_file)
    gtf_list_out = open(stringtie_merge_input_file, "w+")
    for gtf_file in glob.glob(stringtie_dir + "/*.stringtie.gtf"):
        transcript_count = utils.check_gtf_content(gtf_file, "transcript")
        if transcript_count > 0:
            gtf_list_out.write(gtf_file + "\n")
        else:
            logger.warning(
                "Warning, skipping file with no transcripts. Path:\n" + gtf_file
            )
    gtf_list_out.close()

    logger.info("Merging Stringtie results. Writing to the following file:")
    logger.info(stringtie_merge_output_file)
    subprocess.run(
        [
            stringtie_path,
            "--merge",
            "-o",
            stringtie_merge_output_file,
            stringtie_merge_input_file,
        ]
    )


def run_scallop_assemble(scallop_path, stringtie_path, main_output_dir):

    if not scallop_path:
        scallop_path = shutil.which("scallop")
    utils.check_exe(scallop_path)

    if not stringtie_path:
        stringtie_path = shutil.which("stringtie")
    utils.check_exe(stringtie_path)

    scallop_dir = utils.create_dir(main_output_dir, "scallop_output")

    logger.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(scallop_dir, "annotation.gtf")
    if os.path.exists(output_file):
        transcript_count = utils.check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Scallop gtf file exists")
            return
    else:
        logger.info("No gtf file, go on with the analysis")

    stringtie_merge_input_file = os.path.join(scallop_dir, "scallop_assemblies.txt")
    stringtie_merge_output_file = os.path.join(scallop_dir, "annotation.gtf")
    star_dir = os.path.join(main_output_dir, "star_output")
    memory_limit = 40 * 1024**3

    if os.path.exists(star_dir):
        logger.info("Found a Star output dir, will load bam files")

    sorted_bam_files = []
    for bam_file in glob.glob(star_dir + "/*.bam"):
        sorted_bam_files.append(bam_file)

    if not sorted_bam_files:
        raise IndexError(
            "The list of sorted bam files is empty, expected \
            them in Star output dir. Star dir:\n%s"
            % star_dir
        )

    # Don't know why this isn't multiprocessed, probably cos it was
    # fast enough in serial. But consider multiprocessing if
    # the mem usage is low
    for sorted_bam_file in sorted_bam_files:
        sorted_bam_file_name = os.path.basename(sorted_bam_file)
        transcript_file_name = re.sub(".bam", ".scallop.gtf", sorted_bam_file_name)
        transcript_file_path = os.path.join(scallop_dir, transcript_file_name)

        if os.path.exists(transcript_file_path):
            logger.info(
                "Found an existing scallop gtf file, will not overwrite. File found:"
            )
            logger.info(transcript_file_path)
        else:
            logger.info("Running Scallop on: " + sorted_bam_file_name)
            logger.info("Writing output to: " + transcript_file_path)
            scallop_cmd = [
                scallop_path,
                "-i",
                sorted_bam_file,
                "-o",
                transcript_file_path,
                "--min_flank_length",
                "10",
            ]
            if memory_limit is not None:
                scallop_cmd = prlimit_command(scallop_cmd, memory_limit)

            return_value = None
            try:
                return_value = subprocess.check_output(scallop_cmd)
            except subprocess.CalledProcessError as ex:
                logger.error("Issue processing the following region with scallop")
                logger.error("Return value: " + str(return_value))

    # Now need to merge
    logger.info("Creating Stringtie merge input file: " + stringtie_merge_input_file)

    gtf_list_out = open(stringtie_merge_input_file, "w+")
    for gtf_file in glob.glob(scallop_dir + "/*.scallop.gtf"):
        transcript_count = utils.check_gtf_content(gtf_file, "transcript")
        if transcript_count > 0:
            gtf_list_out.write(gtf_file + "\n")
        else:
            logger.warning(
                "Warning, skipping file with no transcripts. Path:\n" + gtf_file
            )
    gtf_list_out.close()

    logger.info("Merging Scallop results. Writing to the following file:")
    logger.info(stringtie_merge_output_file)
    subprocess.run(
        [
            stringtie_path,
            "--merge",
            "-o",
            stringtie_merge_output_file,
            stringtie_merge_input_file,
        ]
    ) 

def model_builder(work_dir):

    star_output_dir = os.path.join(work_dir, "star_output")

    all_junctions_file = os.path.join(star_output_dir, "all_junctions.sj")
    sjf_out = open(all_junctions_file, "w+")

    for sj_tab_file in glob.glob(input_dir + "/*.sj.tab"):
        sjf_in = open(sj_tab_file)
        sjf_lines = sjf_in.readlines()
        for line in sjf_lines:
            elements = line.split("\t")
            strand = "+"

            #    my $slice_name = $eles[0];
            #    my $start = $eles[1];
            #    my $end = $eles[2];
            #    my $strand = $eles[3];

            # If the strand is undefined then skip, Augustus expects a strand
            if elements[3] == "0":
                continue
            elif elements[3] == "2":
                strand = "-"

            junction_length = int(elements[2]) - int(elements[1]) + 1
            if junction_length < 100:
                continue

            if not elements[4] and elements[7] < 10:
                continue

            # For the moment treat multimapping and single
            # mapping things as a combined score
            score = float(elements[6]) + float(elements[7])
            score = str(score)
            output_line = [
                elements[0],
                "RNASEQ",
                "intron",
                elements[1],
                elements[2],
                score,
                strand,
                ".",
                ("src=W;mul=" + score + ";"),
            ]
            sjf_out.write("\t".join(output_line) + "\n")

    sjf_out.close()
    
