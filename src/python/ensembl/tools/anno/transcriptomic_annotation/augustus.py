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