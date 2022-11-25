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
import gc
import glob
import io
import math
import multiprocessing
import os
import pathlib
import random
import re
import shutil
import signal
import subprocess
import tempfile
import sys
import errno
import logging

from typing import List, Union, Tuple

from utils import (
    add_log_file_handler,
    check_exe,
    check_file,
    check_gtf_content,
    create_dir,
    get_seq_region_lengths,
    logger,
    prlimit_command,
    load_results_to_ensembl_db,
    generic_load_records_to_ensembl_db,
    multiprocess_load_records_to_ensembl_db,
    batch_gtf_records,
    run_find_orfs,
    find_orf_phased_region,
    slice_output_to_gtf,
    convert_gff_to_gtf,
    set_attributes,
    create_slice_ids,
    update_gtf_genes,
    read_gtf_genes,
    fasta_to_dict,
    splice_junction_to_gff,
    split_genome,
    multiprocess_generic,
    reverse_complement,
    get_seq_region_names,
    slice_genome,
    subprocess_run_and_log,
    get_sequence,
    seq_region_names,
    list_to_string,
#    coallate_results,
)
def run_subsample_script(subsample_script_path, fastq_file, fastq_file_pair):
    subsample_script_cmd = [
        "python3",
        subsample_script_path,
        "--fastq_file",
        fastq_file,
    ]
    if fastq_file_pair:
        subsample_script_cmd.extend(["--fastq_file_pair", fastq_file_pair])
    subprocess.run(subsample_script_cmd)



def check_for_fastq_subsamples(
    fastq_paired_paths: List[List[pathlib.Path]],
) -> List[List[pathlib.Path]]:
    """
    Replace FASTQ files with their corresponding subsampled versions if they exist.
    Original comment:
    This should probably removed at some point as it is needlessly complicated
    Would be better to just build into the previous step
    Mainly just about making sure that if you have subsamples they're used and if you have pairs they're paired
    Args:
        fastq_paired_paths: list of lists containing the FASTQ paired file paths
    Returns:
        the input list with FASTQ file paths replaced with their corresponding
        subsampled versions
    """
    for index, fastq_files in enumerate(fastq_paired_paths):
        fastq_file = fastq_files[0]
        fastq_file_subsampled = fastq_file.parent / f"{fastq_file.name}.sub"

        if len(fastq_files) == 2:
            fastq_file_pair = fastq_files[1]
            fastq_file_pair_subsampled = (
                fastq_file_pair.parent / f"{fastq_file_pair.name}.sub"
            )
        # len(fastq_files) == 1
        else:
            fastq_file_pair = None
            fastq_file_pair_subsampled = None

        if fastq_file_subsampled.exists():
            logger.info(
                "Found a subsampled file extension, will use that instead of the original file:\n%s"
                % fastq_file_subsampled
            )
            fastq_paired_paths[index] = [fastq_file_subsampled]
        else:
            fastq_paired_paths[index] = [fastq_file]

        if fastq_file_pair_subsampled.exists():
            logger.info(
                "Found a subsampled paired file extension, will use that instead of the original file:\n%s"
                % fastq_file_pair_subsampled
            )
            fastq_paired_paths[index].append(fastq_file_pair_subsampled)
        elif fastq_file_pair:
            fastq_paired_paths[index].append(fastq_file_pair)

        logger.info("Entry at current index:\n%s" % fastq_paired_paths[index])

    return fastq_paired_paths

def run_minimap2_align(
    minimap2_path,
    paftools_path,
    main_output_dir,
    long_read_fastq_dir,
    genome_file: pathlib.Path,
    max_intron_length,
    num_threads: int,
):
    if not minimap2_path:
        minimap2_path = "minimap2"

    check_exe(minimap2_path)

    if not paftools_path:
        paftools_path = "paftools.js"

    check_exe(paftools_path)

    minimap2_output_dir = create_dir(main_output_dir, "minimap2_output")

    output_file = minimap2_output_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Minimap2 GTF file exists skipping analysis")
            return

    genome_file_name = genome_file.name
    genome_file_index = f"{genome_file_name}.mmi"
    minimap2_index_file = minimap2_output_dir / genome_file_index
    minimap2_hints_file = minimap2_output_dir / "minimap2_hints.gff"

    fastq_file_list = []
    for fastq_file in glob.glob(f"{long_read_fastq_dir}/*.fastq"):
        fastq_file_list.append(fastq_file)

    for fastq_file in glob.glob(f"{long_read_fastq_dir}/*.fq"):
        fastq_file_list.append(fastq_file)

    if not fastq_file_list:
        # NOTE: should update this to have a param that says it's okay to be empty if there's an override param
        # This is because it's important that an external user would have an exception if they provided a long read dir
        # but there was nothing in it. Whereas in the Ensembl pipeline there might not be long read data, but we don't
        # know this when the command line is being constructed. This is also true for the short read data. An alternative
        # would be to put in another analysis into the Ensembl pipeline to construct this commandline after the transcriptomic
        # data has been searched for
        return
        # raise IndexError('The list of fastq files is empty. Fastq dir:\n%s' % long_read_fastq_dir)

    if not minimap2_index_file.exists():
        logger.info("Did not find an index file for minimap2. Will create now")
        minimap2_command = [
             minimap2_path,
             "-t",
             str(num_threads),
             "-d",
             minimap2_index_file,
             genome_file,
         ]
        subprocess.run(minimap2_command)

    if not minimap2_index_file.exists():
        raise FileNotFoundError(
            "Failed to create the minimap2 index file at:\n%s" % minimap2_index_file
        )

    logger.info("Running minimap2 on the files in the long read fastq dir")
    for fastq_file_path in fastq_file_list:
        fastq_file_name = fastq_file_path.name
        sam_file = minimap2_output_dir / f"{fastq_file_name}.sam"
        bed_file = minimap2_output_dir / f"{fastq_file_name}.bed"
        bed_file_out = open(bed_file, "w+")
        logger.info("Processing %s" % fastq_file)
        minimap2_command = [
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
        subprocess.run(minimap2_command)
        logger.info("Creating bed file from SAM")
        paftools_command = [paftools_path, "splice2bed", sam_file]
        subprocess.run(paftools_command, stdout=bed_file_out)
        bed_file_out.close()

    bed_to_gtf(minimap2_output_dir)

    logger.info("Completed running minimap2")


def bed_to_gtf(minimap2_output_dir: pathlib.Path):
    gtf_file_path = minimap2_output_dir / "annotation.gtf"
    with open(gtf_file_path, "w+") as gtf_out:
        gene_id = 1
        for bed_file in minimap2_output_dir.glob("*.bed"):
            logger.info("Converting bed to GTF:\n%s" % bed_file)
            with open(bed_file) as bed_in:
                for line in bed_in:
                    line = line.rstrip()
                    elements = line.split("\t")
                    seq_region_name = elements[0]
                    offset = int(elements[1])
                    hit_name = elements[3]
                    strand = elements[5]
                    block_sizes = elements[10].split(",")
                    block_sizes = list(filter(None, block_sizes))
                    block_starts = elements[11].split(",")
                    block_starts = list(filter(None, block_starts))
                    exons = bed_to_exons(block_sizes, block_starts, offset)
                    transcript_line = [
                        seq_region_name,
                        "minimap",
                        "transcript",
                        0,
                        0,
                        ".",
                        strand,
                        ".",
                        f'gene_id "minimap_{gene_id}"; transcript_id "minimap_{gene_id}"',
                    ]
                    transcript_start = None
                    transcript_end = None
                    exon_records = []
                    for index, exon_coords in enumerate(exons, start=1):
                        if (
                            transcript_start is None
                            or exon_coords[0] < transcript_start
                        ):
                            transcript_start = exon_coords[0]

                        if transcript_end is None or exon_coords[1] > transcript_end:
                            transcript_end = exon_coords[1]

                        exon_line = [
                            seq_region_name,
                            "minimap",
                            "exon",
                            str(exon_coords[0]),
                            str(exon_coords[1]),
                            ".",
                            strand,
                            ".",
                            f'gene_id "minimap_{gene_id}"; transcript_id "minimap_{gene_id}"; exon_number "{index}";',
                        ]

                        exon_records.append(exon_line)

                    transcript_line[3] = str(transcript_start)
                    transcript_line[4] = str(transcript_end)

                    gtf_out.write("\t".join(transcript_line) + "\n")
                    for exon_line in exon_records:
                        gtf_out.write("\t".join(exon_line) + "\n")

                    gene_id += 1

def bed_to_gff(input_dir, hints_file):

    gff_out = open(hints_file, "w+")
    exons_dict = {}
    for bed_file in glob.glob(input_dir + "/*.bed"):
        logger.info("Processing file for hints:")
        logger.info(bed_file)
        bed_in = open(bed_file)
        bed_lines = bed_in.readlines()
        for line in bed_lines:
            line = line.rstrip()
            elements = line.split("\t")
            seq_region_name = elements[0]
            offset = int(elements[1])
            hit_name = elements[3]
            strand = elements[5]
            block_sizes = elements[10].split(",")
            block_sizes = list(filter(None, block_sizes))
            block_starts = elements[11].split(",")
            block_starts = list(filter(None, block_starts))
            exons = bed_to_exons(block_sizes, block_starts, offset)
            for i, element in enumerate(exons):
                exon_coords = exons[i]
                exon_key = (
                    seq_region_name
                    + ":"
                    + exon_coords[0]
                    + ":"
                    + exon_coords[1]
                    + ":"
                    + strand
                )
                if exon_key in exons_dict:
                    exons_dict[exon_key][5] += 1
                else:
                    gff_list = [
                        seq_region_name,
                        "CDNA",
                        "exon",
                        exon_coords[0],
                        exon_coords[1],
                        1.0,
                        strand,
                        ".",
                    ]
                    exons_dict[exon_key] = gff_list

    for exon_key, gff_list in exons_dict.items():
        gff_list[5] = str(gff_list[5])
        gff_line = "\t".join(gff_list) + "\tsrc=W;mul=" + gff_list[5] + ";\n"
        gff_out.write(gff_line)

    gff_out.close()

    sorted_hints_out = open((hints_file + ".srt"), "w+")
    subprocess.run(
        ["sort", "-k1,1", "-k7,7", "-k4,4", "-k5,5", hints_file],
        stdout=sorted_hints_out,
    )
    sorted_hints_out.close()

def bed_to_exons(block_sizes, block_starts, offset):
    exons = []
    for i, element in enumerate(block_sizes):
        block_start = offset + int(block_starts[i]) + 1
        block_end = block_start + int(block_sizes[i]) - 1

        if block_end < block_start:
            logger.warning("Warning: block end is less than block start, skipping exon")
            continue

        exon_coords = [str(block_start), str(block_end)]
        exons.append(exon_coords)

    return exons

def check_transcriptomic_output(main_output_dir: pathlib.Path):
    # This will check across the various transcriptomic dirs and check there's actually some data
    transcriptomic_dirs = ["scallop_output", "stringtie_output", "minimap2_output"]
    total_lines = 0
    min_lines = 100_000
    for transcriptomic_dir in transcriptomic_dirs:
        full_file_path = main_output_dir / transcriptomic_dir / "annotation.gtf"
        if not full_file_path.exists():
            logger.warning(
                'Warning, no annotation.gtf found for "%s". This might be fine, e.g. no long read data were provided'
                % transcriptomic_dir
            )
            continue
        num_lines = sum(1 for line in open(full_file_path))
        total_lines += num_lines
        logger.info(
            'For "%s" found a total of %s in the annotation.gtf file'
            % (transcriptomic_dir, num_lines)
        )
    if total_lines == 0:
        raise IOError(
            "Anno was run with transcriptomic mode enabled, but the transcriptomic annotation files are empty"
        )
    elif total_lines <= min_lines:
        raise IOError(
            "Anno was run with transcriptomic mode enabled, but the total number of lines in the output files were less than the min expected value\nFound: %s\nMin allowed: %s"
            % (total_lines, min_lines)
        )
    else:
        logger.info(
            "Found %s total lines across the transcriptomic files. Checks passed"
            % total_lines
        )


# start gene g1
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


def augustus_output_to_gtf(
    augustus_output_dir: pathlib.Path,
    augustus_genome_dir: pathlib.Path,
):
    gtf_file_path = augustus_output_dir / "annotation.gtf"
    with open(gtf_file_path, "w+") as gtf_out:
        record_count = 1
        for gff_file_path in augustus_genome_dir.glob("*.aug"):
            gff_file_name = gff_file_path.name
            match = re.search(r"\.rs(\d+)\.re(\d+)\.", gff_file_name)
            start_offset = int(match.group(1))

            exon_number = 1
            current_exon_hints_total = 0
            current_exon_hints_match = 0
            current_intron_hints_total = 0
            current_intron_hints_match = 0
            current_record = []
            with open(gff_file_path, "r") as gff_in:
                for line in gff_in:
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
                        # fmt: off
                        values[8] = f'gene_id "aug{record_count}"; transcript_id "aug{record_count}";'
                        # fmt: on
                        if values[2] == "exon":
                            values[8] += f' exon_number "{exon_number}";'
                            exon_number += 1
                        values[8] += "\n"
                        current_record.append("\t".join(values))


def run_augustus_predict(
    augustus_path,
    main_output_dir: pathlib.Path,
    masked_genome_file: pathlib.Path,
    num_threads: int,
):
    min_seq_length = 1000

    if not augustus_path:
        augustus_path = (
            "/hps/nobackup2/production/ensembl/jma/src/Augustus/bin/augustus"
        )
    check_exe(augustus_path)

    bam2hints_path = "/homes/fergal/bin/bam2hints"
    bam2wig_path = "/homes/fergal/bin/bam2wig"
    wig2hints_path = "/homes/fergal/bin/wig2hints"

    # Run bam2hints, bam2wig, wig2hints, then combine the hints into a single file
    # Multiprocess with all three steps in the MP as that would be fastest

    augustus_dir = create_dir(main_output_dir, "augustus_output")
    augustus_hints_dir = create_dir(augustus_dir, "hints")
    augustus_genome_dir = create_dir(augustus_dir, "genome_dir")
    augustus_evidence_dir = create_dir(augustus_dir, "evidence")
    augustus_hints_file = augustus_evidence_dir / "augustus_hints.gff"
    star_dir = main_output_dir / "star_output"
    # minimap2_output_dir = main_output_dir / "minimap2_output"

    if star_dir.exists():
        logger.info("Found a STAR output dir, generating hints from any .sj.tab files")
        generate_hints(
            bam2hints_path,
            bam2wig_path,
            wig2hints_path,
            augustus_hints_dir,
            star_dir,
            num_threads,
        )
        with open(augustus_hints_file, "w+") as hints_out:
            for gff_file in augustus_hints_dir.glob("*.bam.hints.gff"):
                with open(gff_file, "r") as gff_in:
                    for line in gff_in:
                        hints_out.write(line)

    seq_region_lengths = get_seq_region_lengths(genome_file, min_seq_length=5000)
    slice_ids = create_slice_ids(
        seq_region_lengths, slice_size=1_000_000, overlap=100_000, min_length=5000
    )

    generic_augustus_cmd = [
        augustus_path,
        "--species=human",
        "--UTR=on",
        (
            "--extrinsicCfgFile="
            + "/hps/nobackup2/production/ensembl/jma/src/Augustus/config/extrinsic/extrinsic.M.RM.E.W.P.cfg"
        ),
    ]

    pool = multiprocessing.Pool(num_threads)
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
    augustus_hints_dir: Union[pathlib.Path, str],
    star_dir: pathlib.Path,
    num_threads: int,
):
    pool = multiprocessing.Pool(num_threads)
    for bam_file in star_dir.glob("*.bam"):
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
    bam2hints_path,
    bam2wig_path,
    wig2hints_path,
    bam_file: pathlib.Path,
    augustus_hints_dir: pathlib.Path,
):
    bam_file_name = bam_file.name
    logger.info("Processing %s for Augustus hints" % bam_file_name)

    bam2hints_file_name = f"{bam_file_name}.hints.gff"
    bam2hints_file_path = augustus_hints_dir / bam2hints_file_name
    bam2hints_cmd = [
        bam2hints_path,
        f"--in={bam_file}",
        f"--out={bam2hints_file_path}",
        "--maxintronlen=100000",
    ]
    logger.info("bam2hints command:%s" % list_to_string(bam2hints_cmd))
    subprocess.run(bam2hints_cmd)

    # bam2wig_cmd = [bam2wig_path,'-D',augustus_hints_dir,bam_file]
    # print("bam2wig command:\n" + ' '.join(bam2wig_cmd))
    # subprocess.run(bam2wig_cmd)

    # wig2hints is odd in that it runs directly off STDIN and then just prints to STDOUT,
    # so the code below is implemented in steps as it's not good practice to use pipes and
    # redirects in a subprocess command
    # wig_file_name = re.sub('.bam','.wig',bam_file_name)
    # wig_file_path = os.path.join(augustus_hints_dir,wig_file_name)
    # wig_hints_file_name = (wig_file_name + '.hints.gff')
    # wig_hints_file_path =  os.path.join(augustus_hints_dir,wig_hints_file_name)
    # print("Writing wig file info to hints file:\n" + wig_hints_file_name)
    # wig2hints_out = open(wig_hints_file_path,'w+')
    # wigcat = subprocess.Popen(('cat',wig_file_path), stdout=subprocess.PIPE)
    # subprocess.run(wig2hints_path, stdin=wigcat.stdout, stdout=wig2hints_out)
    # wig2hints_out.close()


def multiprocess_augustus_id(
    generic_augustus_cmd,
    slice_id: Tuple[str, int, int],
    genome_file: pathlib.Path,
    hints_file,
    output_dir: pathlib.Path,
):
    region = slice_id[0]
    start = slice_id[1]
    end = slice_id[2]
    seq = get_sequence(
        seq_region=region,
        start=start,
        end=end,
        strand=1,
        fasta_file=genome_file,
        output_dir=output_dir,
    )

    region_fasta_file_name = f"{region}.rs{start}.re{end}.fa"
    region_fasta_file_path = output_dir / region_fasta_file_name
    region_augustus_file_path = output_dir / f"{region_fasta_file_name}.aug"

    with open(region_fasta_file_path, "w+") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")

    region_hints_file = create_slice_hints_file(
        region, start, end, hints_file, region_fasta_file_path
    )

    aug_out = open(region_augustus_file_path, "w+")

    augustus_forward_command = generic_augustus_cmd.copy()
    augustus_forward_command.extend(
        [f"--hintsfile={region_hints_file}", "--strand=forward", region_fasta_file_path]
    )
    subprocess.run(augustus_forward, stdout=aug_out)

    augustus_backward_command = generic_augustus_cmd.copy()
    augustus_backward_command.extend(
        [
            f"--hintsfile={region_hints_file}",
            "--strand=backward",
            region_fasta_file_path,
        ]
    )
    subprocess.run(augustus_backward_command, stdout=aug_out)

    aug_out.close()


def create_slice_hints_file(region, start, end, hints_file, region_fasta_file_path):
    """
    Note this is trying to be memory and file efficient at the cost of speed
    So files are only created as needed and the hints are being read line by line as written as needed
    This comes with the downside of being slow, but it's only a very small amount of time relative
    to how slow the step is in total. Given that this step in general eats up a low of memory, saving as much
    as possible here is not a bad thing even if it's adding in an overhead by continuously reading the hints file
    """
    region_hints_file_path = f"{region_fasta_file_path}.gff"
    with open(hints_file) as hints_in, open(region_hints_file_path, "w+") as hints_out:
        for hint_line in hints_in:
            hint_line_values = hint_line.split("\t")
            if not len(hint_line_values) == 9:
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

    return region_hints_file_path


def run_stringtie_assemble(
    stringtie_path, samtools_path, main_output_dir, num_threads: int
):
    if not stringtie_path:
        stringtie_path = shutil.which("stringtie")
    check_exe(stringtie_path)

    if not samtools_path:
        samtools_path = shutil.which("samtools")
    check_exe(samtools_path)

    stringtie_dir = create_dir(main_output_dir, "stringtie_output")

    output_file = stringtie_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("StringTie GTF file exists skipping analysis")
            return

    stringtie_merge_input_file = stringtie_dir / "stringtie_assemblies.txt"
    stringtie_merge_output_file = stringtie_dir / "annotation.gtf"
    star_dir = main_output_dir / "star_output"

    if star_dir.exists():
        logger.info("Found a STAR output dir, will load sam file")

    sorted_bam_files = list(star_dir.glob("*.bam"))

    if not sorted_bam_files:
        raise IndexError(
            "The list of sorted bam files is empty, expected them in STAR output dir:\n%s"
            % star_dir
        )

    # Don't know why this isn't multiprocessed, probably cos it was fast enough in serial. But consider multiprocessing if
    # the mem usage is low
    for sorted_bam_file in sorted_bam_files:
        sorted_bam_file_name = sorted_bam_file.name
        transcript_file_name = re.sub(".bam", ".stringtie.gtf", sorted_bam_file_name)
        transcript_file_path = stringtie_dir / transcript_file_name

        if transcript_file_path.exists():
            logger.info(
                "Found an existing stringtie gtf file, will not overwrite:\n%s"
                % transcript_file_path
            )
        else:
            logger.info(
                "Running Stringtie on: %s, writing output to:\n%s"
                % (sorted_bam_file_name, transcript_file_path)
            )
            stringtie_command = [
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
            subprocess.run(stringtie_command)

    # Now need to merge
    logger.info("Creating Stringtie merge input file: %s" % stringtie_merge_input_file)
    with open(stringtie_merge_input_file, "w+") as gtf_list_out:
        for gtf_file in stringtie_dir.glob("*.stringtie.gtf"):#glob.glob(str(stringtie_dir / "*.stringtie.gtf")):
            transcript_count = check_gtf_content(gtf_file, "transcript")
            if transcript_count > 0:
                gtf_list_out.write(f"{gtf_file}\n")
            else:
                logger.warning(
                    "Warning, skipping file with no transcripts:\n%s" % gtf_file
                )

    logger.info("Merging Stringtie results to:\n%s" % stringtie_merge_output_file)
    stringtie_command = [
         stringtie_path,
         "--merge",
         "-o",
         stringtie_merge_output_file,
         stringtie_merge_input_file,
     ]
    subprocess.run(stringtie_command)


def run_scallop_assemble(scallop_path, stringtie_path, main_output_dir: pathlib.Path):
    if not scallop_path:
        scallop_path = shutil.which("scallop")
    check_exe(scallop_path)

    if not stringtie_path:
        stringtie_path = shutil.which("stringtie")
    check_exe(stringtie_path)

    scallop_dir = create_dir(main_output_dir, "scallop_output")

    output_file = scallop_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Scallop gtf file exists, skipping analysis")
            return

    stringtie_merge_input_file = scallop_dir / "scallop_assemblies.txt"
    stringtie_merge_output_file = scallop_dir / "annotation.gtf"
    star_dir = main_output_dir / "star_output"
    memory_limit = 40 * 1024**3

    if star_dir.exists():
        logger.info("Found a STAR output dir, will load bam files")

    sorted_bam_files = list(star_dir.glob("*.bam"))

    if not sorted_bam_files:
        raise IndexError(
            "Empty list of sorted bam files, expected them in STAR output dir:\n%s"
            % star_dir
        )

    # Don't know why this isn't multiprocessed, probably cos it was fast enough in serial.
    # But consider multiprocessing if the mem usage is low.
    for sorted_bam_file in sorted_bam_files:
        sorted_bam_file_name = sorted_bam_file.name
        transcript_file_name = re.sub(".bam", ".scallop.gtf", sorted_bam_file_name)
        transcript_file_path = scallop_dir / transcript_file_name

        if transcript_file_path.exists():
            logger.info(
                "Found an existing scallop gtf file, will not overwrite:\n%s"
                % transcript_file_path
            )
        else:
            logger.info(
                "Running Scallop on: %s, writing output to:\n%s"
                % (sorted_bam_file_name, transcript_file_path)
            )
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
                logger.error(
                    "Issue processing the following region with scallop\nReturn value: %s"
                    % return_value
                )

    #      subprocess.run([scallop_path,'-i',sorted_bam_file,'-o',transcript_file_path,'--min_flank_length','10'])

    # Now need to merge
    logger.info("Creating Stringtie merge input file: %s" % stringtie_merge_input_file)

    with open(stringtie_merge_input_file, "w+") as gtf_list_out:
        for gtf_file in scallop_dir.glob("*.scallop.gtf"):#glob.glob(str(scallop_dir / "*.scallop.gtf")):
            transcript_count = check_gtf_content(gtf_file, "transcript")
            if transcript_count > 0:
                gtf_list_out.write(f"{gtf_file}\n")
            else:
                logger.warning(
                    "Warning, skipping file with no transcripts:\n%s" % gtf_file
                )

    logger.info("Merging Scallop results to:\n%s" % stringtie_merge_output_file)
    stringtie_command = [
         stringtie_path,
         "--merge",
         "-o",
         stringtie_merge_output_file,
         stringtie_merge_input_file,
     ]
    subprocess.run(stringtie_command)


def run_trimming(
    main_output_dir: Union[pathlib.Path, str],
    short_read_fastq_dir,
    delete_pre_trim_fastq,
    num_threads: int,
):
    trim_galore_path = "trim_galore"
    check_exe(trim_galore_path)

    # TODO
    # update type upstream
    short_read_fastq_dir = pathlib.Path(short_read_fastq_dir)

    trim_dir = create_dir(main_output_dir, "trim_galore_output")

    file_types = ["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"]
    fastq_file_list = []
    for file_type in file_types:
        fastq_file_list.extend(short_read_fastq_dir.glob(file_type))

    fastq_paired_paths = create_paired_paths(fastq_file_list)

    for fastq_paired_files in fastq_paired_paths:
        logger.info(
            "FASTQ file path(s): %s" % list_to_string(fastq_paired_files, separator=",")
        )

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

    pool = multiprocessing.Pool(num_threads)
    for fastq_paired_files in fastq_paired_paths:
        pool.apply_async(
            multiprocess_trim_galore,
            args=(
                generic_trim_galore_cmd,
                fastq_paired_files,
                trim_dir,
            ),
        )
        if delete_pre_trim_fastq:
            for fastq_path in fastq_paired_files:
                logger.info(
                    "Removing original fastq file post trimming:\n%s" % fastq_path
                )
                subprocess.run(["rm", fastq_path])

    pool.close()
    pool.join()

    for trimmed_fastq_path in trim_dir.glob("*.fq.gz"):
        logger.info("Trimmed file path:\n%s" % trimmed_fastq_path)
        sub_patterns = r"|".join(("_val_1.fq", "_val_2.fq", "_trimmed.fq"))
        updated_file_path = re.sub(sub_patterns, ".fq", str(trimmed_fastq_path))
        logger.info("Updated file path:\n%s" % updated_file_path)
        subprocess.run(["mv", trimmed_fastq_path, updated_file_path])

        files_to_delete_list = []
        for file_type in file_types:
            files_to_delete_list.extend(short_read_fastq_dir.glob(file_type))


def multiprocess_trim_galore(
    generic_trim_galore_cmd,
    fastq_paired_files: List[pathlib.Path],
    trim_dir: pathlib.Path,
):
    fastq_file = fastq_files[0]

    if len(fastq_files) == 2:
        fastq_file_pair = fastq_files[1]
    # len(fastq_files) == 1
    else:
        fastq_file_pair = None

    trim_galore_cmd = generic_trim_galore_cmd
    if fastq_file_pair:
        trim_galore_cmd.append("--paired")

    trim_galore_cmd.append(fastq_file)

    if fastq_file_pair:
        trim_galore_cmd.append(fastq_file_pair)

    logger.info(
        "Running Trim Galore with the following command:\n%s"
        % list_to_string(trim_galore_cmd)
    )
    subprocess.run(trim_galore_cmd)


def run_star_align(
    star_path,
    trim_fastq,
    subsample_script_path,
    main_output_dir: pathlib.Path,
    short_read_fastq_dir,
    genome_file: Union[pathlib.Path, str],
    max_reads_per_sample,
    max_total_reads,
    max_intron_length: int,
    num_threads: int,
    main_script_dir: pathlib.Path,
):
    # TODO
    # !!! Need to add in samtools path above instead of just using 'samtools' in command

    # TODO
    # use and pass Path objects directly as arguments to the function
    short_read_fastq_dir = pathlib.Path(short_read_fastq_dir)

    if not star_path:
        star_path = "STAR"

    check_exe(star_path)

    # If trimming has been enabled then switch the path for short_read_fastq_dir from the original location to the trimmed fastq dir
    if trim_fastq:
        short_read_fastq_dir = main_output_dir / "trim_galore_output"

    if (
        subsample_script_path is None
        or not pathlib.Path(subsample_script_path).exists()
    ):
        subsample_script_path = (
            main_script_dir / "support_scripts" / "subsample_fastq.py"
        )

    star_dir = create_dir(main_output_dir, "star_output")

    log_final_out = star_dir / "Log.final.out"
    log_out = star_dir / "Log.out"
    log_progress_out = star_dir / "Log.progress.out"
    if log_final_out.is_file() and log_out.is_file() and log_progress_out.is_file():
        logger.info("Skipping analysis, log files already exist")
        return

    star_tmp_dir = star_dir / "tmp"
    # delete directory if already exists
    shutil.rmtree(star_tmp_dir, ignore_errors=True)

    star_index_file = star_dir / "SAindex"

    fastq_file_list = []
    file_types = ["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"]
    for file_type in file_types:
        fastq_file_list.extend(short_read_fastq_dir.glob(file_type))

    # This works out if the files are paired or not
    fastq_paired_paths = create_paired_paths(fastq_file_list)

    # Subsamples in parallel if there's a value set
    if max_reads_per_sample:
        pool = multiprocessing.Pool(num_threads)
        for fastq_files in fastq_paired_paths:
            fastq_file = fastq_files[0]

            if len(fastq_files) == 2:
                fastq_file_pair = fastq_files[1]
            # len(fastq_files) == 1
            else:
                fastq_file_pair = None

            fastq_file_subsampled = fastq_file.parent / f"{fastq_file.name}.sub"
            fastq_file_pair_subsampled = (
                fastq_file_pair.parent / f"{fastq_file_pair.name}.sub"
            )
            if (
                fastq_file_pair
                and fastq_file_subsampled.exists()
                and fastq_file_pair_subsampled.exists()
            ):
                logger.info(
                    "Found existing .sub files on the fastq path for both members of the pair, will use those instead of subsampling again:\n%s\n%s"
                    % (fastq_file_subsampled, fastq_file_pair_subsampled)
                )
            elif fastq_file_pair:
                pool.apply_async(
                    run_subsample_script,
                    args=(
                        subsample_script_path,
                        fastq_file,
                        fastq_file_pair,
                    ),
                )
            elif fastq_file_subsampled.exists():
                logger.info(
                    "Found an existing .sub file on the fastq path, will use that instead:\n%s"
                    % f"{fastq_file}.sub"
                )
            else:
                pool.apply_async(
                    run_subsample_script,
                    args=(
                        subsample_script_path,
                        fastq_file,
                        fastq_file_pair,
                    ),
                )

        pool.close()
        pool.join()

    fastq_paired_paths = check_for_fastq_subsamples(fastq_paired_paths)

    if not fastq_paired_paths:
        raise IndexError(f"Empty fastq files list. Fastq dir:\n{short_read_fastq_dir}")

    if not star_index_file.is_file():
        logger.info("STAR index file does not exist, generating new one.")
        seq_region_lengths = get_seq_region_lengths(genome_file, 0)
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
                f"{star_dir}/",
                "--genomeDir",
                star_dir,
                "--genomeSAindexNbases",
                str(index_bases),
                "--genomeFastaFiles",
                genome_file,
            ]
        )
    if not star_index_file.is_file():
        raise IOError(f"STAR index file failed to be generated at:\n{star_index_file}")

    logger.info("Running STAR on the files in the fastq dir")
    for fastq_files in fastq_paired_paths:
        # STAR input string of comma separated fastq_files paths
        fastq_files_string = list_to_string(fastq_files, separator=",")
        # The name of the last file in the fastq_files list is used to generate
        # derivative files, in order to mirror the original program implementation.
        # Maybe a better base name to denote all files in fastq_files could be used here.
        representative_fastq_file = fastq_files[-1]
        logger.info("representative_fastq_file: %s" % representative_fastq_file)
        representative_fastq_file_name = representative_fastq_file.name
        check_compression = re.search(r".gz$", representative_fastq_file_name)
        # If there's a tmp dir already, the most likely cause is that STAR failed on the previous input file(s)
        # In this case STAR would effectively break for all the rest of the files, as it won't run if the tmp
        # dir exists. So clean it up and put out a warning. Another approach would be to just name the tmp dir
        # uniquely. Also there should be code checking the return on STAR anyway. So later when this is being
        # cleaned up to have proper tests, need to decide on the best implementation
        if star_tmp_dir.exists():
            logger.error(
                "Found existing tmp dir, implies potential failure on previous file, deleting"
            )
            shutil.rmtree(star_tmp_dir)

        sam_file_path = star_dir / f"{representative_fastq_file_name}.sam"
        sam_temp_file_path = star_dir / f"{sam_file_path.name}.tmp"
        bam_sort_file_path = sam_file_path.with_suffix(".bam")

        if bam_sort_file_path.is_file() and bam_sort_file_path.stat().st_size > 0:
            logger.info(
                "Existing bam file for fastq file found, skipping processing file"
            )
            continue

        logger.info("Processing %s" % fastq_files_string)
        # star_command = [star_path,'--outFilterIntronMotifs','RemoveNoncanonicalUnannotated','--outSAMstrandField','intronMotif','--runThreadN',str(num_threads),'--twopassMode','Basic','--runMode','alignReads','--genomeDir',star_dir,'--readFilesIn',fastq_file,'--outFileNamePrefix',(star_dir + '/'),'--outTmpDir',star_tmp_dir,'--outSAMtype','SAM','--alignIntronMax',str(max_intron_length),'--outSJfilterIntronMaxVsReadN','5000','10000','25000','40000','50000','50000','50000','50000','50000','100000']

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
            fastq_files_string,
            "--outFileNamePrefix",
            f"{star_dir}/",
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
        (star_dir / "Aligned.out.sam").rename(sam_file_path)
        junctions_file_path = star_dir / f"{representative_fastq_file_name}.sj.tab"
        (star_dir / "SJ.out.tab").rename(junctions_file_path)

        logger.info("Converting samfile into sorted bam file:\n%s" % bam_sort_file_path)
        samtools_command = [
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
        subprocess.run(samtools_command)

        os.remove(sam_file_path)

    logger.info("Completed running STAR")    

def create_paired_paths(
    fastq_file_paths: List[pathlib.Path],
) -> List[List[pathlib.Path]]:
    paired_paths = []
    path_dict = {}
    for fastq_file_path in fastq_file_paths:
        match = re.search(r"(.+)_\d+\.(fastq|fq)", str(fastq_file_path))
        if match:
            prefix = match.group(1)
            if prefix not in path_dict:
                path_dict[prefix] = []
            path_dict[prefix].append(fastq_file_path)
        # not match
        else:
            logger.error(
                "Could not find _1 or _2 at the end of the prefix for file. Assuming file is not paired: %s"
                % fastq_file_path
            )
            paired_paths.append([fastq_file_path])

    for pair in path_dict.values():
        paired_paths.append(pair)

    return paired_paths   


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

            # For the moment treat multimapping and single mapping things as a combined score
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

def coallate_results(main_output_dir):

    results_dir = create_dir(main_output_dir, "results")
    output_dirs = [
        "augustus_output",
        "cpg_output",
        "dust_output",
        "eponine_output",
        "red_output",
        "rfam_output",
        "trf_output",
        "trnascan_output",
    ]
    for output_dir in output_dirs:
        match = re.search(r"(.+)_output", output_dir)
        result_type = match.group(1)
        results_path = os.path.join(main_output_dir, output_dir, "annotation.gtf")
        copy_path = os.path.join(results_dir, (result_type + ".gtf"))
        if os.path.exists(results_path):
            cpy_cmd = ["cp", results_path, copy_path]
            subprocess.run(cpy_cmd)
    
