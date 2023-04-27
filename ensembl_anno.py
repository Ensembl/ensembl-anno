# pylint: disable=logging-not-lazy, invalid-name, missing-function-docstring, subprocess-run-check, unused-variable, redefined-outer-name, too-many-arguments, too-many-locals, too-many-branches, too-many-statements, unused-argument, no-else-return, undefined-variable, no-else-continue, no-else-raise, missing-docstring, consider-swap-variables, consider-using-in, too-many-lines, unused-import
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
import sys
import argparse
import errno
import gc
import glob
import io
import json
import logging
import logging.config
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

import repeatmasking_utils
import simple_feature_utils
import utils

with open(os.environ["ENSCODE"] + "/ensembl-anno/config.json", "r") as f:
    config = json.load(f)


def load_results_to_ensembl_db(
    main_script_dir,
    load_to_ensembl_db,
    genome_file,
    main_output_dir,
    db_details,
    num_threads,
):

    db_loading_script = os.path.join(
        main_script_dir, "support_scripts_perl", "load_gtf_ensembl.pl"
    )
    db_loading_dir = utils.create_dir(main_output_dir, "db_loading")

    # Should collapse this into a function
    annotation_results_gtf_file = os.path.join(
        main_output_dir, "annotation_output", "annotation.gtf"
    )
    if os.path.exists(annotation_results_gtf_file):
        logger.info("Loading main geneset to db")
        batch_size = 200
        load_type = "gene"
        analysis_name = "ensembl"
        gtf_records = batch_gtf_records(
            annotation_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logger.error(
            "Did not find the main gene annotation file, so not loading. Path checked:\n"
            + annotation_results_gtf_file
        )

    rfam_results_gtf_file = os.path.join(main_output_dir, "rfam_output", "annotation.gtf")
    if os.path.exists(rfam_results_gtf_file):
        logger.info("Loading Rfam-based sncRNA genes to db")
        batch_size = 500
        load_type = "gene"
        analysis_name = "ncrna"
        gtf_records = batch_gtf_records(
            rfam_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logger.error(
            "Did not find an Rfam annotation file, so not loading. Path checked:\n"
            + rfam_results_gtf_file
        )

    trnascan_results_gtf_file = os.path.join(
        main_output_dir, "trnascan_output", "annotation.gtf"
    )
    if os.path.exists(trnascan_results_gtf_file):
        logger.info("Loading tRNAScan-SE tRNA genes to db")
        batch_size = 500
        load_type = "gene"
        analysis_name = "ncrna"
        gtf_records = batch_gtf_records(
            trnascan_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logger.error(
            "Did not find an tRNAScan-SE annotation file, so not loading. Path checked:\n"
            + trnascan_results_gtf_file
        )

    dust_results_gtf_file = os.path.join(main_output_dir, "dust_output", "annotation.gtf")
    if os.path.exists(dust_results_gtf_file):
        logger.info("Loading Dust repeats to db")
        batch_size = 500
        load_type = "single_line_feature"
        analysis_name = "dust"
        gtf_records = batch_gtf_records(
            dust_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logger.error(
            "Did not find a Dust annotation file, so not loading. Path checked:\n"
            + dust_results_gtf_file
        )

    red_results_gtf_file = os.path.join(main_output_dir, "red_output", "annotation.gtf")
    if os.path.exists(red_results_gtf_file):
        logger.info("Loading Red repeats to db")
        batch_size = 500
        load_type = "single_line_feature"
        analysis_name = "repeatdetector"
        gtf_records = batch_gtf_records(
            red_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logger.error(
            "Did not find a Red annotation file, so not loading. Path checked:\n"
            + red_results_gtf_file
        )

    trf_results_gtf_file = os.path.join(main_output_dir, "trf_output", "annotation.gtf")
    if os.path.exists(trf_results_gtf_file):
        logger.info("Loading TRF repeats to db")
        batch_size = 500
        load_type = "single_line_feature"
        analysis_name = "trf"
        gtf_records = batch_gtf_records(
            trf_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logger.error(
            "Did not find a TRF annotation file, so not loading. Path checked:\n"
            + trf_results_gtf_file
        )

    cpg_results_gtf_file = os.path.join(main_output_dir, "cpg_output", "annotation.gtf")
    if os.path.exists(cpg_results_gtf_file):
        logger.info("Loading CpG islands to db")
        batch_size = 500
        load_type = "single_line_feature"
        analysis_name = "cpg"
        gtf_records = batch_gtf_records(
            cpg_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logger.error(
            "Did not find a CpG annotation file, so not loading. Path checked:\n"
            + cpg_results_gtf_file
        )

    eponine_results_gtf_file = os.path.join(
        main_output_dir, "eponine_output", "annotation.gtf"
    )
    if os.path.exists(eponine_results_gtf_file):
        logger.info("Loading Eponine repeats to db")
        batch_size = 500
        load_type = "single_line_feature"
        analysis_name = "eponine"
        gtf_records = batch_gtf_records(
            eponine_results_gtf_file, batch_size, db_loading_dir, load_type
        )
        generic_load_records_to_ensembl_db(
            load_to_ensembl_db,
            db_loading_script,
            genome_file,
            db_details,
            db_loading_dir,
            load_type,
            analysis_name,
            gtf_records,
            num_threads,
        )
    else:
        logger.error(
            "Did not find an Eponine annotation file, so not loading. Path checked:\n"
            + eponine_results_gtf_file
        )

    logger.info("Finished loading records to db")


def generic_load_records_to_ensembl_db(
    load_to_ensembl_db,
    db_loading_script,
    genome_file,
    db_details,
    db_loading_dir,
    load_type,
    analysis_name,
    gtf_records,
    num_threads,
):

    pool = multiprocessing.Pool(int(num_threads))
    for record_batch in gtf_records:
        pool.apply_async(
            multiprocess_load_records_to_ensembl_db,
            args=(
                load_to_ensembl_db,
                db_loading_script,
                genome_file,
                db_details,
                db_loading_dir,
                load_type,
                analysis_name,
                record_batch,
            ),
        )

    pool.close()
    pool.join()


def multiprocess_load_records_to_ensembl_db(
    load_to_ensembl_db,
    db_loading_script,
    genome_file,
    db_details,
    output_dir,
    load_type,
    analysis_name,
    record_batch,
):

    with tempfile.NamedTemporaryFile(
        mode="w+t", delete=False, dir=output_dir
    ) as gtf_temp_out:
        for line in record_batch:
            gtf_temp_out.writelines(line)
            gtf_temp_file_path = gtf_temp_out.name

    (db_name, db_host, db_port, db_user, db_pass) = db_details.split(",")

    loading_cmd = [
        "perl",
        db_loading_script,
        "-genome_file",
        genome_file,
        "-dbname",
        db_name,
        "-host",
        db_host,
        "-port",
        str(db_port),
        "-user",
        db_user,
        "-pass",
        db_pass,
        "-gtf_file",
        gtf_temp_file_path,
        "-analysis_name",
        analysis_name,
        "-load_type",
        load_type,
    ]

    if load_type == "gene" and analysis_name == "ensembl":
        loading_cmd.extend(
            [
                "-protein_coding_biotype",
                "anno_protein_coding",
                "-non_coding_biotype",
                "anno_lncRNA",
            ]
        )

        if load_to_ensembl_db == "single_transcript_genes":
            loading_cmd.append("-make_single_transcript_genes")

    logger.info(" ".join(loading_cmd))
    subprocess.run(loading_cmd)
    gtf_temp_out.close()
    os.remove(gtf_temp_file_path)  # doesn't seem to be working
    logger.info("Finished: " + gtf_temp_file_path)
    gc.collect()


def batch_gtf_records(input_gtf_file, batch_size, output_dir, record_type):

    records = []
    gtf_in = open(input_gtf_file)
    if record_type == "gene":

        # NOTE that the neverending variations on GTF reading/writing/merging
        # is becoming very messy
        # need to create a set of utility methods outside of this script
        # This one assumes the file has unique ids for the parent features.
        # It then batches them into
        # sets of records based on the batch size passed in
        record_counter = 0
        current_record_batch = []
        current_gene_id = ""
        line = gtf_in.readline()
        while line:
            if re.search(r"^#", line):
                line = gtf_in.readline()
                continue

            eles = line.split("\t")
            if not len(eles) == 9:
                line = gtf_in.readline()
                continue

            match = re.search(r'gene_id "([^"]+)"', line)
            gene_id = match.group(1)

            if not current_gene_id:
                record_counter += 1
                current_gene_id = gene_id

            if not gene_id == current_gene_id:
                record_counter += 1
                if record_counter % batch_size == 0:
                    records.append(current_record_batch)
                    current_record_batch = []
                current_gene_id = gene_id

            current_record_batch.append(line)
            line = gtf_in.readline()

        records.append(current_record_batch)

    elif record_type == "single_line_feature":
        record_counter = 0
        current_record_batch = []
        current_gene_id = ""
        line = gtf_in.readline()
        while line:
            if re.search(r"^#", line):
                line = gtf_in.readline()
                continue

            eles = line.split("\t")
            if not len(eles) == 9:
                line = gtf_in.readline()
                continue

            record_counter += 1

            if record_counter % batch_size == 0:
                records.append(current_record_batch)
                current_record_batch = []

            current_record_batch.append(line)
            line = gtf_in.readline()

        records.append(current_record_batch)

    gtf_in.close()

    return records


def run_find_orfs(genome_file, main_output_dir):

    min_orf_length = 600

    orf_output_dir = utils.create_dir(main_output_dir, "orf_output")
    seq_region_lengths = utils.get_seq_region_lengths(genome_file, 5000)
    for region_name in seq_region_lengths:
        region_length = seq_region_lengths[region_name]
        seq = utils.get_sequence(
            region_name, 1, region_length, 1, genome_file, orf_output_dir
        )
        for phase in range(0, 6):
            find_orf_phased_region(
                region_name, seq, phase, min_orf_length, orf_output_dir
            )


def find_orf_phased_region(region_name, seq, phase, min_orf_length, orf_output_dir):

    current_index = phase
    orf_counter = 1
    if phase > 2:
        seq = reverse_complement(seq)
        current_index = current_index % 3

    orf_file_path = os.path.join(
        orf_output_dir, (region_name + ".phase" + str(phase) + ".orf.fa")
    )
    orf_out = open(orf_file_path, "w+")

    while current_index < len(seq):
        codon = seq[current_index : current_index + 3]
        if codon == "ATG":
            orf_seq = codon
            for j in range(current_index + 3, len(seq), 3):
                next_codon = seq[j : j + 3]
                if next_codon == "TAA" or next_codon == "TAG" or next_codon == "TGA":
                    orf_seq += next_codon
                    if len(orf_seq) >= min_orf_length:
                        orf_out.write(
                            ">"
                            + region_name
                            + "_phase"
                            + str(phase)
                            + "_orf"
                            + str(orf_counter)
                            + "\n"
                        )
                        orf_out.write(orf_seq + "\n")
                        orf_counter += 1
                        orf_seq = ""
                        break

                # If there's another met in phase, then put i to the start of
                # the codon after j so that only the longest ORF is found
                if next_codon == "ATG":
                    current_index = j + 3
                orf_seq += next_codon
        current_index += 3
    orf_out.close()





def prlimit_command(command_list, virtual_memory_limit):

    """
    Prepend memory limiting arguments to a command list to be run with subprocess.

    This method uses the `prlimit` program to set the memory limit.

    The `virtual_memory_limit` size is in bytes.

    prlimit arguments:
    -v, --as[=limits]
           Address space limit.
    """
    return ["prlimit", f"-v{virtual_memory_limit}"] + command_list




def convert_gff_to_gtf(gff_file):
    gtf_string = ""
    file_in = open(gff_file)
    line = file_in.readline()
    while line:
        #    match = re.search(r"genBlastG",line)
        #    if match:
        results = line.split()
        if not len(results) == 9:
            line = file_in.readline()
            continue
        if results[2] == "coding_exon":
            results[2] = "exon"
        attributes = set_attributes(results[8], results[2])
        results[8] = attributes
        converted_line = "\t".join(results)
        gtf_string += converted_line + "\n"
        line = file_in.readline()
    file_in.close()
    return gtf_string


def set_attributes(attributes, feature_type):

    converted_attributes = ""
    split_attributes = attributes.split(";")
    if feature_type == "transcript":
        match = re.search(r"Name\=(.+)$", split_attributes[1])
        name = match.group(1)
        converted_attributes = 'gene_id "' + name + '"; transcript_id "' + name + '";'
    elif feature_type == "exon":
        match = re.search(r"\-E(\d+);Parent\=(.+)\-R\d+\-\d+\-", attributes)
        exon_rank = match.group(1)
        name = match.group(2)
        converted_attributes = (
            'gene_id "'
            + name
            + '"; transcript_id "'
            + name
            + '"; exon_number "'
            + exon_rank
            + '";'
        )

    return converted_attributes



def bed_to_gtf(minimap2_output_dir):

    gtf_file_path = os.path.join(minimap2_output_dir, "annotation.gtf")
    gtf_out = open(gtf_file_path, "w+")
    exons_dict = {}
    gene_id = 1
    for bed_file in glob.glob(minimap2_output_dir + "/*.bed"):
        logger.info("Converting bed to GTF:")
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
            transcript_line = [
                seq_region_name,
                "minimap",
                "transcript",
                0,
                0,
                ".",
                strand,
                ".",
                'gene_id "minimap_'
                + str(gene_id)
                + '"; '
                + 'transcript_id "minimap_'
                + str(gene_id)
                + '"',
            ]
            transcript_start = None
            transcript_end = None
            exon_records = []
            for i, exon_coords in enumerate(exons):
                if transcript_start is None or exon_coords[0] < transcript_start:
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
                    'gene_id "minimap_'
                    + str(gene_id)
                    + '"; '
                    + 'transcript_id "minimap_'
                    + str(gene_id)
                    + '"; exon_number "'
                    + str(i + 1)
                    + '";',
                ]

                exon_records.append(exon_line)

            transcript_line[3] = str(transcript_start)
            transcript_line[4] = str(transcript_end)

            gtf_out.write("\t".join(transcript_line) + "\n")
            for exon_line in exon_records:
                gtf_out.write("\t".join(exon_line) + "\n")

            gene_id += 1

    gtf_out.close()


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




def splice_junction_to_gff(input_dir, hints_file):

    sjf_out = open(hints_file, "w+")

    for sj_tab_file in glob.glob(input_dir + "/*.sj.tab"):
        sjf_in = open(sj_tab_file)
        sjf_lines = sjf_in.readlines()
        for line in sjf_lines:
            elements = line.split("\t")
            strand = "+"
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




def split_genome(genome_file, target_dir, min_seq_length):
    # This is the lazy initial way of just splitting into a dir
    # of files based on the toplevel sequence with a min sequence length filter
    # There are a couple of obvious improvements:
    # 1) Instead of making files for all seqs, just process N seqs
    #    parallel, where N = num_threads. Then you could clean up the seq file
    #    after each seq finishes, thus avoiding potentially
    #    having thousands of file in a dir
    # 2) Split the seq into even slices and process these in parallel
    #    (which the same cleanup as in 1). For sequences smaller than the
    #    target slice size, bundle them up together into a single file. Vastly
    #    more complex, partially implemented in the splice_genome method
    #    Allows for more consistency with parallelisation (since there should
    #    be no large outliers). But require a mapping strategy for the
    #    coords and sequence names and all the hints will need to be adjusted
    current_header = ""
    current_seq = ""

    file_in = open(genome_file)
    line = file_in.readline()
    while line:
        match = re.search(r">(.+)$", line)
        if match and current_header:
            if len(current_seq) > min_seq_length:
                file_out_name = os.path.join(target_dir, (current_header + ".split.fa"))
                if not os.path.exists(file_out_name):
                    file_out = open(file_out_name, "w+")
                    file_out.write(">" + current_header + "\n" + current_seq + "\n")
                    file_out.close()

                else:
                    print(
                        "Found an existing split file, so will not overwrite. File found:"
                    )
                    print(file_out_name)

            current_seq = ""
            current_header = match.group(1)
        elif match:
            current_header = match.group(1)
        else:
            current_seq += line.rstrip()

        line = file_in.readline()

    if len(current_seq) > min_seq_length:
        file_out_name = os.path.join(target_dir, (current_header + ".split.fa"))
        if not os.path.exists(file_out_name):
            file_out = open(file_out_name, "w+")
            file_out.write(">" + current_header + "\n" + current_seq + "\n")
            file_out.close()

        else:
            logger.info(
                "Found an existing split file, so will not overwrite. File found:"
            )
            logger.info(file_out_name)

    file_in.close()


def run_finalise_geneset(
    main_script_dir,
    main_output_dir,
    genome_file,
    seq_region_names,
    validation_type,
    diamond_validation_db,
    num_threads,
):

    if validation_type is None:
        logger.info("Setting validation type to relaxed")
    else:
        logger.info("Setting validation type to " + validation_type)

    final_annotation_dir = utils.create_dir(main_output_dir, "annotation_output")
    region_annotation_dir = utils.create_dir(final_annotation_dir, "initial_region_gtfs")
    final_region_annotation_dir = utils.create_dir(
        final_annotation_dir, "final_region_gtfs"
    )
    utr_region_annotation_dir = utils.create_dir(final_annotation_dir, "utr_region_gtfs")
    validation_dir = utils.create_dir(final_annotation_dir, "cds_validation")
    seq_region_lengths = utils.get_seq_region_lengths(genome_file, 0)

    logger.info("Skip analysis if the gtf file already exists")
    output_file = os.path.join(final_annotation_dir, "annotation.gtf")
    if os.path.exists(output_file):
        logger.info("final_annotation_dir exists")
        transcript_count = utils.check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Final gtf file exists")
            return
    else:
        logger.info("No gtf file, go on with the analysis")
    # This used to be a list of output dirs and a loop which was neat,
    # I'm coverting to a list of conditions as
    # it's more straightforward with the renaming
    # and having to merge scallop and stringtie
    protein_annotation_raw = os.path.join(
        main_output_dir, "genblast_output", "annotation.gtf"
    )
    minimap2_annotation_raw = os.path.join(
        main_output_dir, "minimap2_output", "annotation.gtf"
    )
    stringtie_annotation_raw = os.path.join(
        main_output_dir, "stringtie_output", "annotation.gtf"
    )
    scallop_annotation_raw = os.path.join(
        main_output_dir, "scallop_output", "annotation.gtf"
    )
    busco_annotation_raw = os.path.join(main_output_dir, "busco_output", "annotation.gtf")

    transcript_selector_script = os.path.join(
        main_script_dir, "support_scripts_perl", "select_best_transcripts.pl"
    )
    finalise_geneset_script = os.path.join(
        main_script_dir, "support_scripts_perl", "finalise_geneset.pl"
    )
    clean_geneset_script = os.path.join(
        main_script_dir, "support_scripts_perl", "clean_geneset.pl"
    )
    clean_utrs_script = os.path.join(
        main_script_dir, "support_scripts_perl", "clean_utrs_and_lncRNAs.pl"
    )
    gtf_to_seq_script = os.path.join(
        main_script_dir, "support_scripts_perl", "gtf_to_seq.pl"
    )

    transcriptomic_annotation_raw = os.path.join(
        final_annotation_dir, "transcriptomic_raw.gtf"
    )
    file_out = open(transcriptomic_annotation_raw, "w+")
    for transcriptomic_file in [
        minimap2_annotation_raw,
        scallop_annotation_raw,
        stringtie_annotation_raw,
    ]:

        if not os.path.exists(transcriptomic_file):
            logger.info(
                "No annotation.gtf file found in " + transcriptomic_file + ", skipping"
            )
            continue

        file_in = open(transcriptomic_file)
        line = file_in.readline()
        while line:
            print(line.rstrip(), file=file_out)
            line = file_in.readline()
        file_in.close()
    file_out.close()

    # Copy the raw files into the annotation dir, this is not needed
    # as such, but collecting them in one place and relabelling is
    # helpful for a user
    if os.path.exists(busco_annotation_raw):
        subprocess.run(
            [
                "cp",
                busco_annotation_raw,
                os.path.join(final_annotation_dir, "busco_raw.gtf"),
            ]
        )

    if os.path.exists(protein_annotation_raw):
        subprocess.run(
            [
                "cp",
                protein_annotation_raw,
                os.path.join(final_annotation_dir, "protein_raw.gtf"),
            ]
        )

    gtf_files = ["transcriptomic_raw.gtf", "protein_raw.gtf", "busco_raw.gtf"]
    generic_select_cmd = [
        "perl",
        transcript_selector_script,
        "-genome_file",
        genome_file,
    ]
    pool = multiprocessing.Pool(int(num_threads))
    for seq_region_name in seq_region_names:
        # The selection script needs different params depending on
        # whether the seqs are from transcriptomic data or not
        region_details = (
            seq_region_name + ".rs1" + ".re" + str(seq_region_lengths[seq_region_name])
        )
        transcriptomic_region_gtf_path = os.path.join(
            region_annotation_dir, (region_details + ".trans.gtf")
        )
        busco_region_gtf_path = os.path.join(
            region_annotation_dir, (region_details + ".busco.gtf")
        )
        protein_region_gtf_path = os.path.join(
            region_annotation_dir, (region_details + ".protein.gtf")
        )

        if os.path.exists(transcriptomic_annotation_raw):
            logger.info("Finalising transcriptomic data for: " + seq_region_name)
            transcriptomic_annotation_select = re.sub(
                "_raw.gtf", "_sel.gtf", transcriptomic_annotation_raw
            )
            cmd = generic_select_cmd.copy()
            cmd.extend(
                [
                    "-region_details",
                    region_details,
                    "-input_gtf_file",
                    transcriptomic_annotation_raw,
                    "-output_gtf_file",
                    transcriptomic_region_gtf_path,
                    "-cds_search",
                    "-final_biotype",
                    "transcriptomic",
                ]
            )
            pool.apply_async(multiprocess_finalise_geneset, args=(cmd,))

        if os.path.exists(busco_annotation_raw):
            logger.info("Finalising BUSCO data for: " + seq_region_name)
            busco_annotation_select = re.sub("_raw.gtf", "_sel.gtf", busco_annotation_raw)
            cmd = generic_select_cmd.copy()
            cmd.extend(
                [
                    "-region_details",
                    region_details,
                    "-input_gtf_file",
                    busco_annotation_raw,
                    "-output_gtf_file",
                    busco_region_gtf_path,
                    "-all_cds_exons",
                    "-final_biotype",
                    "busco",
                ]
            )
            pool.apply_async(multiprocess_finalise_geneset, args=(cmd,))

        if os.path.exists(protein_annotation_raw):
            logger.info("Finalising protein data for: " + seq_region_name)
            protein_annotation_select = re.sub(
                "_raw.gtf", "_sel.gtf", protein_annotation_raw
            )
            cmd = generic_select_cmd.copy()
            cmd.extend(
                [
                    "-region_details",
                    region_details,
                    "-input_gtf_file",
                    protein_annotation_raw,
                    "-output_gtf_file",
                    protein_region_gtf_path,
                    "-clean_transcripts",
                    "-all_cds_exons",
                    "-final_biotype",
                    "protein",
                ]
            )
            pool.apply_async(multiprocess_finalise_geneset, args=(cmd,))

    pool.close()
    pool.join()

    # At this point we will have the region files for all the,
    merge_finalise_output_files(
        final_annotation_dir,
        region_annotation_dir,
        ".trans.gtf",
        "transcriptomic",
    )
    merge_finalise_output_files(
        final_annotation_dir, region_annotation_dir, ".busco.gtf", "busco"
    )
    merge_finalise_output_files(
        final_annotation_dir, region_annotation_dir, ".protein.gtf", "protein"
    )

    # Create a single GTF file with all the selected transcripts
    # now that they have proper ids
    fully_merged_gtf_path = os.path.join(
        final_annotation_dir, "all_selected_transcripts.gtf"
    )
    fully_merged_gtf_out = open(fully_merged_gtf_path, "w+")

    merge_gtf_cmd = ["cat"]
    merge_gtf_cmd.extend(glob.glob(final_annotation_dir + "/*_sel.gtf"))
    subprocess.run(merge_gtf_cmd, stdout=fully_merged_gtf_out)
    fully_merged_gtf_out.close()

    # Now collapse the gene set
    generic_finalise_cmd = [
        "perl",
        finalise_geneset_script,
        "-genome_file",
        genome_file,
    ]

    pool = multiprocessing.Pool(int(num_threads))
    for seq_region_name in seq_region_names:
        region_details = (
            seq_region_name + ".rs1" + ".re" + str(seq_region_lengths[seq_region_name])
        )
        final_region_gtf_path = os.path.join(
            final_region_annotation_dir, (region_details + ".final.gtf")
        )

        cmd = generic_finalise_cmd.copy()
        cmd.extend(
            [
                "-region_details",
                region_details,
                "-input_gtf_file",
                fully_merged_gtf_path,
                "-output_gtf_file",
                final_region_gtf_path,
            ]
        )
        pool.apply_async(multiprocess_finalise_geneset, args=(cmd,))

    pool.close()
    pool.join()

    merge_finalise_output_files(
        final_annotation_dir,
        final_region_annotation_dir,
        ".final.gtf",
        "prevalidation",
    )
    merged_gtf_file = os.path.join(final_annotation_dir, ("prevalidation_sel.gtf"))
    merged_cdna_file = os.path.join(final_annotation_dir, ("prevalidation_sel.cdna.fa"))
    merged_amino_acid_file = os.path.join(
        final_annotation_dir, ("prevalidation_sel.prot.fa")
    )
    updated_gtf_lines = validate_coding_transcripts(
        merged_cdna_file,
        merged_amino_acid_file,
        validation_dir,
        validation_type,
        diamond_validation_db,
        merged_gtf_file,
        num_threads,
    )
    postvalidation_gtf_file = os.path.join(final_annotation_dir, ("postvalidation.gtf"))
    file_out = open(postvalidation_gtf_file, "w+")
    for line in updated_gtf_lines:
        file_out.write(line)
    file_out.close()

    cleaned_initial_gtf_file = os.path.join(final_annotation_dir, ("cleaned_pre_utr.gtf"))
    cleaned_utr_gtf_file = os.path.join(final_annotation_dir, ("annotation.gtf"))

    logger.info("Cleaning initial set")
    cleaning_cmd = [
        "perl",
        clean_geneset_script,
        "-genome_file",
        genome_file,
        "-gtf_file",
        postvalidation_gtf_file,
        "-output_gtf_file",
        cleaned_initial_gtf_file,
    ]
    logger.info(" ".join(cleaning_cmd))
    subprocess.run(cleaning_cmd)

    # Clean UTRs
    generic_clean_utrs_cmd = [
        "perl",
        clean_utrs_script,
        "-genome_file",
        genome_file,
        "-input_gtf_file",
        cleaned_initial_gtf_file,
    ]
    pool = multiprocessing.Pool(int(num_threads))
    for seq_region_name in seq_region_names:
        region_details = (
            seq_region_name + ".rs1" + ".re" + str(seq_region_lengths[seq_region_name])
        )
        utr_region_gtf_path = os.path.join(
            utr_region_annotation_dir, (region_details + ".utr.gtf")
        )

        cmd = generic_clean_utrs_cmd.copy()
        cmd.extend(
            [
                "-region_details",
                region_details,
                "-input_gtf_file",
                cleaned_initial_gtf_file,
                "-output_gtf_file",
                utr_region_gtf_path,
            ]
        )
        pool.apply_async(multiprocess_generic, args=(cmd,))
    pool.close()
    pool.join()

    merge_finalise_output_files(
        final_annotation_dir,
        utr_region_annotation_dir,
        ".utr.gtf",
        "annotation",
    )
    subprocess.run(
        [
            "mv",
            os.path.join(final_annotation_dir, "annotation_sel.gtf"),
            cleaned_utr_gtf_file,
        ]
    )

    logger.info("Dumping transcript and translation sequences")
    dumping_cmd = [
        "perl",
        gtf_to_seq_script,
        "-genome_file",
        genome_file,
        "-gtf_file",
        cleaned_utr_gtf_file,
    ]
    logger.info(" ".join(dumping_cmd))
    subprocess.run(dumping_cmd)

    logger.info("Finished creating gene set")


def validate_coding_transcripts(
    cdna_file,
    amino_acid_file,
    validation_dir,
    validation_type,
    diamond_validation_db,
    gtf_file,
    num_threads,
):

    logger.info("Running CDS validation with RNAsamba and CPC2")
    rnasamba_weights = config["rnasamba"]["weights"]
    rnasamba_output_path = os.path.join(validation_dir, "rnasamba.tsv.txt")
    cpc2_output_path = os.path.join(validation_dir, "cpc2.tsv")
    rnasamba_volume = validation_dir + "/:/app:rw"
    rnasamba_cmd = [
        "singularity",
        "exec",
        "--bind",
        rnasamba_volume,
        config["rnasamba"]["software"],
        "rnasamba",
        "classify",
        rnasamba_output_path,
        cdna_file,
        rnasamba_weights,
    ]
    logger.info(" ".join(rnasamba_cmd))
    subprocess.run(rnasamba_cmd)
    cpc2_volume = validation_dir + "/:/app:rw"
    cpc2_cmd = [
        "singularity",
        "exec",
        "--bind",
        cpc2_volume,
        config["cpc2"]["software"],
        "python3",
        "/CPC2_standalone-1.0.1/bin/CPC2.py",
        "-i",
        cdna_file,
        "--ORF",
        "-o",
        cpc2_output_path,
    ]
    logger.info(" ".join(cpc2_cmd))
    subprocess.run(cpc2_cmd)
    cpc2_output_path = cpc2_output_path + ".txt"

    check_file(rnasamba_output_path)
    check_file(cpc2_output_path)

    logger.info("diamond validation")
    diamond_results = None
    if diamond_validation_db is not None:
        diamond_output_dir = utils.create_dir(validation_dir, "diamond_output")
        diamond_validation(
            diamond_validation_db,
            amino_acid_file,
            diamond_output_dir,
            num_threads,
        )
        diamond_results = read_diamond_results(diamond_output_dir)

    logger.info("read results")
    rnasamba_results = read_rnasamba_results(rnasamba_output_path)
    cpc2_results = read_cpc2_results(cpc2_output_path)
    combined_results = combine_results(rnasamba_results, cpc2_results, diamond_results)
    logger.info("read gtf genes")
    parsed_gtf_genes = read_gtf_genes(gtf_file)
    updated_gtf_lines = update_gtf_genes(
        parsed_gtf_genes, combined_results, validation_type
    )

    return updated_gtf_lines


def diamond_validation(
    diamond_validation_db, amino_acid_file, diamond_output_dir, num_threads
):

    batched_protein_files = split_protein_file(amino_acid_file, diamond_output_dir, 100)

    pool = multiprocessing.Pool(int(num_threads))
    for batched_protein_file in batched_protein_files:
        pool.apply_async(
            multiprocess_diamond,
            args=(
                batched_protein_file,
                diamond_output_dir,
                diamond_validation_db,
            ),
        )
    pool.close()
    pool.join()


def multiprocess_diamond(
    batched_protein_file,
    diamond_output_dir,
    diamond_validation_db,
):

    batch_num = os.path.splitext(batched_protein_file)[0]
    batch_dir = os.path.dirname(batched_protein_file)
    diamond_output_file = batched_protein_file + ".dmdout"
    logger.info("Running diamond on " + batched_protein_file + ":")

    diamond_cmd = [
        "diamond",
        "blastp",
        "--query",
        batched_protein_file,
        "--db",
        diamond_validation_db,
        "--out",
        diamond_output_file,
    ]

    logger.info(" ".join(diamond_cmd))
    subprocess.run(diamond_cmd)
    subprocess.run(["mv", diamond_output_file, diamond_output_dir])


def update_gtf_genes(parsed_gtf_genes, combined_results, validation_type):

    output_lines = []

    for gene_id in parsed_gtf_genes.keys():
        transcript_ids = parsed_gtf_genes[gene_id].keys()
        for transcript_id in transcript_ids:
            transcript_line = parsed_gtf_genes[gene_id][transcript_id]["transcript"]
            single_cds_exon_transcript = 0
            translation_match = re.search(
                r'; translation_coords "([^"]+)";', transcript_line
            )
            if translation_match:
                translation_coords = translation_match.group(1)
                translation_coords_list = translation_coords.split(":")
                # If the start exon coords of both exons are the same,
                # then it's the same exon and thus a single exon cds
                if translation_coords_list[0] == translation_coords_list[3]:
                    single_cds_exon_transcript = 1

            exon_lines = parsed_gtf_genes[gene_id][transcript_id]["exons"]
            validation_results = combined_results[transcript_id]
            rnasamba_coding_probability = float(validation_results[0])
            rnasamba_coding_potential = validation_results[1]
            cpc2_coding_probability = float(validation_results[2])
            cpc2_coding_potential = validation_results[3]
            transcript_length = int(validation_results[4])
            peptide_length = int(validation_results[5])
            diamond_e_value = None
            if len(validation_results) == 7:
                diamond_e_value = validation_results[6]

            avg_coding_probability = (
                rnasamba_coding_probability + cpc2_coding_probability
            ) / 2
            max_coding_probability = max(
                rnasamba_coding_probability, cpc2_coding_probability
            )

            match = re.search(r'; biotype "([^"]+)";', transcript_line)
            biotype = match.group(1)
            if biotype == "busco" or biotype == "protein":
                transcript_line = re.sub(
                    '; biotype "' + biotype + '";',
                    '; biotype "protein_coding";',
                    transcript_line,
                )
                output_lines.append(transcript_line)
                output_lines.extend(exon_lines)
                continue

            min_single_exon_pep_length = 100
            min_multi_exon_pep_length = 75
            min_single_source_probability = 0.8
            min_single_exon_probability = 0.9

            # Note that the below looks at validating things
            # under different levels of strictness
            # There are a few different continue statements,
            # where transcripts will be skipped resulting
            # in a smaller post validation file. It mainly
            # removes single coding exon genes with no real
            # support or for multi-exon lncRNAs that are less than 200bp long
            if single_cds_exon_transcript == 1 and validation_type == "relaxed":
                if diamond_e_value is not None:
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    rnasamba_coding_potential == "coding"
                    and cpc2_coding_potential == "coding"
                    and peptide_length >= min_single_exon_pep_length
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    (
                        rnasamba_coding_potential == "coding"
                        or cpc2_coding_potential == "coding"
                    )
                    and peptide_length >= min_single_exon_pep_length
                    and max_coding_probability >= min_single_source_probability
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                else:
                    continue
            elif single_cds_exon_transcript == 1 and validation_type == "moderate":
                if (
                    diamond_e_value is not None
                    and peptide_length >= min_single_exon_pep_length
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    (
                        rnasamba_coding_potential == "coding"
                        and cpc2_coding_potential == "coding"
                    )
                    and peptide_length >= min_single_exon_pep_length
                    and avg_coding_probability >= min_single_exon_probability
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                else:
                    continue
            else:
                if diamond_e_value is not None:
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    rnasamba_coding_potential == "coding"
                    and cpc2_coding_potential == "coding"
                    and peptide_length >= min_multi_exon_pep_length
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    (
                        rnasamba_coding_potential == "coding"
                        or cpc2_coding_potential == "coding"
                    )
                    and peptide_length >= min_multi_exon_pep_length
                    and max_coding_probability >= min_single_source_probability
                ):
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif transcript_length >= 200:
                    transcript_line = re.sub(
                        '; biotype "' + biotype + '";',
                        '; biotype "lncRNA";',
                        transcript_line,
                    )
                    transcript_line = re.sub(
                        ' translation_coords "[^"]+";', "", transcript_line
                    )
                else:
                    continue

            output_lines.append(transcript_line)
            output_lines.extend(exon_lines)

    return output_lines


def read_rnasamba_results(file_path):

    results = []

    file_in = open(file_path)
    line = file_in.readline()
    while line:
        line = line.rstrip()
        match = re.search(r"^sequence_name", line)
        if match:
            line = file_in.readline()
            continue

        eles = line.split("\t")
        if not len(eles) == 3:
            line = file_in.readline()
            continue

        transcript_id = eles[0]
        coding_probability = eles[1]
        coding_potential = eles[2]
        results.append([transcript_id, coding_probability, coding_potential])
        line = file_in.readline()
    file_in.close()

    return results


def read_cpc2_results(file_path):

    results = []

    file_in = open(file_path)
    line = file_in.readline()
    while line:
        line = line.rstrip()
        match = re.search(r"^#ID", line)
        if match:
            line = file_in.readline()
            continue

        eles = line.split("\t")
        if not len(eles) == 9:
            line = file_in.readline()
            continue

        transcript_id = eles[0]
        transcript_length = eles[1]
        peptide_length = eles[2]
        coding_probability = eles[7]
        coding_potential = eles[8]
        results.append(
            [
                transcript_id,
                coding_probability,
                coding_potential,
                transcript_length,
                peptide_length,
            ]
        )
        line = file_in.readline()
    file_in.close()

    return results


def read_diamond_results(diamond_output_dir):

    results = []
    diamond_files = glob.glob(diamond_output_dir + "/*.dmdout")
    for file_path in diamond_files:
        file_in = open(file_path)
        line = file_in.readline()
        while line:
            line = line.rstrip()

            eles = line.split("\t")
            if not len(eles) == 12:
                line = file_in.readline()
                continue

            transcript_id = eles[0]
            e_value = eles[10]
            results.append([transcript_id, e_value])
            line = file_in.readline()
    file_in.close()

    return results


def combine_results(rnasamba_results, cpc2_results, diamond_results):

    transcript_ids = {}

    for result in rnasamba_results:
        transcript_id = result[0]
        coding_probability = result[1]
        coding_potential = result[2]

        if transcript_id not in transcript_ids:
            transcript_ids[transcript_id] = [
                coding_probability,
                coding_potential,
            ]

    for result in cpc2_results:
        transcript_id = result[0]
        coding_probability = result[1]
        coding_potential = result[2]
        transcript_length = result[3]
        peptide_length = result[4]
        transcript_ids[transcript_id].extend(
            [
                coding_probability,
                coding_potential,
                transcript_length,
                peptide_length,
            ]
        )

    if diamond_results is not None:
        for result in diamond_results:
            transcript_id = result[0]
            e_value = result[1]
            # There seems to be an issue where there are a small number
            # of sequences that don't make it into the cpc2/rnasamba output
            # Should code in a system for this, but it would be good to
            # understand why it happens to begin with. Seems to be the same
            # number of missing seqs in both, so maybe a shared cut-off
            if transcript_id in transcript_ids:
                transcript_ids[transcript_id].extend([e_value])

    return transcript_ids


def read_gtf_genes(gtf_file):

    gtf_genes = {}
    gtf_in = open(gtf_file)
    line = gtf_in.readline()
    while line:
        eles = line.split("\t")
        if not len(eles) == 9:
            line = gtf_in.readline()
            continue

        match = re.search(r'gene_id "([^"]+)".+transcript_id "([^"]+)"', line)

        if not match:
            line = gtf_in.readline()
            continue

        gene_id = match.group(1)
        transcript_id = match.group(2)
        feature_type = eles[2]
        if gene_id not in gtf_genes:
            gtf_genes[gene_id] = {}
        if feature_type == "transcript":
            gtf_genes[gene_id][transcript_id] = {}
            gtf_genes[gene_id][transcript_id]["transcript"] = line
            gtf_genes[gene_id][transcript_id]["exons"] = []
        elif feature_type == "exon":
            gtf_genes[gene_id][transcript_id]["exons"].append(line)
        line = gtf_in.readline()
    gtf_in.close()

    return gtf_genes


def merge_finalise_output_files(
    final_annotation_dir, region_annotation_dir, extension, id_label
):

    gtf_files = glob.glob(region_annotation_dir + "/*" + extension)

    merged_gtf_file = os.path.join(final_annotation_dir, (id_label + "_sel.gtf"))
    merged_cdna_file = os.path.join(final_annotation_dir, (id_label + "_sel.cdna.fa"))
    merged_amino_acid_file = os.path.join(
        final_annotation_dir, (id_label + "_sel.prot.fa")
    )

    # The below is not great, it's a bit messy because there might be
    # some cases where there aren't translations. So it's not as
    # straightforward as reading the records across all three files
    # in parallel. The solution is to just load the seqs into
    # memory and index them on the current header, which should
    # correspond to a transcript/gene id in the GTF. When writing the
    # results into the single merged files the ids will be updated to
    # be unique and consistent across the header,  three file types

    gene_id_counter = 0
    transcript_id_counter = 0
    gtf_out = open(merged_gtf_file, "w+")
    cdna_out = open(merged_cdna_file, "w+")
    amino_acid_out = open(merged_amino_acid_file, "w+")
    for gtf_file in gtf_files:
        logger.info("GTF file: " + gtf_file)
        cdna_seq_index = {}
        amino_acid_seq_index = {}
        cdna_file = gtf_file + ".cdna"
        amino_acid_file = gtf_file + ".prot"
        cdna_in = open(cdna_file)
        amino_acid_in = open(amino_acid_file)
        cdna_seq_index = fasta_to_dict(cdna_in.readlines())
        amino_acid_seq_index = fasta_to_dict(amino_acid_in.readlines())
        cdna_in.close()
        amino_acid_in.close()

        current_gene_id = ""
        gtf_in = open(gtf_file)
        line = gtf_in.readline()
        while line:
            if re.search(r"^#", line):
                line = gtf_in.readline()
                continue

            eles = line.split("\t")
            if not len(eles) == 9:
                line = gtf_in.readline()
                continue

            match = re.search(r'gene_id "([^"]+)".+transcript_id "([^"]+)"', line)
            if match and eles[2] == "transcript":
                transcript_id_counter += 1

            gene_id = match.group(1)
            transcript_id = match.group(2)

            if not current_gene_id:
                gene_id_counter += 1
                current_gene_id = gene_id

            if not gene_id == current_gene_id:
                gene_id_counter += 1
                current_gene_id = gene_id

            new_gene_id = id_label + "_" + str(gene_id_counter)
            new_transcript_id = id_label + "_" + str(transcript_id_counter)
            line = re.sub(
                'gene_id "' + gene_id + '"',
                ('gene_id "' + new_gene_id + '"'),
                line,
            )
            line = re.sub(
                'transcript_id "' + transcript_id + '"',
                ('transcript_id "' + new_transcript_id + '"'),
                line,
            )
            gtf_out.write(line)
            line = gtf_in.readline()

            if eles[2] == "transcript":
                new_header = ">" + new_transcript_id + "\n"
                cdna_out.write(new_header + cdna_seq_index[transcript_id])

                if transcript_id in amino_acid_seq_index:
                    amino_acid_out.write(new_header + amino_acid_seq_index[transcript_id])

    gtf_out.close()
    cdna_out.close()
    amino_acid_out.close()


def fasta_to_dict(fasta_list):

    index = {}
    it = iter(fasta_list)
    for header in it:
        match = re.search(r">(.+)\n$", header)
        header = match.group(1)
        seq = next(it)
        index[header] = seq
    return index


# def merge_gtf_files(file_paths,id_label):

#  gtf_file_path = os.path.join(output_dir,'annotation.gtf')
#  gtf_out = open(gtf_file_path,'w+')
#  for gtf_file_path in gtf_files:
#    gtf_file_name = os.path.basename(gtf_file_path)
#    match = re.search(r'\.rs(\d+)\.re(\d+)\.',gtf_file_name)
#    start_offset = int(match.group(1))
#    gtf_in = open(gtf_file_path,'r')
#    line = gtf_in.readline()
#    while line:
#      values = line.split("\t")
#      if len(values) == 9 and (values[2] in feature_types):
#        values[3] = str(int(values[3]) + (start_offset - 1))
#        values[4] = str(int(values[4]) + (start_offset - 1))
#        gtf_out.write("\t".join(values))
#        line = gtf_in.readline()
#    gtf_in.close()
#  gtf_out.close()


def multiprocess_generic(cmd):
    print(" ".join(cmd))
    subprocess.run(cmd)


def multiprocess_finalise_geneset(cmd):

    print(" ".join(cmd))
    subprocess.run(cmd)


def seq_region_names(genome_file):
    region_list = []

    file_in = open(genome_file)
    line = file_in.readline()
    while line:
        match = re.search(r">([^\s]+)", line)
        if match:
            region_name = match.group(1)
            if region_name == "MT":
                logger.info("Skipping region named MT")
                line = file_in.readline()
                continue
            else:
                region_list.append(match.group(1))
        line = file_in.readline()

    return region_list


def slice_genome(genome_file, target_dir, target_slice_size):
    # The below is sort of tested
    # Without the
    target_seq_length = 50000000
    min_seq_length = 1000
    current_header = ""
    current_seq = ""
    seq_dict = {}
    for line in seq:  # pylint:disable=used-before-assignment
        match = re.search(r">(.+)$", line)
        if match and current_header:
            seq_dict[current_header] = current_seq
            current_seq = ""
            current_header = match.group(1)
        elif match:
            current_header = match.group(1)
        else:
            current_seq += line.rstrip()

    seq_dict[current_header] = current_seq

    seq_buffer = 0
    file_number = 0
    file_name = "genome_file_" + str(file_number)

    for header in seq_dict:
        seq_iterator = 0
        seq = seq_dict[header]

        while len(seq) > target_seq_length:
            file_out = open(os.path.join(target_dir, file_name), "w+")
            subseq = seq[0:target_seq_length]
            file_out.write(
                ">" + header + "_sli" + str(seq_iterator) + "\n" + subseq + "\n"
            )
            file_out.close()
            seq = seq[target_seq_length:]
            seq_iterator += 1
            file_number += 1
            file_name = "genome_file_" + str(file_number)

        if len(seq) >= min_seq_length:
            file_name = "genome_file_" + str(file_number)
            file_out = open(os.path.join(file_name), "w+")
            file_out.write(">" + header + "_sli" + str(seq_iterator) + "\n" + seq + "\n")
            file_out.close()
            file_number += 1
            file_name = "genome_file_" + str(file_number)




def check_file(file_path):

    if not os.path.exists(file_path):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file_path)


def coallate_results(main_output_dir):

    results_dir = utils.create_dir(main_output_dir, "results")
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


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output_dir",
        type=str,
        help="Path where the output and temp files will write to. \
        Uses current dir by default",
    )
    parser.add_argument(
        "--genome_file",
        type=str,
        help="Path to the fasta genome file",
        required=True,
    )
    parser.add_argument(
        "--num_threads", type=int, default=1, help="Number of threads to use"
    )
    parser.add_argument(
        "--run_masking",
        action="store_true",
        help="Run Red to find repeats and softmask the genome. \
        Otherwise provide a softmasked genome",
    )
    parser.add_argument(
        "--red_path",
        type=str,
        help="Path to Red executable. See http://toolsmith.ens.utulsa.edu",
    )
    parser.add_argument(
        "--genblast_path",
        type=str,
        help="Path to GenBlast executable. See http://genome.sfu.ca/\
        genblast/download.html",
    )
    parser.add_argument(
        "--convert2blastmask_path",
        type=str,
        help="Path to convert2blastmask executable",
    )
    parser.add_argument(
        "--makeblastdb_path", type=str, help="Path to makeblastdb executable"
    )
    parser.add_argument(
        "--run_genblast",
        action="store_true",
        help="Run GenBlast to align protein sequences",
    )
    parser.add_argument(
        "--genblast_timeout",
        type=int,
        help="GenBlast timeout in seconds",
        default=10800,
    )
    parser.add_argument(
        "--run_busco",
        action="store_true",
        help="Run GenBlast to align BUSCO protein sequences",
    )
    parser.add_argument(
        "--protein_file",
        type=str,
        help="Path to a fasta file with protein sequences",
    )
    parser.add_argument(
        "--busco_protein_file",
        type=str,
        help="Path to a fasta file with BUSCO protein sequences",
    )
    parser.add_argument(
        "--rfam_accessions_file",
        type=str,
        help="Path to a file with Rfam CM accessions, one accession\
        per line, to use with cmsearch",
    )
    parser.add_argument(
        "--run_star",
        action="store_true",
        help="Run Star for short read alignment",
    )
    parser.add_argument(
        "--star_path", type=str, help="Path to Star for short read alignment"
    )
    parser.add_argument(
        "--max_reads_per_sample",
        nargs="?",
        type=int,
        default=0,
        help="The maximum number of reads to use per sample. Default=0 (unlimited)",
    )
    parser.add_argument(
        "--max_total_reads",
        nargs="?",
        type=int,
        default=0,
        help="The maximum total number of reads. Default=0 (unlimited)",
    )
    parser.add_argument(
        "--short_read_fastq_dir",
        help="Path to short read fastq dir for running with Star",
    )
    parser.add_argument(
        "--max_intron_length",
        nargs="?",
        type=int,
        default=100000,
        help="The maximum intron size for alignments. Default=100000",
    )
    parser.add_argument(
        "--run_minimap2",
        action="store_true",
        help="Run minimap2 for long read alignment",
    )
    parser.add_argument(
        "--minimap2_path",
        type=str,
        help="Path to minimap2 for long read alignment",
    )
    parser.add_argument(
        "--paftools_path",
        type=str,
        help="Path to paftools for SAM to BED conversion",
    )
    parser.add_argument(
        "--long_read_fastq_dir",
        type=str,
        help="Path to long read fastq dir for running with minimap2",
    )
    parser.add_argument(
        "--run_augustus",
        action="store_true",
        help="Run Augustus with hints for gene/transcript prediction",
    )
    parser.add_argument("--augustus_path", type=str, help="Path to Augustus")
    parser.add_argument(
        "--run_stringtie",
        action="store_true",
        help="Run Stringtie on the results from the STAR alignments",
    )
    parser.add_argument(
        "--run_scallop",
        action="store_true",
        help="Run Scallop on the results from the STAR alignments",
    )
    parser.add_argument("--stringtie_path", type=str, help="Path to Stringtie")
    parser.add_argument("--scallop_path", type=str, help="Path to Scallop")
    parser.add_argument(
        "--subsample_script_path",
        type=str,
        help="Path to ensembl-anno subsampling script",
    )
    parser.add_argument("--samtools_path", type=str, help="Path to subsampling script")
    parser.add_argument(
        "--finalise_geneset",
        action="store_true",
        help="Used to finalise the gene set from the various GTF files generated",
    )
    parser.add_argument(
        "--db_details",
        type=str,
        help="A comma separated string of dbname,host,port,user,pass",
    )
    parser.add_argument(
        "--run_cmsearch",
        action="store_true",
        help="Search for sncRNA structures using Rfam and cmsearch",
    )
    parser.add_argument(
        "--run_trf", action="store_true", help="Run TRF to find tandem repeats"
    )
    parser.add_argument("--trf_path", type=str, help="Path to TRF")
    parser.add_argument(
        "--run_dust",
        action="store_true",
        help="Run Dust to find low complexity regions",
    )
    parser.add_argument("--dust_path", type=str, help="Path to Dust")
    parser.add_argument(
        "--run_repeatmasker",
        action="store_true",
        help="Run RepeatMasker to find repeat regions",
    )
    parser.add_argument(
        "--repeatmasker_path", action="store_true", help="Path to RepeatMasker"
    )
    parser.add_argument(
        "--run_trnascan",
        action="store_true",
        help="Run tRNAscan-SE to find tRNAs",
    )
    parser.add_argument("--trnascan_path", type=str, help="Path to tRNAscan-SE")
    parser.add_argument(
        "--trnascan_filter_path",
        type=str,
        help="Path to tRNAscan-SE high confidence filter",
    )
    parser.add_argument(
        "--run_cpg", action="store_true", help="Run cpg_lh to find CpG islands"
    )
    parser.add_argument("--cpg_path", type=str, help="Path to cpg_lh")
    parser.add_argument(
        "--run_eponine",
        action="store_true",
        help="Run Eponine to find transcription start sites",
    )
    parser.add_argument("--eponine_path", type=str, help="Path to Eponine jar file")
    parser.add_argument("--java_path", type=str, help="Path to Java for use with Eponine")
    parser.add_argument(
        "--run_full_annotation",
        action="store_true",
        help="Run a full annotation, will automatically check \
        for input data and run tools based on that",
    )
    parser.add_argument("--run_repeats", action="store_true", help="Run Red, Dust, TRF")
    parser.add_argument(
        "--run_simple_features", action="store_true", help="Run CpG, Eponine"
    )
    parser.add_argument(
        "--run_sncrnas", action="store_true", help="Run Rfam, tRNAscan-SE"
    )
    parser.add_argument(
        "--run_transcriptomic",
        action="store_true",
        help="Run STAR, Stringtie2, Scallop, minimap2 (if short_read_fastq_dir \
        and/or long_read_fastq_dir are provided)",
    )
    parser.add_argument(
        "--run_proteins",
        action="store_true",
        help="Run GenBlast if protein_file and/or busco_protein_file",
    )
    parser.add_argument(
        "--diamond_validation_db",
        type=str,
        help="Use a Diamond db with blastp mode to help validate cds sequences",
    )
    parser.add_argument(
        "--validation_type",
        type=str,
        help='The strength of evidence needed to validate and ORF as \
        protein coding, can be "relaxed" or "moderate"',
    )
    parser.add_argument(
        "--load_to_ensembl_db",
        action="store_true",
        help="Load results to an Ensembl db, must also provide the db_details flag",
    )
    parser.add_argument(
        "--trim_fastq",
        action="store_true",
        help="Trim the short read files using Trim Galore",
    )
    parser.add_argument(
        "--delete_pre_trim_fastq",
        action="store_true",
        help="Delete the original fastq files after trimming",
    )
    parser.add_argument(
        "--repeatmasker_library",
        type=str,
        help="Specify library for repeatmasker",
    )
    parser.add_argument(
        "--repeatmasker_species",
        type=str,
        help="Specify species for repeatmasker (default homo)",
    )
    args = parser.parse_args()

    work_dir = args.output_dir
    genome_file = args.genome_file
    num_threads = args.num_threads
    # masked_genome_file = genome_file  # This will be updated later if Red is run
    run_masking = args.run_masking
    red_path = args.red_path
    genblast_path = args.genblast_path
    convert2blastmask_path = args.convert2blastmask_path
    makeblastdb_path = args.makeblastdb_path
    run_genblast = args.run_genblast
    genblast_timeout = args.genblast_timeout
    run_busco = args.run_busco
    protein_file = args.protein_file
    busco_protein_file = args.busco_protein_file
    rfam_accessions_file = args.rfam_accessions_file
    run_star = args.run_star
    star_path = args.star_path
    short_read_fastq_dir = args.short_read_fastq_dir
    max_intron_length = args.max_intron_length
    max_reads_per_sample = args.max_reads_per_sample
    max_total_reads = args.max_total_reads
    run_minimap2 = args.run_minimap2
    minimap2_path = args.minimap2_path
    paftools_path = args.paftools_path
    long_read_fastq_dir = args.long_read_fastq_dir
    run_augustus = args.run_augustus
    augustus_path = args.augustus_path
    run_stringtie = args.run_stringtie
    run_scallop = args.run_scallop
    stringtie_path = args.stringtie_path
    scallop_path = args.scallop_path
    subsample_script_path = args.subsample_script_path
    samtools_path = args.samtools_path
    finalise_geneset = args.finalise_geneset
    db_details = args.db_details
    run_cmsearch = args.run_cmsearch
    run_trf = args.run_trf
    trf_path = args.trf_path
    run_dust = args.run_dust
    dust_path = args.dust_path
    run_trnascan = args.run_trnascan
    trnascan_path = args.trnascan_path
    trnascan_filter_path = args.trnascan_filter_path
    run_cpg = args.run_cpg
    cpg_path = args.cpg_path
    run_eponine = args.run_eponine
    eponine_path = args.eponine_path
    java_path = args.java_path
    run_repeatmasker = args.run_repeatmasker
    repeatmasker_path = args.repeatmasker_path
    run_full_annotation = args.run_full_annotation
    run_repeats = args.run_repeats
    run_simple_features = args.run_simple_features
    run_sncrnas = args.run_sncrnas
    run_transcriptomic = args.run_transcriptomic
    run_proteins = args.run_proteins
    diamond_validation_db = args.diamond_validation_db
    validation_type = args.validation_type
    load_to_ensembl_db = args.load_to_ensembl_db
    trim_fastq = args.trim_fastq
    delete_pre_trim_fastq = args.delete_pre_trim_fastq
    library = args.repeatmasker_library
    species = args.repeatmasker_species

    main_script_dir = os.path.dirname(os.path.realpath(__file__))
    # work_dir=glob.glob(work_dir)
    if not os.path.exists(genome_file):
        raise IOError("File does not exist: %s" % genome_file)

    if not work_dir:
        work_dir = os.getcwd()
        # work_dir=glob.glob(work_dir)

    # set up logger
    log_file_path = pathlib.Path(work_dir) / "ensembl_anno.log"
    loginipath = pathlib.Path(os.environ["ENSCODE"] + "/ensembl-anno/logging.conf")
    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": log_file_path},
        disable_existing_loggers=False,
    )
    logger = logging.getLogger()
    logger.propagate = False

    logger.info("work directory: %s" % work_dir)
    if not os.path.exists(work_dir):
        logger.info("Work dir does not exist, will create")
        utils.create_dir(work_dir, None)

    if num_threads == 1:
        logger.info("Thread count is set to the default value 1; this might be slow.")

    if os.path.exists(
        os.path.join(work_dir, "red_output", "mask_output")
    ) or os.path.join(work_dir, "red_output", "mask_output").endswith(".msk"):
        red_genome_file = [
            f
            for f in os.listdir(os.path.join(work_dir, "red_output", "mask_output"))
            if f.endswith(".msk")
        ]
        logger.info("red_genome_file %s", red_genome_file)
        masked_genome_file = os.path.join(
            work_dir, "red_output", "mask_output", red_genome_file[0]
        )
    else:
        masked_genome_file = genome_file
    logger.info("Masked genome file %s", masked_genome_file)
    # If the run_full_annotation flag is set then we want to set
    # a standardised set of analyses
    if run_full_annotation:
        run_repeats = True
        run_simple_features = True
        run_sncrnas = True
        run_transcriptomic = True
        run_proteins = True
        finalise_geneset = True

    # These are subsets of the analyses that can be run, group by type
    if run_repeats:
        run_masking = True
        run_dust = True
        run_trf = True

    if run_simple_features:
        run_cpg = True
        run_eponine = True

    if run_sncrnas:
        if rfam_accessions_file:
            run_cmsearch = True
        run_trnascan = True

    if run_transcriptomic:
        if short_read_fastq_dir:
            run_star = True
            run_scallop = True
            run_stringtie = True
        if long_read_fastq_dir:
            run_minimap2 = True

    if run_proteins:
        if protein_file:
            run_genblast = True
        if busco_protein_file:
            run_busco = True

    # Collect a list of seq region names, most useful for multiprocessing regions
    seq_region_names = seq_region_names(genome_file)
    for i in seq_region_names:
        logger.info(i)

    #################################
    # Repeat analyses
    #################################
    if run_masking:
        logger.info("Running masking via Red")
        masked_genome_file = repeatmasking_utils.run_red(red_path, work_dir, genome_file)
        logger.info("Masked genome file: " + masked_genome_file)

    else:
        logger.info("Not running masking, presuming the genome file is softmasked")

    if run_dust:
        logger.info("Annotating low complexity regions")
        logger.info("run_dust genome file %s", genome_file)
        repeatmasking_utils.run_dust_regions(
            genome_file, dust_path, work_dir, num_threads
        )

    if run_trf:
        logger.info("Annotating tandem repeats")
        logger.info("run_trf genome file %s", genome_file)
        repeatmasking_utils.run_trf_repeats(genome_file, trf_path, work_dir, num_threads)

    if run_repeatmasker:
        logger.info("Annotating repeats with RepeatMasker")
        logger.info("run_repeatmasker genome file %s", genome_file)
        repeatmasking_utils.run_repeatmasker_regions(
            genome_file,
            repeatmasker_path,
            library,
            species,
            work_dir,
            num_threads,
        )

    #################################
    # Simple feature analyses
    #################################
    if run_cpg:
        logger.info("Annotating CpG islands")
        logger.info("run_cpg genome file %s", genome_file)
        simple_feature_utils.run_cpg_regions(genome_file, cpg_path, work_dir, num_threads)

    if run_eponine:
        logger.info("Running Eponine to find transcription start sites")
        logger.info("run_eponine genome file %s", genome_file)
        simple_feature_utils.run_eponine_regions(
            genome_file, java_path, eponine_path, work_dir, num_threads
        )

    #################################
    # sncRNA analyses
    #################################
    # Search Rfam with cmsearch
    if run_cmsearch:
        logger.info("Annotating sncRNAs")
        logger.info("run_cmsearch genome file %s", genome_file)
        run_cmsearch_regions(
            genome_file,
            None,
            None,
            None,
            rfam_accessions_file,
            work_dir,
            num_threads,
        )

    if run_trnascan:
        logger.info("Annotating tRNAs")
        logger.info("run_trnascan genome file %s", genome_file)
        run_trnascan_regions(
            genome_file,
            trnascan_path,
            trnascan_filter_path,
            work_dir,
            num_threads,
        )

    #################################
    # Transcriptomic analyses
    #################################
    if trim_fastq:
        run_trimming(work_dir, short_read_fastq_dir, delete_pre_trim_fastq, num_threads)

    # Run STAR
    if run_star:
        logger.info("Running Star")
        logger.info("run_star genome file %s", genome_file)
        run_star_align(
            star_path,
            trim_fastq,
            subsample_script_path,
            work_dir,
            short_read_fastq_dir,
            genome_file,
            max_reads_per_sample,
            max_total_reads,
            max_intron_length,
            num_threads,
        )

    # Run Scallop
    if run_scallop:
        logger.info("Running Scallop")
        run_scallop_assemble(scallop_path, stringtie_path, work_dir)

    # Run Stringtie
    if run_stringtie:
        logger.info("Running Stringtie")
        logger.info("run_stringtie genome file %s", genome_file)
        run_stringtie_assemble(
            stringtie_path, samtools_path, work_dir, genome_file, num_threads
        )

    # Run minimap2
    if run_minimap2:
        logger.info("Running minimap2")
        logger.info("run_minimap2 genome file %s", genome_file)
        run_minimap2_align(
            minimap2_path,
            paftools_path,
            work_dir,
            long_read_fastq_dir,
            genome_file,
            max_intron_length,
            num_threads,
        )

    if run_transcriptomic:
        check_transcriptomic_output(work_dir)

    #################################
    # Protein analyses
    #################################
    # Run GenBlast
    if run_genblast:
        logger.info("Running GenBlast")
        logger.info("run_genblast genome file %s", masked_genome_file)
        run_genblast_align(
            genblast_path,
            convert2blastmask_path,
            makeblastdb_path,
            os.path.join(work_dir, "genblast_output"),
            protein_file,
            masked_genome_file,
            max_intron_length,
            num_threads,
            genblast_timeout,
        )

    # Run GenBlast on BUSCO set, gives higher priority when creating thei
    # final genes in cases where transcriptomic data are missing or fragmented
    if run_busco:
        logger.info("Running GenBlast of BUSCO proteins")
        logger.info("run_busco genome file %s", masked_genome_file)
        run_genblast_align(
            genblast_path,
            convert2blastmask_path,
            makeblastdb_path,
            os.path.join(work_dir, "busco_output"),
            busco_protein_file,
            masked_genome_file,
            max_intron_length,
            num_threads,
            genblast_timeout,
        )

    #################################
    # Finalisation analyses
    #################################
    # Do some magic
    if finalise_geneset:
        logger.info("Finalise geneset")
        logger.info("finalise geneset genome file %s", genome_file)
        run_finalise_geneset(
            main_script_dir,
            work_dir,
            genome_file,
            seq_region_names,
            validation_type,
            diamond_validation_db,
            num_threads,
        )

    #################################
    # Other analyses
    #################################
    # Run Augustus
    if run_augustus:
        logger.info("Running Augustus")
        run_augustus_predict(augustus_path, work_dir, masked_genome_file, num_threads)

    if load_to_ensembl_db:
        load_results_to_ensembl_db(
            main_script_dir,
            load_to_ensembl_db,
            genome_file,
            work_dir,
            db_details,
            num_threads,
        )

#  coallate_results(work_dir)

#  run_find_orfs(genome_file,work_dir)
