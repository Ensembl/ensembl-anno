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

from typing import Union

# project imports
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
)

from repeatmasking_utils import (
    run_repeatmasker_regions,
    multiprocess_repeatmasker,
    create_repeatmasker_gtf,
    run_dust_regions,
    multiprocess_dust,
    create_dust_gtf,
    run_trf_repeats,
    multiprocess_trf,
    create_trf_gtf,
    run_red,
    create_red_gtf,
    )

from simple_features_utils import (
    run_eponine_regions,
    multiprocess_eponine,
    create_eponine_gtf,
    run_cpg_regions,
    multiprocess_cpg,
    create_cpg_gtf,
    )

from sncRNA_utils import (
    run_trnascan_regions,
    multiprocess_trnascan,
    create_trnascan_gtf,
    run_cmsearch_regions,
    multiprocess_cmsearch,
    get_rfam_seed_descriptions,
    extract_rfam_metrics,
    parse_rfam_tblout,
    remove_rfam_overlap,
    filter_rfam_results,
    create_rfam_gtf,
    check_rnafold_structure,
)


def slice_output_to_gtf(
    output_dir: pathlib.Path,
    extension,
    unique_ids,
    feature_id_label,
    new_id_prefix,
):
    if not extension:
        extension = ".gtf"

    # Note that this does not make unique ids at the moment
    # In many cases this is fine because the ids are unique by seq region, but in cases like batching it can cause problems
    # So will add in a helper method to make ids unique

    # This holds keys of the current slice details with the gene id to form unique keys. Each time a new key is added
    # the overall gene counter is incremented and the value of the key is set to the new gene id. Any subsequent
    # lines with the same region/gene id key will then just get the new id without incrementing the counter
    gene_id_index = {}
    gene_transcript_id_index = {}
    gene_counter = 1

    # Similar to the gene id index, this will have a key that is based on the slice details, gene id and transcript id. If there
    # is no existing entry, the transcript key will be added and the transcript counter is incremented. If there is a key then
    # the transcript id will be replaced with the new transcript id (which is based on the new gene id and transcript counter)
    # Example key KS8000.rs1.re1000000.gene_1.transcript_1 =
    transcript_id_count_index = {}

    feature_counter = 1

    feature_types = ["exon", "transcript", "repeat", "simple_feature"]
    gtf_output_file_path = output_dir / "annotation.gtf"
    with open(gtf_output_file_path, "w+") as gtf_out:
        for gtf_input_file_path in output_dir.glob("*{extension}"):
            if gtf_input_file_path.stat().st_size == 0:
                logger.info("File is empty, will skip:\n%s" % gtf_input_file_path)
                continue

            gtf_file_name = gtf_input_file_path.name
            match = re.search(r"\.rs(\d+)\.re(\d+)\.", gtf_file_name)
            start_offset = int(match.group(1))
            with open(gtf_input_file_path, "r") as gtf_in:
                for line in gtf_in:
                    values = line.split("\t")
                    if len(values) == 9 and (values[2] in feature_types):
                        values[3] = str(int(values[3]) + (start_offset - 1))
                        values[4] = str(int(values[4]) + (start_offset - 1))
                        if unique_ids:
                            # Maybe make a unique id based on the feature type
                            # Basically region/feature id should be unique at this point, so could use region_id and current_id is key, value is the unique id that is incremented
                            attribs = values[8]

                            # This bit assigns unique gene/transcript ids if the line contains gene_id/transcript_id
                            match_gene_type = re.search(
                                r'(gene_id +"([^"]+)").+(transcript_id +"([^"]+)")',
                                line,
                            )
                            if match_gene_type:
                                full_gene_id_string = match_gene_type.group(1)
                                current_gene_id = match_gene_type.group(2)
                                full_transcript_id_string = match_gene_type.group(3)
                                current_transcript_id = match_gene_type.group(4)
                                gene_id_key = f"{gtf_file_name}.{current_gene_id}"
                                transcript_id_key = (
                                    f"{gene_id_key}.{current_transcript_id}"
                                )
                                if gene_id_key not in gene_id_index:
                                    new_gene_id = f"gene{gene_counter}"
                                    gene_id_index[gene_id_key] = new_gene_id
                                    attribs = re.sub(
                                        full_gene_id_string,
                                        f'gene_id "{new_gene_id}"',
                                        attribs,
                                    )
                                    transcript_id_count_index[gene_id_key] = 1
                                    gene_counter += 1
                                else:
                                    new_gene_id = gene_id_index[gene_id_key]
                                    attribs = re.sub(
                                        full_gene_id_string,
                                        f'gene_id "{new_gene_id}"',
                                        attribs,
                                    )
                                if transcript_id_key not in gene_transcript_id_index:
                                    new_transcript_id = (
                                        gene_id_index[gene_id_key]
                                        + ".t"
                                        + str(transcript_id_count_index[gene_id_key])
                                    )
                                    gene_transcript_id_index[
                                        transcript_id_key
                                    ] = new_transcript_id
                                    attribs = re.sub(
                                        full_transcript_id_string,
                                        f'transcript_id "{new_transcript_id}"',
                                        attribs,
                                    )
                                    transcript_id_count_index[gene_id_key] += 1
                                else:
                                    new_transcript_id = gene_transcript_id_index[
                                        transcript_id_key
                                    ]
                                    attribs = re.sub(
                                        full_transcript_id_string,
                                        f'transcript_id "{new_transcript_id}"',
                                        attribs,
                                    )
                                values[8] = attribs

                            # If you don't match a gene line, try a feature line
                            else:
                                match_feature_type = re.search(
                                    r"(" + feature_id_label + ' +"([^"]+)")', line
                                )
                                if match_feature_type:
                                    full_feature_id_string = match_feature_type.group(1)
                                    current_feature_id = match_feature_type.group(2)
                                    new_feature_id = f"{new_id_prefix}{feature_counter}"
                                    attribs = re.sub(
                                        full_feature_id_string,
                                        f'{feature_id_label} "{new_feature_id}"',
                                        attribs,
                                    )
                                    feature_counter += 1
                                    values[8] = attribs

                        gtf_out.write("\t".join(values))
                    else:
                        logger.info(
                            "Feature type not recognised, will skip: %s" % values[2]
                        )




def convert_gff_to_gtf(gff_file):
    gtf_string = ""
    with open(gff_file) as file_in:
        for line in file_in:
            # match = re.search(r"genBlastG",line)
            # if match:
            results = line.split()
            if not len(results) == 9:
                continue
            if results[2] == "coding_exon":
                results[2] = "exon"
            attributes = set_attributes(results[8], results[2])
            results[8] = attributes
            converted_line = "\t".join(results)
            gtf_string += f"{converted_line}\n"
    return gtf_string


def set_attributes(attributes, feature_type):
    converted_attributes = ""
    split_attributes = attributes.split(";")
    if feature_type == "transcript":
        match = re.search(r"Name\=(.+)$", split_attributes[1])
        name = match.group(1)
        converted_attributes = f'gene_id "{name}"; transcript_id "{name}";'
    elif feature_type == "exon":
        match = re.search(r"\-E(\d+);Parent\=(.+)\-R\d+\-\d+\-", attributes)
        exon_rank = match.group(1)
        name = match.group(2)
        converted_attributes = (
            f'gene_id "{name}"; transcript_id "{name}"; exon_number "{exon_rank}";'
        )

    return converted_attributes


# Example genBlast output
# 1       genBlastG       transcript      131128674       131137049       252.729 -       .       ID=259447-R1-1-A1;Name=259447;PID=84.65;Coverage=94.22;Note=PID:84.65-Cover:94.22
# 1       genBlastG       coding_exon     131137031       131137049       .       -       .       ID=259447-R1-1-A1-E1;Parent=259447-R1-1-A1
# 1       genBlastG       coding_exon     131136260       131136333       .       -       .       ID=259447-R1-1-A1-E2;Parent=259447-R1-1-A1
# 1       genBlastG       coding_exon     131128674       131130245       .       -       .       ID=259447-R1-1-A1-E3;Parent=259447-R1-1-A1
##sequence-region       1_group1        1       4534
# 1       genBlastG       transcript      161503457       161503804       30.94   +       .       ID=259453-R1-1-A1;Name=259453;PID=39.46;Coverage=64.97;Note=PID:39.46-Cover:64.97
# 1       genBlastG       coding_exon     161503457       161503804       .       +       .       ID=259453-R1-1-A1-E1;Parent=259453-R1-1-A1
##sequence-region       5_group1        1       4684
# 5       genBlastG       transcript      69461063        69461741        86.16   +       .       ID=259454-R1-1-A1;Name=259454;PID=82.02;Coverage=91.67;Note=PID:82.02-Cover:91.67
# 5       genBlastG       coding_exon     69461063        69461081        .       +       .       ID=259454-R1-1-A1-E1;Parent=259454-R1-1-A1
# 5       genBlastG       coding_exon     69461131        69461741        .       +       .       ID=259454-R1-1-A1-E2;Parent=259454-R1-1-A1

from  protein_utils {
    run_genblast_align,
    multiprocess_genblast,
    generate_genblast_gtf,
    split_protein_file,
    run_convert2blastmask,
    run_makeblastdb,
    }



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

from transcriptomic_utils inport {
    run_trimming,
    multiprocess_trim_galore,
    create_paired_paths,
    run_star_align,
    run_subsample_script,
    check_for_fastq_subsamples,
    run_minimap2_align,
    check_transcriptomic_output,
    augustus_output_to_gtf,
    run_augustus_predict,
    generate_hints,
    multiprocess_augustus_hints,
    multiprocess_augustus_id,
    create_slice_hints_file,
    run_stringtie_assemble,
    run_scallop_assemble,
    model_builder,
}


def create_slice_ids(
    seq_region_lengths,
    slice_size: int = 1_000_000,
    overlap: int = 0,
    min_length: int = 0,
):
    slice_ids = []
    for region, region_length in seq_region_lengths.items():
        if region_length < min_length:
            continue

        if region_length <= slice_size:
            slice_ids.append([region, 1, region_length])
            continue

        start = 1
        end = start + slice_size - 1
        while end < region_length:
            start = start - overlap
            if start < 1:
                start = 1

            end = start + slice_size - 1
            if end > region_length:
                end = region_length
            if (end - start + 1) >= min_length:
                slice_ids.append([region, start, end])
            start = end + 1

    return slice_ids



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
                # If the start exon coords of both exons are the same, then it's the same exon and thus a single exon cds
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
                    '; biotype "{biotype}";',
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

            # Note that the below looks at validating things under different levels of strictness
            # There are a few different continue statements, where transcripts will be skipped resulting
            # in a smaller post validation file. It mainly removes single coding exon genes with no real
            # support or for multi-exon lncRNAs that are less than 200bp long
            if single_cds_exon_transcript == 1 and validation_type == "relaxed":
                if diamond_e_value is not None:
                    transcript_line = re.sub(
                        '; biotype "{biotype}";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    rnasamba_coding_potential == "coding"
                    and cpc2_coding_potential == "coding"
                    and peptide_length >= min_single_exon_pep_length
                ):
                    transcript_line = re.sub(
                        '; biotype "{biotype}";',
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
                        '; biotype "{biotype}";',
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
                        '; biotype "{biotype}";',
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
                        '; biotype "{biotype}";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                else:
                    continue
            else:
                if diamond_e_value is not None:
                    transcript_line = re.sub(
                        '; biotype "{biotype}";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif (
                    rnasamba_coding_potential == "coding"
                    and cpc2_coding_potential == "coding"
                    and peptide_length >= min_multi_exon_pep_length
                ):
                    transcript_line = re.sub(
                        '; biotype "{biotype}";',
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
                        '; biotype "{biotype}";',
                        '; biotype "protein_coding";',
                        transcript_line,
                    )
                elif transcript_length >= 200:
                    transcript_line = re.sub(
                        '; biotype "{biotype}";', '; biotype "lncRNA";', transcript_line
                    )
                    transcript_line = re.sub(
                        ' translation_coords "[^"]+";', "", transcript_line
                    )
                else:
                    continue

            output_lines.append(transcript_line)
            output_lines.extend(exon_lines)

    return output_lines

from fianlisation_utils import {
    run_finalise_geneset,
    validate_coding_transcripts,
    diamond_validation,
    multiprocess_diamond,
    read_rnasamba_results,
    read_cpc2_results,
    read_diamond_results,
    combine_results,
    merge_finalise_output_files,
}



def read_gtf_genes(gtf_file):
    gtf_genes = {}
    with open(gtf_file) as gtf_in:
        for line in gtf_in:
            eles = line.split("\t")
            if not len(eles) == 9:
                continue

            match = re.search(r'gene_id "([^"]+)".+transcript_id "([^"]+)"', line)

            if not match:
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

    return gtf_genes



def fasta_to_dict(fasta_list):
    index = {}
    it = iter(fasta_list)
    for header in it:
        match = re.search(r">(.+)\n$", header)
        header = match.group(1)
        seq = next(it)
        index[header] = seq
    return index


def subprocess_run_and_log(command):
    logger.info("subprocess_run_and_log command: %s" % " ".join(command))
    subprocess.run(command)


def get_sequence(
    seq_region,
    start: int,
    end: int,
    strand: int,
    fasta_file,
    output_dir: Union[pathlib.Path, str],
):
    start -= 1
    bedtools_path = "bedtools"

    # This creates a tempfile and writes the bed info to it based on whatever information
    # has been passed in about the sequence. Then runs bedtools getfasta. The fasta file
    # should have a faidx. This can be created with the create_faidx static method prior
    # to fetching sequence
    with tempfile.NamedTemporaryFile(
        mode="w+t", delete=False, dir=output_dir
    ) as bed_temp_file:
        bed_temp_file.write(f"{seq_region}\t{start}\t{end}")
        bed_temp_file.close()

    bedtools_command = [
        bedtools_path,
        "getfasta",
        "-fi",
        fasta_file,
        "-bed",
        bed_temp_file.name,
    ]
    bedtools_output = subprocess.Popen(bedtools_command, stdout=subprocess.PIPE)
    for idx, line in enumerate(
        io.TextIOWrapper(bedtools_output.stdout, encoding="utf-8")
    ):
        if idx == 1:
            if strand == 1:
                sequence = line.rstrip()
            else:
                sequence = reverse_complement(line.rstrip())

    os.remove(bed_temp_file.name)
    return sequence


def reverse_complement(sequence: str) -> str:
    rev_matrix = str.maketrans("atgcATGC", "tacgTACG")
    return sequence.translate(rev_matrix)[::-1]


def get_seq_region_names(genome_file: Union[pathlib.Path, str]):
    region_list = []
    with open(genome_file) as file_in:
        for line in file_in:
            match = re.search(r">([^\s]+)", line)
            if match:
                region_name = match.group(1)
                if region_name == "MT":
                    logger.info('Skipping region named "MT"')
                    continue
                else:
                    region_list.append(match.group(1))

    return region_list


def main():
    """
    main function
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output_dir",
        type=str,
        help="Path where the output and temp files will write to. Uses current dir by default",
    )
    parser.add_argument(
        "--genome_file", type=str, required=True, help="Path to the fasta genome file"
    )
    parser.add_argument(
        "--num_threads", type=int, default=1, help="Number of threads to use"
    )
    parser.add_argument(
        "--run_masking",
        action="store_true",
        help="Run Red to find repeats and softmask the genome. Otherwise provide a softmasked genome",
    )
    parser.add_argument(
        "--red_path",
        type=str,
        help="Path to Red executable. See http://toolsmith.ens.utulsa.edu",
    )
    parser.add_argument(
        "--genblast_path",
        type=str,
        help="Path to GenBlast executable. See http://genome.sfu.ca/genblast/download.html",
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
        default=10800,
        help="GenBlast timeout in seconds",
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
        help="Path to a file with Rfam CM accessions, one accession per line, to use with cmsearch",
    )
    parser.add_argument(
        "--run_star", action="store_true", help="Run STAR for short read alignment"
    )
    parser.add_argument(
        "--star_path", type=str, help="Path to STAR for short read alignment"
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
        help="Path to short read fastq dir for running with STAR",
    )
    parser.add_argument(
        "--max_intron_length",
        nargs="?",
        type=int,
        default=100_000,
        help="The maximum intron size for alignments. Default=100,000",
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
        "--run_trnascan", action="store_true", help="Run tRNAscan-SE to find tRNAs"
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
    parser.add_argument(
        "--java_path", type=str, help="Path to Java for use with Eponine"
    )
    parser.add_argument(
        "--run_full_annotation",
        action="store_true",
        help="Run a full annotation, will automatically check for input data and run tools based on that",
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
        help="Run STAR, Stringtie2, Scallop, minimap2 (if short_read_fastq_dir and/or long_read_fastq_dir are provided)",
    )
    parser.add_argument(
        "--run_proteins",
        action="store_true",
        help="Run GenBlast if protein_file and/or busco_protein_file have also been set",
    )
    parser.add_argument(
        "--diamond_validation_db",
        type=str,
        help="Use a Diamond db with blastp mode to help validate cds sequences",
    )
    parser.add_argument(
        "--validation_type",
        type=str,
        help='The strength of evidence needed to validate an ORF as protein coding, can be "relaxed" or "moderate"',
    )
    parser.add_argument(
        "--load_to_ensembl_db",
        # TODO
        # resolve boolean vs string argument type
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
        "--repeatmasker_library", type=str, help="Specify library for repeatmasker"
    )
    parser.add_argument(
        "--repeatmasker_species",
        type=str,
        # TODO
        # does a default value have to been defined here?
        # run_repeatmasker_regions function doesn't assume an existing value for species
        # default="homo",
        help="Specify species for repeatmasker (default homo)",
    )
    args = parser.parse_args()

    work_dir = args.output_dir
    genome_file = args.genome_file
    num_threads = args.num_threads
    masked_genome_file = genome_file  # This will be updated later if Red is run
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

    main_script_dir = pathlib.Path(os.path.realpath(__file__)).parent

    genome_file = pathlib.Path(genome_file)
    if not genome_file.exists():
        raise FileNotFoundError("File does not exist: %s" % genome_file)

    if not work_dir:
        work_dir = pathlib.Path.cwd()
    work_dir = pathlib.Path(work_dir)

    # save log to a file
    log_file_path = work_dir / "ensembl_anno.log"
    add_log_file_handler(logger, log_file_path)

    logger.info("work directory: %s" % work_dir)

    if not work_dir.exists():
        logger.info("Work dir does not exist, creating:")
        create_dir(work_dir)

    if num_threads == 1:
        logger.warning(
            "Thread count is set to the default value 1; this might be slow."
        )

    # If the run_full_annotation flag is set then we want to set a standardised set of analyses
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
    seq_region_names = get_seq_region_names(genome_file)
    logger.debug("seq region names:\n%s" % seq_region_names)

    # Repeat analyses
    ############################################################################
    if run_masking:
        logger.info("Running masking via Red")
        masked_genome_file = run_red(red_path, genome_file, work_dir)
        logger.info("Masked genome file: %s" % masked_genome_file)
    else:
        logger.info("Not running masking, presuming the genome file is softmasked")

    if run_dust:
        logger.info("Annotating low complexity regions")
        run_dust_regions(genome_file, dust_path, work_dir, num_threads)

    if run_trf:
        logger.info("Annotating tandem repeats")
        run_trf_repeats(genome_file, trf_path, work_dir, num_threads)

    if run_repeatmasker:
        logger.info("Annotating repeats with RepeatMasker")
        run_repeatmasker_regions(
            genome_file, repeatmasker_path, library, species, work_dir, num_threads
        )
    ############################################################################

    # Simple feature analyses
    ############################################################################
    if run_cpg:
        logger.info("Annotating CpG islands")
        run_cpg_regions(genome_file, cpg_path, work_dir, num_threads)

    if run_eponine:
        logger.info("Running Eponine to find transcription start sites")
        run_eponine_regions(genome_file, java_path, eponine_path, work_dir, num_threads)
    ############################################################################

    # sncRNA analyses
    ############################################################################
    # Search Rfam with cmsearch
    if run_cmsearch:
        logger.info("Annotating sncRNAs")
        run_cmsearch_regions(
            genome_file, None, None, None, rfam_accessions_file, work_dir, num_threads
        )

    if run_trnascan:
        logger.info("Annotating tRNAs")
        run_trnascan_regions(
            genome_file, trnascan_path, trnascan_filter_path, work_dir, num_threads
        )
    ############################################################################

    # Transcriptomic analyses
    ############################################################################
    if trim_fastq:
        run_trimming(work_dir, short_read_fastq_dir, delete_pre_trim_fastq, num_threads)

    # Run STAR
    if run_star:
        logger.info("Running STAR")
        run_star_align(
            star_path=star_path,
            trim_fastq=trim_fastq,
            subsample_script_path=subsample_script_path,
            main_output_dir=work_dir,
            short_read_fastq_dir=short_read_fastq_dir,
            genome_file=genome_file,
            max_reads_per_sample=max_reads_per_sample,
            max_total_reads=max_total_reads,
            max_intron_length=max_intron_length,
            num_threads=num_threads,
            main_script_dir=main_script_dir,
        )

    # Run Scallop
    if run_scallop:
        logger.info("Running Scallop")
        run_scallop_assemble(scallop_path, stringtie_path, work_dir)

    # Run Stringtie
    if run_stringtie:
        logger.info("Running StringTie")
        run_stringtie_assemble(stringtie_path, samtools_path, work_dir, num_threads)

    # Run minimap2
    if run_minimap2:
        logger.info("Running minimap2")
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
    ############################################################################

    # Protein analyses
    ############################################################################
    # Run GenBlast
    if run_genblast:
        logger.info("Running GenBlast")
        run_genblast_align(
            genblast_path,
            convert2blastmask_path,
            makeblastdb_path,
            work_dir / "genblast_output",
            protein_file,
            masked_genome_file,
            max_intron_length,
            num_threads,
            genblast_timeout,
            main_script_dir,
        )

    # Run GenBlast on BUSCO set, gives higher priority when creating the final genes in cases where transcriptomic data are missing or fragmented
    if run_busco:
        logger.info("Running GenBlast of BUSCO proteins")
        run_genblast_align(
            genblast_path,
            convert2blastmask_path,
            makeblastdb_path,
            work_dir / "busco_output",
            busco_protein_file,
            masked_genome_file,
            max_intron_length,
            num_threads,
            genblast_timeout,
            main_script_dir,
        )
    ############################################################################

    # Finalisation analyses
    ############################################################################
    # Do some magic
    if finalise_geneset:
        logger.info("Finalise geneset")
        run_finalise_geneset(
            main_script_dir,
            work_dir,
            genome_file,
            seq_region_names,
            validation_type,
            diamond_validation_db,
            num_threads,
        )
    ############################################################################

    # Other analyses
    ############################################################################
    # Run Augustus
    if run_augustus:
        logger.info("Running Augustus")
        run_augustus_predict(augustus_path, work_dir, masked_genome_file, num_threads)
    ############################################################################

    if load_to_ensembl_db:
        load_results_to_ensembl_db(
            main_script_dir,
            load_to_ensembl_db,
            genome_file,
            work_dir,
            db_details,
            num_threads,
        )


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.info("Interrupted with CTRL-C, exiting...")
