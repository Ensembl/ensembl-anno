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
"""
Set of functions used to collect the results of the subpipelines and generate the final output
This logic will be revised in the nextflow pipeline
"""
__all__ = ["run_finalise_geneset"]

import logging
import logging.config
import multiprocessing
import os
from pathlib import Path
import re
import subprocess
from typing import List
import shutil

from ensembl.tools.anno.utils._utils import (
    check_file,
    get_seq_region_length,
    fasta_to_dict,
    check_exe,
    create_dir,
    check_gtf_content,
    split_protein_file
)

logger = logging.getLogger(__name__)

def run_finalise_geneset(
    main_script_dir: Path,
    main_output_dir: Path,
    genome_file: Path,
    seq_region_names: List[str],
    diamond_validation_db :Path,
    validation_type: str = "relaxed",
    num_threads: int = 1,
    cpc2_bin: Path = Path("/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/test_cpc2.sif"),
    rnasamba_bin: Path = Path("/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/rnasamba_latest.sif"),
    rnasamba_weights: Path= Path("/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/rnasamba_data/full_length_weights.hdf5"),

):
    """Collect results 

        :param main_script_dir:
        :type main_script_dir: Path
        :param main_output_dir: _description_
        :type main_output_dir: Path
        :param genome_file: Path for genome fasta file
        :type genome_file: Path
        :param seq_region_names: list of seq region names
        :type seq_region_names: List[str]
        :param diamond_validation_db: Path for diamond validation db
        :type diamond_validation_db: Path
        :param validation_type: type of validation
        :type validation_type: str, default "relaxed"
        :param num_threads: Num of threads. Defaults to 1.
        :type num_threads: int, default 1 
        :param cpc2_bin: CPC2 software path
        :type cpc2_bin: Path, default "/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/test_cpc2.sif"
        :param rnasamba_bin: RNASamba software path
        :type rnasamba_bin: Path, default "/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/rnasamba_latest.sif"
        :param rnasamba_weights: RNASamba weights path
        :type rnasamba_weights: Path, default "/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/rnasamba_data/full_length_weights.hdf5"

    """


    final_annotation_dir = create_dir(main_output_dir, "annotation_output")
    region_annotation_dir = create_dir(final_annotation_dir, "initial_region_gtfs")
    final_region_annotation_dir = create_dir(
        final_annotation_dir, "final_region_gtfs"
    )
    utr_region_annotation_dir = create_dir(final_annotation_dir, "utr_region_gtfs")
    validation_dir = create_dir(final_annotation_dir, "cds_validation")
    seq_region_lengths = get_seq_region_length(genome_file, 0)
    output_file = final_annotation_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Final gtf file exists, skipping analysis")
            return

    protein_annotation_raw = main_output_dir / "genblast_output" / "annotation.gtf"
    minimap2_annotation_raw = main_output_dir / "minimap2_output" / "annotation.gtf"
    stringtie_annotation_raw = main_output_dir / "stringtie_output" / "annotation.gtf"
    scallop_annotation_raw = main_output_dir / "scallop_output" / "annotation.gtf"
    busco_annotation_raw = main_output_dir / "busco_output" / "annotation.gtf"
    transcript_selector_script = main_output_dir / "support_scripts_perl" / "select_best_transcripts.pl"
    finalise_geneset_script = main_output_dir / "support_scripts_perl" / "finalise_geneset.pl"
    clean_geneset_script = main_output_dir / "support_scripts_perl" / "clean_geneset.pl"
    clean_utrs_script = main_output_dir / "support_scripts_perl" / "clean_utrs_and_lncRNAs.pl"
    gtf_to_seq_script = main_output_dir / "support_scripts_perl" / "gtf_to_seq.pl"

    transcriptomic_annotation_raw = main_output_dir / "transcriptomic_raw.gtf"

    with transcriptomic_annotation_raw.open("w") as file_out:
        for transcriptomic_file in [
        minimap2_annotation_raw,
        scallop_annotation_raw,
        stringtie_annotation_raw,
        ]:
            if not transcriptomic_file.exists():
                logger.info("No annotation.gtf file found in %s, skipping", transcriptomic_file)
                continue

            with transcriptomic_file.open("r") as file_in:
                for line in file_in:
                    file_out.write(line)

    # Copy the raw files into the annotation dir
    copy_raw_files(busco_annotation_raw, final_annotation_dir / "busco_raw.gtf")
    copy_raw_files(protein_annotation_raw, final_annotation_dir / "protein_raw.gtf")

    #select best transcript
    generic_select_cmd = [
        "perl",
        str(transcript_selector_script),
        "-genome_file",
        str(genome_file),
    ]
    pool = multiprocessing.Pool(int(num_threads)) # pylint: disable=consider-using-with
    for seq_region_name in seq_region_names:
        # The selection script needs different params depending on
        # whether the seqs are from transcriptomic data or not
        region_details = f"{seq_region_name}.rs1.re{seq_region_lengths[seq_region_name]}"
        transcriptomic_region_gtf_path = region_annotation_dir / f"{region_details}.trans.gtf"
        busco_region_gtf_path = region_annotation_dir / f"{region_details}.busco.gtf"
        protein_region_gtf_path = region_annotation_dir / f"{region_details}.protein.gtf"
        if transcriptomic_annotation_raw:
            cmd = generic_select_cmd.copy()
            cmd.extend(
                [
                    "-region_details",
                    region_details,
                    "-input_gtf_file",
                    str(transcriptomic_annotation_raw),
                    "-output_gtf_file",
                    str(transcriptomic_region_gtf_path),
                    "-cds_search",
                    "-final_biotype",
                    "transcriptomic",
                ]
            )
            pool.apply_async(run_command, args=(cmd,))
        if busco_annotation_raw:
            cmd = generic_select_cmd.copy()
            cmd.extend(
                [
                    "-region_details",
                    region_details,
                    "-input_gtf_file",
                    str(busco_annotation_raw),
                    "-output_gtf_file",
                    str(busco_region_gtf_path),
                    "-all_cds_exons",
                    "-final_biotype",
                    "busco",
                ]
            )
            pool.apply_async(run_command, args=(cmd,))
        if protein_annotation_raw:
            cmd = generic_select_cmd.copy()
            cmd.extend(
                [
                    "-region_details",
                    region_details,
                    "-input_gtf_file",
                    str(protein_annotation_raw),
                    "-output_gtf_file",
                    str(protein_region_gtf_path),
                    "-clean_transcripts",
                    "-all_cds_exons",
                    "-final_biotype",
                    "protein",
                ]
            )
            pool.apply_async(run_command, args=(cmd,))

    pool.close()
    pool.join()

    # At this point we will have the region files for each seq region and we can merge them
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
    fully_merged_gtf_path = final_annotation_dir / "all_selected_transcripts.gtf"
    # Collect all *_sel.gtf files
    gtf_files = list(final_annotation_dir.glob("*_sel.gtf"))

    # Ensure there are files to merge
    if not gtf_files:
        raise FileNotFoundError("No *_sel.gtf files found in %s",final_annotation_dir)


    try:
        with fully_merged_gtf_path.open("w+") as fully_merged_gtf_out:
            merge_gtf_cmd = ["cat"] + [str(gtf_file) for gtf_file in gtf_files]
            subprocess.run(merge_gtf_cmd, stdout=fully_merged_gtf_out, check=True)
            logger.info("Merged GTF files into %s", fully_merged_gtf_path)
    except Exception as e:
        print("An error occurred while merging GTF files: %s",e)
        
    # Now collapse the gene set
    generic_finalise_cmd = [
        "perl",
        str(finalise_geneset_script),
        "-genome_file",
        str(genome_file),
    ]

    pool = multiprocessing.Pool(int(num_threads)) # pylint:disable=consider-using-with
    for seq_region_name in seq_region_names:
        region_details = f"{seq_region_name}.rs1.re{seq_region_lengths[seq_region_name]}"
        final_region_gtf_path = final_region_annotation_dir / f"{region_details}.final.gtf"
        cmd = generic_finalise_cmd.copy()
        cmd.extend(
            [
                "-region_details",
                region_details,
                "-input_gtf_file",
                str(fully_merged_gtf_path),
                "-output_gtf_file",
                str(final_region_gtf_path),
            ]
        )
        pool.apply_async(run_command, args=(cmd,))

    pool.close()
    pool.join()

    merge_finalise_output_files(
        final_annotation_dir,
        final_region_annotation_dir,
        ".final.gtf",
        "prevalidation",
    )
    merged_gtf_file = final_annotation_dir / "prevalidation_sel.gtf"
    merged_cdna_file = final_annotation_dir / "prevalidation_sel.cdna.fa"
    merged_amino_acid_file = final_annotation_dir /"prevalidation_sel.prot.fa"
    #validation
    updated_gtf_lines = validate_coding_transcripts(
        merged_cdna_file,
        merged_amino_acid_file,
        validation_dir,
        validation_type,
        diamond_validation_db,
        merged_gtf_file,
        num_threads,
        cpc2_bin,
        rnasamba_bin,
        rnasamba_weights
    )
    
    ####
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
    cdna_file:Path,
    amino_acid_file:Path,
    validation_dir:Path,
    validation_type:str,
    diamond_validation_db:Path,
    gtf_file:Path,
    num_threads:int,
    cpc2_bin: Path = Path("/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/test_cpc2.sif"),
    rnasamba_bin: Path = Path("/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/rnasamba_latest.sif"),
    rnasamba_weights: Path= Path("/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/rnasamba_data/full_length_weights.hdf5"),

):
    """Run validation on protein coding transcript 
    using CPC2, RNASamba, Diamond

    Args:
        cdna_file (Path): Input cdna file
        amino_acid_file (Path): Input protein file
        validation_dir (Path): Validation output directory
        validation_type (str): _description_
        diamond_validation_db (Path): _description_
        gtf_file (Path): _description_
        num_threads (int): number of threads

    Returns:
        _type_: _description_
    """

    logger.info("Running CDS validation with RNAsamba and CPC2")
    #rnasamba_weights = config["rnasamba"]["weights"]

    

    
    

   

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
            # of sequences that don't make it into the cpc2/rnasamba oqutput
            # Should code in a system for this, but it would be good to
            # understand why it happens to begin with. Seems to be the same
            # number of missing seqs in both, so maybe a shared cut-off
            if transcript_id in transcript_ids:
                transcript_ids[transcript_id].extend([e_value])

    return transcript_ids

def merge_finalise_output_files(
    final_annotation_dir: Path, region_annotation_dir: Path, extension: str, id_label: str
)-> None:
    """Merge the resulting gene models from all the seq regions for the specified extension

    Args:
        final_annotation_dir (Path): Output directory
        region_annotation_dir (Path): Input directory
        extension (str): File extension
        id_label (str): Type of data or process to merge
    """
    gtf_files = list(region_annotation_dir.glob(f"*{extension}"))
    merged_gtf_file = final_annotation_dir / f"{id_label}_sel.gtf"
    merged_cdna_file = final_annotation_dir / f"{id_label}_sel.cdna.fa"
    merged_amino_acid_file = final_annotation_dir / f"{id_label}_sel.prot.fa"

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
    
    with merged_gtf_file.open("w") as gtf_out, \
        merged_cdna_file.open("w") as cdna_out, \
        merged_amino_acid_file.open("w") as amino_acid_out:

        for gtf_file in gtf_files:
            logger.info("GTF file: %s", gtf_file)
            cdna_seq_index = {}
            amino_acid_seq_index = {}
            cdna_file = gtf_file.with_suffix(".cdna")
            amino_acid_file = gtf_file.with_suffix(".prot")
            if cdna_file.exists():
                with cdna_file.open() as cdna_in:
                    cdna_seq_index = fasta_to_dict(cdna_in.readlines())
            
            if amino_acid_file.exists():
                with amino_acid_file.open() as amino_acid_in:
                    amino_acid_seq_index = fasta_to_dict(amino_acid_in.readlines())
            current_gene_id = ""
            with gtf_file.open() as gtf_in:
                for line in gtf_in:
                    if re.search(r"^#", line):
                    #if line.startswith("#"):
                        continue
                    eles = line.split("\t")
                    if len(eles) != 9:
                        continue

                    match = re.search(r'gene_id "([^"]+)".+transcript_id "([^"]+)"', line)
                    if match and eles[2] == "transcript":
                        transcript_id_counter += 1
                    gene_id = match.group(1)
                    transcript_id = match.group(2)

                    if not current_gene_id:
                        gene_id_counter += 1
                        current_gene_id = gene_id

                    if gene_id != current_gene_id:
                        gene_id_counter += 1
                        current_gene_id = gene_id
                    new_gene_id = f"{id_label}_{gene_id_counter}"
                    new_transcript_id = f"{id_label}_{transcript_id_counter}"

                    line = re.sub(f'gene_id "{gene_id}"', f'gene_id "{new_gene_id}"', line)
                    line = re.sub(f'transcript_id "{transcript_id}"', f'transcript_id "{new_transcript_id}"', line)
                    gtf_out.write(line)
                    if eles[2] == "transcript":
                        new_header = f">{new_transcript_id}\n"
                        cdna_out.write(new_header + cdna_seq_index[transcript_id])
                        if transcript_id in amino_acid_seq_index:
                            amino_acid_out.write(new_header + amino_acid_seq_index[transcript_id])


def multiprocess_finalise_geneset(cmd):

    print(" ".join(cmd))
    subprocess.run(cmd)
    


def run_command(cmd)->None:
    """Submit command in multiprocessing action 

    Args:
        cmd (List[str]): command to execute
    """
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info("Command executed successfully: %s", cmd)
        logger.info("Output: %s", result.stdout)
    except subprocess.CalledProcessError as e:
        logger.error("Command failed: %s", cmd)
        logger.error("Error: %s", e.stderr)

def copy_raw_files(raw_file: Path, destination_file :Path)-> None:
    """Copy file in a new destination

    Args:
        raw_file (Path): File to copy
        destination_file (Path): New destination
    """
    if raw_file.exists():
        try:
            shutil.copy(raw_file, destination_file)
            logger.info("Copied %s to %s", raw_file, destination_file)
        except Exception as e:
            logger.error("Failed to copy %s to %s: %s", raw_file, destination_file, e)
    else:
        logger.info("No file found at %s, skipping", raw_file)



    
    
    
    
    
    
    
    
    
    
    
    
    
    

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Genblast arguments")
    parser.add_argument(
        "--masked_genome_file", required=True, help="Masked genome file path"
    )
    parser.add_argument("--output_dir", required=True, help="Output directory path")
    parser.add_argument("--protein_file", required=True, help="Path for the protein dataset")
    parser.add_argument(
        "--genblast_timeout_secs", type=int, default=10800, help="Genblast timeout period"
    )
    parser.add_argument(
        "--max_intron_length", type=int, required=True, help="Maximum intron length"
    )
    parser.add_argument(
        "--genblast_bin",
        default="genblast",
        help="Genblast executable path",
    )
    parser.add_argument(
        "--convert2blastmask_bin",
        default="convert2blastmask",
        help="convert2blastmask executable path",
    )
    parser.add_argument(
        "--makeblastdb_bin",
        default="makeblastdb",
        help="makeblastdb executable path",
    )
    parser.add_argument("--num_threads", type=int, default=1, help="Number of threads")
    parser.add_argument(
        "--protein_set",
        required=True,
        choices=["uniprot", "orthodb"],
        help="Protein set [uniprot, orthodb]",
    )
    return parser.parse_args()

def main():
    """Genblast's entry-point."""
    args = parse_args()

    log_file_path = create_dir(args.output_dir, "log") / "genblast.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"

    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )

    run_genblast(
        Path(args.masked_genome_file),
        Path(args.output_dir),
        Path(args.protein_file),
        args.max_intron_length,
        args.genblast_timeout_secs,
        Path(args.genblast_bin),
        Path(args.convert2blastmask_bin),
        Path(args.makeblastdb_bin),
        args.num_threads,
        args.protein_set,
    )

if __name__ == "__main__":
    main()
