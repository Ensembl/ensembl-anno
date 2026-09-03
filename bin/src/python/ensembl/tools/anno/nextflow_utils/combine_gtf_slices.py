import argparse
import os
from typing import List
from pathlib import Path
import re

def slice_output_to_gtf(  # pylint: disable=too-many-branches, too-many-statements, too-many-locals
    output_gtf: Path,
    sliced_gtf_list: List,
    feature_id_label: str = "",
    unique_ids: bool = True,
    tool: str = ""
) -> None:
    """
    Collect all the gtf files per file extension and
    merge them in a final gtf file  assigning unique ids.

    This holds keys of the current slice details with
    the gene id to form unique keys.
    Each time a new key is added the overall gene
    counter is incremented
    and the value of the key is set to the new gene id.
    Any subsequent
    lines with the same region/gene id key will then just get
    the new id without incrementing the counter.

    Args:
    output_dir : Output directory.
    feature_id_label : Feature identifier.
    new_id_prefix : New feature identifier.
    unique_ids : If True assign unique ids for the same
    feature type.
    """
    feature_types = ["exon", "transcript", "repeat", "simple_feature"]
    new_id_prefix = ""
    if tool == "repeatmasker":
        print("new_id_prefix set")
        new_id_prefix == "repeatmask"
    # Initialise gene and feature counter
    gene_counter = 1
    feature_counter = 1
    # Initialise dictionaries that will store the list of
    # gene and transcript indexes
    gene_id_collection = {}
    gene_transcript_id_collection = {}
    transcript_id_count_gene_id = {}  # one gene can have multiple transcripts
    with open(output_gtf, "w+", encoding="utf8") as output_file:
        for input_file in sliced_gtf_list:
            if os.stat(input_file).st_size == 0:
                print("File is empty, will skip %s", input_file)
                continue

            # I've replaced the original search logic here with logic to strip the start coordinates out of the 
            # filename. This is pretty dirty - I've done this as a temporary measure. In the long term maybe we
            # could look at gtf parsing and operations across all of our repositories and build something more
            # robust.

            # This is what was originally here. It doesn't work because the filenames are now slightly different
            # match = re.search(r"\.rs(\d+)\.re(\d+)\.", input_file)
            # assert match is not None
            # start_offset = int(match.group(1))

            # Here is my temporary fix - we should plan to replace this!!!
            # Under the new naming scheme the filename is chr:start-end'.1.bed'. Therefore...
            start_offset = input_file.split(':')[1].split("-")[0]
            try: 
                start_offset = int(start_offset)
            except:
                raise ValueError(f"Filenames are not as expected - the third field in the filename should be the " \
                "slice start coordinate but instead it is {start_offset}. The filename was {input_file}. This shouldn't " \
                "happen, exiting.")
            
            with open(input_file, "r", encoding="utf8") as gtf_in:
                for line in gtf_in:
                    values = line.split("\t")
                    print(str(len(values)))
                    print(str((values[2] in feature_types)))
                    if len(values) == 9 and (values[2] in feature_types) and unique_ids:
                        # each slice start from 0 so we need to add the
                        # offset to get the real coordinates
                        values[3] = str(int(values[3]) + (start_offset - 1))
                        values[4] = str(int(values[4]) + (start_offset - 1))
                        # Unique id based on gene and transcript ids
                        attribs = values[8]
                        # Get the slice details, gene id and transcript id.
                        match_gene_type = re.search(
                            r'(gene_id +"([^"]+)").+(transcript_id +"([^"]+)")',
                            line,
                        )
                        # gene_id "1"; transcript_id "1"; biotype "tRNA_pseudogene";
                        if match_gene_type:
                            full_gene_id_string = match_gene_type.group(1)
                            current_gene_id = match_gene_type.group(2)
                            full_transcript_id_string = match_gene_type.group(3)
                            current_transcript_id = match_gene_type.group(4)
                            # Example key KS8000.rs1.re1000000.1
                            gene_id_slice = input_file.name + "." + str(current_gene_id)
                            # Example key KS8000.rs1.re1000000.1.transcript.1
                            transcript_id_slice = gene_id_slice + "." + str(current_transcript_id)
                            # If there is no existing entry, the gene key will be added
                            # and the gene counter is incremented.
                            # gene_id "gene1"; transcript_id "gene1.t1"
                            if gene_id_slice not in gene_id_collection:
                                new_gene_id = f"gene{gene_counter}"
                                gene_id_collection[gene_id_slice] = new_gene_id
                                transcript_id_count_gene_id[gene_id_slice] = 1
                                gene_counter += 1
                            else:
                                # If there is a key then the gene
                                # id will be replaced with the new gene id
                                new_gene_id = gene_id_collection[gene_id_slice]
                            attribs = re.sub(
                                full_gene_id_string,
                                'gene_id "' + new_gene_id + '"',
                                attribs,
                            )
                            # If there is no existing entry, the transcript key will be added
                            # and the transcript counter is incremented.
                            if transcript_id_slice not in gene_transcript_id_collection:
                                new_transcript_id = (
                                    gene_id_collection[gene_id_slice]
                                    + ".t"
                                    + str(
                                        transcript_id_count_gene_id[gene_id_slice]
                                    )  # pylint:disable=line-too-long
                                )
                                gene_transcript_id_collection[transcript_id_slice] = (
                                    new_transcript_id  # pylint:disable=line-too-long
                                )
                                transcript_id_count_gene_id[gene_id_slice] += 1
                            else:
                                # If a transcript of the same set is already present,
                                # the new id with the incremented counter is added
                                new_transcript_id = gene_transcript_id_collection[
                                    transcript_id_slice
                                ]  # pylint:disable=line-too-long
                            attribs = re.sub(
                                full_transcript_id_string,
                                'transcript_id "' + new_transcript_id + '"',
                                attribs,
                            )
                            values[8] = attribs
                            output_file.write("\t".join(values))
                            # Unique id based on the feature type
                        else:
                            if new_id_prefix == "repeatmask":
                                match_feature_type = re.search(
                                    r"("
                                    + feature_id_label
                                    + "\d+)",  # pylint:disable=line-too-long, anomalous-backslash-in-string
                                    line,
                                )
                            else:
                                match_feature_type = re.search(
                                    r"("
                                    + feature_id_label
                                    + ' +"([^"]+)")',  # "\d+)",# pylint:disable=line-too-long
                                    line,
                                )
                            if match_feature_type:
                                full_feature_id = match_feature_type.group(1)
                                new_feature_id = new_id_prefix + str(feature_counter)
                                attribs = re.sub(
                                    full_feature_id,
                                    feature_id_label + ' "' + new_feature_id + '"',
                                    attribs,
                                )
                                feature_counter += 1
                                values[8] = attribs
                                output_file.write("\t".join(values))
                    else:
                        print(
                            "Feature type not recognised, will skip. Feature type: %s",
                            values[2],
                        )

def parse_args():
    parser = argparse.ArgumentParser(description="Arguments for script to check contents of transcriptomic gtfs")
    parser.add_argument("--sliced_gtfs", nargs = '+', help="Path to input to convert to gtf")
    parser.add_argument("--output_gtf", help="Path to output logfile recording status of each gtf")
    parser.add_argument("--tool")
    args = parser.parse_args()
    return args
    

if __name__ == "__main__":
    args = parse_args()
    slice_output_to_gtf(output_gtf = args.output_gtf, 
                        sliced_gtf_list=args.sliced_gtfs,
                        tool = args.tool)