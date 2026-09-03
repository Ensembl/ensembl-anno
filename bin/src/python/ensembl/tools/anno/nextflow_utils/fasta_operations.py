import argparse
import re
import math
from os import PathLike
from typing import Dict, List

def get_seq_region_length(genome_file: PathLike, min_seq_length: int = 0) -> Dict:
    """
    Split the genome file according to the header and store in a dictionary
    all the sequences whose length is greater than min_seq_length.
    Args:
        genome_file: Genome file path.
        min_seq_length: maximum slice length.
    Return: Dictionary of sequence headers and the corresponding sequence length
    """
    current_header = ""
    current_seq = ""

    seq_region_to_length = {}
    with open(genome_file, "r", encoding="utf8") as file_in:
        for line in file_in:
            match = re.search(r">(.+)$", line)
            if match and current_header:
                if len(current_seq) > min_seq_length:
                    seq_region_to_length[current_header] = len(current_seq)

                current_seq = ""
                current_header = match.group(1)
            elif match:
                current_header = match.group(1)
            else:
                current_seq += line.rstrip()

        if len(current_seq) > min_seq_length:
            seq_region_to_length[current_header] = len(current_seq)
    return seq_region_to_length


def choose_index_length(genome_file, min_seq_length=0, maximum_index_bases=14):
    seq_region_to_length = get_seq_region_length(genome_file, min_seq_length)
    genome_size = sum(seq_region_to_length.values())

    # This calculates the base-2 logarithm of the genome_size. The logarithm of the genome size is
    # a measure of how many bits are needed to represent the genome size in binary.
    #
    # The choice of 14 as the default maximum value for maximum_index_bases is likely based on 
    # empirical observations and optimization considerations. Too large of a seed length can 
    # lead to increased memory usage and potentially slower indexing, while a seed length that 
    # is too small might affect alignment accuracy.
    index_bases = min(maximum_index_bases, math.floor((math.log(genome_size, 2) / 2) - 1))
    # print to stdout - this is what is picked up by Nextflow
    print(index_bases) 
    return index_bases


def get_slice_id(
    seq_region_to_length: Dict,
    slice_size: int = 1000000,
    overlap: int = 0,
    min_length: int = 0,
) -> List:
    """
    Get list of ids for a genomic slice
    Arg:
    seq_region_to_length: Dictionary with the sequence headers
    as keys and the sequence lengths as values
    slice_size: Size of the slice
    overlap: Overlap length between two slices
    min_length: Min length of the slice
    Return: List of IDs for each genomic slice

    Future TODO (not in nextflowification PR) - clean up start coordinate logic to return a zero based coordinate directly
    for bedfiles. I think the current logic is correct but it's unneccessarily complicated.
    """

    slice_ids_per_region = []
    for seq_region in seq_region_to_length:
        seq_region_length = int(seq_region_to_length[seq_region])
        if seq_region_length < min_length:
            continue

        if seq_region_length <= slice_size:
            slice_ids_per_region.append([seq_region, 1, seq_region_length])
            continue

        start = 1
        end = start + slice_size - 1
        while end < seq_region_length:
            start = start - overlap
            start = max(start, 1)
            end = start + slice_size - 1
            end = min(end, seq_region_length)
            if (end - start + 1) >= min_length:
                slice_ids_per_region.append([seq_region, start, end])
            start = end + 1
    return slice_ids_per_region

def write_fasta_intervals_to_bed(args):

    seq_region_to_length_dict = get_seq_region_length(args.genome_file, 
                                                      args.min_seq_length)
    slice_ids = get_slice_id(seq_region_to_length_dict, 
                             slice_size= args.slice_size,
                             overlap = args.overlap,
                             min_length = args.min_seq_length) # nb. pass different value on cmdline

    for slice_id in slice_ids:
        seq_region = slice_id[0]
        start = slice_id[1]
        end = slice_id[2]
    
        start -= 1 

        bed_file = f"{seq_region}:{start}-{end}.1.bed"
        with open(bed_file, 'w') as bed_handle:
            bed_handle.write(f"{seq_region}\t{start}\t{end}\n")



def parse_args():
    parser = argparse.ArgumentParser(description="Arguments for script to calculate index base value to pass to STAR in genomeGenerate mode")
    parser.add_argument("--genome_file", help="Genome file path")
    parser.add_argument("--calculate_genomeSAindexNbases", action='store_true', help="calculate star index param")
    parser.add_argument("--splitFasta", action='store_true')
    parser.add_argument("--min_seq_length", required=False, type=int, default = 0, help="maximum chromosome or slice size to consider when calculating genome size")
    parser.add_argument("--slice_size", required=False, type=int, default=1000000)
    parser.add_argument("--overlap", required=False, type=int, default=0)
    parser.add_argument("--maximum_index_bases", required=False, type=int, default=14, help="maximum index bases value")
    args = parser.parse_args()
    return args
    
if __name__ == "__main__":
    args = parse_args()
    if args.calculate_genomeSAindexNbases:
        choose_index_length(genome_file=args.genome_file, 
                            min_seq_length=args.min_seq_length, 
                            maximum_index_bases=args.maximum_index_bases)

    if args.splitFasta:
        write_fasta_intervals_to_bed(args)
    
