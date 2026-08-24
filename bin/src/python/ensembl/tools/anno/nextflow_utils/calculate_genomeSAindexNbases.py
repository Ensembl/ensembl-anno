import argparse
import re
import math
from os import PathLike
from typing import Dict

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


def parse_args():
    parser = argparse.ArgumentParser(description="Arguments for script to calculate index base value to pass to STAR in genomeGenerate mode")
    parser.add_argument("--genome_file", help="Genome file path")
    parser.add_argument("--min_seq_length", required=False, default = 0, help="maximum chromosome or slice size to consider when calculating genome size")
    parser.add_argument("--maximum_index_bases", required=False, default=14, help="maximum index bases value")
    args = parser.parse_args()
    return args
    
if __name__ == "__main__":
    args = parse_args()
    choose_index_length(genome_file=args.genome_file, 
                        min_seq_length=args.min_seq_length, 
                        maximum_index_bases=args.maximum_index_bases)
    
