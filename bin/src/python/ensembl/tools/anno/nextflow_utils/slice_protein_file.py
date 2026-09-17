import argparse
import random
import re
from pathlib import Path
from typing import List
import os



def split_protein_file(protein_dataset: Path,  protein_source,  batch_size: int = 200) -> List:
    """
    The protein dataset file is split by a number of sequence
    equals to the batch_size
    in batch files stored in 10 output directories.
    protein_dataset : Path for the protein dataset.
    protein_source : Output directory path.
    batch_size : Size of the batch, it needs to be equals to the
    number of threads
    to parallelise the sequence processing for each file.
    """
    batched_protein_files = []
    os.mkdir(protein_source)

    for i in range(0, 10):
        os.mkdir(f"{protein_source}/bin_{i}")
    with open(protein_dataset, "r", encoding="utf8") as file_in:
        seq_count = 0
        batch_count = 0
        current_record = ""
        initial_seq = True
        for line in file_in:
            match = re.search(r">(.+)$", line)
            # match header and is not first sequence, if the number
            # of stored sequences in each file equals
            # the number of batch_size, a new file will be created
            # and the current_record reset
            if match and not initial_seq and seq_count % batch_size == 0:
                bin_num = random.randint(0, 9)
                batch_file = f"{protein_source}/bin_{str(bin_num)}/{str(batch_count)}.fa"
                with open(batch_file, "w+") as file_out:
                    file_out.write(current_record)
                batch_count += 1
                seq_count += 1
                current_record = line
                batched_protein_files.append(batch_file)
            # match header and is the first sequence
            elif match:
                current_record += line
                initial_seq = False
                seq_count += 1
            # other lines
            else:
                current_record += line

        if current_record:
            bin_num = random.randint(0, 9)
            batch_file = f"{protein_source}/bin_{str(bin_num)}/{str(batch_count)}.fa"
            with open(batch_file, "w+") as file_out:
                file_out.write(current_record)
            batched_protein_files.append(batch_file)
    return batched_protein_files

def parse_args():
    parser = argparse.ArgumentParser(description="Arguments for script to check contents of transcriptomic gtfs")
    parser.add_argument("--proteins", help="Path to input protein file")
    parser.add_argument("--protein_source", help="name of protein source")
    args = parser.parse_args()
    return args
    

if __name__ == "__main__":
    args = parse_args()
    split_protein_file(protein_dataset=args.proteins, 
                       protein_source=args.protein_source)
