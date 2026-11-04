import argparse
import os
import pandas as pd


def process_current_record(current_record, imgt_gca_clade, output_dir):

    # Throw out short sequences (10 amino acids = 30 nucleotides)
    # 10 amino acids is an arbitrary number
    name = current_record['name']
    fasta_header = current_record['fasta_header']
    sequence = current_record['sequence']
    if len(sequence) < 10:
        return None

    clade = imgt_gca_clade[imgt_gca_clade['Original-labels'] == name]['Clade'].to_list()
    if len(clade) != 1:
        if len(clade) ==0:
            return None
        else:
            raise ValueError(f'{str(len(clade))} entries found for {name}, expecting 1 entry')

    lines_to_write = f'{fasta_header}{sequence}\n'
    file_to_write_to = f'{clade[0]}_imgt_reference_sequences.fa'
    with open(os.path.join(output_dir, file_to_write_to), 'a') as out_handle:
        out_handle.writelines(lines_to_write)


def process_imgt_sequences(imgt_sequences, imgt_gca_clade_map, output_dir):

    imgt_gca_clade = pd.read_csv(imgt_gca_clade_map)

    with open(imgt_sequences, 'r') as imgt_sequence_handle:

        current_record = {
            'name': '',
            'fasta_header': '',
            'sequence': ''
        }


        for line in imgt_sequence_handle:
            if line[0] == '>':
                process_current_record(current_record, imgt_gca_clade, output_dir)
                current_record = {
                    'name': '',
                    'fasta_header': '',
                    'sequence': ''
                }

                fields = line.split("|")
                ab_name = fields[2]

                current_record['name'] = ab_name
                current_record['fasta_header'] = line

            else:
                sequence_so_far = current_record['sequence']
                current_record['sequence'] = f'{sequence_so_far}{line.strip()}'



def parse_args():
    parser = argparse.ArgumentParser(description="Arguments for script filter IMGT sequences and split by clade")
    parser.add_argument("--imgt_sequences", help="Path to IMGT sequence fasta")
    parser.add_argument("--imgt_gca_clade_map", help="Path to csv mapping GCA, clade and IMGT species names")
    parser.add_argument("--output_dir", help='where to write output gtfs')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    process_imgt_sequences(args.imgt_sequences, args.imgt_gca_clade_map, args.output_dir)