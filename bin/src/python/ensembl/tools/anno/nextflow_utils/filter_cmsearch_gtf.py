import argparse


def parse_rnafold_predictions(rnafold_predictions):
    id_list = []
    with open(rnafold_predictions, "r") as handle:
        line_number = 0

        for line in handle:
            if line_number % 3 == 0:
                if not line[0] == ">":
                    raise ValueError("RNAfold output file format not as expected - exiting")
                id = line.split(" ")[0].split(">")[1]
                id_list.append(id)
            line_number += 1

    return id_list


def extract_id_from_gtf_line(line):
    ninth_field = line.split("\t")[8]
    fields_in_ninth_field = ninth_field.split(";")
    repeat_id_field = [x for x in fields_in_ninth_field if "gene_id" in x]
    if len(repeat_id_field) != 1:
        raise ValueError(
            f"Expected one repeat_id but got a different number. \
                         Repeat id field: {','.join(repeat_id_field)}, line: {line}"
        )
    repeat_id = repeat_id_field[0].split(" ")[1].strip('"')

    # The purpose of this try/except is just to ensure that the repeat id is a number and give a clear
    # error message if not
    try:
        integer_repeat_id = int(repeat_id)
    except:
        raise ValueError(f"Could not coerce repeat id {repeat_id} to int - exiting. Full field {ninth_field}")

    return repeat_id


def filter_gtf(input_gtf, output_gtf, rnafold_predictions):
    ncRNA_ids_with_rnafold_hits = parse_rnafold_predictions(rnafold_predictions)

    with open(input_gtf, "r") as input_gtf_handle:
        with open(output_gtf, "w") as output_gtf_handle:
            for line in input_gtf_handle:
                ncRNA_id = extract_id_from_gtf_line(line)
                if ncRNA_id in ncRNA_ids_with_rnafold_hits:
                    output_gtf_handle.write(line)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Arguments for script to filter cmsearch hits not validated by RNAfold"
    )
    parser.add_argument("--input_gtf", help="Input gtf with all cmsearch hits")
    parser.add_argument("--output_gtf", help="Output filtered gtf")
    parser.add_argument("--rnafold_predictions", help="File containing RNAfold free energy predictions")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    filter_gtf(
        input_gtf=args.input_gtf, output_gtf=args.output_gtf, rnafold_predictions=args.rnafold_predictions
    )
