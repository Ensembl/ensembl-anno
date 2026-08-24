import argparse
from typing import List

def beds_to_gtf(bedfile_list, gtf_path) -> None:  # pylint:disable = too-many-locals
    """
    Convert bed file into gtf file
    Args:
        output_dir : Working directory path.
    """
    with open(gtf_path, "w+", encoding="utf8") as gtf_out:
        gene_id = 1
        for bed_path in bedfile_list:
            print("Converting bed to GTF: %s", str(bed_path))
            with open(bed_path, "r", encoding="utf8") as bed_in:
                for line in bed_in:
                    elements = line.rstrip().split("\t")
                    seq_region_name = elements[0]
                    offset = int(elements[1])
                    strand = elements[5]
                    # sizes of individual block of exons
                    block_sizes = [size for size in elements[10].split(",") if size]
                    block_starts = [size for size in elements[11].split(",") if size]
                    exons = bed_block_to_exons(block_sizes, block_starts, offset)
                    transcript_start = None
                    transcript_end = None
                    exon_records = []
                    for i, exon_coords in enumerate(exons):
                        if transcript_start is None or exon_coords[0] < transcript_start:
                            transcript_start = exon_coords[0]

                        if transcript_end is None or exon_coords[1] > transcript_end:
                            transcript_end = exon_coords[1]

                        exon_line = (
                            f"{seq_region_name}\tminimap\texon\t{exon_coords[0]}\t"
                            f"{exon_coords[1]}\t.\t{strand}\t.\t"
                            f'gene_id "minimap_{gene_id}"; transcript_id "minimap_{gene_id}"; '
                            f'exon_number "{i+ 1}";\n'
                        )
                        exon_records.append(exon_line)
                    transcript_line = (
                        f"{seq_region_name}\tminimap\ttranscript\t{transcript_start}\t"
                        f"{transcript_end}\t.\t{strand}\t.\t"
                        f'gene_id "minimap_{gene_id}"; transcript_id "minimap_{gene_id}"\n'
                    )
                    gtf_out.write(transcript_line)
                    for exon_line in exon_records:
                        gtf_out.write(exon_line)
                    gene_id += 1


def bed_block_to_exons(block_sizes: List, block_starts: List, offset: int) -> List:
    """
    Extract exon size and start from exon feature block
    Args:
        block_sizes : Block feature size.
        block_starts : Block feature starts.
        offset : Feature offset.

    Returns:
        List of exon coordinates
    """
    exons = []
    for i, _ in enumerate(block_sizes):
        block_start = offset + int(block_starts[i]) + 1
        block_end = block_start + int(block_sizes[i]) - 1
        if block_end < block_start:
            print("Warning: block end is less than block start, skipping exon")
            continue
        exon_coords = [str(block_start), str(block_end)]
        exons.append(exon_coords)
    return exons


def parse_args():
    parser = argparse.ArgumentParser(description="Arguments for script to convert paftools bedfiles to gtfs")
    parser.add_argument("--bedfile_list", nargs='+', help="Path to input bedfile")
    parser.add_argument("--gtf_path", help="Path to output gtf")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    beds_to_gtf(bedfile_list=args.bedfile_list, 
                gtf_path=args.gtf_path)
