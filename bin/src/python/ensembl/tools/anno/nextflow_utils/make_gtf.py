import argparse
from pathlib import Path
import re
from typing import Union

def create_red_gtf(repeat_coords_file: Path, output_file: Path):
    """
    Create Red gtf file from masked genome file

    Args:
        repeat_coords_file: Coordinates for repeats.
        output_file : GTF file with the final results.
    """
    with (
        open(repeat_coords_file, "r", encoding="utf8") as red_in,
        open(output_file, "w+", encoding="utf8") as red_out,
    ):
        for repeat_id, line in enumerate(red_in, start=1):
            result_match = re.search(r"^\>(.+)\:(\d+)\-(\d+)", line)
            if result_match:
                region_name = result_match.group(1)
                # Note that Red is 0-based, so add 1
                start = int(result_match.group(2)) + 1
                end = int(result_match.group(3)) + 1
                gtf_line = (
                    f"{region_name}\tRed\trepeat\t{start}\t"
                    f'{end}\t.\t+\t.\trepeat_id "{repeat_id}";\n'  # pylint:disable=line-too-long
                )
                red_out.write(gtf_line)

def create_dust_gtf(
    input_file: Path,
    output_gtf: Path,
    region_name: str,
) -> None:
    """
    All the genomic slices are collected in a single gtf output
    Args:
        input_file : GTF file with final results.
        output_gtf : GTF file with the results per region.
        region_name :Coordinates of genomic slice.
    """
    with (
        open(input_file, "r", encoding="utf8") as dust_in,
        open(output_gtf, "w+", encoding="utf8") as dust_out,
    ):
        repeat_count = 1
        for line in dust_in:
            result_match = re.search(r"(\d+)\ - (\d+)", line)
            if result_match:
                start = int(result_match.group(1)) + 1
                end = int(result_match.group(2)) + 1
                gtf_line = (
                    f"{region_name}\tDust\trepeat\t{start}\t"
                    f'{end}\t.\t+\t.\trepeat_id "{repeat_count}";\n'  # pylint:disable=line-too-long
                )
                dust_out.write(gtf_line)
                repeat_count += 1

# Function to find the repeat class based on the mappings
def get_repeat_type(repeat_type: str) -> str:
    """Get the repeat type based on the provided repeat_type string.

    Args:
        repeat_type (str): The repeat type string to match against the mappings.

    Returns:
        str: The corresponding repeat type description if a match is found,
        otherwise "Unknown".
    """
    mappings = {
        r"^Low_Comp": "Low complexity regions",
        r"^LINE": "Type I Transposons/LINE",
        r"^SINE": "Type I Transposons/SINE",
        r"^DNA": "Type II Transposons",
        r"^LTR": "LTRs",
        r"^Other": "Other repeats",
        r"^Satelli": "Satellite repeats",
        r"^Simple": "Simple repeats",
        r"^Tandem": "Tandem repeats",
        r"^TRF": "Tandem repeats",
        r"^Waterman": "Waterman",
        r"^Recon": "Recon",
        r"^Tet_repeat": "Tetraodon repeats",
        r"^MaskRegion": "Mask region",
        r"^dust": "Dust",
        r"^Unknown": "Unknown",
        r"RNA$": "RNA repeats",
    }
    for pattern, description in mappings.items():
        if re.match(pattern, repeat_type):
            return description
    return "Unknown"  # Default if no match is found


def create_repeatmasker_gtf(  # pylint: disable=too-many-locals
    input_file: Path,
    output_gtf: Path,
    region_name: str,
) -> None:
    """

    All the genomic slices are collected in a single gtf output with the following format:
    SW    perc perc perc query    position in query matching repeat       position in repeat
    score div. del. ins. sequence begin end (left)  repeat   class/family begin end  (left)  ID
    Args:
        input_file : GTF file with final results.
        output_gtf_path : GTF file with results per region.
        region_name : Coordinates of genomic slice.
    """
    with (
        open(input_file, "r", encoding="utf8") as repeatmasker_in,
        open(output_gtf, "w+", encoding="utf8") as repeatmasker_out,
    ):
        repeat_count = 1
        for line in repeatmasker_in:
            result_match = re.search(r"^\s*\d+\s+", line)
            if result_match:
                results = line.split()
                if results[-1] == "*":
                    results.pop()
                if len(results) != 15:
                    continue
                score = results[0]
                start = results[5]
                end = results[6]
                strand = results[8]
                repeat_name = results[9]
                repeat_class = results[10]
                repeat_type = get_repeat_type(results[10])
                if strand == "+":
                    repeat_start = results[11]
                    repeat_end = results[12]
                else:
                    repeat_start = results[13]
                    repeat_end = results[12]
                    strand = "-"
                gtf_line = (
                    f"{region_name}\tRepeatMasker\trepeat\t{start}\t{end}\t.\t"
                    f"{strand}\t.\trepeat_id {repeat_count}; "
                    f'repeat_name "{repeat_name}"; repeat_class "{repeat_class}"; '
                    f'repeat_type "{repeat_type}"; repeat_start "{repeat_start}"; '
                    f'repeat_end "{repeat_end}"; score "{score}";\n'
                )
                repeatmasker_out.write(gtf_line)
                repeat_count += 1


def create_trf_gtf(  # pylint:disable=too-many-locals, too-many-branches
    input_file: Path,
    output_gtf: Path,
    region_name: str,
) -> None:
    """

    TRF output format:
    cols 1+2:  Indices of the repeat relative to the start of the sequence
    col 3:     Period size of the repeat
    col 4:     Number of copies aligned with the consensus pattern
    col 5:     Size of consensus pattern (may differ slightly from the period size)
    col 6:     Percent of matches between adjacent copies overall
    col 7:     Percent of indels between adjacent copies overall
    col 8:     Alignment score
    cols 9-12: Percent composition for each of the four nucleotides
    col 13:    Entropy measure based on percent composition
    col 14:    Consensus sequence
    col 15:    Repeat sequence
    Args:
       input_file : GTF file with final results.
       output_gtf : GTF file with results per region.
       region_name : Coordinates of genomic slice.
    """
    with (
        open(input_file, "r", encoding="utf8") as trf_in,
        open(output_gtf, "w+", encoding="utf8") as trf_out,
    ):
        repeat_count = 1
        for line in trf_in:
            result_match = re.search(r"^\d+", line)
            if result_match:
                results = line.split()
                if len(results) != 15:
                    continue
                start = results[0]
                end = results[1]
                period = float(results[2])
                copy_number = float(results[3])
                percent_matches = float(results[5])
                score = float(results[7])
                repeat_consensus = results[13]
                if (  # pylint: disable=too-many-boolean-expressions
                    score < 50 and percent_matches >= 80 and copy_number > 2 and period < 10
                ) or (copy_number >= 2 and percent_matches >= 70 and score >= 50):
                    gtf_line = (
                        f"{region_name}\tTRF\trepeat\t{start}\t{end}\t.\t+\t.\t"
                        f'repeat_id "{repeat_count}"; score "{score}"; '
                        f'repeat_consensus "{repeat_consensus}";\n'
                    )
                    trf_out.write(gtf_line)
                    repeat_count += 1


def create_cpg_gtf(  # pylint:disable=too-many-arguments, too-many-locals, too-many-branches, too-many-positional-arguments
    input_file: Path,
    output_gtf: Path,
    region_name: str,
    cpg_min_length: int = 400,
    cpg_min_gc_content: int = 50,
    cpg_min_oe: float = 0.6,
) -> None:
    """
    Read the fasta file and save the content in gtf format
    All the genomic slices are collected in a single gtf output
    Args:
        input_file : GTF file with final results.
        output_gtf : GTF file with the results per region.
        region_name :Coordinates of genomic slice.
        cpg_dir : Output dir.
        cpg_min_length : Min length of CpG islands
        cpg_min_gc_content : Min GC frequency percentage
        cpg_min_oe :  Min ratio of the observed to expected number of CpG (CpGo/e)
    """
    with (
        open(input_file, "r", encoding="utf8") as cpg_in,
        open(output_gtf, "w+", encoding="utf8") as cpg_out,
    ):
        feature_count = 1
        for line in cpg_in:
            result_match = re.search(r"^" + region_name, line)
            if result_match:
                results = line.split()
                start = int(results[1])
                end = int(results[2])
                length = end - start + 1
                score = float(results[3])
                gc_content = float(results[6])
                oe_score_str = results[7]
                oe_score: Union[float, int]
                if oe_score_str in ("-", "inf"):
                    oe_score = 0
                else:
                    oe_score = float(oe_score_str)
                if (
                    int(length) >= int(cpg_min_length)
                    and gc_content >= int(cpg_min_gc_content)
                    and oe_score >= float(cpg_min_oe)
                ):
                    gtf_line = (
                        f"{region_name}\tCpG\tsimple_feature\t{start}\t"
                        f'{end}\t.\t+\t.\tfeature_id "{feature_count}"; score "{score}";\n'
                    )
                    cpg_out.write(gtf_line)

def create_eponine_gtf(
    input_file: Path,
    output_gtf: Path,
    region_name: str,
) -> None:
    """
    Read the fasta file and save the content in gtf format
    All the genomic slices are collected in a single gtf output
    Args:
        input_file: GTF file with final results.
        output_gtf: GTF file with the results per region.
        region_name: Coordinates of genomic slice.
    """
    with (
        open(input_file, "r", encoding="utf8") as eponine_in,
        open(output_gtf, "w+", encoding="utf8") as eponine_out,
    ):
        feature_count = 1
        for line in eponine_in:
            result_match = re.search(r"^" + region_name, line)
            if result_match:
                results = line.split()
                start = int(results[3])
                end = int(results[4])
                score = float(results[5])
                strand = results[6]
                print(results)
                # There's a one base offset on the reverse strand
                if strand == "-":
                    start -= 1
                    end -= 1

                gtf_line = (
                    f"{region_name}\tEponine\tsimple_feature\t"
                    f"{start}\t{end}\t.\t{strand}\t.\t"
                    f'feature_id "{feature_count}"; score "{score}";\n'
                )
                eponine_out.write(gtf_line)
                feature_count += 1

def parse_args():
    parser = argparse.ArgumentParser(description="Arguments for script to check contents of transcriptomic gtfs")
    parser.add_argument("--input_file", help="Path to input to convert to gtf")
    parser.add_argument("--output_gtf", help="Path to output logfile recording status of each gtf")
    parser.add_argument("--region_name", default = None, help= "Optional region name field")
    parser.add_argument("--red", action='store_true', help="convert red output to gtf")
    parser.add_argument("--dust", action='store_true', help="convert red output to gtf")
    parser.add_argument("--repeatmasker", action='store_true', help="convert red output to gtf")
    parser.add_argument("--trf", action='store_true', help="convert trf output to gtf")
    parser.add_argument("--cpg", action='store_true', help="convert cpg output to gtf")
    parser.add_argument("--eponine", action='store_true', help="convert eponine output to gtf")

    args = parser.parse_args()
    return args
    

if __name__ == "__main__":
    args = parse_args()
    if args.red:
        create_red_gtf(args.input_file, args.output_gtf)

    if args.dust:
        create_dust_gtf(args.input_file, args.output_gtf, args.region_name)

    if args.repeatmasker:
        create_repeatmasker_gtf(args.input_file, args.output_gtf, args.region_name)
    
    if args.trf:
        create_trf_gtf(args.input_file, args.output_gtf, args.region_name)

    if args.cpg:
        create_cpg_gtf(args.input_file, args.output_gtf, args.region_name)

    if args.eponine:
        create_eponine_gtf(args.input_file, args.output_gtf, args.region_name)
    