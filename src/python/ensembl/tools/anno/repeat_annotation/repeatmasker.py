# See the NOTICE file distributed with this work for additional information #pylint: disable=missing-module-docstring
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
RepeatMasker is a program that screens DNA sequences for interspersed
repeats and low complexity DNA sequences.

References
----------
:cite:`repeatmasker`
"""

__all__ = ["run_repeatmasker"]

import argparse
import logging
import logging.config
import multiprocessing
import re
import subprocess
from os import PathLike
from pathlib import Path
from collections import defaultdict
from ensembl.tools.anno.utils._utils import (
    check_exe,
    check_gtf_content,
    create_dir,
    get_seq_region_length,
    get_sequence,
    get_slice_id,
    slice_output_to_gtf,
)

logger = logging.getLogger("__name__")


def run_repeatmasker(  # pylint:disable=too-many-arguments, too-many-positional-arguments, too-many-locals
    genome_file: PathLike,
    output_dir: Path,
    repeatmasker_bin: Path = Path("RepeatMasker"),
    library: str = "",
    repeatmasker_engine: str = "rmblast",
    species: str = "",
    num_threads: int = 1,
    bedtools_bin: str = "/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/bin/bedtools",  # pylint:disable=line-too-long
) -> None:
    """
    Executes RepeatMasker on the genome slices and stores the
    final annotation.gtf in repeatmasker_output

        :param genome_file: Genome file path.
        :type genome_file: PathLike
        :param output_dir: Output directory path.
        :type output_dir: Path
        :param repeatmasker_bin: RepeatMasker executable path.
        :type repeatmasker_bin: Path, default RepeatMasker
        :param library: Custom repeat library.
        :type library: str
        :param repeatmasker_engine: RepeatMasker engine.
        :type repeatmasker_engine: str, default rmblast
        :param species: Species name.
        :type species: str
        :param num_threads: Number of threads.
        :type num_threads: int, default 1
        :param bedtools_bin: Bedtools executable path.
        :type bedtools_bin: str, default "/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/bin/bedtools"#pylint:disable=line-too-long

        :return: None
        :rtype: None
    """

    # Use default path if user didn't supply one
    repeatmasker_bin = repeatmasker_bin or Path("RepeatMasker")

    check_exe(repeatmasker_bin)
    repeatmasker_dir = create_dir(output_dir, "repeatmasker_output")

    output_file = repeatmasker_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "repeat")
        if transcript_count > 0:
            logger.info("Repeatmasker gtf file exists")
            return
    logger.info("Creating list of genomic slices")
    seq_region_to_length = get_seq_region_length(genome_file, 5000)
    slice_ids_per_region = get_slice_id(
        seq_region_to_length, slice_size=1000000, overlap=0, min_length=5000
    )  # pylint:disable=line-too-long
    repeatmasker_cmd = [
        str(repeatmasker_bin),
        "-nolow",  # does not display simple repeats or low_complexity DNA in the annotation
        "-engine",
        repeatmasker_engine,
        "-dir",
        str(repeatmasker_dir),
    ]
    if not library:
        if not species:
            species = "homo"
        repeatmasker_cmd.extend(["-species", species])
    else:
        repeatmasker_cmd.extend(["-lib", library])
    logger.info("Running RepeatMasker %s", repeatmasker_cmd)
    pool = multiprocessing.Pool(num_threads)  # pylint: disable=consider-using-with
    for slice_id in slice_ids_per_region:
        pool.apply_async(
            _multiprocess_repeatmasker,
            args=(
                repeatmasker_cmd,
                slice_id,
                genome_file,
                repeatmasker_dir,
                bedtools_bin,
            ),
        )
    pool.close()
    pool.join()
    slice_output_to_gtf(repeatmasker_dir, "repeat_id", "repeatmask", True, ".rm.gtf")
    slice_output_to_masked_fasta(
        output_dir=repeatmasker_dir,
        output_fasta="genome.masked.fa",
    )
    slice_output_to_repeatmasker_out(
        output_dir=repeatmasker_dir,
        output_file="annotation.rm.tsv",
        provenance_file="annotation.provenance.tsv",
    )
    # for gtf_file in repeatmasker_dir.glob("*.rm.gtf"):
    #    gtf_file.unlink()


def _multiprocess_repeatmasker(  # pylint: disable=too-many-locals
    repeatmasker_cmd: list[str],
    slice_id: list[str],
    genome_file: Path,
    repeatmasker_dir: Path,
    bedtools_bin: str,
) -> None:
    """
    Run Repeatmasker on genomic slice

    Args:
        repeatmasker_cmd: RepeatMasker command to execute.
        slice_id: Slice ID to run RepeatMasker on.
        genome_file : Genome file path.
        repeatmasker_dir : RepeatMasker output directory path.
        bedtools_bin : Bedtools executable path.
    """

    region_name, start, end = slice_id
    logger.info(
        "Processing slice to find repeats with RepeatMasker: %s:%s:%s",
        region_name,
        start,
        end,
    )
    seq = get_sequence(
        region_name, int(start), int(end), 1, genome_file, repeatmasker_dir, bedtools_bin
    )  # pylint:disable=line-too-long
    slice_file_name = f"{region_name}.rs{start}.re{end}"
    region_file = repeatmasker_dir / f"{slice_file_name}.fa"
    with open(region_file, "w+", encoding="utf8") as region_fasta_out:
        region_fasta_out.write(f">{region_name}\n{seq}\n")
    region_results_file = Path(f"{region_file}.rm.gtf")
    output_file = Path(f"{region_file}.out")
    masked_file = Path(f"{region_file}.masked")
    tbl_file = Path(f"{region_file}.tbl")
    log_file = Path(f"{region_file}.log")
    cat_file = Path(f"{region_file}.cat")
    repeatmasker_cmd = repeatmasker_cmd.copy()
    repeatmasker_cmd.append(str(region_file))
    logger.info(repeatmasker_cmd)
    subprocess.run(repeatmasker_cmd, check=True)
    _create_repeatmasker_gtf(output_file, region_results_file, region_name)
    # output_file.unlink()
    # region_file.unlink()
    # masked_file.unlink(missing_ok=True)
    tbl_file.unlink(missing_ok=True)
    log_file.unlink(missing_ok=True)
    cat_file.unlink(missing_ok=True)


def slice_output_to_masked_fasta(
    output_dir: Path,
    output_fasta: str = "genome.masked.fa",
) -> None:
    """
    Reconstruct a repeat-masked genome FASTA from RepeatMasker
    slice outputs.

    Expected filenames:
        chromosome.rs1.re1000000.fa.masked
        chromosome.rs1000001.re2000000.fa.masked

    Handles overlapping slices by retaining only the
    non-overlapping sequence from each subsequent slice.

    Streams sequences to avoid storing whole chromosomes
    in memory.
    """

    slices_by_region = defaultdict(list)

    # Collect masked FASTA files
    for input_file in output_dir.glob("*.fa.masked"):

        match = re.match(
            r"(.+)\.rs(\d+)\.re(\d+)\.fa\.masked$",
            input_file.name,
        )

        if match is None:
            continue

        seq_region = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))

        slices_by_region[seq_region].append((start, end, input_file))

    output_path = output_dir / output_fasta

    with open(output_path, "w", encoding="utf8") as output:

        for seq_region, slices in sorted(slices_by_region.items()):

            # Sort slices by genomic start
            slices.sort(key=lambda x: (x[0], x[1]))

            output.write(f">{seq_region}\n")

            previous_end = 0

            for start, end, input_file in slices:

                # Number of bases already written for this region
                skip_bases = max(0, previous_end - start + 1)

                written = 0
                skipped = 0

                with open(input_file, "r", encoding="utf8") as fasta:

                    for line in fasta:

                        if line.startswith(">"):
                            continue

                        sequence = line.strip()

                        # Skip overlapping bases
                        if skipped < skip_bases:
                            skip = min(
                                len(sequence),
                                skip_bases - skipped,
                            )

                            sequence = sequence[skip:]
                            skipped += skip

                        if sequence:
                            output.write(sequence + "\n")
                            written += len(sequence)

                previous_end = max(previous_end, end)

            output.write("\n")


import re
from pathlib import Path
from collections import defaultdict


def slice_output_to_repeatmasker_out(
    output_dir: Path,
    output_file: str = "annotation.rm.tsv",
    provenance_file: str = "annotation.provenance.tsv",
) -> None:
    """
    Merge per-slice RepeatMasker .out files into a genome-coordinate
    annotation table, assigning globally unique repeat IDs.

    Expected input filenames:
        chr1.rs1.re1000000.fa.out
        chr1.rs1000001.re2000000.fa.out

    Outputs:
        annotation.rm.tsv
            Merged repeat annotations with genome coordinates,
            repeat names, classes/families and unique IDs.

        annotation.provenance.tsv
            Mapping between new repeat IDs, original RepeatMasker
            IDs and source slice filenames.

    Coordinates:
        RepeatMasker query coordinates are 1-based inclusive.
        The merged output retains 1-based inclusive coordinates.

    Notes:
        - Repeat fragments sharing the same original ID within
          a slice receive the same new ID.
        - Original IDs from different slices are treated separately.
        - Overlapping annotations from overlapping slices are
          not deduplicated.
        - Repeats crossing slice boundaries are not merged.
        - The query-left field is set to (0), because full
          chromosome lengths are not available here.
        - Memory usage is low because files are streamed.
    """

    output_dir = Path(output_dir)

    slices_by_region = defaultdict(list)

    # ---------------------------------------------------------
    # Collect RepeatMasker output files
    # ---------------------------------------------------------

    for input_file in output_dir.glob("*.fa.out"):

        match = re.fullmatch(
            r"(.+)\.rs(\d+)\.re(\d+)\.fa\.out",
            input_file.name,
        )

        if match is None:
            continue

        seq_region = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))

        if start < 1 or end < start:
            raise ValueError(f"Invalid slice coordinates in {input_file.name}")

        slices_by_region[seq_region].append((start, end, input_file))

    if not slices_by_region:
        raise ValueError(f"No RepeatMasker .out slices found in {output_dir}")

    # ---------------------------------------------------------
    # Initialise globally unique repeat ID mapping
    # ---------------------------------------------------------

    repeat_counter = 1
    repeat_id_mapping = {}

    annotation_path = output_dir / output_file
    provenance_path = output_dir / provenance_file

    # ---------------------------------------------------------
    # Open annotation and provenance outputs
    # ---------------------------------------------------------

    with (
        open(annotation_path, "w", encoding="utf8") as output,
        open(provenance_path, "w", encoding="utf8") as provenance,
    ):

        # Public annotation header
        output.write(
            "SW_score\tperc_div\tperc_del\tperc_ins\t"
            "query\tq_begin\tq_end\tq_left\tstrand\t"
            "repeat_name\tclass_family\t"
            "r_begin\tr_end\tr_left\trepeat_id\toverlap\n"
        )

        # Internal provenance header
        provenance.write("repeat_id\toriginal_repeat_id\tsource_slice\n")

        # -----------------------------------------------------
        # Process regions in deterministic order
        # -----------------------------------------------------

        for seq_region, slices in sorted(slices_by_region.items()):

            slices.sort(key=lambda x: (x[0], x[1]))

            previous_end = 0

            for start, end, input_file in slices:

                # Check for missing slices between intervals
                if start > previous_end + 1:
                    raise ValueError(
                        f"Missing RepeatMasker slice in "
                        f"{seq_region}: "
                        f"{previous_end + 1}-{start - 1}"
                    )

                offset = start - 1

                # -------------------------------------------------
                # Read RepeatMasker output
                # -------------------------------------------------

                with open(
                    input_file,
                    "r",
                    encoding="utf8",
                ) as rm_in:

                    for line_number, line in enumerate(rm_in, start=1):

                        fields = line.strip().split()

                        # Skip blank lines and standard headers
                        if not fields:
                            continue

                        if not fields[0].isdigit():
                            continue

                        if len(fields) < 15:
                            raise ValueError(
                                f"Malformed RepeatMasker record "
                                f"in {input_file.name}, "
                                f"line {line_number}: {line.strip()}"
                            )

                        # -----------------------------------------
                        # Extract native RepeatMasker fields
                        # -----------------------------------------

                        sw_score = fields[0]
                        perc_div = fields[1]
                        perc_del = fields[2]
                        perc_ins = fields[3]

                        q_begin_local = int(fields[5])
                        q_end_local = int(fields[6])

                        # Validate coordinates against slice length
                        slice_length = end - start + 1

                        if (
                            q_begin_local < 1
                            or q_end_local < q_begin_local
                            or q_end_local > slice_length
                        ):
                            raise ValueError(
                                f"Invalid repeat coordinates in "
                                f"{input_file.name}, "
                                f"line {line_number}: "
                                f"{q_begin_local}-{q_end_local}"
                            )

                        # Convert to genome coordinates
                        q_begin = q_begin_local + offset
                        q_end = q_end_local + offset

                        strand = fields[8]
                        repeat_name = fields[9]
                        class_family = fields[10]

                        # Preserve repeat consensus coordinates
                        r_begin = fields[11]
                        r_end = fields[12]
                        r_left = fields[13]

                        original_repeat_id = fields[14]

                        overlap = fields[15] if len(fields) > 15 else "."

                        # -----------------------------------------
                        # Assign globally unique repeat ID
                        # -----------------------------------------

                        repeat_key = (
                            seq_region,
                            start,
                            end,
                            original_repeat_id,
                        )

                        if repeat_key not in repeat_id_mapping:

                            repeat_id_mapping[repeat_key] = repeat_counter

                            # Write provenance once per new ID
                            provenance.write(
                                f"{repeat_counter}\t"
                                f"{original_repeat_id}\t"
                                f"{input_file.name}\n"
                            )

                            repeat_counter += 1

                        repeat_id = str(repeat_id_mapping[repeat_key])

                        # -----------------------------------------
                        # Write merged annotation
                        # -----------------------------------------

                        output.write(
                            "\t".join(
                                [
                                    sw_score,
                                    perc_div,
                                    perc_del,
                                    perc_ins,
                                    seq_region,
                                    str(q_begin),
                                    str(q_end),
                                    "(0)",
                                    strand,
                                    repeat_name,
                                    class_family,
                                    r_begin,
                                    r_end,
                                    r_left,
                                    repeat_id,
                                    overlap,
                                ]
                            )
                            + "\n"
                        )

                previous_end = max(previous_end, end)


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


def _create_repeatmasker_gtf(  # pylint: disable=too-many-locals
    output_file: Path,
    region_results_file: Path,
    region_name: str,
) -> None:
    """
    Read the fasta file and save the content in gtf format

    All the genomic slices are collected in a single gtf output with the following format:
    SW    perc perc perc query    position in query matching repeat       position in repeat
    score div. del. ins. sequence begin end (left)  repeat   class/family begin end  (left)  ID
    Args:
        output_file : GTF file with final results.
        region_results_file_path : GTF file with results per region.
        region_name : Coordinates of genomic slice.
    """
    with (
        open(output_file, "r", encoding="utf8") as repeatmasker_in,
        open(region_results_file, "w+", encoding="utf8") as repeatmasker_out,
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
                    f"{strand}\t.\trepeat_id{repeat_count}; "
                    f'repeat_name "{repeat_name}"; repeat_class "{repeat_class}"; '
                    f'repeat_type "{repeat_type}"; repeat_start "{repeat_start}"; '
                    f'repeat_end "{repeat_end}"; score "{score}";\n'
                )
                print(gtf_line)
                repeatmasker_out.write(gtf_line)
                repeat_count += 1


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="RepeatMasker's arguments")
    parser.add_argument("--genome_file", required=True, help="Genome file path")
    parser.add_argument("--output_dir", required=True, help="Output directory path")
    parser.add_argument(
        "--repeatmasker_bin",
        default="RepeatMasker",
        help="RepeatMasker executable path",
    )
    parser.add_argument(
        "--library",
        default="",
        help="Custom repeat library",
    )
    parser.add_argument(
        "--repeatmasker_engine",
        default="rmblast",
        help="RepeatMasker engine",
    )
    parser.add_argument(
        "--species",
        default="homo",
        help="Species name (used if no library is provided)",
    )
    parser.add_argument(
        "--num_threads",
        type=int,
        default=1,
        help="Number of threads",
    )
    parser.add_argument(
        "--bedtools_bin",
        default="/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/bin/bedtools",
        help="Bedtools executable path",
    )

    return parser.parse_args()


def main():
    """RepeatMasker's entry-point."""
    args = parse_args()

    log_file_path = create_dir(args.output_dir, "log") / "repeatmasking.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"

    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )

    run_repeatmasker(
        args.genome_file,
        args.output_dir,
        args.repeatmasker_bin,
        args.library,
        args.repeatmasker_engine,
        args.species,
        args.num_threads,
        args.bedtools_bin,
    )


if __name__ == "__main__":
    main()
