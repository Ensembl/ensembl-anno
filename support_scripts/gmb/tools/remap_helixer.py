"""Rename the sequence IDs of a GFF3 (e.g. Helixer output) using an NCBI assembly report.

Upstream normalisation step: every GMB input must use the genome FASTA's sequence names, and
``gmb-preflight`` fails a track whose names are absent from the genome. Maps the GenBank
accession column to the Assigned-Molecule column (e.g. ``CM001196.1`` -> ``1``), including
``##sequence-region`` headers. Coordinates are unchanged, so this is only valid when each
accession is the whole assigned molecule (checked: the mapping must be one-to-one).
Sequences absent from the report are an error unless ``--allow-unmapped`` is given.

    python tools/remap_helixer.py --input helixer.gff3 \
        --assembly-report GCA_xxx_assembly_report.txt --output helixer_remapped.gff3
"""

import argparse
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Remap Helixer GFF3 sequence IDs using an NCBI assembly report."
    )
    parser.add_argument("--input", required=True, help="Input Helixer GFF3 file")
    parser.add_argument("--assembly-report", required=True, help="NCBI assembly report TXT file")
    parser.add_argument("--output", required=True, help="Output remapped GFF3 file")
    parser.add_argument(
        "--allow-unmapped",
        action="store_true",
        help="Keep records on sequences absent from the report under their original name "
        "(default: fail, because GMB would reject them anyway).",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    mapping = {}
    assigned_to_genbank = {}
    print(f"Parsing mapping from {args.assembly_report}...")
    with open(args.assembly_report) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            # 6: Assigned-Molecule (1, 2, 3...)
            # 4: GenBank-Accn (CM076438.1...)
            # Check indices based on:
            # Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn
            # 0             1               2                   3                               4

            if len(parts) >= 5:
                genbank = parts[4]
                assigned = parts[2]
                if assigned != "na":
                    mapping[genbank] = assigned
                    assigned_to_genbank.setdefault(assigned, []).append(genbank)

    collisions = {
        assigned: genbanks
        for assigned, genbanks in assigned_to_genbank.items()
        if len(genbanks) > 1
    }
    if collisions:
        examples = ", ".join(
            f"{assigned}<-{','.join(genbanks[:3])}"
            for assigned, genbanks in list(collisions.items())[:5]
        )
        raise SystemExit(
            "Assembly report mapping is not one-to-one; simple seqname remapping "
            "would collapse multiple sequence records without coordinate offsets "
            f"({examples}). Use accession seqnames end-to-end, or build true "
            "pseudomolecules and transform coordinates."
        )

    print(f"Loaded {len(mapping)} mappings: {mapping}")

    print(f"Remapping {args.input}...")
    remapped_count = 0
    unmapped = set()
    with open(args.input) as infile, open(args.output, "w") as outfile:
        for line in infile:
            if line.startswith("##sequence-region"):
                fields = line.split()
                if len(fields) >= 2 and fields[1] in mapping:
                    fields[1] = mapping[fields[1]]
                    line = " ".join(fields) + "\n"
                outfile.write(line)
                continue
            if line.startswith("#"):
                outfile.write(line)
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                outfile.write(line)
                continue

            seq_id = parts[0]
            if seq_id in mapping:
                parts[0] = mapping[seq_id]
                outfile.write("\t".join(parts) + "\n")
                remapped_count += 1
            else:
                unmapped.add(seq_id)
                outfile.write(line)

    if unmapped and not args.allow_unmapped:
        sys.exit(
            f"ERROR: {len(unmapped)} sequence(s) are not in the assembly report "
            f"(e.g. {sorted(unmapped)[:3]}); {args.output} is incomplete. Re-run with "
            "--allow-unmapped to keep them under their original names."
        )
    print(
        f"Created {args.output}: {remapped_count} feature lines remapped, "
        f"{len(unmapped)} unmapped sequence(s)."
    )


if __name__ == "__main__":
    main()
