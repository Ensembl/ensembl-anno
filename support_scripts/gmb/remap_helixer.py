import argparse
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Remap Helixer GFF3 sequence IDs using an NCBI assembly report."
    )
    parser.add_argument("--input", required=True, help="Input Helixer GFF3 file")
    parser.add_argument("--assembly-report", required=True, help="NCBI assembly report TXT file")
    parser.add_argument("--output", required=True, help="Output remapped GFF3 file")
    return parser.parse_args()


def main():
    args = parse_args()

    mapping = {}
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

    print(f"Loaded {len(mapping)} mappings: {mapping}")

    print(f"Remapping {args.input}...")
    remapped_count = 0
    with open(args.input, "r") as infile, open(args.output, "w") as outfile:
        for line in infile:
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
                # If not in mapping, keep original or warn?
                # Usually keep original if it's a scaffold not in the main chromosome list
                outfile.write(line)

    print(f"Created {args.output} with {remapped_count} overlapping lines remapped.")


if __name__ == "__main__":
    main()
