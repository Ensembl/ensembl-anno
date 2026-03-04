import pandas as pd
import csv
import sys

# URL for assembly report
REPORT_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/759/435/GCA_002759435.3_Cand_auris_B8441_V3/GCA_002759435.3_Cand_auris_B8441_V3_assembly_report.txt"
HELIXER_FILE = "GCA_002759435.3_Cand_auris_B8441_V3_genomic.gff3"
OUTPUT_FILE = "helixer_remapped.gff3"

def main():
    print("Downloading assembly report...")
    # We can use pandas to read the table directly if it's clean, but it has comments.
    # Let's use requests or just standard lib for robustness if we had it, but here we can just curl it first or use pandas with comment char.
    
    try:
        # Skip comment lines
        report = pd.read_csv(REPORT_URL, sep='\t', comment='#', header=None)
        # The columns in the file without comments might be tricky. 
        # Let's inspect the file structure again. 
        # The header line starts with # Sequence-Name, so pandas comment='#' might skip the header.
        # It's safer to download and parse manually.
        pass
    except:
        pass

    # Let's just download it using shell command in the workflow, but here inside python:
    import urllib.request
    urllib.request.urlretrieve(REPORT_URL, "assembly_report.txt")
    
    mapping = {}
    print("Parsing mapping...")
    with open("assembly_report.txt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            # 6: Assigned-Molecule (1, 2, 3...)
            # 4: GenBank-Accn (CM076438.1...)
            # Check indices based on:
            # Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn
            # 0             1               2                   3                               4
            
            if len(parts) >= 5:
                genbank = parts[4]
                assigned = parts[2]
                if assigned != 'na':
                    mapping[genbank] = assigned
                    
    print(f"Loaded {len(mapping)} mappings: {mapping}")
    
    print(f"Remapping {HELIXER_FILE}...")
    remapped_count = 0
    with open(HELIXER_FILE, 'r') as infile, open(OUTPUT_FILE, 'w') as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                outfile.write(line)
                continue
                
            seq_id = parts[0]
            if seq_id in mapping:
                parts[0] = mapping[seq_id]
                outfile.write('\t'.join(parts) + '\n')
                remapped_count += 1
            else:
                # If not in mapping, keep original or warn?
                # Usually keep original if it's a scaffold not in the main chromosome list
                outfile.write(line)
                
    print(f"Created {OUTPUT_FILE} with {remapped_count} overlapping lines remapped.")

if __name__ == "__main__":
    main()
