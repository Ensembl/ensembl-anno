import argparse
import re

def write_rfam_selected_models_file(rfam_accession_file, 
                                    rfam_cm_db, 
                                    rfam_selected_models_file):

    with open(rfam_accession_file, encoding="utf-8") as rfam_accessions_in:
        rfam_accessions = rfam_accessions_in.read().splitlines()
    
    with open(rfam_cm_db, "r", encoding="utf-8") as rfam_cm_in:
        rfam_data = rfam_cm_in.read()

    rfam_models = rfam_data.split("//\n")  # get a chunck containing a model
    with open(rfam_selected_models_file, "w+", encoding="utf-8") as rfam_cm_out:
        for model in rfam_models:
            # The Rfam.cm file has INFERNAL and HMMR models, both are needed at this point
            # Later we just want the INFERNAL ones for looking at thresholds
            match = re.search(r"(RF\d+)", model)
            if match:
                model_accession = match.group(1)
                if model_accession in rfam_accessions:
                    rfam_cm_out.write(model + "//\n")
     


def parse_args():
    parser = argparse.ArgumentParser(description="Arguments for script to write a file of rfam models to use with cmsearch.")
    parser.add_argument("--rfam_accession_file", help="Path to rfam accession file")
    parser.add_argument("--rfam_cm_db", help="Path to rfam cm db")
    parser.add_argument("--rfam_selected_models_file", help="Path to output rfam models file")

    args = parser.parse_args()
    return args
    

if __name__ == "__main__":
    args = parse_args()

    write_rfam_selected_models_file(args.rfam_accession_file,
                                    args.rfam_cm_db,
                                    args.rfam_selected_model_file)



