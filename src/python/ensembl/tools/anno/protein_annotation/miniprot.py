# See the NOTICE file distributed with this work for additional information
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
Miniprot aligns protein sequences to genomic DNA to recover spliced, protein-to-genome
alignments that are useful for gene annotation.

One of Miniprot’s key strengths is its ability to model introns (splicing) and
frameshifts while using an affine gap scoring scheme, allowing it to map proteins
directly onto whole genomes and infer protein-coding gene structures in newly
assembled species. This makes it especially valuable in comparative genomics
workflows where curated proteins from related organisms are used to transfer gene
models and functional context to a new genome.

Miniprot is designed for speed at genome scale, combining a seed–chain–extend style
mapping strategy with modern acceleration techniques (e.g., k-mer sketching and
vectorized/SIMD dynamic programming) to achieve high throughput while maintaining
accuracy. It is commonly run as a standalone command-line tool and is frequently
embedded in genome annotation pipelines, producing alignments in formats such as
PAF and/or GFF-based outputs.

Reference:
Li, H. (2023). Protein-to-genome alignment with miniprot. Bioinformatics, 39(1),
btad014. doi:10.1093/bioinformatics/btad014
"""

__all__ = ["run_miniprot"]

import logging
import logging.config
import os, sys
from pathlib import Path
import numpy as np
import re
import subprocess
import argparse


from ensembl.tools.anno.utils._utils import (
    check_exe,
    create_dir,
    check_gtf_content,
)

logger = logging.getLogger(__name__)


def run_miniprot(#pylint:disable=dangerous-default-value
    masked_genome: Path,
    output_dir: Path,
    protein_dataset: Path,
    miniprot_bin: Path = Path("miniprot"),
    num_threads: int = 1,
    top_n: int = 1,
    outs: float = 1.0,
    protein_set: str = "uniprot",
) -> None:
    

    check_exe(miniprot_bin)
    
    if protein_set not in {"uniprot", "orthodb"}:
        raise ValueError("protein_set must be either 'uniprot' or 'orthodb'")
    if protein_set == "uniprot":
        miniprot_dir = create_dir(output_dir, "miniprot_output/uniprot_output")
    elif protein_set == "orthodb":
        miniprot_dir = create_dir(output_dir, "miniprot_output/orthodb_output")
    	
    initial_output_file = miniprot_dir / "initial_annotation.gff"
    output_file = miniprot_dir / "annotation.gtf"
    if output_file.exists():
        transcript_count = check_gtf_content(output_file, "transcript")
        if transcript_count > 0:
            logger.info("Miniprot gtf file exists, skipping analysis")
            return

    if not masked_genome.exists():
        raise IOError(f"Masked genome file does not exist: {masked_genome}")
    if not protein_dataset.exists():
        raise IOError(f"Protein file does not exist: {protein_dataset}")
    
    miniprot_index_file = Path(str(masked_genome) + ".mpi")
    if not miniprot_index_file.exists():
        run_miniprot_index(miniprot_bin, masked_genome, miniprot_index_file, num_threads)
    else:
        logger.info ("Found an existing miniprot index, so will skip indexing")
	
    miniprot_cmd = [
        str(miniprot_bin),
        "-t", str(num_threads),
        "-N", str(top_n), #get exactly one alignment per protein, primary alignment
        f"--outs={outs}", #Keep only chains equal to best
        "--gff",
        str(miniprot_index_file),
        str(protein_dataset),
    ]
    logger.info("Running miniprot mapping")
    logger.info(" ".join(miniprot_cmd))
    
    with open(initial_output_file, 'w') as process_output_file:
        subprocess.run(miniprot_cmd, stdout=process_output_file,check=True)
    logger.info("Completed running miniprot")
    logger.info("Creating standardised GFF")
    generate_miniprot_gtf(miniprot_dir) 

def generate_miniprot_gtf(miniprot_dir:str)-> None:
    logger.info("generate_miniprot_gff")
    file_out_name: str = os.path.join(miniprot_dir, "annotation.gtf")
    for root, dirs, files in os.walk(miniprot_dir):
        for miniprot_file in files:
            miniprot_file=os.path.join(root, miniprot_file)
            if miniprot_file.endswith(".gff"):
                convert_miniprot_gff_to_gtf(input_file=miniprot_file,output_file=file_out_name)

def convert_miniprot_gff_to_gtf(input_file: str | Path, output_file: str | Path) -> None:
    """
    params  input_file: input filename
    params  output_file: output gtf filename
    """
    input_file = Path(input_file)
    output_file = Path(output_file)

    blocks =  open(input_file, 'r').read().split('\n#') # read the gff file and split it in blocks, where each block contain mapping of a single transcript, #PAF alignment blocks

    file_out=open(output_file, 'w+')

    for block in blocks:
        nblock_lines = [x for x in block.split("\n") if x != ""]
        if not nblock_lines:
            continue
        header_line = nblock_lines[0]
        nblock = [x for x in block.split("\n") if x!=''][1:] # split the block by newline and remove empty lines
        nblock = np.array([x.split("\t") for x in nblock if x!='']) # split each line in the nblock by tab
        nblock = nblock.astype("object") #avoid fixed-width string truncation
        if "fs:i:" in header_line: ##   FILTER: keep only fs:i:0 and st:i:0
            m_fs = re.search(r"fs:i:(\d+)", header_line)
            if m_fs and int(m_fs.group(1)) != 0:
                continue  # skip this block
        
        if "st:i:" in header_line:
            m_st = re.search(r"st:i:(\d+)", header_line) #   FILTER: keep only fs:i:0 and st:i:0
            if m_st and int(m_st.group(1)) != 0:
                continue  # skip this block

        if nblock.shape[0] != 0:
            nrows, ncols= nblock.shape
            nblock[0,2] =    nblock[0,2].replace('mRNA', 'transcript') # get the first line and substitute mRNA to transcript in column 2
            nblock[1:(nrows),2] = [x.replace('CDS', 'exon') for x in nblock[1:(nrows), 2]] # from second line onward substitute CDS to exon in column2

            target_info = [x.replace("Target=", "") for x in (re.split(';|\\s', nblock[0,8])) if "Target" in x][0] #retrive the protein Id from the Target feature
            gi_ti = 'gene_id "%s"; transcript_id "%s";'%(target_info, target_info)

            for i in range(0, nrows):
                if i ==0:
                    nblock[i,8] = gi_ti
                else:
                    nblock[i,8] = '%s exon_number "%s";'%(gi_ti, i)

            if nblock[(nrows-1),2] == "stop_codon":
                if nblock[0,6] == "-":
                    nblock[(nrows-2),3] = nblock[(nrows-1),3] # adjust the CDS position as per stop_codon line
                    nblock = nblock[:-1] # remove the last line containing stop_codon info

                if nblock[0,6] == "+":
                    nblock[(nrows-2),4] = nblock[(nrows-1),4] # adjust the CDS position as per stop_codon line
                    nblock = nblock[:-1] # remove the last line containing stop_codon info
            else:
                pass

            for ele in nblock:
                file_out.write("%s\n"% "\t".join(ele))
    file_out.close()

def run_miniprot_index(
    miniprot_bin:str | Path, 
    masked_genome: str | Path, 
    miniprot_index_file: str | Path, 
    num_threads: int,
) -> None:
    """
    Executes Miniprot indexing on masked genome
    Args:
    	   masked_genome : Masked genome file path.
	       miniprot_bin  : Software path.
	       num_threads   : number of threads
    	   Command line options:
	       -d FILE      save index to FILE
           -t INT       number of threads [4]
    """
    miniprot_bin = Path(miniprot_bin)
    masked_genome = Path(masked_genome)
    miniprot_index_file = Path(miniprot_index_file)

    miniprot_index_cmd = [
        str(miniprot_bin),
        "-t" + str(num_threads),
        "-d",
        str(miniprot_index_file),
        str(masked_genome),
    ]
    logger.info(" ".join(miniprot_index_cmd))
    subprocess.run(miniprot_index_cmd, check=True)
    logger.info("Completed running miniprot indexing")

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Miniprot arguments")
    parser.add_argument(
        "--masked_genome_file", required=True, help="Masked genome file path"
    )
    parser.add_argument("--output_dir", required=True, help="Output directory path")
    parser.add_argument("--protein_file", required=True, help="Path for the protein dataset")
    parser.add_argument(
        "--miniprot_bin",
        default="miniprot",
        help="Miniprot executable path",
    )
    parser.add_argument("--num_threads", type=int, default=1, help="Number of threads")
    parser.add_argument(
        "--protein_set",
        required=True,
        choices=["uniprot", "orthodb"],
        help="Protein set [uniprot, orthodb]",
    )
    return parser.parse_args()

def main():
    """Miniprot's entry-point."""
    args = parse_args()

    log_file_path = create_dir(args.output_dir, "log") / "miniprot.log"
    loginipath = Path(__file__).parents[6] / "conf" / "logging.conf"

    logging.config.fileConfig(
        loginipath,
        defaults={"logfilename": str(log_file_path)},
        disable_existing_loggers=False,
    )

    run_miniprot(
        Path(args.masked_genome_file),
        Path(args.output_dir),
        Path(args.protein_file),
        Path(args.miniprot_bin),
        args.num_threads,
        args.protein_set,
    )

if __name__ == "__main__":
    main()
