import sys

sys.path.insert(0, ".")
import pandas as pd

from compare_annotations import load_consensus_genes, load_gff, plot_comparison_locus
from subset_utils import load_assembly_mapping, remap_df_seqnames

# Load mapping
mapping = load_assembly_mapping("assembly_report.txt")

# Load reference
ref_exons, ref_cds, ref_genes = load_gff(
    "GCA_002759435.3_Cand_auris_B8441_V3_genomic.gff.gz", "GenBank"
)
ref_exons = remap_df_seqnames(ref_exons, mapping)
ref_cds = remap_df_seqnames(ref_cds, mapping)

# Load consensus
cons_genes, cons_exons, cons_cds, cons_mrna = load_consensus_genes("output/consensus.gff3")

# Load evidence
sc_ex, _, _ = load_gff("scallop_annotation.gtf", "Scallop")
st_ex, _, _ = load_gff("stringtie_annotation.gtf", "StringTie")
hx_ex, hx_cds, _ = load_gff("helixer_remapped.gff3", "Helixer")
od_ex, _, _ = load_gff("orthodb_annotation.gtf", "OrthoDB")
up_ex, _, _ = load_gff("uniprot_annotation.gtf", "UniProt")

tracks = {
    "GenBank": (ref_exons, ref_cds),
    "Consensus": (cons_exons, cons_cds),
    "Scallop": (sc_ex, pd.DataFrame()),
    "StringTie": (st_ex, pd.DataFrame()),
    "Helixer": (hx_ex, hx_cds),
    "OrthoDB": (od_ex, pd.DataFrame()),
    "UniProt": (up_ex, pd.DataFrame()),
}

# 5 diverse structural mismatches
loci = [
    ("2", 193814, 206264, "gene-B9J08_02051", "C_AURIS_4405"),  # 12.5kb chr2
    ("1", 1567469, 1577135, "gene-B9J08_00736", "C_AURIS_1807"),  # 9.7kb chr1
    ("4", 1057545, 1065273, "gene-B9J08_04266", "C_AURIS_6893"),  # 7.7kb chr4
    ("6", 793502, 798395, "gene-B9J08_05184", "C_AURIS_8945"),  # 4.9kb chr6
    ("3", 278665, 285262, "gene-B9J08_03151", "C_AURIS_5783"),  # 6.6kb chr3
]

for i, (chrom, s, e, ref_id, cons_id) in enumerate(loci):
    outf = f"output/validation/qc/structural_mismatch_extra_{i+1}_{chrom}_{s}.png"
    plot_comparison_locus(chrom, s, e, "Structural_Mismatch", ref_id, cons_id, tracks, outf)
    print(f"  Saved {outf}")

print("Done")
