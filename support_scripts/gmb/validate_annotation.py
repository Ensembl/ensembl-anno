#!/usr/bin/env python3
"""
Validate Annotation
===================
Compare a predicted gene annotation (GFF3) against a reference annotation
and report sensitivity, precision, and gene counts.

Usage:
    python validate_annotation.py \\
        --prediction consensus_genes.gff3 \\
        --reference Candida_Auris.gff3
"""

import argparse
import pyranges as pr


def validate(pred_path, ref_path):
    """Run validation and print results."""
    print(f"Loading prediction: {pred_path}")
    pred_gr = pr.read_gff3(pred_path)
    pred_genes = pred_gr[pred_gr.Feature == 'gene']

    print(f"Loading reference:  {ref_path}")
    ref_gr = pr.read_gff3(ref_path)
    ref_genes = ref_gr[ref_gr.Feature == 'gene']

    total_ref = ref_genes.df['ID'].nunique()
    total_pred = len(pred_genes)

    # --- Stranded ---
    valid_pred = pred_genes[pred_genes.Strand.isin(['+', '-'])]
    total_valid = valid_pred.df['ID'].nunique()

    ovl_sens = ref_genes.overlap(valid_pred, strandedness='same')
    found = ovl_sens.df['ID'].nunique() if not ovl_sens.df.empty else 0
    sensitivity = found / total_ref * 100 if total_ref else 0

    ovl_prec = valid_pred.overlap(ref_genes, strandedness='same')
    confirmed = ovl_prec.df['ID'].nunique() if not ovl_prec.df.empty else 0
    precision = confirmed / total_valid * 100 if total_valid else 0

    # --- Unstranded ---
    ovl_un = ref_genes.overlap(pred_genes, strandedness=False)
    found_un = ovl_un.df['ID'].nunique() if not ovl_un.df.empty else 0
    sens_un = found_un / total_ref * 100 if total_ref else 0

    # --- Report ---
    print()
    print('=' * 60)
    print('VALIDATION RESULTS')
    print('=' * 60)
    print(f"  Reference genes:           {total_ref}")
    print(f"  Predicted genes (total):   {total_pred}")
    print(f"  Predicted genes (stranded):{total_valid}")
    print()
    print('  Stranded:')
    print(f"    Sensitivity: {found} / {total_ref} ({sensitivity:.2f}%)")
    print(f"    Precision:   {confirmed} / {total_valid} ({precision:.2f}%)")
    print()
    print('  Unstranded:')
    print(f"    Sensitivity: {found_un} / {total_ref} ({sens_un:.2f}%)")
    print('=' * 60)


def main():
    parser = argparse.ArgumentParser(
        description='Validate a predicted annotation against a reference')
    parser.add_argument('--prediction', required=True,
                        help='Predicted GFF3 file')
    parser.add_argument('--reference', required=True,
                        help='Reference GFF3 file')
    args = parser.parse_args()

    validate(args.prediction, args.reference)


if __name__ == '__main__':
    main()
