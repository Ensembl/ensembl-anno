#!/usr/bin/env bash
# Create pre-subsetted fixture files for the fast z_tritici integration test.
#
# Run from support_scripts/gmb/:
#   bash tests/create_fixtures.sh
#
# The fixtures cover the first 500 kb of chromosome 1 of the bundled
# Zymoseptoria tritici dataset.  This region contains ~47 K OrthoDB lines,
# ~700 Scallop/StringTie lines, and ~1500 Helixer lines, allowing the
# pipeline to run in ~8 seconds instead of the ~25 minutes needed for the
# full chromosome 1.
#
# After running this script, generate golden fixtures with:
#   python tests/generate_golden_fixtures.py

set -euo pipefail

REGION_COORD=500000
SRC="z_tritici"
OUTDIR="tests/fixtures/z_tritici_region1"

mkdir -p "$OUTDIR"

echo "Subsetting to chr 1, coords 1..${REGION_COORD} …"

# GTF files: tab-separated, field 1 = seqname, field 5 = end coordinate
awk -v c="$REGION_COORD" '$1=="1" && $5<=c' "$SRC/orthodb_geneset.gtf"  > "$OUTDIR/orthodb_geneset.gtf"
awk -v c="$REGION_COORD" '$1=="1" && $5<=c' "$SRC/scallop_geneset.gtf"  > "$OUTDIR/scallop_geneset.gtf"
awk -v c="$REGION_COORD" '$1=="1" && $5<=c' "$SRC/stringtie_geneset.gtf" > "$OUTDIR/stringtie_geneset.gtf"
awk -v c="$REGION_COORD" '$1=="1" && $5<=c' "$SRC/uniprot_geneset.gtf"  > "$OUTDIR/uniprot_geneset.gtf"

# GFF3 file: preserve the ##gff-version header
{
  grep "^#" "$SRC/helixer_remapped.gff3" | head -1
  awk -v c="$REGION_COORD" '$1=="1" && $5<=c' "$SRC/helixer_remapped.gff3"
} > "$OUTDIR/helixer_remapped.gff3"

# Genome FASTA: extract only the first 500 kb of chr 1
python3 - <<'PYEOF'
import sys

coord_limit = 500000
outdir = "tests/fixtures/z_tritici_region1"
with open("z_tritici/zymoseptoria_tritici.fa") as fin, \
     open(f"{outdir}/genome.fa", "w") as fout:
    capture = False
    seq_lines = []
    for line in fin:
        if line.startswith(">1"):
            capture = True
            fout.write(line)
        elif line.startswith(">") and capture:
            break
        elif capture:
            seq_lines.append(line.rstrip())
    seq = "".join(seq_lines)[:coord_limit]
    for i in range(0, len(seq), 60):
        fout.write(seq[i:i+60] + "\n")
print(f"Wrote {coord_limit:,} bp genome for chr 1")
PYEOF

# Copy the assembly report (not required for region fixtures but kept for reference)
cp "$SRC/GCF_000219625.1_MYCGR_v2.0_assembly_report.txt" "$OUTDIR/"

echo ""
echo "Line counts in $OUTDIR:"
for f in "$OUTDIR"/*.gtf "$OUTDIR"/*.gff3 "$OUTDIR"/*.fa; do
    printf "  %-40s %d lines\n" "$(basename $f)" "$(wc -l < "$f")"
done

echo ""
echo "Done.  Generate golden fixtures with:"
echo "  python tests/generate_golden_fixtures.py"
