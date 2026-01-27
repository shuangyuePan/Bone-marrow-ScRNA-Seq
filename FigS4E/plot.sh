#!/bin/bash
set -euo pipefail

# Set input files relative to the script location or project root
# Adjust these paths as needed based on where the actual data files are located
ATAC_NC="./Control_merged.rpkm.bw"
ATAC_KO="./KO_merged.rpkm.bw"

# Output configuration
OUT_DIR="./results"
mkdir -p "$OUT_DIR"
ATAC_DIFF="$OUT_DIR/ATAC_diff.bw"

# Generate log2 ratio bigWig if not exists
if [ ! -s "$ATAC_DIFF" ]; then
  bigwigCompare -b1 "$ATAC_KO" -b2 "$ATAC_NC" --operation log2 --pseudocount 1e-6 -p 8 -o "$ATAC_DIFF"
fi

# Download and prepare TSS reference if needed
REF_TSS="$OUT_DIR/hg38_TSS.bed"
if [ ! -s "$REF_TSS" ]; then
  curl -L "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz" -o "$OUT_DIR/refGene.txt.gz"
  gunzip -f "$OUT_DIR/refGene.txt.gz"
  awk 'BEGIN{OFS="\t"} {
    chrom=$3; strand=$4; txStart=$5; txEnd=$6; gene=$13;
    if (substr(chrom,1,3)!="chr") chrom="chr"chrom;
    if (strand=="+") tss=txStart; else tss=txEnd;
    print chrom, tss, tss+1, gene, 0, strand;
  }' "$OUT_DIR/refGene.txt" > "$REF_TSS"
fi

# Compute matrix centered on TSS
# Note: Using all TSS from refGene as per standard practice for this type of plot,
# or user can substitute with a specific TSS bed file if filtering is required.
MATRIX_FILE="$OUT_DIR/TSS_ATAC_log2_matrix.gz"

if [ ! -s "$MATRIX_FILE" ]; then
  computeMatrix reference-point \
    --referencePoint TSS \
    -S "$ATAC_DIFF" \
    -R "$REF_TSS" \
    -a 2000 -b 2000 \
    --skipZeros \
    -p 8 \
    --binSize 20 \
    -o "$MATRIX_FILE"
fi

# Plot Heatmap and Profile
plotHeatmap -m "$MATRIX_FILE" \
  --heatmapHeight 10 \
  --heatmapWidth 4 \
  --colorMap 'seismic' \
  --zMin -3 --zMax 3 \
  --samplesLabel 'ATAC log2(KO/NC)' \
  -o "FigS4E_ATAC_TSS_heatmap.pdf" \
  --plotFileFormat pdf

echo "FigS4E generated: FigS4E_ATAC_TSS_heatmap.pdf"
