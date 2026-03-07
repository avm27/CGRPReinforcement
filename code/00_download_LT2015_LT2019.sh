#!/bin/bash
## ============================================================
## Download Script: LT2015 + LT2019 Fastq Files from SRA
## For use with: 01_preprocessing_LT2019_LT2015_merged.sh
## ============================================================
##
## Downloads 9 single-end fastq files from the SRA trace server.
## Both the LT2015 (Tuesta et al. 2015) and LT2019 (Tuesta et al.
## 2019) datasets are processed together in the merged pipeline.
## All files are placed in the same fastq directory.
##
## DEPENDENCIES: wget
##
## USAGE:
##   bash 00_download_LT2015_LT2019.sh
## ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

outdir="$REPO_ROOT/data/raw/LT19_PBN2_LT15_Merged/fastq"
mkdir -p "$outdir"

echo "Downloading LT2015 + LT2019 SRA fastq files to: $outdir"
echo ""

## LT2015 — Tuesta et al. 2015
declare -A lt2015=(
    ["SRR1647854"]="https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR1647854"
    ["SRR1647855"]="https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR1647855"
    ["SRR1647856"]="https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR1647856"
    ["SRR1647857"]="https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR1647857"
    ["SRR1647858"]="https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR1647858"
    ["SRR1647859"]="https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR1647859"
)

## LT2019 — Tuesta et al. 2019
declare -A lt2019=(
    ["SRR6294441"]="https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR6294441"
    ["SRR6294443"]="https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR6294443"
    ["SRR6294444"]="https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR6294444"
)

echo "=== LT2015 Samples ==="
for acc in "${!lt2015[@]}"; do
    url="${lt2015[$acc]}"
    outfile="$outdir/${acc}.fastq.gz"
    if [ -f "$outfile" ]; then
        echo "[SKIP] Already exists: $outfile"
    else
        echo "[DOWNLOAD] $acc"
        wget -q --show-progress -O "$outfile" "$url"
        echo "[OK] $acc -> $outfile"
    fi
done

echo ""
echo "=== LT2019 Samples ==="
for acc in "${!lt2019[@]}"; do
    url="${lt2019[$acc]}"
    outfile="$outdir/${acc}.fastq.gz"
    if [ -f "$outfile" ]; then
        echo "[SKIP] Already exists: $outfile"
    else
        echo "[DOWNLOAD] $acc"
        wget -q --show-progress -O "$outfile" "$url"
        echo "[OK] $acc -> $outfile"
    fi
done

echo ""
echo "LT2015 + LT2019 download complete. Files in: $outdir"
