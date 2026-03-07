#!/bin/bash
## ============================================================
## Script 03: Generate Merged Gene Count Matrix
## Project:   CGRPReinforcement
## Output:    gene_count_matrix_LT19_PBN2_LT15_Merged.csv
## ============================================================
##
## This script collects StringTie GTF outputs from both
## processed datasets and generates a single merged gene count
## matrix using prepDE.py.
##
## Run AFTER both preprocessing scripts have completed:
##   01_preprocessing_LT2019_LT2015_merged.sh
##   02_preprocessing_PBN2.sh
##
## INPUT — StringTie GTF directories (relative to repo root):
##   data/raw/CGRPReinforcement_LT19LT15/StringTie/
##   data/raw/CGRPReinforcement_PBN2/StringTie/
##
## OUTPUT:
##   data/gene_count_matrix_LT19_PBN2_LT15_Merged.csv
##
##
## FASTQ NAMING CONVENTIONS HANDLED:
##   SRA single-end  : <sampleName>.fastq.gz
##                     (strips last 9 chars: .fastq.gz)
##   PBN2 paired-end: <sampleName>_R1_001.fastq.gz
##                     (strips last 16 chars: _R1_001.fastq.gz)
##
## DEPENDENCIES: python, prepDE.py
##
## USAGE:
##   bash 03_generate_gene_count_matrix.sh
## ============================================================

set -euo pipefail

## ---- Experiment identifier ----------------------------------
experimentname="LT19_PBN2_LT15_Merged"

## ---- Resolve repository root (2 levels up from code/) ------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

## ---- Python / prepDE path -----------------------------------
python="${python:-$(command -v python3)}"
stringtiescript="${stringtiescript:-$(command -v prepDE.py 2>/dev/null || echo '')}"

## ---- Output paths -------------------------------------------
outdir="$REPO_ROOT/data"
gene_count_matrix="$outdir/gene_count_matrix_${experimentname}.csv"

## ---- Temporary sample list (deleted after use) --------------
sample_lst="$(mktemp /tmp/sample_lst_${experimentname}_XXXXXX.txt)"

## ---- Dataset definitions ------------------------------------
## Each entry: "fastqPath|stringtiePath|convention"
##   convention = "sra"     → strips 9 chars  (.fastq.gz)
##   convention = "paired"  → strips 16 chars (_R1_001.fastq.gz)
##                            only R1 files are used to get sample names

declare -a datasets=(
    "$REPO_ROOT/data/raw/CGRPReinforcement_LT19LT15/fastq|$REPO_ROOT/data/raw/CGRPReinforcement_LT19LT15/StringTie|sra"
    "$REPO_ROOT/data/raw/CGRPReinforcement_PBN2/fastq|$REPO_ROOT/data/raw/CGRPReinforcement_PBN2/StringTie|paired"
)

## ---- Sanity checks ------------------------------------------
echo "============================================================"
echo "  Generating Merged Gene Count Matrix"
echo "  Experiment: $experimentname"
echo "  Start: $(date)"
echo "============================================================"

if [ -z "$stringtiescript" ] || [ ! -f "$stringtiescript" ]; then
    echo "[ERROR] prepDE.py not found. Set the stringtiescript env var:"
    echo "        stringtiescript=/path/to/prepDE.py bash $0"
    exit 1
fi

mkdir -p "$outdir"

## ---- Build sample list from all three datasets --------------
echo ""
echo "Building sample list..."

for entry in "${datasets[@]}"; do
    fastqPath="$(echo "$entry"  | cut -d'|' -f1)"
    stringtiedir="$(echo "$entry" | cut -d'|' -f2)"
    convention="$(echo "$entry"  | cut -d'|' -f3)"

    if [ ! -d "$fastqPath" ]; then
        echo "[ERROR] fastq directory not found: $fastqPath"
        echo "        Ensure preprocessing has been run for all datasets."
        rm -f "$sample_lst"
        exit 1
    fi

    if [ ! -d "$stringtiedir" ]; then
        echo "[ERROR] StringTie directory not found: $stringtiedir"
        echo "        Ensure preprocessing has been run for all datasets."
        rm -f "$sample_lst"
        exit 1
    fi

    echo "  Adding samples from: $fastqPath ($convention)"

    if [ "$convention" = "sra" ]; then
        ## Single-end SRA files: <sampleName>.fastq.gz
        ## Strip last 9 characters (.fastq.gz)
        find "$fastqPath" -maxdepth 1 -type f -name "*.fastq.gz" \
            | while read F; do basename "$F" | rev | cut -c 10- | rev; done \
            | sort | uniq \
            | while read sampleName; do
                gtf="$stringtiedir/${sampleName}.gtf"
                if [ ! -f "$gtf" ]; then
                    echo "[WARNING] GTF not found, skipping: $gtf" >&2
                else
                    echo -e "${sampleName}\t${gtf}"
                fi
            done >> "$sample_lst"

    elif [ "$convention" = "paired" ]; then
        ## Paired-end files: <sampleName>_R1_001.fastq.gz
        ## Use R1 files only to extract sample names
        ## Strip last 16 characters (_R1_001.fastq.gz)
        find "$fastqPath" -maxdepth 1 -type f -name "*_R1_001.fastq.gz" \
            | while read F; do basename "$F" | rev | cut -c 17- | rev; done \
            | sort | uniq \
            | while read sampleName; do
                gtf="$stringtiedir/${sampleName}.gtf"
                if [ ! -f "$gtf" ]; then
                    echo "[WARNING] GTF not found, skipping: $gtf" >&2
                else
                    echo -e "${sampleName}\t${gtf}"
                fi
            done >> "$sample_lst"
    fi
done

## ---- Verify sample list is not empty ------------------------
sample_count=$(wc -l < "$sample_lst")
if [ "$sample_count" -eq 0 ]; then
    echo "[ERROR] Sample list is empty. Check that preprocessing completed"
    echo "        successfully for all three datasets and that GTF files exist."
    rm -f "$sample_lst"
    exit 1
fi

echo ""
echo "  Samples found: $sample_count"
echo ""
cat "$sample_lst"
echo ""

## ---- Run prepDE.py ------------------------------------------
echo "Generating gene count matrix..."
"$python" "$stringtiescript" \
    -i "$sample_lst" \
    -g "$gene_count_matrix"

## ---- Clean up temporary sample list ------------------------
rm -f "$sample_lst"
echo "Sample list removed."

echo ""
echo "Gene count matrix written to:"
echo "  $gene_count_matrix"
echo ""
echo "============================================================"
echo "  $experimentname — Gene Count Matrix Complete"
echo "  $(date)"
echo "============================================================"
