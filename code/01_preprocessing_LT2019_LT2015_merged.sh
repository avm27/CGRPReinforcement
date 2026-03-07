#!/bin/bash
## ============================================================
## Script 01: LT2019 + LT2015 Merged RNA-Seq Pre-Processing
## Project:   CGRPReinforcement
## Data:      Tuesta et al. 2015 + 2019 — Public SRA Data (merged)
## Read Type: Single-end
## Genome:    mm10 (NCBI RefSeq)
## ============================================================
##
## DEPENDENCIES:
##   trim_galore, STAR, samtools, stringtie, pigz
##
## DATA DOWNLOAD:
##   Run the companion download script first:
##     bash 00_download_LT2015_LT2019.sh
##
##   LT2015 SRA accessions:
##     SRR1647854, SRR1647855, SRR1647856
##     SRR1647857, SRR1647858, SRR1647859
##
##   LT2019 SRA accessions:
##     SRR6294441, SRR6294443, SRR6294444
##
##   All 9 fastq files should be placed together in the fastq/
##   directory for this experiment. Both datasets share the same
##   single-end processing parameters.
##
## USAGE:
##   bash 01_preprocessing_LT2019_LT2015_merged.sh
## ============================================================

set -euo pipefail

## ---- Experiment identifiers ---------------------------------
experimentname="CGRPReinforcement_LT19LT15"
experimentnamefull="CGRPReinforcement — LT2019 + LT2015 Datasets"

## ---- Resolve repository root (2 levels up from code/) ------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

## ---- Reference genome paths (edit to match your system) ----
mm10genomedir="${mm10genomedir:-$HOME/genome/mm10/ncbi_STAR}"
mm10gtf="${mm10gtf:-$HOME/genome/mm10/mm10.ncbiRefSeq.gtf}"
mm10fai="${mm10fai:-$HOME/genome/mm10/mm10.fa.fai}"

## ---- Python / prepDE path -----------------------------------
python="${python:-$(command -v python3)}"
stringtiescript="${stringtiescript:-$(command -v prepDE.py 2>/dev/null || echo '')}"

## ---- Project directory layout (relative to repo root) -------
projPath="$REPO_ROOT/data/raw/$experimentname"
fastqPath="$projPath/fastq"
trimdir="$projPath/trimmed_fastqs"
starPath="$projPath/STAR"
bamdir="$projPath/bams"
stringtiedir="$projPath/StringTie"
FPKMdir="$projPath/FPKM"
logsdir="$projPath/Logs"

## ---- Resources ----------------------------------------------
cores="10"

## ---- Sanity checks ------------------------------------------
echo "============================================================"
echo "  $experimentnamefull"
echo "  Linux Pre-Processing — Single-End Pipeline (Merged)"
echo "  Start: $(date)"
echo "============================================================"

for ref in "$mm10genomedir" "$mm10gtf" "$mm10fai"; do
    if [ ! -e "$ref" ]; then
        echo "[ERROR] Reference file/directory not found: $ref"
        exit 1
    fi
done

if [ ! -d "$fastqPath" ] || [ -z "$(find "$fastqPath" -maxdepth 1 -name '*.fastq.gz' 2>/dev/null)" ]; then
    echo "[ERROR] No .fastq.gz files found in: $fastqPath"
    echo "        Run 00_download_LT2015_LT2019.sh first."
    exit 1
fi

## ---- Create output directories ------------------------------
mkdir -p "$trimdir" "$starPath" "$bamdir" "$stringtiedir" "$FPKMdir" "$logsdir"

## ---- Main processing loop -----------------------------------
## Input files are expected as: <sampleName>.fastq.gz
## The loop strips the .fastq.gz suffix (last 9 characters).
## Both LT2015 and LT2019 fastqs are processed identically here
## since both are single-end libraries with the same parameters.

for sampleName in $(find "$fastqPath" -maxdepth 1 -type f -name "*.fastq.gz" \
    | while read F; do basename "$F" | rev | cut -c 10- | rev; done \
    | sort | uniq)
do
    echo ""
    echo "------------------------------------------------------------"
    echo "  Processing sample: ${sampleName}"
    echo "------------------------------------------------------------"

    ## -- Step 1: Adapter/quality trimming (single-end) ---------
    echo "[$(date +%T)] Trimming: ${sampleName}"
    trim_galore \
        --cores "$cores" \
        --output_dir "$trimdir" \
        --basename "${sampleName}" \
        --dont_gzip \
        "$fastqPath/${sampleName}.fastq.gz"

    mv "$trimdir"/*trimming_report.txt "$logsdir"
    echo "[$(date +%T)] Trimming complete: ${sampleName}"

    ## -- Step 2: STAR alignment (single-end) -------------------
    echo "[$(date +%T)] Aligning to mm10: ${sampleName}"
    STAR \
        --runThreadN "$cores" \
        --runMode alignReads \
        --genomeDir "$mm10genomedir" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outFilterType BySJout \
        --outSAMattributes NH HI AS NM MD XS \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNoverLmax 0.3 \
        --outFileNamePrefix "$starPath/${sampleName}" \
        --quantMode TranscriptomeSAM GeneCounts \
        --readFilesIn "$trimdir/${sampleName}_trimmed.fq" \
        --outTmpDir "$starPath/${sampleName}_tmp"
    echo "[$(date +%T)] Alignment complete: ${sampleName}"

    ## -- Step 3: Filter uniquely mapped reads ------------------
    echo "[$(date +%T)] Filtering unique alignments: ${sampleName}"
    samtools view "$starPath/${sampleName}Aligned.sortedByCoord.out.bam" \
        | grep -w 'NH:i:1' \
        | samtools view -bt "$mm10fai" > "$bamdir/${sampleName}.bam"

    ## -- Step 4: Index BAM ------------------------------------
    echo "[$(date +%T)] Indexing BAM: ${sampleName}"
    samtools index "$bamdir/${sampleName}.bam"

    ## -- Step 5: StringTie quantification ---------------------
    echo "[$(date +%T)] StringTie quantification: ${sampleName}"
    stringtie \
        -p "$cores" \
        -G "$mm10gtf" \
        -e \
        -o "$stringtiedir/${sampleName}.gtf" \
        -A "$FPKMdir/${sampleName}_FPKM.txt" \
        -l "${sampleName}" \
        "$bamdir/${sampleName}.bam"

    ## -- Step 6: Compress trimmed FASTQ and move logs ---------
    pigz "$trimdir/${sampleName}_trimmed.fq"
    mv "$starPath"/*Log.final.out "$logsdir" 2>/dev/null || true
    mv "$starPath"/*Log.out       "$logsdir" 2>/dev/null || true

    ## -- Step 7: Remove large intermediate files --------------
    rm -f "$starPath/${sampleName}Aligned.sortedByCoord.out.bam"
    rm -f "$starPath/${sampleName}Aligned.toTranscriptome.out.bam"
    rm -f "$bamdir/${sampleName}.bam"
    rm -f "$bamdir/${sampleName}.bam.bai"
    rm -f "$starPath"/*ReadsPerGene.out.tab
    rm -f "$starPath"/*SJ.out.tab
    rm -f "$starPath"/*Log.progress.out

    echo "[$(date +%T)] Sample ${sampleName} complete."
done

wait

echo ""
echo "============================================================"
echo "  $experimentnamefull"
echo "  Pre-Processing Complete — $(date)"
echo ""
echo "  Next step: once ALL datasets are processed, run:"
echo "    bash 03_generate_gene_count_matrix.sh"
echo "============================================================"
