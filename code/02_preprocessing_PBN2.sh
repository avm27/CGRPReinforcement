#!/bin/bash
## ============================================================
## Script 02: PBN2 RNA-Seq Pre-Processing
## Project:   CGRPReinforcement
## Data:      Lab-generated data — GEO Accession: GSE324173
## Read Type: Paired-end  *** NOTE: DIFFERENT FROM LT2015/LT2019 ***
## Genome:    mm10 (NCBI RefSeq)
## ============================================================
##
## IMPORTANT — KEY DIFFERENCES FROM SCRIPTS 01 & 02:
##   - Paired-end reads (R1 + R2) vs. single-end
##   - trim_galore uses --paired flag; output files are
##     *_val_1.fq.gz and *_val_2.fq.gz (NOT *_trimmed.fq)
##   - Trimmed files are kept gzip-compressed (--gzip)
##   - STAR uses --readFilesCommand zcat for compressed input
##   - Input filename convention: <sampleName>_R1_001.fastq.gz
##     The loop strips the last 16 characters (_R1_001.fastq.gz)
##   - Trimming and alignment are run as separate loops to
##     allow all samples to finish trimming before alignment begins
##
## DEPENDENCIES (see check_versions.sh):
##   trim_galore, STAR, samtools, stringtie, pigz
##
## DATA DOWNLOAD:
##   Files are available from GEO: GSE324173
##   Download fastq files and place them in the fastq/ directory
##   for this experiment (see directory structure in README.md).
##   Expected naming convention: <sampleName>_R1_001.fastq.gz
##                                <sampleName>_R2_001.fastq.gz
##
## USAGE:
##   bash 02_preprocessing_PBN2.sh
## ============================================================

set -euo pipefail

## ---- Experiment identifiers ---------------------------------
experimentname="CGRPReinforcement_PBN2"
experimentnamefull="CGRPReinforcement -- PBN2"

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
cores_trim="8"  
cores_align="10"

## ---- Sanity checks ------------------------------------------
echo "============================================================"
echo "  $experimentnamefull"
echo "  Linux Pre-Processing — Paired-End Pipeline"
echo "  GEO Accession: GSE324173"
echo "  Start: $(date)"
echo "============================================================"

for ref in "$mm10genomedir" "$mm10gtf" "$mm10fai"; do
    if [ ! -e "$ref" ]; then
        echo "[ERROR] Reference file/directory not found: $ref"
        exit 1
    fi
done

if [ ! -d "$fastqPath" ] || \
   [ -z "$(find "$fastqPath" -maxdepth 1 -name '*_R1_001.fastq.gz' 2>/dev/null)" ]; then
    echo "[ERROR] No *_R1_001.fastq.gz files found in: $fastqPath"
    echo "        Download files from GEO accession GSE324173 first."
    exit 1
fi

## ---- Create output directories ------------------------------
mkdir -p "$trimdir" "$starPath" "$bamdir" "$stringtiedir" "$FPKMdir" "$logsdir"

## ============================================================
## LOOP 1: Paired-end trimming (all samples trimmed first)
## Input convention: <sampleName>_R1_001.fastq.gz
## Loop strips last 16 chars to isolate sampleName
## ============================================================
echo ""
echo "--- PHASE 1: Trimming (paired-end) ---"

for sampleName in $(find "$fastqPath" -maxdepth 1 -type f -name "*_R1_001.fastq.gz" \
    | while read F; do basename "$F" | rev | cut -c 17- | rev; done \
    | sort | uniq)
do
    echo "[$(date +%T)] Trimming (paired-end): ${sampleName}"
    trim_galore \
        --cores "$cores_trim" \
        --output_dir "$trimdir" \
        --basename "${sampleName}" \
        --gzip \
        --paired \
        "$fastqPath/${sampleName}_R1_001.fastq.gz" \
        "$fastqPath/${sampleName}_R2_001.fastq.gz"

    mv "$trimdir"/*trimming_report.txt "$logsdir"
    echo "[$(date +%T)] Trimming complete: ${sampleName}"
done

## ============================================================
## LOOP 2: Alignment, filtering, quantification
## Input: trimmed paired files from Loop 1
##   *_val_1.fq.gz (R1) and *_val_2.fq.gz (R2)
## ============================================================
echo ""
echo "--- PHASE 2: Alignment and Quantification (paired-end) ---"

for sampleName in $(find "$fastqPath" -maxdepth 1 -type f -name "*_R1_001.fastq.gz" \
    | while read F; do basename "$F" | rev | cut -c 17- | rev; done \
    | sort | uniq)
do
    echo ""
    echo "------------------------------------------------------------"
    echo "  Processing sample: ${sampleName}"
    echo "------------------------------------------------------------"

    ## -- Step 1: STAR alignment (paired-end, compressed input) -
    echo "[$(date +%T)] Aligning to mm10: ${sampleName}"
    STAR \
        --runThreadN "$cores_align" \
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
        --readFilesCommand zcat \
        --readFilesIn \
            "$trimdir/${sampleName}_val_1.fq.gz" \
            "$trimdir/${sampleName}_val_2.fq.gz" \
        --outTmpDir "$starPath/${sampleName}_tmp"

    mv "$starPath"/*Log.out          "$logsdir" 2>/dev/null || true
    mv "$starPath"/*Log.final.out    "$logsdir" 2>/dev/null || true
    mv "$starPath"/*Log.progress.out "$logsdir" 2>/dev/null || true
    echo "[$(date +%T)] Alignment complete: ${sampleName}"

    ## -- Step 2: Filter uniquely mapped reads ------------------
    echo "[$(date +%T)] Filtering unique alignments: ${sampleName}"
    samtools view "$starPath/${sampleName}Aligned.sortedByCoord.out.bam" \
        | grep -w 'NH:i:1' \
        | samtools view -bt "$mm10fai" > "$bamdir/${sampleName}.bam"

    ## -- Step 3: Index BAM ------------------------------------
    echo "[$(date +%T)] Indexing BAM: ${sampleName}"
    samtools index "$bamdir/${sampleName}.bam"

    ## -- Step 4: StringTie quantification ---------------------
    echo "[$(date +%T)] StringTie quantification: ${sampleName}"
    stringtie \
        -p "$cores_align" \
        -G "$mm10gtf" \
        -e \
        -o "$stringtiedir/${sampleName}.gtf" \
        -A "$FPKMdir/${sampleName}_FPKM.txt" \
        -l "${sampleName}" \
        "$bamdir/${sampleName}.bam"

    ## -- Step 5: Remove large intermediate files --------------
    rm -f "$starPath/${sampleName}Aligned.sortedByCoord.out.bam"
    rm -f "$starPath/${sampleName}Aligned.toTranscriptome.out.bam"
    rm -f "$bamdir/${sampleName}.bam"
    rm -f "$bamdir/${sampleName}.bam.bai"
    rm -f "$starPath"/*ReadsPerGene.out.tab
    rm -f "$starPath"/*SJ.out.tab

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
