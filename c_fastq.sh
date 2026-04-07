#!/bin/bash

#####################################
# Concatenate multi-lane FASTQ files
# Output: SAMPLENAME_R1/R2.fastq.gz
# Usage: bash concatenate_fastq.sh /path/to/fastq/dir
#####################################

set -euo pipefail

# Input directory from argument, or current directory
INPUT_DIR="${1:-.}"
OUTPUT_DIR="${INPUT_DIR}/concatenated_fastq"
LOG_FILE="${OUTPUT_DIR}/concatenate_fastq.log"

# Colors
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'

info()    { echo -e "${BLUE}[INFO]${NC}    $1"; }
success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
warn()    { echo -e "${YELLOW}[WARN]${NC}    $1"; }
error()   { echo -e "${RED}[ERROR]${NC}   $1"; }

echo "========================================"
echo "   FASTQ Multi-lane Concatenation"
echo "   Input: $INPUT_DIR"
echo "========================================"
echo ""

# Check input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    error "Input directory not found: $INPUT_DIR"
    exit 1
fi

# Check FASTQ files exist
if ! ls "${INPUT_DIR}"/*_L00[0-9]_R1*.fastq.gz 1>/dev/null 2>&1; then
    error "No lane-split FASTQ files found in: $INPUT_DIR"
    error "Expected pattern: SAMPLE_L001_R1.fastq.gz"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
success "Output directory: $OUTPUT_DIR"
echo ""

# Start log
echo "Started: $(date)"          > "$LOG_FILE"
echo "Input:   $INPUT_DIR"      >> "$LOG_FILE"
echo "Output:  $OUTPUT_DIR"     >> "$LOG_FILE"
echo ""                         >> "$LOG_FILE"

# Extract unique sample base names (everything before _L00X)
mapfile -t SAMPLES < <(
    ls "${INPUT_DIR}"/*_L00[0-9]_R1*.fastq.gz \
    | xargs -n1 basename \
    | sed 's/_L00[0-9]_R1.*//' \
    | sort -u
)

TOTAL=${#SAMPLES[@]}
info "Found $TOTAL samples to process"
echo ""

SUCCESSFUL=0
FAILED=0

for sample in "${SAMPLES[@]}"; do

    # Simplified name: strip barcode/flowcell suffix (_XXXXXXXX_BARCODE)
    SIMPLE=$(echo "$sample" | sed 's/_[A-Z0-9]\{8,\}_.*//')

    R1_FILES=($(ls "${INPUT_DIR}/${sample}_L00"[0-9]"_R1"*.fastq.gz 2>/dev/null | sort))
    R2_FILES=($(ls "${INPUT_DIR}/${sample}_L00"[0-9]"_R2"*.fastq.gz 2>/dev/null | sort))

    N_R1=${#R1_FILES[@]}
    N_R2=${#R2_FILES[@]}

    # Validate
    if [ "$N_R1" -eq 0 ]; then
        warn "No R1 files found for $sample — skipping"
        echo "SKIP $sample: no R1 files" >> "$LOG_FILE"
        ((FAILED++)); continue
    fi
    if [ "$N_R1" -ne "$N_R2" ]; then
        error "R1/R2 count mismatch for $sample ($N_R1 vs $N_R2) — skipping"
        echo "SKIP $sample: R1/R2 mismatch" >> "$LOG_FILE"
        ((FAILED++)); continue
    fi

    R1_OUT="${OUTPUT_DIR}/${SIMPLE}_R1.fastq.gz"
    R2_OUT="${OUTPUT_DIR}/${SIMPLE}_R2.fastq.gz"

    if [ -f "$R1_OUT" ]; then
        warn "Output already exists for $SIMPLE — skipping"
        echo "SKIP $sample: output exists" >> "$LOG_FILE"
        ((FAILED++)); continue
    fi

    info "Processing $SIMPLE ($N_R1 lanes)..."

    cat "${R1_FILES[@]}" > "$R1_OUT"
    cat "${R2_FILES[@]}" > "$R2_OUT"

    R1_SIZE=$(du -h "$R1_OUT" | cut -f1)
    R2_SIZE=$(du -h "$R2_OUT" | cut -f1)

    success "$sample → $SIMPLE  [R1: $R1_SIZE | R2: $R2_SIZE]"
    echo "OK $sample → $SIMPLE ($N_R1 lanes) R1: $R1_SIZE R2: $R2_SIZE" >> "$LOG_FILE"
    ((SUCCESSFUL++))

done

echo ""
echo "========================================"
echo ""
success "Done! $SUCCESSFUL succeeded, $FAILED failed/skipped"
echo ""
info "Output files:"
ls -lh "${OUTPUT_DIR}"/*.fastq.gz 2>/dev/null | awk '{print "  "$9, "("$5")"}'
echo ""
echo "Finished: $(date)" >> "$LOG_FILE"
echo "Succeeded: $SUCCESSFUL  Failed: $FAILED" >> "$LOG_FILE"
success "Log: $LOG_FILE"