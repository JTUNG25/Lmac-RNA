#!/bin/bash

#####################################
# Concatenate multi-lane FASTQ files
# With SIMPLIFIED OUTPUT NAMES
# FIXED: Better error reporting
# Usage: bash concatenate_fastq.sh
#####################################

# Don't exit on first error - we want to see what fails
# set -e

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'  # No Color

# Configuration
INPUT_DIR="."
OUTPUT_DIR="concatenated_fastq"
LOG_FILE="concatenate_fastq.log"

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Header
echo "========================================"
echo "   FASTQ Multi-lane Concatenation"
echo "   (Simplified output names)"
echo "========================================"
echo ""

# Check if any FASTQ files exist
if ! ls ${INPUT_DIR}/*_R1*.fastq.gz 1> /dev/null 2>&1; then
    print_error "No FASTQ files found in ${INPUT_DIR}/"
    print_error "Make sure your files end with _R1*.fastq.gz or _R2*.fastq.gz"
    exit 1
fi

# Create output directory
if [ -d "$OUTPUT_DIR" ]; then
    print_warning "Output directory '$OUTPUT_DIR' already exists"
    read -p "Do you want to overwrite it? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_status "Exiting without overwriting"
        exit 0
    fi
    rm -rf "$OUTPUT_DIR"
fi

mkdir -p "$OUTPUT_DIR"
print_success "Created output directory: $OUTPUT_DIR"
echo ""

# Initialize log
echo "Concatenation started at $(date)" > "$LOG_FILE"
echo "Input directory: $INPUT_DIR" >> "$LOG_FILE"
echo "Output directory: $OUTPUT_DIR" >> "$LOG_FILE"
echo "Output naming: Simplified (R2-5-2_R1.fastq.gz format)" >> "$LOG_FILE"
echo "" >> "$LOG_FILE"

# Extract unique sample names (everything before _L00X_R1)
print_status "Scanning for samples..."
declare -a SAMPLES
SAMPLES=($(ls ${INPUT_DIR}/*_R1*.fastq.gz 2>/dev/null | xargs -n1 basename | sed 's/_L00[0-9]_R1.*//' | sort -u))

TOTAL_SAMPLES=${#SAMPLES[@]}
print_success "Found $TOTAL_SAMPLES unique samples"
echo ""

# Check if we found any samples
if [ $TOTAL_SAMPLES -eq 0 ]; then
    print_error "No samples found. Check your filename format."
    exit 1
fi

# Process each sample
SUCCESSFUL=0
FAILED=0

print_status "Starting concatenation..."
echo "========================================"
echo ""

for sample in "${SAMPLES[@]}"; do
    # Find all R1 and R2 files for this sample
    R1_FILES=($(ls ${INPUT_DIR}/${sample}_L00[0-9]_R1*.fastq.gz 2>/dev/null | sort))
    R2_FILES=($(ls ${INPUT_DIR}/${sample}_L00[0-9]_R2*.fastq.gz 2>/dev/null | sort))
    
    R1_COUNT=${#R1_FILES[@]}
    R2_COUNT=${#R2_FILES[@]}
    
    # Validate that we have matching R1 and R2 files
    if [ $R1_COUNT -ne $R2_COUNT ]; then
        print_error "Sample $sample: R1 files ($R1_COUNT) don't match R2 files ($R2_COUNT)"
        echo "Sample $sample: ERROR - R1/R2 mismatch" >> "$LOG_FILE"
        ((FAILED++))
        continue
    fi
    
    if [ $R1_COUNT -eq 0 ]; then
        print_warning "Sample $sample: No files found"
        echo "Sample $sample: WARNING - No files found" >> "$LOG_FILE"
        ((FAILED++))
        continue
    fi
    
    # EXTRACT SIMPLIFIED SAMPLE NAME
    SIMPLIFIED_NAME=$(echo "$sample" | sed 's/_[A-Z0-9].*//')

    # Concatenate R1 and R2 files using SIMPLIFIED names
    R1_OUTPUT="${OUTPUT_DIR}/${SIMPLIFIED_NAME}_R1.fastq.gz"
    R2_OUTPUT="${OUTPUT_DIR}/${SIMPLIFIED_NAME}_R2.fastq.gz"
    
    # Check if output files already exist (might cause issues)
    if [ -f "$R1_OUTPUT" ]; then
        print_warning "Output $R1_OUTPUT already exists, skipping to avoid overwrite"
        echo "Sample $sample: WARNING - Output file already exists" >> "$LOG_FILE"
        ((FAILED++))
        continue
    fi
    
    # Try to concatenate and check for errors
    if cat ${R1_FILES[@]} > "$R1_OUTPUT" 2>&1; then
        print_status "✓ Concatenated R1 for $SIMPLIFIED_NAME"
    else
        print_error "Failed to concatenate R1 for $sample"
        echo "Sample $sample: ERROR - Failed to concatenate R1" >> "$LOG_FILE"
        ((FAILED++))
        continue
    fi
    
    if cat ${R2_FILES[@]} > "$R2_OUTPUT" 2>&1; then
        print_status "✓ Concatenated R2 for $SIMPLIFIED_NAME"
    else
        print_error "Failed to concatenate R2 for $sample"
        echo "Sample $sample: ERROR - Failed to concatenate R2" >> "$LOG_FILE"
        ((FAILED++))
        continue
    fi
    
    # Check file sizes to verify concatenation
    R1_SIZE=$(du -h "$R1_OUTPUT" | cut -f1)
    R2_SIZE=$(du -h "$R2_OUTPUT" | cut -f1)
    
    print_success "$sample → $SIMPLIFIED_NAME"
    echo "  → R1: ${SIMPLIFIED_NAME}_R1.fastq.gz ($R1_SIZE from $R1_COUNT files)"
    echo "  → R2: ${SIMPLIFIED_NAME}_R2.fastq.gz ($R2_SIZE from $R2_COUNT files)"
    
    # Log the operation
    echo "Sample $sample → $SIMPLIFIED_NAME: OK ($R1_COUNT lanes) - R1: $R1_SIZE, R2: $R2_SIZE" >> "$LOG_FILE"
    ((SUCCESSFUL++))
done

echo ""
echo "========================================"
echo ""

# Summary
print_status "Concatenation complete!"
print_success "Successfully processed: $SUCCESSFUL samples"
if [ $FAILED -gt 0 ]; then
    print_warning "Failed/Skipped: $FAILED samples"
fi
echo ""

# Final statistics
TOTAL_FILES=$(ls ${OUTPUT_DIR}/*.fastq.gz 2>/dev/null | wc -l)
TOTAL_SIZE=$(du -sh ${OUTPUT_DIR} | cut -f1)

print_status "Output directory: $OUTPUT_DIR"
print_status "Total files created: $TOTAL_FILES"
print_status "Total size: $TOTAL_SIZE"
echo ""

# List all output files
echo "Output files:"
ls -lh ${OUTPUT_DIR}/*.fastq.gz 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'
echo ""

# Log completion
echo "" >> "$LOG_FILE"
echo "Concatenation completed at $(date)" >> "$LOG_FILE"
echo "Successfully processed: $SUCCESSFUL samples" >> "$LOG_FILE"
echo "Failed/Skipped: $FAILED samples" >> "$LOG_FILE"
echo "Output directory size: $TOTAL_SIZE" >> "$LOG_FILE"

print_success "Log file saved: $LOG_FILE"
echo ""
print_status "Ready for Galaxy upload!"
echo ""