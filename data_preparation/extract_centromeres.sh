#!/bin/bash
#PBS -N extract_centromeres_per_file
#PBS -l select=1:ncpus=1:mem=16gb:scratch_local=400gb
#PBS -l walltime=01:00:00

# ==========================================
# 1. INPUT ARGUMENTS
# ==========================================
if [ "$#" -lt 3 ]; then
    echo "Error: Missing arguments."
    echo "Usage: qsub $0 <reference.fasta> <centromeres.bed> <output_dir>"
    exit 1
fi

REFERENCE_FASTA=$1
CENTROMERE_BED=$2
FINAL_OUT_DIR=$3
INDIVIDUAL_DIR="${FINAL_OUT_DIR}/per_chromosome"

# ==========================================
# 2. ENVIRONMENT SETUP (USER DEFINED)
# ==========================================
# >>> ADD YOUR ENVIRONMENT SETUP HERE <<<
# Example for Conda users:
# source /path/to/conda/etc/profile.d/conda.sh
# conda activate your_environment_name

# Validation: Check if required tools are accessible
if ! command -v bedtools &> /dev/null || ! command -v samtools &> /dev/null; then
    echo "Error: bedtools or samtools not found in PATH."
    echo "Please edit the ENVIRONMENT SETUP section in this script."
    exit 1
fi

# ==========================================
# 3. SCRATCH SETUP AND DATA TRANSFER
# ==========================================
echo "Setting up scratch directory..."
cd "$SCRATCHDIR" || exit 1

cp "$REFERENCE_FASTA" "input_assembly.fasta"
cp "$CENTROMERE_BED" "centromeres.bed"

# Index the main reference
samtools faidx "input_assembly.fasta"
mkdir -p "individual_tmp"

# ==========================================
# 4. SUBSET EXTRACTION
# ==========================================
echo "Extracting regions from BED file..."

# 4a. Create a single combined FASTA file
bedtools getfasta -fi "input_assembly.fasta" -bed "centromeres.bed" -fo "reference.centromeres.fasta"
samtools faidx "reference.centromeres.fasta"

# 4b. Create individual FASTA files for each row in the BED file
echo "Generating individual files..."
while read -r chrom start end rest; do
    # Skip empty lines, comments, or track headers
    [[ -z "$chrom" || "$chrom" == "track"* || "$chrom" == "#"* ]] && continue

    FILENAME="${chrom}_${start}_${end}.fasta"

    # Extract specific region
    echo -e "${chrom}\t${start}\t${end}" > temp_region.bed
    bedtools getfasta -fi "input_assembly.fasta" -bed temp_region.bed -fo "individual_tmp/${FILENAME}"

    # Index the small file
    samtools faidx "individual_tmp/${FILENAME}"
done < "centromeres.bed"

# ==========================================
# 5. TRANSFER RESULTS
# ==========================================
echo "Saving results to: $FINAL_OUT_DIR"
mkdir -p "$INDIVIDUAL_DIR"

# Copy combined files
cp "reference.centromeres.fasta" "$FINAL_OUT_DIR/"
cp "reference.centromeres.fasta.fai" "$FINAL_OUT_DIR/"

# Copy individual files
cp individual_tmp/*.fasta "$INDIVIDUAL_DIR/"
cp individual_tmp/*.fai "$INDIVIDUAL_DIR/"

echo "Process complete."