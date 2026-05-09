#!/bin/bash
#PBS -N extract_centromeres_per_file
#PBS -l select=1:ncpus=1:mem=16gb:scratch_local=400gb
#PBS -l walltime=01:00:00
#PBS -j oe

# ==============================================================================
# DOCUMENTATION:
# ==============================================================================
# Extracts centromeric regions from a reference FASTA using a BED file.
# Produces a single combined FASTA and individual per-region FASTA files.
#
# USAGE:
#    1. Copy config.example.sh to config.sh and fill in your paths.
#    2. Submit with: qsub extract_centromeres.sh
#
# NOTE: PBS directives are parsed before the shell runs and cannot use
#       variables, so any PBS-level paths must be set directly in the header.
# ==============================================================================

# --- LOAD USER CONFIGURATION ---
SCRIPT_DIR="$(dirname "$0")"
CONFIG="$SCRIPT_DIR/config.sh"
if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: config.sh not found. Copy config.example.sh to config.sh and fill in your paths."
    exit 1
fi
# shellcheck source=config.sh
source "$CONFIG"

# ==============================================================================
# 1. VALIDATE CONFIGURATION
# ==============================================================================
if [ -z "$REF" ] || [ -z "$BED" ] || [ -z "$OUT" ]; then
    echo "ERROR: Missing required variables REF, BED, or OUT in config.sh."
    exit 1
fi

# ==============================================================================
# 2. ENVIRONMENT SETUP
# ==============================================================================
# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$ENV_PATH" || { echo "ERROR: Failed to activate conda environment: $ENV_PATH"; exit 1; }

if ! command -v bedtools &> /dev/null || ! command -v samtools &> /dev/null; then
    echo "ERROR: bedtools or samtools not found in PATH."
    exit 1
fi

# ==============================================================================
# 3. SCRATCH SETUP AND DATA TRANSFER
# ==============================================================================
echo "Setting up scratch directory..."
cd "$SCRATCHDIR" || exit 1

cp "$REF" input_assembly.fasta
cp "$BED" centromeres.bed

samtools faidx input_assembly.fasta
mkdir -p individual_tmp

# ==============================================================================
# 4. SUBSET EXTRACTION
# ==============================================================================
echo "Extracting regions from BED file..."

# 4a. Create a single combined FASTA file
bedtools getfasta -fi input_assembly.fasta -bed centromeres.bed -fo reference.centromeres.fasta
samtools faidx reference.centromeres.fasta

# 4b. Create individual FASTA files for each row in the BED file
# Extra BED columns beyond chrom/start/end are captured in $rest and ignored
echo "Generating individual files..."
while read -r chrom start end rest; do
    # Skip empty lines, comments, or track headers
    [[ -z "$chrom" || "$chrom" == "track"* || "$chrom" == "#"* ]] && continue

    FILENAME="${chrom}_${start}_${end}.fasta"

    # Extract specific region into a temporary single-region BED
    echo -e "${chrom}\t${start}\t${end}" > temp_region.bed
    bedtools getfasta -fi input_assembly.fasta -bed temp_region.bed -fo "individual_tmp/${FILENAME}"
    samtools faidx "individual_tmp/${FILENAME}"

    # Clean up the temporary single-region BED immediately after use
    rm -f temp_region.bed
done < centromeres.bed

# ==============================================================================
# 5. TRANSFER RESULTS
# ==============================================================================
echo "Saving results to: $OUT"
mkdir -p "$OUT/per_chromosome"

cp reference.centromeres.fasta "$OUT/" || exit 1
cp reference.centromeres.fasta.fai "$OUT/" || exit 1
cp individual_tmp/*.fasta "$OUT/per_chromosome/" || exit 1
cp individual_tmp/*.fai "$OUT/per_chromosome/" || exit 1

# Clean up all intermediate scratch files
rm -f input_assembly.fasta input_assembly.fasta.fai \
       centromeres.bed reference.centromeres.fasta reference.centromeres.fasta.fai
rm -rf individual_tmp

echo "Process complete."