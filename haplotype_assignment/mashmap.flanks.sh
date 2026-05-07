#!/bin/bash
#PBS -N mashmap_align_flanks
#PBS -l select=1:ncpus=16:mem=64gb:scratch_local=200gb
#PBS -l walltime=02:00:00
#PBS -j oe

# ==========================================
# 1. INPUT ARGUMENTS & USAGE EXAMPLE
# ==========================================
# Example Run Command:
# qsub script.sh -v CHR="chrX",QRY="PAN027.paternal",H1="PAN028.h1",H2="PAN028.h2",IN="/path/flanks",OUT="/path/results"
#
# Arguments provided via -v (PBS variables):
# CHR - Chromosome name (e.g., chrX)
# QRY - Prefix for the query flank files
# H1  - Prefix for the Haplotype 1 flank files
# H2  - Prefix for the Haplotype 2 flank files
# IN  - Input directory containing the .fasta flank files
# OUT - Base output directory

# Check for required variables
if [ -z "$CHR" ] || [ -z "$QRY" ] || [ -z "$H1" ] || [ -z "$H2" ] || [ -z "$IN" ] || [ -z "$OUT" ]; then
    echo "Error: Missing required variables."
    echo "Usage: qsub $0 -v CHR=\"chrX\",QRY=\"prefix\",H1=\"h1\",H2=\"h2\",IN=\"/in\",OUT=\"/out\""
    exit 1
fi

# Map variables for internal use
CHR_QUERY=$QRY
CHR_H1=$1
CHR_H2=$H2
INPUT_DIR=$IN
FINAL_OUT_DIR="$OUT/$CHR"

# Create output directory
mkdir -p "$FINAL_OUT_DIR"

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory $INPUT_DIR does not exist."
    exit 1
fi

# ==========================================
# 2. ENVIRONMENT SETUP (USER DEFINED)
# ==========================================
# >>> ADD YOUR ENVIRONMENT SETUP HERE <<<
# Example:
# source /path/to/conda/etc/profile.d/conda.sh
# conda activate bioinf_env

# Validation: Check if mashmap is accessible
if ! command -v mashmap &> /dev/null; then
    echo "Error: mashmap not found in PATH."
    echo "Please edit the ENVIRONMENT SETUP section in this script."
    exit 1
fi

# ==========================================
# 3. SCRATCH SETUP
# ==========================================
echo "Setting up scratch and moving to $SCRATCHDIR"
cd "$SCRATCHDIR" || exit 1

# ==========================================
# 4. MAPPING
# ==========================================

# --- LEFT FLANK ---
echo "Mapping LEFT flank..."
# Prepare targets with H1/H2 prefixes to distinguish them in output
sed "s/>/>H1_L_/" "${INPUT_DIR}/${CHR_H1}_left_flank.fasta" > targets_L.fasta
sed "s/>/>H2_L_/" "${INPUT_DIR}/${CHR_H2}_left_flank.fasta" >> targets_L.fasta

mashmap -r targets_L.fasta \
        -q "${INPUT_DIR}/${CHR_QUERY}_left_flank.fasta" \
        -t 16 -n 2 --pi 95 -o "${CHR_QUERY}_left.out"

# --- RIGHT FLANK ---
echo "Mapping RIGHT flank..."
# Prepare targets with H1/H2 prefixes
sed "s/>/>H1_R_/" "${INPUT_DIR}/${CHR_H1}_right_flank.fasta" > targets_R.fasta
sed "s/>/>H2_R_/" "${INPUT_DIR}/${CHR_H2}_right_flank.fasta" >> targets_R.fasta

mashmap -r targets_R.fasta \
        -q "${INPUT_DIR}/${CHR_QUERY}_right_flank.fasta" \
        -t 16 -n 2 --pi 95 -o "${CHR_QUERY}_right.out"

# ==========================================
# 5. POST-PROCESSING (BEST HITS)
# ==========================================
echo "Selecting best hits..."
FINAL_BEST_FILE="${CHR_QUERY}_best_hits_combined.out"

# Select best hit for Left flank (sorting by identity/score column)
if [ -f "${CHR_QUERY}_left.out" ]; then
    sort -k11,11n "${CHR_QUERY}_left.out" | tail -n 1 > "$FINAL_BEST_FILE"
fi

# Append best hit for Right flank
if [ -f "${CHR_QUERY}_right.out" ]; then
    sort -k11,11n "${CHR_QUERY}_right.out" | tail -n 1 >> "$FINAL_BEST_FILE"
fi

# ==========================================
# 6. TRANSFER RESULTS
# ==========================================
echo "Saving results to $FINAL_OUT_DIR"
cp "$FINAL_BEST_FILE" "$FINAL_OUT_DIR/"
cp "${CHR_QUERY}_left.out" "$FINAL_OUT_DIR/"
cp "${CHR_QUERY}_right.out" "$FINAL_OUT_DIR/"

echo "Process complete. Results saved in: $FINAL_OUT_DIR"