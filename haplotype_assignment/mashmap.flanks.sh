#!/bin/bash
#PBS -N mashmap_align_flanks
#PBS -l select=1:ncpus=16:mem=64gb:scratch_local=200gb
#PBS -l walltime=02:00:00
#PBS -j oe

# ==========================================
# 1. INPUT ARGUMENTS
# ==========================================
# Check for the 6 required arguments
if [ "$#" -lt 6 ]; then
    echo "Error: Missing arguments."
    echo "Usage: qsub $0 <chr_name> <query_prefix> <h1_prefix> <h2_prefix> <input_dir> <output_dir>"
    echo ""
    echo "Example:"
    echo "  qsub $0 chrX PAN027.paternal PAN028.h1 PAN028.h2 /path/to/flanks /path/to/results"
    exit 1
fi

CHR=$1           # e.g., chrX
CHR_QUERY=$2     # e.g., PAN027.chrX.paternal
CHR_H1=$3        # e.g., PAN028.chrX.haplotype1
CHR_H2=$4        # e.g., PAN028.chrX.haplotype2
INPUT_DIR=$5     # Where your .fasta flank files are located
OUTPUT_DIR=$6    # Where you want the results saved

# Construct final output path (includes chromosome subfolder)
FINAL_OUT_DIR="$OUTPUT_DIR/$CHR"

# Create output directory if it doesn't exist
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
# source /path/to/conda/etc/profile.d/conda.sh
# conda activate bioinf

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