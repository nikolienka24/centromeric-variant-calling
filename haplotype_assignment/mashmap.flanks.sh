#!/bin/bash
#PBS -N mashmap_align_flanks
#PBS -l select=1:ncpus=16:mem=64gb:scratch_local=200gb
#PBS -l walltime=02:00:00
#PBS -j oe

# ==============================================================================
# DOCUMENTATION:
# ==============================================================================
# Aligns left and right centromere flanking regions of a query sequence against
# two parental haplotypes using MashMap to determine haplotype identity.
#
# USAGE:
#    1. Copy config.example.sh to config.sh and fill in your paths.
#    2. Submit with: qsub mashmap.flanks.sh
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
if [ -z "$CHR" ] || [ -z "$QRY" ] || [ -z "$H1" ] || [ -z "$H2" ] || [ -z "$IN" ] || [ -z "$OUT" ]; then
    echo "ERROR: Missing required variables CHR, QRY, H1, H2, IN, or OUT in config.sh."
    exit 1
fi

if [ ! -d "$IN" ]; then
    echo "ERROR: Input directory $IN does not exist."
    exit 1
fi

mkdir -p "$OUT/$CHR"

# ==============================================================================
# 2. ENVIRONMENT SETUP
# ==============================================================================
# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$ENV_PATH" || { echo "ERROR: Failed to activate conda environment: $ENV_PATH"; exit 1; }

if ! command -v mashmap &> /dev/null; then
    echo "ERROR: mashmap not found in PATH."
    exit 1
fi

# ==============================================================================
# 3. SCRATCH SETUP
# ==============================================================================
echo "Setting up scratch and moving to $SCRATCHDIR"
cd "$SCRATCHDIR" || exit 1

# ==============================================================================
# 4. MAPPING
# ==============================================================================

# --- LEFT FLANK ---
echo "Mapping LEFT flank..."
# Prepare targets with H1/H2 prefixes to distinguish them in output
sed "s/>/>H1_L_/" "${IN}/${H1}_left_flank.fasta" > targets_L.fasta
sed "s/>/>H2_L_/" "${IN}/${H2}_left_flank.fasta" >> targets_L.fasta

mashmap -r targets_L.fasta \
        -q "${IN}/${QRY}_left_flank.fasta" \
        -t 16 -n 2 --pi 95 -o "${QRY}_left.out"

# --- RIGHT FLANK ---
echo "Mapping RIGHT flank..."
sed "s/>/>H1_R_/" "${IN}/${H1}_right_flank.fasta" > targets_R.fasta
sed "s/>/>H2_R_/" "${IN}/${H2}_right_flank.fasta" >> targets_R.fasta

mashmap -r targets_R.fasta \
        -q "${IN}/${QRY}_right_flank.fasta" \
        -t 16 -n 2 --pi 95 -o "${QRY}_right.out"

# ==============================================================================
# 5. POST-PROCESSING (BEST HITS)
# ==============================================================================
echo "Selecting best hits..."
FINAL_BEST_FILE="${QRY}_best_hits_combined.out"

# Select best hit for left flank (sorting by identity/score column)
if [ -f "${QRY}_left.out" ]; then
    sort -k11,11n "${QRY}_left.out" | tail -n 1 > "$FINAL_BEST_FILE"
fi

# Append best hit for right flank
if [ -f "${QRY}_right.out" ]; then
    sort -k11,11n "${QRY}_right.out" | tail -n 1 >> "$FINAL_BEST_FILE"
fi

# ==============================================================================
# 6. TRANSFER RESULTS
# ==============================================================================
echo "Saving results to $OUT/$CHR"
cp "$FINAL_BEST_FILE" "$OUT/$CHR/" || exit 1
cp "${QRY}_left.out" "$OUT/$CHR/" || exit 1
cp "${QRY}_right.out" "$OUT/$CHR/" || exit 1

# Clean up all intermediate scratch files
rm -f targets_L.fasta targets_R.fasta \
       "${QRY}_left.out" "${QRY}_right.out" \
       "$FINAL_BEST_FILE"

echo "Process complete. Results saved in: $OUT/$CHR"