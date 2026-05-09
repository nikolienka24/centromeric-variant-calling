#!/bin/bash
#PBS -N mashmap_simple_fasta
#PBS -l select=1:ncpus=16:mem=64gb:scratch_local=200gb
#PBS -l walltime=04:00:00
#PBS -j oe

# ==============================================================================
# DOCUMENTATION:
# ==============================================================================
# Aligns a query centromere FASTA against two parental haplotypes using MashMap
# and selects the best-matching haplotype based on sequence identity.
#
# USAGE:
#    1. Copy config.example.sh to config.sh and fill in your paths.
#    2. Submit with: qsub mashmap.centromeres.sh
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
if [ -z "$QRY" ] || [ -z "$H1" ] || [ -z "$H2" ] || [ -z "$OUT" ]; then
    echo "ERROR: Missing required variables QRY, H1, H2, or OUT in config.sh."
    exit 1
fi

mkdir -p "$OUT"

# FIX: strip both .fasta and .fa extensions to handle either naming convention
QUERY_NAME=$(basename "$QRY" .fasta)
QUERY_NAME=$(basename "$QUERY_NAME" .fa)
FINAL_OUT="$OUT/mashmap_best_${QUERY_NAME}.out"

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
# 3. PREPARE TARGETS (SCRATCH)
# ==============================================================================
echo "Moving to scratch: $SCRATCHDIR"
cd "$SCRATCHDIR" || exit 1

# Combine haplotypes into one reference file
# Headers are prefixed to distinguish haplotypes in the results
sed "s/>/>HAP1_/" "$H1" > targets.fasta
sed "s/>/>HAP2_/" "$H2" >> targets.fasta

# ==============================================================================
# 4. RUN MASHMAP
# ==============================================================================
echo "Running MashMap alignment..."

# --noSplit: Treat the query as one unit
# -n 2: Get hits for both haplotypes
mashmap -r targets.fasta \
        -q "$QRY" \
        -t 16 \
        --noSplit \
        -n 2 \
        --pi 90 \
        -o mashmap_raw.out

# ==============================================================================
# 5. SELECT BEST HIT & SAVE
# ==============================================================================
if [ -s mashmap_raw.out ]; then
    # Sort by column 10 (Identity %) numerically, highest first
    sort -k10,10nr mashmap_raw.out | head -n 1 > best_match.out

    WINNER=$(awk '{print $6}' best_match.out)
    ID_VAL=$(awk '{print $10}' best_match.out)
    echo "Alignment complete. Winner: $WINNER ($ID_VAL% identity)"

    cp best_match.out "$FINAL_OUT"
else
    echo "ERROR: No alignment found between query and haplotypes."
    exit 1
fi

# Clean up all intermediate scratch files
rm -f targets.fasta mashmap_raw.out best_match.out

echo "Result saved to: $FINAL_OUT"