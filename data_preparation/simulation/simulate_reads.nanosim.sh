#!/bin/bash
#PBS -N NanoSim_Simulation
#PBS -l select=1:ncpus=4:mem=256gb:scratch_local=1000gb
#PBS -l walltime=48:00:00
#PBS -j oe

# ==============================================================================
# DOCUMENTATION:
# ==============================================================================
# Simulates Oxford Nanopore long reads from a reference FASTA using NanoSim,
# based on a pre-trained error and read-length model.
#
# USAGE:
#    1. Copy config.example.sh to config.sh and fill in your paths.
#    2. Submit with: qsub simulate_reads.nanosim.sh
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
if [ -z "$REF" ] || [ -z "$MODEL" ] || [ -z "$OUT" ]; then
    echo "ERROR: Missing required variables REF, MODEL, or OUT in config.sh."
    exit 1
fi

# Default coverage to 60 if not set in config
COV=${COV:-60}

# ==============================================================================
# 2. ENVIRONMENT SETUP
# ==============================================================================
# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$ENV_PATH" || { echo "ERROR: Failed to activate conda environment: $ENV_PATH"; exit 1; }

if ! command -v simulator.py &> /dev/null; then
    echo "ERROR: simulator.py (NanoSim) not found in PATH."
    exit 1
fi

# ==============================================================================
# 3. SIMULATION EXECUTION
# ==============================================================================
mkdir -p "$OUT"

echo "---------------------------------------------------"
echo "Starting NanoSim Simulation"
echo "Reference: $REF"
echo "Model:     $MODEL"
echo "Coverage:  $COV"
echo "Output:    $OUT"
echo "---------------------------------------------------"

# NanoSim Command Breakdown:
# -rg: Reference genome
# -c:  Trained model prefix (characterizes error profile and read length)
# -o:  Output file prefix
# -x:  Coverage depth
# -t:  Number of threads (aligned with PBS ncpus)
# --fastq: Generate output in FASTQ format
# -hp: Enable homopolymer-aware simulation
# -k:  k-mer size for homopolymer simulation
simulator.py genome \
    -rg "$REF" \
    -c "$MODEL" \
    -o "$OUT/simulated_reads" \
    -x "$COV" \
    -t 4 \
    --fastq \
    -hp -k 5

echo "---------------------------------------------------"
echo "Simulation complete. Results located in: $OUT"
echo "---------------------------------------------------"