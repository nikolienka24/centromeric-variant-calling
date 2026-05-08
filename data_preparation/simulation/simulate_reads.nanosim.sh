#!/bin/bash
#PBS -N NanoSim_Simulation
#PBS -l select=1:ncpus=4:mem=256gb:scratch_local=1000gb
#PBS -l walltime=48:00:00
#PBS -j oe

# ==========================================
# 1. INPUT ARGUMENTS & USAGE EXAMPLE
# ==========================================
# Example Run Command:
# qsub simulation_script.sh -v REF="genom.fa",MODEL="/path/to/model",OUT="/path/to/out",COV=30
#
# Arguments provided via -v (PBS variables):
# REF   - Path to the reference FASTA file
# MODEL - Prefix path to the trained NanoSim model
# OUT   - Directory where simulated reads will be saved
# COV   - Desired coverage (e.g., 30, 60)

if [ -z "$REF" ] || [ -z "$MODEL" ] || [ -z "$OUT" ]; then
    echo "Error: Missing required variables REF, MODEL, or OUT."
    echo "Usage: qsub $0 -v REF=\"ref.fa\",MODEL=\"/model/prefix\",OUT=\"/out/dir\",[COV=60]"
    exit 1
fi

# Set default coverage to 60 if not provided
COVERAGE=${COV:-60}

# ==========================================
# 2. ENVIRONMENT SETUP (USER DEFINED)
# ==========================================
# >>> ADD YOUR ENVIRONMENT SETUP HERE <<<
CONDA_BASE="/cvmfs/software.metacentrum.cz/conda/envs/miniforge3-25.3.1-0"
ENV_PATH="/storage/praha5-elixir/projects/bioinf-fi/polakova/apps/miniconda3/envs/bioinf"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$ENV_PATH"

# Validation: Check if the NanoSim simulator is accessible
if ! command -v simulator.py &> /dev/null; then
    echo "Error: simulator.py (NanoSim) not found in PATH."
    echo "Please check your ENVIRONMENT SETUP section."
    exit 1
fi

# ==========================================
# 3. SIMULATION EXECUTION
# ==========================================
mkdir -p "$OUT"

echo "---------------------------------------------------"
echo "Starting NanoSim Simulation"
echo "Reference: $REF"
echo "Model:     $MODEL"
echo "Coverage:  $COVERAGE"
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
    -x "$COVERAGE" \
    -t 4 \
    --fastq \
    -hp -k 5

echo "---------------------------------------------------"
echo "Simulation complete."
echo "Results located in: $OUT"
echo "---------------------------------------------------"