#!/bin/bash
#PBS -N nanoplot_analysis
#PBS -l select=1:ncpus=2:mem=64gb
#PBS -l walltime=01:00:00
#PBS -j oe

# ==========================================
# 1. INPUT ARGUMENTS & USAGE EXAMPLE
# ==========================================
# Example Run Command:
# qsub script.sh -v IN="data.fastq",OUT="/path/to/results"
#
# Arguments provided via -v (PBS variables):
# IN  - Path to the input FASTQ file (Long-read data)
# OUT - Directory where the NanoPlot reports/plots will be saved

if [ -z "$IN" ] || [ -z "$OUT" ]; then
    echo "Error: Missing required variables IN or OUT."
    echo "Usage: qsub $0 -v IN=\"data.fastq\",OUT=\"/output/dir\""
    exit 1
fi

mkdir -p "$OUT"

# ==========================================
# 2. ENVIRONMENT SETUP
# ==========================================
# >>> UPDATE THESE PATHS TO MATCH YOUR CLUSTER SETUP <<<
CONDA_BASE="/path/to/your/conda/base"
ENV_PATH="/path/to/your/conda/env"

# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$ENV_PATH" || { echo "ERROR: Failed to activate conda environment: $ENV_PATH"; exit 1; }

if ! command -v NanoPlot &> /dev/null; then
    echo "Error: NanoPlot not found in PATH."
    exit 1
fi

# ==========================================
# 3. ANALYSIS
# ==========================================
echo "---------------------------------------------------"
echo "Starting NanoPlot QC Analysis"
echo "Input FASTQ: $IN"
echo "Output Dir:  $OUT"
echo "---------------------------------------------------"

NanoPlot --fastq "$IN" -o "$OUT"

echo "---------------------------------------------------"
echo "Analysis complete. Reports generated in $OUT"
echo "---------------------------------------------------"