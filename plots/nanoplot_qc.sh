#!/bin/bash
#PBS -N nanoplot_analysis
#PBS -l select=1:ncpus=2:mem=64gb
#PBS -l walltime=01:00:00
#PBS -j oe

# ==========================================
# 1. INPUT ARGUMENTS
# ==========================================
# Check if at least 2 arguments are provided
if [ "$#" -lt 2 ]; then
    echo "Error: Missing arguments."
    echo "Usage: qsub $0 <input.fastq> <output_dir>"
    exit 1
fi

INPUT_FASTQ=$1
OUT_DIR=$2

mkdir -p "$OUT_DIR"

# ==========================================
# 2. ENVIRONMENT SETUP (USER DEFINED)
# ==========================================
# >>> ADD YOUR ENVIRONMENT SETUP HERE <<<
# Example:
# source /path/to/conda/etc/profile.d/conda.sh
# conda activate nanoplot_env

# Validation: Check if NanoPlot is accessible
if ! command -v NanoPlot &> /dev/null; then
    echo "Error: NanoPlot not found in PATH."
    echo "Please edit the ENVIRONMENT SETUP section in this script."
    exit 1
fi

# ==========================================
# 3. ANALYSIS
# ==========================================
echo "Starting NanoPlot analysis..."
echo "Input: $INPUT_FASTQ"
echo "Output directory: $OUT_DIR"

NanoPlot --fastq "$INPUT_FASTQ" -o "$OUT_DIR"

echo "Analysis complete."