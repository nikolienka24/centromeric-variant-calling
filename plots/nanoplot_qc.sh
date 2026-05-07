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

INPUT_FASTQ=$IN
OUT_DIR=$OUT

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
echo "---------------------------------------------------"
echo "Starting NanoPlot QC Analysis"
echo "Input FASTQ: $INPUT_FASTQ"
echo "Output Dir:  $OUT_DIR"
echo "---------------------------------------------------"

# NanoPlot generates high-quality plots for long-read sequencing (Nanopore/PacBio)
# providing metrics like read length N50, quality scores, and yield.
NanoPlot --fastq "$INPUT_FASTQ" -o "$OUT_DIR"

echo "---------------------------------------------------"
echo "Analysis complete. Reports generated in $OUT_DIR"
echo "---------------------------------------------------"