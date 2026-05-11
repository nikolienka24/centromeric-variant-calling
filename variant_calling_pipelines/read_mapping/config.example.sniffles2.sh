#!/usr/bin/env bash
# ==============================================================================
# config.example.sh — Example configuration for Sniffles2 variant calling
# ==============================================================================
# Copy this file to config.sh and fill in your actual paths.
# ==============================================================================

# ------------------------------------------------------------------------------
# SAMPLE IDENTIFICATION
# ------------------------------------------------------------------------------
SAMPLE_NAME="PAN027.chr7.maternal"
REF_ID="chr7_mat"

# ------------------------------------------------------------------------------
# INPUT FILES
# ------------------------------------------------------------------------------
REF="/storage/brno2/home/user/references/PAN011_hap2/PAN011.chr7.hap2.fasta"
READS_BAM="/storage/brno2/home/user/reads/PAN027.chr7.maternal.bam"
BAM_DIR="/storage/brno2/home/user/alignments/sniffles2"

# ------------------------------------------------------------------------------
# FILTER BED FILES
# ------------------------------------------------------------------------------
BED_OFFSETS="/storage/brno2/home/user/beds/offsets.bed"
BED_REGIONS="/storage/brno2/home/user/beds/centromeric_regions.bed"
BED_PROBLEMATIC="/storage/brno2/home/user/beds/problematic_regions.bed"
BED_GAPS="/storage/brno2/home/user/beds/gaps.bed"

# ------------------------------------------------------------------------------
# OUTPUT
# ------------------------------------------------------------------------------
OUT="/storage/brno2/home/user/results/sniffles2"

# ------------------------------------------------------------------------------
# CONDA
# ------------------------------------------------------------------------------
CONDA_BASE="/storage/brno2/home/user/miniconda3"
CONDA_ENV="sniffles-env"