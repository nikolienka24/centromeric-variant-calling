# --- CONDA SETUP ---
CONDA_BASE="/path/to/your/conda/base"
ENV_PATH="/path/to/your/conda/env"

# --- RUN CONFIGURATION ---
# Chromosome name (e.g., chrX)
CHR="chrX"

# Prefix for the query flank files (e.g., PAN027.paternal)
QRY="PAN027.paternal"

# Prefix for haplotype 1 flank files (e.g., PAN028.h1)
H1="PAN028.h1"

# Prefix for haplotype 2 flank files (e.g., PAN028.h2)
H2="PAN028.h2"

# --- INPUT / OUTPUT ---
# Input directory containing the .fasta flank files
IN="/path/to/your/flanks/directory"

# Base output directory — results go into <OUT>/<CHR>/
OUT="/path/to/your/output/folder"