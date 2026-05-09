# Root project directory
PROJ_DIR="/path/to/your/project"

# --- CONDA SETUP ---
CONDA_BASE="/path/to/your/conda/base"
ENV_PATH="/path/to/your/conda/env"

# --- RUN CONFIGURATION ---
# Chromosome identifier for output naming
CHR_ID="chrX"

# Reference sequence ID (Gen1)
REF_ID="PAN010"

# Query sequence ID (Gen2)
QRY_ID="PAN027"

# Final destination directory for results
OUT="$PROJ_DIR/path/to/output/folder"

# --- INPUT FASTA FILES ---
# Full paths to the reference and query FASTA files
FASTA_REF="/path/to/your/reference.fasta"
FASTA_QRY="/path/to/your/query.fasta"

# --- INPUT BED FILES ---
# BED file containing genomic offsets (column 1 = sequence name, column 2 = offset)
BED_OFFSETS="$PROJ_DIR/path/to/offsets.bed"

# BED file containing genomic gap regions to filter out
BED_GAPS="$PROJ_DIR/path/to/gaps.bed"

# BED file containing problematic regions to filter out (used for both ref and query)
BED_PROBLEMATIC="$PROJ_DIR/path/to/problematic.bed"