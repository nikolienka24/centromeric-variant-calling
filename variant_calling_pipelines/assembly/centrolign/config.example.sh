# Root project directory
PROJ_DIR="/path/to/your/project/directory"

# Conda setup
CONDA_BASE="/path/to/your/conda/base"
ENV_PATH="/path/to/your/conda/env"

# --- INPUT FASTA FILES ---
FASTA_REF="$PROJ_DIR/path/to/reference.fasta"
FASTA_MUT="$PROJ_DIR/path/to/query.fasta"

# --- INPUT BED FILES ---
# BED file containing genomic offsets (column 1 = sequence name, column 2 = offset)
BED_OFFSETS="$PROJ_DIR/path/to/offsets.bed"

# BED file containing genomic gap regions to filter out
BED_GAPS="$PROJ_DIR/path/to/gaps.bed"

# BED file containing problematic regions to filter out (used for both ref and mut)
BED_PROBLEMATIC="$PROJ_DIR/path/to/problematic.bed"

# --- OUTPUT ---
# Folder where all results will be written
OUTPUT_FOLDER="$PROJ_DIR/path/to/output/folder"

# --- TOOLS ---
CENTROLIGN="$PROJ_DIR/path/to/centrolign/bin/centrolign"
PYTHON_SCRIPT="$PROJ_DIR/path/to/cigar_to_tsv.only_mutations.py"