# Root project directory
PROJ_DIR="/path/to/your/project"

# --- RUN CONFIGURATION ---
# Chromosome name (e.g., chr13_maternal)
CHR="chr13"

# Label for the reference genome (e.g., chm13.chr13)
REF_ID="chm13.chr13"

# Label for sample 1 (Generation 1 relative)
S1_ID="PAN010.chr13.haplotype1"

# Label for sample 2 (The primary query/proband)
S2_ID="PAN027.chr13.maternal"

# Label for sample 3 (Generation 2 relative)
S3_ID="PAN028.chr13.maternal"

# Final destination directory for results
OUT="$PROJ_DIR/path/to/output/folder"

# --- INPUT FASTA FILES ---
# Full paths to the input FASTA files
FASTA_REF="/path/to/your/chm13_chromosome.fasta"
FASTA_S1="/path/to/your/sample1.fasta"
FASTA_S2="/path/to/your/sample2.fasta"
FASTA_S3="/path/to/your/sample3.fasta"

# --- GUIDE TREE ---
# Path to a pre-computed guide tree in Newick format
GUIDE_TREE="/path/to/your/guide_tree.nwk"

# --- TOOLS ---
CENTROLIGN_BIN="$PROJ_DIR/path/to/centrolign/bin/centrolign"
VAR_MAT_BIN="$PROJ_DIR/path/to/centrolign/make_var_mat"
PYTHON_CONVERT="$PROJ_DIR/path/to/convert_matrix.py"