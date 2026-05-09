# Root project directory
PROJ_DIR="/path/to/your/project"

# --- RUN CONFIGURATION ---
# Unique ID for the output file (e.g., PAN027.chr3.maternal_2)
CHR_ID="PAN027.chr3.maternal"

# Exact name of sequence 1 in the BED/FASTA (e.g., PAN027.chr3.maternal)
S1="PAN027.chr3.maternal"

# Exact name of sequence 2 in the BED/FASTA (e.g., PAN010.chr3.haplotype2)
S2="PAN010.chr3.haplotype2"

# Final destination directory for results
OUT="$PROJ_DIR/path/to/output/folder"

# --- INPUT FASTA FILES ---
# Full paths to the input FASTA files
FASTA_REF="/path/to/your/sequence1.fasta"
FASTA_MUT="/path/to/your/sequence2.fasta"

# --- INPUT BED FILES ---
# BED file containing genomic offsets (column 1 = sequence name, column 2 = offset)
BED_OFFSETS="$PROJ_DIR/path/to/offsets.bed"

# BED file containing genomic gap regions to filter out
BED_GAPS="$PROJ_DIR/path/to/gaps.bed"

# BED file containing problematic regions to filter out
BED_PROBLEMATIC="$PROJ_DIR/path/to/problematic.bed"