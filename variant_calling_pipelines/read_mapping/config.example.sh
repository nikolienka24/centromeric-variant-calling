# Root project directory
PROJECT_DIR="/path/to/your/project"

# --- CONDA SETUP ---
CONDA_BASE="/path/to/your/conda/base"
CONDA_ENV="/path/to/your/conda/env"

# --- RUN CONFIGURATION ---
# Full sample name used for output file naming (e.g., PAN027.chr22.maternal)
SAMPLE_NAME="PAN027.chr22.maternal"

# Reference sequence ID — must match column 1 of BED_OFFSETS
REF_ID="PAN011.chr22.haplotype2"

# Final destination directory for results
OUT="$PROJECT_DIR/path/to/output/folder"

# --- INPUT FILES ---
# Reference FASTA
REF="/path/to/your/reference.fasta"

# ONT reads BAM file
READS_BAM="/path/to/your/reads.bam"

# BED file to restrict variant calling to centromeric regions
BED_REGIONS="/path/to/your/centromeres.bed"

# BED file containing genomic offsets (column 1 = sequence name, column 2 = offset)
BED_OFFSETS="$PROJECT_DIR/path/to/offsets.bed"

# BED file containing problematic regions to filter out
BED_PROBLEMATIC="$PROJECT_DIR/path/to/problematic.bed"

# BED file containing genomic gap regions to filter out
BED_GAPS="$PROJECT_DIR/path/to/gaps.bed"

# Pre-computed BAM directory (if alignment already exists, it will be reused)
BAM_DIR="/path/to/your/precomputed_bam_dir"

# --- TOOLS ---
# DeepVariant Singularity container image
CONTAINER_IMG="/path/to/your/deepvariant.sif"