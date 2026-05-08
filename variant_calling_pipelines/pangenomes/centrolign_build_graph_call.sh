#!/usr/bin/env bash
#PBS -N Centrolign_build_graph
#PBS -l select=1:ncpus=1:mem=1500gb:scratch_local=300gb
#PBS -l walltime=08:00:00
#PBS -j oe

# ==========================================
# 1. INPUT ARGUMENTS & USAGE EXAMPLE
# ==========================================
# Example Run Command:
# qsub script.sh -v CHR="chr13",REF_ID="chm13",S1_ID="P010",S2_ID="P027",S3_ID="P028",OUT="/path/to/results"
#
# Arguments provided via -v (PBS variables):
# CHR    - Chromosome name (e.g., chr13_maternal)
# REF_ID - Label for the reference genome (e.g., chm13.chr13)
# S1_ID  - Label for sample 1 (Generation 1 relative)
# S2_ID  - Label for sample 2 (The primary query/proband)
# S3_ID  - Label for sample 3 (Generation 2 relative)
# OUT    - Path to the final output directory

if [ -z "$CHR" ] || [ -z "$REF_ID" ] || [ -z "$S1_ID" ] || [ -z "$S2_ID" ] || [ -z "$S3_ID" ] || [ -z "$OUT" ]; then
    echo "Error: Missing required variables CHR, REF_ID, S1_ID, S2_ID, S3_ID, or OUT."
    echo "Usage: qsub $0 -v CHR=\"c\",REF_ID=\"r\",S1_ID=\"s1\",S2_ID=\"s2\",S3_ID=\"s3\",OUT=\"/path\""
    exit 1
fi

# Map variables for clarity
OUT_DIR=$OUT
mkdir -p "$OUT_DIR"

# ==========================================
# 2. ENVIRONMENT & PATH SETUP
# ==========================================
PROJ_DIR="/storage/praha5-elixir/projects/bioinf-fi/polakova/BP"

# Tool paths
CENTROLIGN_BIN="$PROJ_DIR/centrolign/bin/centrolign"
VAR_MAT_BIN="$PROJ_DIR/centrolign/centrolign/make_var_mat"
PYTHON_CONVERT="$PROJ_DIR/__scripts/pangenomes/convert_matrix.py"

# ==========================================
# 2. ENVIRONMENT SETUP (USER DEFINED)
# ==========================================
# >>> ADD YOUR ENVIRONMENT SETUP HERE <<<
PROJ_DIR="/storage/praha5-elixir/projects/bioinf-fi/polakova/BP"
CONDA_BASE="/cvmfs/software.metacentrum.cz/conda/envs/miniforge3-25.3.1-0"
CONDA_ENV="/storage/praha1/home/nikolpolakovaa/.conda/envs/mashtree"

source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

# ==========================================
# 3. PREPARE INPUTS (SCRATCH)
# ==========================================
echo "Moving to scratch: $SCRATCHDIR"
cd "$SCRATCHDIR" || exit 1

# Define source file paths based on naming convention
REF_FA_SRC="$PROJ_DIR/__data/CHM13_per_chromosome/${CHR}.fasta"
S1_FA_SRC="$PROJ_DIR/__data/pedigree/extracted_centromeres_v1.1/${S1_ID}.fasta"
S2_FA_SRC="$PROJ_DIR/__data/pedigree/extracted_centromeres_v1.1/${S2_ID}.fasta"
S3_FA_SRC="$PROJ_DIR/__data/pedigree/extracted_centromeres_v1.1/${S3_ID}.fasta"

# Re-header and localize files
echo "Preparing localized FASTA files..."
sed "s/^>.*/>$REF_ID/" "$REF_FA_SRC" > "${REF_ID}.fasta"
sed "s/^>.*/>$S1_ID/" "$S1_FA_SRC" > "${S1_ID}.fasta"
sed "s/^>.*/>$S2_ID/" "$S2_FA_SRC" > "${S2_ID}.fasta"
sed "s/^>.*/>$S3_ID/" "$S3_FA_SRC" > "${S3_ID}.fasta"

cat "${REF_ID}.fasta" "${S1_ID}.fasta" "${S2_ID}.fasta" "${S3_ID}.fasta" > all_sequences.fasta

# ==========================================
# 4. RUN GUIDE TREE & CENTROLIGN
# ==========================================
echo "Building guide tree..."
mashtree --mindepth 0 --numcpus 1 "${REF_ID}.fasta" "${S1_ID}.fasta" "${S2_ID}.fasta" "${S3_ID}.fasta" > raw.nwk
sed 's/\.fasta//g' raw.nwk > guide_tree.nwk

echo "Running Centrolign pangenome construction for $CHR..."
"$CENTROLIGN_BIN" -T guide_tree.nwk -S "subproblems_${CHR}" all_sequences.fasta > pangenome.gfa
CENTRO_EXIT=$?

# ==========================================
# 5. GENERATE MATRIX & SAVE RESULTS
# ==========================================
if [ $CENTRO_EXIT -eq 0 ]; then
    echo "Extracting variant matrix..."
    "$VAR_MAT_BIN" --base --indels --mnvs pangenome.gfa > matrix.tsv

    echo "Running Python variant comparisons..."
    python3 "$PYTHON_CONVERT" matrix.tsv variants_final.1gen.tsv "$S2_ID" "$S1_ID"
    python3 "$PYTHON_CONVERT" matrix.tsv variants_final.2gen.tsv "$S2_ID" "$S3_ID"

    # Copy results to the user-defined OUT directory
    cp pangenome.gfa "$OUT_DIR/${CHR}_pangenome.gfa"
    cp guide_tree.nwk "$OUT_DIR/${CHR}_guide_tree.nwk"
    cp matrix.tsv "$OUT_DIR/${CHR}_matrix.tsv"
    cp variants_final.1gen.tsv "$OUT_DIR/${CHR}_variants_final.1gen.tsv"
    cp variants_final.2gen.tsv "$OUT_DIR/${CHR}_variants_final.2gen.tsv"

    echo "Pipeline complete. Results saved to: $OUT_DIR"
else
    echo "Error: Centrolign failed with exit code $CENTRO_EXIT"
    exit 1
fi