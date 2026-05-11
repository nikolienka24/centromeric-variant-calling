#!/usr/bin/env bash
#PBS -N Centrolign_build_graph
#PBS -l select=1:ncpus=1:mem=1500gb:scratch_local=300gb
#PBS -l walltime=08:00:00
#PBS -j oe

# ==============================================================================
# DOCUMENTATION:
# ==============================================================================
# Builds a pangenome graph from three pedigree samples and a CHM13 reference
# using Centrolign, then extracts variant matrices for two generation comparisons.
#
# USAGE:
#    1. Copy config.example.sh to config.sh and fill in your paths.
#    2. Submit with: qsub build_graph_call.sh
#
# NOTE: PBS directives are parsed before the shell runs and cannot use
#       variables, so any PBS-level paths must be set directly in the header.
# ==============================================================================

# --- LOAD USER CONFIGURATION ---
SCRIPT_DIR="$(dirname "$0")"
CONFIG="$SCRIPT_DIR/config.sh"
if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: config.sh not found. Copy config.example.sh to config.sh and fill in your paths."
    exit 1
fi
# shellcheck source=config.sh
source "$CONFIG"

# ==============================================================================
# 1. VALIDATE CONFIGURATION
# ==============================================================================
if [ -z "$CHR" ] || [ -z "$REF_ID" ] || [ -z "$S1_ID" ] || [ -z "$S2_ID" ] || [ -z "$S3_ID" ] || [ -z "$OUT" ] || [ -z "$GUIDE_TREE" ]; then
    echo "ERROR: Missing required variables CHR, REF_ID, S1_ID, S2_ID, S3_ID, OUT, or GUIDE_TREE in config.sh."
    exit 1
fi

mkdir -p "$OUT"

# ==============================================================================
# 2. PREPARE INPUTS (SCRATCH)
# ==============================================================================
echo "Moving to scratch: $SCRATCHDIR"
cd "$SCRATCHDIR" || exit 1

# Re-header and localize FASTA files to scratch;
# headers are replaced with the sample labels so Centrolign
# uses the correct sequence names in the graph
echo "Preparing localized FASTA files..."
sed "s/^>.*/>$REF_ID/" "$FASTA_REF" > "${REF_ID}.fasta"
sed "s/^>.*/>$S1_ID/" "$FASTA_S1" > "${S1_ID}.fasta"
sed "s/^>.*/>$S2_ID/" "$FASTA_S2" > "${S2_ID}.fasta"
sed "s/^>.*/>$S3_ID/" "$FASTA_S3" > "${S3_ID}.fasta"

cat "${REF_ID}.fasta" "${S1_ID}.fasta" "${S2_ID}.fasta" "${S3_ID}.fasta" > all_sequences.fasta

# ==============================================================================
# 3. RUN CENTROLIGN
# ==============================================================================
echo "Copying guide tree..."
cp "$GUIDE_TREE" guide_tree.nwk || { echo "ERROR: Failed to copy guide tree from $GUIDE_TREE"; exit 1; }

echo "Running Centrolign pangenome construction for $CHR..."
"$CENTROLIGN_BIN" -T guide_tree.nwk -S "subproblems_${CHR}" all_sequences.fasta > pangenome.gfa
CENTRO_EXIT=$?

if [ $CENTRO_EXIT -ne 0 ]; then
    echo "ERROR: Centrolign failed with exit code $CENTRO_EXIT"
    exit 1
fi

# Check that Centrolign produced a non-empty GFA file
if [[ ! -s pangenome.gfa ]]; then
    echo "ERROR: Centrolign produced an empty GFA file."
    exit 1
fi

# ==============================================================================
# 4. GENERATE MATRIX & SAVE RESULTS
# ==============================================================================
echo "Extracting variant matrix..."
"$VAR_MAT_BIN" --base --indels --mnvs pangenome.gfa > matrix.tsv

# Check that make_var_mat produced a non-empty matrix
if [[ ! -s matrix.tsv ]]; then
    echo "ERROR: make_var_mat produced an empty matrix file."
    exit 1
fi

echo "Running Python variant comparisons..."
python3 "$PYTHON_CONVERT" matrix.tsv variants_final.1gen.tsv "$S2_ID" "$S1_ID"
python3 "$PYTHON_CONVERT" matrix.tsv variants_final.2gen.tsv "$S2_ID" "$S3_ID"

# Copy results to the output directory
cp pangenome.gfa "$OUT/${CHR}_pangenome.gfa" || exit 1
cp guide_tree.nwk "$OUT/${CHR}_guide_tree.nwk" || exit 1
cp matrix.tsv "$OUT/${CHR}_matrix.tsv" || exit 1
cp variants_final.1gen.tsv "$OUT/${CHR}_variants_final.1gen.tsv" || exit 1
cp variants_final.2gen.tsv "$OUT/${CHR}_variants_final.2gen.tsv" || exit 1

# Clean up all intermediate scratch files
rm -f "${REF_ID}.fasta" "${S1_ID}.fasta" "${S2_ID}.fasta" "${S3_ID}.fasta" \
       all_sequences.fasta guide_tree.nwk \
       pangenome.gfa matrix.tsv \
       variants_final.1gen.tsv variants_final.2gen.tsv
rm -rf "subproblems_${CHR}"

echo "Pipeline complete. Results saved to: $OUT"