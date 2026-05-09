#!/usr/bin/env bash
#PBS -N centrolign_pairwise_alignment
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=1:mem=64gb:scratch_local=200gb

# ==============================================================================
# DOCUMENTATION:
# ==============================================================================
# 1. OFFSET EXTRACTION:
#    Extracts genomic offsets from BED files by matching the FASTA filename
#    against the 1st column of the BED.
#
# 2. TRIPLE FILTERING STRATEGY (bedtools):
#    Stage 1: Filter against PROBLEMATIC (BED_PROBLEMATIC).
#    Stage 3: Filter against Genomic GAPS (BED_GAPS).
#
# 3. COORDINATE MAPPING:
#    TSV Column 3 = Reference Position.
#    TSV Column 4 = Mutation Position.
#
# USAGE:
#    1. Copy config.example.sh to config.sh and fill in your paths.
#    2. Set CHR_ID below to the chromosome/sample you want to process.
#    3. Submit with: qsub align_call.sh
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

# --- PATH CONFIGURATION ---
# CHR_ID identifies the sample/chromosome being processed
CHR_ID="PAN027.chr7.maternal"

# Derived output paths (based on OUTPUT_FOLDER from config)
OUTPUT_CIGAR="$OUTPUT_FOLDER/$CHR_ID.cigar.txt"
OUTPUT_TSV="$OUTPUT_FOLDER/$CHR_ID.results.tsv"
OUTPUT_TSV_FILTERED="$OUTPUT_FOLDER/$CHR_ID.results.filtered.tsv"

# --- ENVIRONMENT SETUP ---
export LD_LIBRARY_PATH="$PROJ_DIR/centrolign/lib:$LD_LIBRARY_PATH"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$ENV_PATH" || { echo "ERROR: Failed to activate conda environment: $ENV_PATH"; exit 1; }

mkdir -p "$OUTPUT_FOLDER"
cd "$SCRATCHDIR" || exit 1

# --- DYNAMIC ID SETUP ---
ID_REF=$(basename "$FASTA_REF" .fasta)
ID_MUT=$(basename "$FASTA_MUT" .fasta)

echo "Processing Reference: $ID_REF and Mutation: $ID_MUT"

# --- OFFSET EXTRACTION ---
echo "Extracting offsets from BED files..."
OFF_REF=$(awk -v target="$ID_REF" '$1 == target {print $2; exit}' "$BED_OFFSETS")
OFF_MUT=$(awk -v target="$ID_MUT" '$1 == target {print $2; exit}' "$BED_OFFSETS")

if [[ -z "$OFF_REF" || -z "$OFF_MUT" ]]; then
    echo "ERROR: Could not find offsets for $ID_REF or $ID_MUT!"
    exit 1
fi

echo "Offsets detected: REF=$OFF_REF, MUT=$OFF_MUT"

# --- ALIGNMENT ---
echo "Running Centrolign..."
cat "$FASTA_REF" "$FASTA_MUT" > joined.fasta
"$CENTROLIGN" joined.fasta > "$OUTPUT_CIGAR" 2> "$OUTPUT_FOLDER/$CHR_ID.engine.log"

if [[ ! -s "$OUTPUT_CIGAR" ]]; then
    echo "ERROR: Centrolign produced an empty CIGAR file. Check $OUTPUT_FOLDER/$CHR_ID.engine.log"
    exit 1
fi

# --- CONVERSION TO TSV ---
echo "Converting CIGAR to TSV..."
python3 "$PYTHON_SCRIPT" \
    --ref "$FASTA_REF" \
    --mut "$FASTA_MUT" \
    --cigar "$OUTPUT_CIGAR" \
    --output "$OUTPUT_TSV" \
    --chrom "$CHR_ID" \
    --off_ref "$OFF_REF" \
    --off_mut "$OFF_MUT"

EXIT_CODE=$?

# --- MULTI-STAGE FILTERING ---
if [ $EXIT_CODE -eq 0 ]; then
    echo "Starting variant filtering pipeline..."

    # Stage 1: Reference problematic regions
    if [[ -f "$BED_PROBLEMATIC" ]]; then
        echo "Filtering Stage 1: Reference problematic regions..."
        awk -v bname="$ID_REF" 'NR > 1 {print bname"\t"$3"\t"$3"\t"$0}' "$OUTPUT_TSV" > stage1.bed
        bedtools intersect -a stage1.bed -b "$BED_PROBLEMATIC" -v > stage1_clean.bed
        awk 'BEGIN{FS=OFS="\t"} {for (i=4; i<=NF; i++) printf $i (i==NF?ORS:OFS)}' stage1_clean.bed > stage1.tsv
    else
        tail -n +2 "$OUTPUT_TSV" > stage1.tsv
    fi

    # Stage 2: Query problematic regions
    if [[ -f "$BED_PROBLEMATIC" ]]; then
        echo "Filtering Stage 2: Mutation problematic regions..."
        awk -v bname="$ID_MUT" 'BEGIN{FS=OFS="\t"} {print bname"\t"$4"\t"$4"\t"$0}' stage1.tsv > stage2.bed
        bedtools intersect -a stage2.bed -b "$BED_PROBLEMATIC" -v > stage2_clean.bed
        awk 'BEGIN{FS=OFS="\t"} {for (i=4; i<=NF; i++) printf $i (i==NF?ORS:OFS)}' stage2_clean.bed > stage2.tsv
    else
        cat stage1.tsv > stage2.tsv
    fi

    # Stage 3: Genomic Gaps
    if [[ -f "$BED_GAPS" ]]; then
        echo "Filtering Stage 3: Genomic gaps..."
        awk 'BEGIN{FS=OFS="\t"} {print $1"\t"$3"\t"$3"\t"$0}' stage2.tsv > stage3.bed
        bedtools intersect -a stage3.bed -b "$BED_GAPS" -v > stage3_clean.bed

        head -n 1 "$OUTPUT_TSV" > "$OUTPUT_TSV_FILTERED"
        awk 'BEGIN{FS=OFS="\t"} {for (i=4; i<=NF; i++) printf $i (i==NF?ORS:OFS)}' stage3_clean.bed >> "$OUTPUT_TSV_FILTERED"
    else
        head -n 1 "$OUTPUT_TSV" > "$OUTPUT_TSV_FILTERED"
        cat stage2.tsv >> "$OUTPUT_TSV_FILTERED"
    fi

    # Clean up all intermediate files
    rm -f joined.fasta \
           stage1.bed stage1_clean.bed stage1.tsv \
           stage2.bed stage2_clean.bed stage2.tsv \
           stage3.bed stage3_clean.bed

    echo "Filtering complete. Output: $OUTPUT_TSV_FILTERED"
    echo "Job completed successfully."
else
    echo "ERROR: Python conversion step failed with exit code $EXIT_CODE. Filtering skipped."
    exit $EXIT_CODE
fi