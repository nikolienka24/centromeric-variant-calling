#!/usr/bin/env bash
#PBS -N chr7_pedigree_alignment
#PBS -q default@pbs-m1.metacentrum.cz
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=1:mem=64gb:scratch_local=200gb:cl_adan=True
#PBS -j oe
#PBS -o /storage/praha5-elixir/projects/bioinf-fi/polakova/BP/logs/chr7_pedigree_alignment.log
#PBS -m abe

# ==============================================================================
# DOCUMENTATION:
# ==============================================================================
# 1. OFFSET EXTRACTION:
#    Extracts genomic offsets from BED files by matching the FASTA filename
#    against the 4th column of the BED.
#
# 2. TRIPLE FILTERING STRATEGY (bedtools):
#    Stage 1: Filter Ref_Pos against PROBLEMATIC_REF.
#    Stage 2: Filter Mut_Pos against PROBLEMATIC_MUT.
#    Stage 3: Filter against Genomic GAPS (BED_GAPS).
#
# 3. COORDINATE MAPPING:
#    TSV Column 3 = Reference Position.
#    TSV Column 4 = Mutation Position.
# ==============================================================================

# --- PATH CONFIGURATION ---
PROJ_DIR="/storage/praha5-elixir/projects/bioinf-fi/polakova/BP"
CHR_ID="PAN027.chr7.maternal"

# Input FASTA Files
FASTA_REF="$PROJ_DIR/__data/pedigree/extracted_centromeres_v1.1/PAN010.chr7.haplotype2.fasta"
FASTA_MUT="$PROJ_DIR/__data/pedigree/extracted_centromeres_v1.1/PAN027.chr7.maternal.fasta"

# Input BED Files (Offsets & Filters)
BED_REF="$PROJ_DIR/__data/pedigree/annotations_v1.1/centromeres_bed/PAN010_hap2_HiFi_element_final_hap2.polished.cenSat.active_hor_merged_clean.bed"
BED_MUT="$PROJ_DIR/__data/pedigree/annotations_v1.1/centromeres_bed/PAN027_mat_HiFi_element_final_mat.polished.cenSat.active_hor_merged_clean.bed"
BED_GAPS="$PROJ_DIR/__data/pedigree/annotations_v1.1/gaps_all.bed"

# Problematic regions BED
PROBLEMATIC_REF="$PROJ_DIR/__data/pedigree/annotations_v1.1/problematic.PAN010.bed"
PROBLEMATIC_MUT="$PROJ_DIR/__data/pedigree/annotations_v1.1/problematic.PAN027.bed"

# Output Configuration
OUTPUT_FOLDER="$PROJ_DIR/__results/centrolign_align/pedigree_2generations.correction/$CHR_ID"
OUTPUT_CIGAR="$OUTPUT_FOLDER/$CHR_ID.cigar.txt"
OUTPUT_TSV="$OUTPUT_FOLDER/$CHR_ID.results.tsv"
OUTPUT_TSV_FILTERED="$OUTPUT_FOLDER/$CHR_ID.results.filtered.tsv"

# Tools & Environment
CENTROLIGN="$PROJ_DIR/centrolign/bin/centrolign"
PYTHON_SCRIPT="$PROJ_DIR/__scripts/assembly/cigar_to_tsv.only_mutations.py"
CONDA_BASE="/cvmfs/software.metacentrum.cz/conda/envs/miniforge3-25.3.1-0"
ENV_PATH="/storage/praha5-elixir/projects/bioinf-fi/polakova/apps/miniconda3/envs/bioinf"

# --- ENVIRONMENT SETUP ---
export LD_LIBRARY_PATH="$PROJ_DIR/centrolign/lib:$LD_LIBRARY_PATH"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$ENV_PATH"

mkdir -p "$OUTPUT_FOLDER"
cd "$SCRATCHDIR" || exit 1

# --- DYNAMIC ID SETUP ---
ID_REF=$(basename "$FASTA_REF" .fasta)
ID_MUT=$(basename "$FASTA_MUT" .fasta)

echo "Processing Reference: $ID_REF and Mutation: $ID_MUT"

# --- OFFSET EXTRACTION ---
echo "Extracting offsets from BED files..."
OFF_REF=$(awk -v target="$ID_REF" '$4 == target {print $2; exit}' "$BED_REF")
OFF_MUT=$(awk -v target="$ID_MUT" '$4 == target {print $2; exit}' "$BED_MUT")

if [[ -z "$OFF_REF" || -z "$OFF_MUT" ]]; then
    echo "ERROR: Could not find offsets for $ID_REF or $ID_MUT!"
    exit 1
fi

echo "Offsets detected: REF=$OFF_REF, MUT=$OFF_MUT"

# --- ALIGNMENT ---
echo "Running Centrolign..."
cat "$FASTA_REF" "$FASTA_MUT" > joined.fasta
"$CENTROLIGN" joined.fasta > "$OUTPUT_CIGAR" 2> "$OUTPUT_FOLDER/$CHR_ID.engine.log"

# --- CONVERSION TO TSV (Updated with Argparse) ---
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
    if [[ -f "$PROBLEMATIC_REF" ]]; then
        echo "Filtering Stage 1: Reference problematic regions..."
        REF_PROB_NAME=$(awk 'NR==1 {print $1}' "$PROBLEMATIC_REF")
        awk -v bname="$REF_PROB_NAME" 'NR > 1 {print bname"\t"$3"\t"$3"\t"$0}' "$OUTPUT_TSV" > stage1.bed
        bedtools intersect -a stage1.bed -b "$PROBLEMATIC_REF" -v > stage1_clean.bed
        awk 'BEGIN{FS=OFS="\t"} {for (i=4; i<=NF; i++) printf $i (i==NF?ORS:OFS)}' stage1_clean.bed > stage1.tsv
    else
        tail -n +2 "$OUTPUT_TSV" > stage1.tsv
    fi

    # Stage 2: Mutation problematic regions
    if [[ -f "$PROBLEMATIC_MUT" ]]; then
        echo "Filtering Stage 2: Mutation problematic regions..."
        MUT_PROB_NAME=$(awk 'NR==1 {print $1}' "$PROBLEMATIC_MUT")
        awk -v bname="$MUT_PROB_NAME" 'BEGIN{FS=OFS="\t"} {print bname"\t"$4"\t"$4"\t"$0}' stage1.tsv > stage2.bed
        bedtools intersect -a stage2.bed -b "$PROBLEMATIC_MUT" -v > stage2_clean.bed
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

    rm -f stage*.bed stage1_clean.bed stage2_clean.bed stage3_clean.bed stage1.tsv stage2.tsv
    echo "Filtering complete. Output: $OUTPUT_TSV_FILTERED"
fi

echo "Job completed successfully."