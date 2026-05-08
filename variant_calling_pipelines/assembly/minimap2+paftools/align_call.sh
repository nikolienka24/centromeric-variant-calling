#!/usr/bin/env bash
#PBS -N minimap2_paftools_full_pipeline
#PBS -q elixircz@pbs-m1.metacentrum.cz
#PBS -l select=1:ncpus=8:mem=16gb:scratch_local=200gb
#PBS -l walltime=00:45:00
#PBS -j oe

# ==============================================================================
# 1. INPUT ARGUMENTS & USAGE EXAMPLE
# ==============================================================================
# Example Run Command:
# qsub script.sh -v CHR_ID="chrX",REF_ID="PAN010",QRY_ID="PAN027",OUT="/path/to/results"
#
# Variables required via -v:
# CHR_ID - Chromosome identifier for output naming
# REF_ID - Reference sequence ID (Gen1)
# QRY_ID - Query sequence ID (Gen2)
# OUT    - Final destination directory

if [ -z "$CHR_ID" ] || [ -z "$REF_ID" ] || [ -z "$QRY_ID" ] || [ -z "$OUT" ]; then
    echo "Error: Missing required variables CHR_ID, REF_ID, QRY_ID, or OUT."
    exit 1
fi

OUT_DIR=$OUT
mkdir -p "$OUT_DIR"

# ==============================================================================
# 2. ENVIRONMENT SETUP
# ==============================================================================
PROJ_DIR="/storage/praha5-elixir/projects/bioinf-fi/polakova/BP"
CONDA_BASE="/cvmfs/software.metacentrum.cz/conda/envs/miniforge3-25.3.1-0"
ENV_PATH="/storage/praha5-elixir/projects/bioinf-fi/polakova/apps/miniconda3/envs/bioinf"

source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$ENV_PATH"

# ==============================================================================
# 3. PREPARE INPUTS & OFFSETS
# ==============================================================================
echo "Moving to scratch: $SCRATCHDIR"
cd "$SCRATCHDIR" || exit 1

# Define annotation and data paths
BED_OFFSETS="$PROJ_DIR/__data/pedigree/annotations_v1.1/centromeres_bed/offsets.bed"
BED_GAPS="$PROJ_DIR/__data/pedigree/annotations_v1.1/gaps_all.bed"
BED_PROBLEMATIC="$PROJ_DIR/__data/pedigree/annotations_v1.1/problematic_all.bed"

INPUT_GEN1="$PROJ_DIR/__data/pedigree/extracted_centromeres_v1.1/${REF_ID}.fasta"
INPUT_GEN2="$PROJ_DIR/__data/pedigree/extracted_centromeres_v1.1/${QRY_ID}.fasta"

# Clean FASTA headers to prevent paftools.js TypeError (truncates to first word)
awk '{print $1}' "$INPUT_GEN1" > ref_clean.fa
awk '{print $1}' "$INPUT_GEN2" > query_clean.fa
CLEAN_REF_ID=$(grep ">" ref_clean.fa | head -n 1 | sed 's/>//')

# Retrieve genomic offsets for both sequences from the BED file
REF_OFFSET=$(awk -v id="$REF_ID" '$1 == id {print $2; exit}' "$BED_OFFSETS" | tr -d '\r')
QRY_OFFSET=$(awk -v id="$QRY_ID" '$1 == id {print $2; exit}' "$BED_OFFSETS" | tr -d '\r')

# Default to 0 if offset is not found
REF_OFFSET=${REF_OFFSET:-0}
QRY_OFFSET=${QRY_OFFSET:-0}

echo "Offsets detected: Reference=$REF_OFFSET, Query=$QRY_OFFSET"

# ==============================================================================
# 4. RUN ALIGNMENT & VARIANT CALLING
# ==============================================================================
echo "Running minimap2 alignment and paftools variant calling..."
minimap2 -t 8 -cx asm5 --cs ref_clean.fa query_clean.fa > alignment.paf 2> alignment.log

paftools.js call -f ref_clean.fa alignment.paf > variants_raw.vcf

# --- OFFSET CORRECTION (REF CHROM/POS + QSTART in INFO field) ---
# Corrects the VCF chromosome name, position, and the QSTART tag within the INFO field
awk -v r_off="$REF_OFFSET" -v q_off="$QRY_OFFSET" -v r_name="$CLEAN_REF_ID" 'BEGIN {OFS="\t"}
    /^#/ {print $0; next}
    {
        $1 = r_name;
        $2 = $2 + r_off;

        # Use regex to find and update QSTART within the INFO column ($8)
        if (match($8, /QSTART=[0-9]+/)) {
            start = RSTART + 7;
            len = RLENGTH - 7;
            val = substr($8, start, len) + q_off;
            $8 = substr($8, 1, start-1) val substr($8, start+len);
        }
        print $0
    }' variants_raw.vcf > variants_corrected.vcf

# ==============================================================================
# 5. DOUBLE-SIDED FILTERING (BEDTOOLS)
# ==============================================================================
echo "Filtering variants based on Reference and Query problematic regions..."

# A. Reference-side filtering (standard VCF columns 1 and 2)
grep -E "$REF_ID" "$BED_PROBLEMATIC" > r_probs.bed
grep -E "$REF_ID" "$BED_GAPS" > r_gaps.bed

bedtools intersect -header -v -a variants_corrected.vcf -b r_probs.bed | \
bedtools intersect -header -v -a stdin -b r_gaps.bed > tmp_ref_filtered.vcf

# B. Query-side filtering (extracting coordinates from INFO:QSTART)
# Create a temporary BED stream: chromosome = QRY_ID, position = QSTART
grep -v "^#" tmp_ref_filtered.vcf | awk 'BEGIN {OFS="\t"} {
    if (match($8, /QSTART=[0-9]+/)) {
        qpos = substr($8, RSTART+7, RLENGTH-7);
        print "'"$QRY_ID"'", qpos, qpos+1, $0
    }
}' > query_coords_tmp.bed

grep -E "$QRY_ID" "$BED_PROBLEMATIC" > q_probs.bed
grep -E "$QRY_ID" "$BED_GAPS" > q_gaps.bed

# Filter out rows where the query position falls into a problematic zone
bedtools intersect -v -a query_coords_tmp.bed -b q_probs.bed | \
bedtools intersect -v -a stdin -b q_gaps.bed | \
cut -f4- > filtered_body.txt

# Reconstruct VCF: combine header with filtered variants
grep "^#" tmp_ref_filtered.vcf > variants_final.vcf
cat filtered_body.txt >> variants_final.vcf

# ==============================================================================
# 6. EXPORT RESULTS
# ==============================================================================
cp alignment.paf "$OUT_DIR/${CHR_ID}.paf"
cp variants_raw.vcf "$OUT_DIR/${CHR_ID}_raw.vcf"
cp variants_final.vcf "$OUT_DIR/${CHR_ID}.filtered.vcf"
cp alignment.log "$OUT_DIR/"

echo "Job finished successfully. Results exported to $OUT_DIR"