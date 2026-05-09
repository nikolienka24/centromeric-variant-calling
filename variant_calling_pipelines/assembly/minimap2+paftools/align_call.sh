#!/usr/bin/env bash
#PBS -N minimap2_paftools_full_pipeline
#PBS -l select=1:ncpus=8:mem=16gb:scratch_local=200gb
#PBS -l walltime=00:45:00
#PBS -j oe

# ==============================================================================
# DOCUMENTATION:
# ==============================================================================
# 1. OFFSET EXTRACTION:
#    Extracts genomic offsets from BED files by matching the FASTA filename
#    against the 1st column of the BED.
#
# 2. DOUBLE-SIDED FILTERING STRATEGY (bedtools):
#    Stage 1: Filter Ref_Pos against PROBLEMATIC and GAPS (BED_PROBLEMATIC, BED_GAPS).
#    Stage 2: Filter Qry_Pos against PROBLEMATIC and GAPS (BED_PROBLEMATIC, BED_GAPS).
#
# 3. COORDINATE MAPPING:
#    VCF Column 2  = Reference Position.
#    INFO:QSTART   = Query Position.
#
# USAGE:
#    1. Copy config.example.sh to config.sh and fill in your paths.
#    2. Submit with: qsub align_call.sh
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
if [ -z "$CHR_ID" ] || [ -z "$REF_ID" ] || [ -z "$QRY_ID" ] || [ -z "$OUT" ]; then
    echo "ERROR: Missing required variables CHR_ID, REF_ID, QRY_ID, or OUT in config.sh."
    exit 1
fi

mkdir -p "$OUT"

# ==============================================================================
# 2. ENVIRONMENT SETUP
# ==============================================================================
# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$ENV_PATH" || { echo "ERROR: Failed to activate conda environment: $ENV_PATH"; exit 1; }

# ==============================================================================
# 3. PREPARE OFFSETS
# ==============================================================================
echo "Moving to scratch: $SCRATCHDIR"
cd "$SCRATCHDIR" || exit 1

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
minimap2 -t 8 -cx asm5 --cs "$FASTA_REF" "$FASTA_QRY" > alignment.paf 2> alignment.log

if [[ ! -s alignment.paf ]]; then
    echo "ERROR: minimap2 produced an empty PAF file. Check alignment.log"
    exit 1
fi

paftools.js call -f "$FASTA_REF" alignment.paf > variants_raw.vcf

# --- OFFSET CORRECTION (REF POS + QSTART in INFO field) ---
# Adds the genomic offset to the reference position and the QSTART tag in the INFO field
awk -v r_off="$REF_OFFSET" -v q_off="$QRY_OFFSET" 'BEGIN {OFS="\t"}
    /^#/ {print $0; next}
    {
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
grep -F "$REF_ID" "$BED_PROBLEMATIC" > r_probs.bed
grep -F "$REF_ID" "$BED_GAPS" > r_gaps.bed

bedtools intersect -header -v -a variants_corrected.vcf -b r_probs.bed | \
bedtools intersect -header -v -a stdin -b r_gaps.bed > tmp_ref_filtered.vcf

# B. Query-side filtering (extracting coordinates from INFO:QSTART)
# Create a temporary BED stream: chromosome = QRY_ID, position = QSTART
grep -v "^#" tmp_ref_filtered.vcf | awk -v qry_id="$QRY_ID" 'BEGIN {OFS="\t"} {
    if (match($8, /QSTART=[0-9]+/)) {
        qpos = substr($8, RSTART+7, RLENGTH-7);
        print qry_id, qpos, qpos+1, $0
    }
}' > query_coords_tmp.bed

grep -F "$QRY_ID" "$BED_PROBLEMATIC" > q_probs.bed
grep -F "$QRY_ID" "$BED_GAPS" > q_gaps.bed

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
cp alignment.paf "$OUT/${CHR_ID}.paf"
cp variants_raw.vcf "$OUT/${CHR_ID}_raw.vcf"
cp variants_final.vcf "$OUT/${CHR_ID}.filtered.vcf"
cp alignment.log "$OUT/"

# Clean up all intermediate scratch files
rm -f alignment.paf variants_raw.vcf variants_corrected.vcf \
       r_probs.bed r_gaps.bed q_probs.bed q_gaps.bed \
       query_coords_tmp.bed tmp_ref_filtered.vcf filtered_body.txt

echo "Job finished successfully. Results exported to $OUT"