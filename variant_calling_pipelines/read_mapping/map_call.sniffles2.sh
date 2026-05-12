#!/usr/bin/env bash
#PBS -N sniffles2_variant_calling
#PBS -l select=1:ncpus=16:mem=64gb:scratch_local=100gb
#PBS -l walltime=02:00:00
#PBS -j oe

# ==============================================================================
# DOCUMENTATION:
# ==============================================================================
# Maps ONT reads to a parental reference and calls structural variants using
# Sniffles2. If a pre-computed BAM already exists it is reused, skipping
# the alignment step.
#
# 1. OFFSET EXTRACTION:
#    Extracts genomic offset from BED file by matching REF_ID against column 1.
#
# 2. FILTERING STRATEGY:
#    Stage 1: bcftools restricts calls to centromeric regions (BED_REGIONS)
#             and applies quality filters.
#    Stage 2: bedtools removes variants in problematic regions (BED_PROBLEMATIC).
#    Stage 3: bedtools removes variants in genomic gaps (BED_GAPS).
#
# USAGE:
#    1. Copy config.example.sh to config.sh and fill in your paths.
#    2. Submit with: qsub sniffles_call.sh
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
if [ -z "$SAMPLE_NAME" ] || [ -z "$OUT" ]; then
    echo "ERROR: Missing required variables SAMPLE_NAME or OUT in config.sh."
    exit 1
fi

mkdir -p "$OUT/$SAMPLE_NAME"

for TOOL in minimap2 samtools bcftools bedtools sniffles; do
    if ! command -v "$TOOL" &> /dev/null; then
        echo "ERROR: $TOOL not found in PATH."
        exit 1
    fi
done

# ==============================================================================
# 2. ENVIRONMENT SETUP
# ==============================================================================
# shellcheck source=/dev/null
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV" || { echo "ERROR: Failed to activate conda environment: $CONDA_ENV"; exit 1; }

# ==============================================================================
# 3. PREPARE INPUTS & OFFSETS (SCRATCH)
# ==============================================================================
echo "Moving to scratch: $SCRATCHDIR"
cd "$SCRATCHDIR" || exit 1
mkdir -p tmp && export TMPDIR="$SCRATCHDIR/tmp"

# Localize reference to scratch
cp "$REF" ./ref.fasta
samtools faidx ./ref.fasta

# Extract genomic offset for the reference sequence (optional)
REF_OFFSET=0
if [[ -n "$BED_OFFSETS" && -f "$BED_OFFSETS" ]]; then
    REF_OFFSET=$(awk -v id="$REF_ID" '$1 == id {print $2; exit}' "$BED_OFFSETS" | tr -d '\r')
    REF_OFFSET=${REF_OFFSET:-0}
    echo "Offset detected: REF=$REF_OFFSET"
else
    echo "BED_OFFSETS not set. Offset defaulting to 0."
fi

# ==============================================================================
# 4. STEP 1: MAPPING (IF BAM NOT PRESENT)
# ==============================================================================
echo "Checking for existing alignment..."
if [ -f "$READS_BAM" ]; then
    echo "BAM exists. Localizing..."
    cp "$READS_BAM" ./mapped_sorted.bam || { echo "ERROR: Failed to copy BAM."; exit 1; }

    # Re-index if index is missing or older than the BAM
    if [ ! -f "${BAM_DIR}.bai" ] || \
       [ "$READS_BAM" -nt "${BAM_DIR}.bai" ]; then
        echo "BAM index missing or outdated. Indexing on scratch..."
        samtools index ./mapped_sorted.bam || { echo "ERROR: Failed to index BAM."; exit 1; }
    else
        cp "${BAM_DIR}.bai" ./mapped_sorted.bam.bai || { echo "ERROR: Failed to copy BAM index."; exit 1; }
    fi
else
    echo "Running Minimap2 alignment..."

    # Convert BAM to FASTQ first (avoids pipe issues with coordinate-sorted BAM)
    echo "Converting BAM to FASTQ..."
    samtools fastq -@ 8 "$READS_BAM" > reads.fastq \
        || { echo "ERROR: samtools fastq failed."; exit 1; }

    echo "Running Minimap2..."
    minimap2 -ax map-ont -t 16 ref.fasta reads.fastq | \
        samtools sort -@ 8 -o mapped_sorted.bam
    PIPE=(${PIPESTATUS[@]})
    if [ ${PIPE[0]} -ne 0 ] || [ ${PIPE[1]} -ne 0 ]; then
        echo "ERROR: Alignment pipeline failed (minimap2=${PIPE[0]}, sort=${PIPE[1]})."
        exit 1
    fi

    rm -f reads.fastq
    samtools index mapped_sorted.bam || { echo "ERROR: Failed to index BAM."; exit 1; }
fi

# ==============================================================================
# 5. STEP 2: RUN SNIFFLES2
# ==============================================================================
echo "Running Sniffles2 SV calling..."

sniffles \
    --input mapped_sorted.bam \
    --vcf raw.vcf.gz \
    --reference ref.fasta \
    --threads 16 \
    --minsvlen 50 \
    --mapq 20 \
    --output-rnames

SNIFFLES_EXIT=$?

if [ $SNIFFLES_EXIT -ne 0 ]; then
    echo "ERROR: Sniffles2 execution failed with exit code $SNIFFLES_EXIT."
    exit 1
fi

# ==============================================================================
# 6. FILTER & EXPORT RESULTS
# ==============================================================================
echo "Filtering VCF..."
tabix -f -p vcf raw.vcf.gz

# Apply quality filters and restrict to centromeric regions
bcftools view \
    -i 'FILTER="PASS" && QUAL>=20 && INFO/SUPPORT>=10 && INFO/VAF>=0.8' \
    -R "$BED_REGIONS" \
    raw.vcf.gz -O v -o bcftools_filtered.vcf || { echo "ERROR: bcftools filtering failed."; exit 1; }

# Warn if filter produced no variants
if [ ! -s bcftools_filtered.vcf ]; then
    echo "WARNING: bcftools produced an empty VCF. Check filter thresholds or BED_REGIONS."
fi

# --- OPTIONAL: OFFSET CORRECTION ---
if [[ -n "$BED_OFFSETS" && -f "$BED_OFFSETS" ]]; then
    echo "Applying offset correction (REF_OFFSET=$REF_OFFSET)..."
    awk -v r_off="$REF_OFFSET" 'BEGIN {OFS="\t"}
        /^#/ {print $0; next}
        { $2 = $2 + r_off; print $0 }
    ' bcftools_filtered.vcf > offset_corrected.vcf
else
    echo "BED_OFFSETS not set or not found. Skipping offset correction."
    cp bcftools_filtered.vcf offset_corrected.vcf
fi

# --- OPTIONAL STAGE 2: Filter problematic regions ---
if [[ -n "$BED_PROBLEMATIC" && -f "$BED_PROBLEMATIC" ]]; then
    echo "Filtering problematic regions..."
    grep -F "$REF_ID" "$BED_PROBLEMATIC" > ref_probs.bed
    bedtools intersect -header -v -a offset_corrected.vcf -b ref_probs.bed > stage2.vcf
else
    echo "BED_PROBLEMATIC not set or not found. Skipping problematic region filter."
    cp offset_corrected.vcf stage2.vcf
fi

# --- OPTIONAL STAGE 3: Filter genomic gaps ---
if [[ -n "$BED_GAPS" && -f "$BED_GAPS" ]]; then
    echo "Filtering genomic gaps..."
    grep -F "$REF_ID" "$BED_GAPS" > ref_gaps.bed
    bedtools intersect -header -v -a stage2.vcf -b ref_gaps.bed > "${SAMPLE_NAME}.filtered.vcf"
else
    echo "BED_GAPS not set or not found. Skipping gap filter."
    cp stage2.vcf "${SAMPLE_NAME}.filtered.vcf"
fi

# Export results
cp raw.vcf.gz "$OUT/$SAMPLE_NAME/${SAMPLE_NAME}.raw.vcf.gz" || exit 1
cp raw.vcf.gz.tbi "$OUT/$SAMPLE_NAME/${SAMPLE_NAME}.raw.vcf.gz.tbi" || exit 1
cp "${SAMPLE_NAME}.filtered.vcf" "$OUT/$SAMPLE_NAME/" || exit 1

# Clean up all intermediate scratch files
rm -f ref.fasta ref.fasta.fai reads.fastq mapped_sorted.bam mapped_sorted.bam.bai \
       raw.vcf.gz raw.vcf.gz.tbi \
       bcftools_filtered.vcf offset_corrected.vcf \
       ref_probs.bed ref_gaps.bed \
       stage2.vcf "${SAMPLE_NAME}.filtered.vcf"
rm -rf tmp

echo "Results saved to: $OUT/$SAMPLE_NAME"