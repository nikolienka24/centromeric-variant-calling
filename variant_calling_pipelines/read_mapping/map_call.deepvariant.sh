#!/usr/bin/env bash
#PBS -N deepvariant_variant_calling
#PBS -l select=1:ncpus=16:mem=64gb:scratch_local=200gb
#PBS -l walltime=02:00:00
#PBS -j oe

# ==============================================================================
# DOCUMENTATION:
# ==============================================================================
# Maps ONT reads to a parental reference and calls variants using DeepVariant.
# If a pre-computed BAM already exists it is reused, skipping the alignment step.
#
# 1. OFFSET EXTRACTION:
#    Extracts genomic offset from BED file by matching REF_ID against column 1.
#
# 2. FILTERING STRATEGY:
#    Stage 1: bcftools restricts calls to centromeric regions (BED_REGIONS).
#    Stage 2: bedtools removes variants in problematic regions (BED_PROBLEMATIC).
#    Stage 3: bedtools removes variants in genomic gaps (BED_GAPS).
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
if [ -z "$SAMPLE_NAME" ] || [ -z "$OUT" ]; then
    echo "ERROR: Missing required variables SAMPLE_NAME or OUT in config.sh."
    exit 1
fi

mkdir -p "$OUT/$SAMPLE_NAME"

# Validate required tools and container
for TOOL in minimap2 samtools bcftools bedtools singularity; do
    if ! command -v "$TOOL" &> /dev/null; then
        echo "ERROR: $TOOL not found in PATH."
        exit 1
    fi
done

if [ ! -f "$CONTAINER_IMG" ]; then
    echo "ERROR: DeepVariant Singularity image not found at $CONTAINER_IMG"
    exit 1
fi

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

# Extract genomic offset for the reference sequence
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
FINAL_BAM_NAME=$READS_BAM

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
    samtools fastq "$READS_BAM" | \
        minimap2 -ax map-ont -t 16 ref.fasta - | \
        samtools sort -@ 8 -o mapped_sorted.bam
    # Check exit codes of all pipe stages
    if [ ${PIPESTATUS[0]} -ne 0 ] || [ ${PIPESTATUS[1]} -ne 0 ] || [ ${PIPESTATUS[2]} -ne 0 ]; then
        echo "ERROR: Alignment pipeline failed (fastq=${PIPESTATUS[0]}, minimap2=${PIPESTATUS[1]}, sort=${PIPESTATUS[2]})."
        exit 1
    fi
    samtools index mapped_sorted.bam || { echo "ERROR: Failed to index BAM."; exit 1; }
fi

# ==============================================================================
# 5. STEP 2: RUN DEEPVARIANT
# ==============================================================================
echo "Running DeepVariant variant calling..."

singularity exec --bind "$SCRATCHDIR:/data" "$CONTAINER_IMG" \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=ONT_R104 \
    --ref="/data/ref.fasta" \
    --reads="/data/mapped_sorted.bam" \
    --output_vcf="/data/raw.vcf.gz" \
    --num_shards=16 \
    --logging_dir="/data"

DV_EXIT=$?

if [ $DV_EXIT -ne 0 ]; then
    echo "ERROR: DeepVariant execution failed with exit code $DV_EXIT."
    exit 1
fi

# ==============================================================================
# 6. FILTER & EXPORT RESULTS
# ==============================================================================
echo "Filtering VCF..."
tabix -p vcf raw.vcf.gz

# Apply quality filters and restrict to centromeric regions
bcftools view \
    -i 'FILTER="PASS" && QUAL>=20 && FORMAT/DP>=20 && (FORMAT/AD[0:1])/(FORMAT/DP)>=0.8' \
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
rm -f ref.fasta ref.fasta.fai mapped_sorted.bam mapped_sorted.bam.bai \
       raw.vcf.gz raw.vcf.gz.tbi \
       bcftools_filtered.vcf offset_corrected.vcf \
       ref_probs.bed ref_gaps.bed \
       stage2.vcf "${SAMPLE_NAME}.filtered.vcf"
rm -rf tmp

echo "Results saved to: $OUT/$SAMPLE_NAME"