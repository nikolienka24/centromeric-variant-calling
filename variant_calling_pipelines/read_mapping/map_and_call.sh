#!/usr/bin/env bash
#PBS -N deepvariant_variant_calling
#PBS -l select=1:ncpus=16:mem=64gb:scratch_local=200gb
#PBS -l walltime=02:00:00
#PBS -j oe

# ==========================================
# 1. INPUT ARGUMENTS & USAGE EXAMPLE
# ==========================================
# Example Run Command:
# qsub script.sh -v SAMPLE="PAN027",CHR="chr22",PARENT="paternal",OUT="/path/to/results"
#
# Arguments provided via -v (PBS variables):
# SAMPLE - Sample ID (e.g., PAN027)
# CHR    - Chromosome (e.g., chr22)
# PARENT - Parent haplotype (paternal/maternal)
# OUT    - Path to the final output directory

if [ -z "$SAMPLE" ] || [ -z "$CHR" ] || [ -z "$PARENT" ] || [ -z "$OUT" ]; then
    echo "Error: Missing required variables SAMPLE, CHR, PARENT, or OUT."
    echo "Usage: qsub $0 -v SAMPLE=\"S\",CHR=\"C\",PARENT=\"P\",OUT=\"/path\""
    exit 1
fi

# Map variables for clarity
SAMPLE_NAME="${SAMPLE}.${CHR}.${PARENT}"
OUT_DIR="$OUT/$SAMPLE_NAME"
mkdir -p "$OUT_DIR"

# ==========================================
# 2. ENVIRONMENT SETUP (USER DEFINED)
# ==========================================
# >>> ADD YOUR ENVIRONMENT SETUP HERE <<<
PROJECT_DIR="/storage/praha5-elixir/projects/bioinf-fi/polakova/BP"
CONDA_BASE="/cvmfs/software.metacentrum.cz/conda/envs/miniforge3-25.3.1-0"
CONDA_ENV="$PROJECT_DIR/apps/miniconda3/envs/bioinf"

source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

# Tool and Container paths
CONTAINER_IMG="${PROJECT_DIR}/__scripts/deepvariant/deepvariant_1.10.0.sif"

# Validation: Check if required tools/images are accessible
for TOOL in minimap2 samtools bcftools singularity; do
    if ! command -v $TOOL &> /dev/null; then
        echo "Error: $TOOL not found in PATH."
        exit 1
    fi
done

if [ ! -f "$CONTAINER_IMG" ]; then
    echo "Error: DeepVariant Singularity image not found at $CONTAINER_IMG"
    exit 1
fi

# ==========================================
# 3. PREPARE INPUTS (SCRATCH)
# ==========================================
echo "Moving to scratch: $SCRATCHDIR"
cd "$SCRATCHDIR" || exit 1
mkdir -p tmp && export TMPDIR="$SCRATCHDIR/tmp"

# Dynamic Input Mapping (Update these paths to match your folder structure)
REF="$PROJECT_DIR/__data/pedigree/assembliesv1.1/PAN011_hap2_HiFi_element_final_XY_hap2.polished.fasta"
READS_BAM="$PROJECT_DIR/__data/pedigree/map_ont/PAN_realign.${SAMPLE_NAME}.bam"
BED_FILTER="$PROJECT_DIR/__data/pedigree/annotations_v1.1/centromeres_bed/PAN011_hap2_HiFi_element_final_XY_hap2.polished.cenSat.active_hor_merged_clean.bed"

# Localize Reference
cp "$REF" ./ref.fasta
samtools faidx ./ref.fasta

# ==========================================
# 4. STEP 1: MAPPING (IF BAM NOT PRESENT)
# ==========================================
FINAL_BAM_NAME="PAN_realign.${SAMPLE_NAME}.to.PAN011_hap2.bam"

echo "Checking for existing alignment..."
if [ -f "$PROJECT_DIR/__data/pedigree/map_ont_parent_to_GP/$FINAL_BAM_NAME" ]; then
    echo "BAM exists. Localizing..."
    cp "$PROJECT_DIR/__data/pedigree/map_ont_parent_to_GP/$FINAL_BAM_NAME" ./mapped_sorted.bam
    cp "$PROJECT_DIR/__data/pedigree/map_ont_parent_to_GP/${FINAL_BAM_NAME}.bai" ./mapped_sorted.bam.bai
else
    echo "Running Minimap2 alignment..."
    samtools fastq "$READS_BAM" | minimap2 -ax map-ont -t 16 ref.fasta - | samtools sort -@ 8 -o mapped_sorted.bam
    samtools index mapped_sorted.bam
fi

# ==========================================
# 5. STEP 2: RUN DEEPVARIANT
# ==========================================
echo "Running DeepVariant calling..."

singularity exec --bind "$SCRATCHDIR:/data" "$CONTAINER_IMG" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=ONT_R104 \
  --ref="/data/ref.fasta" \
  --reads="/data/mapped_sorted.bam" \
  --output_vcf="/data/raw.vcf.gz" \
  --num_shards=16 \
  --logging_dir="/data"

if [ $? -eq 0 ]; then
    echo "DeepVariant successful. Filtering VCF..."

    tabix -p vcf raw.vcf.gz

    # Apply quality and region filters
    bcftools view -i 'FILTER="PASS" && QUAL>=20 && FORMAT/DP>=20 && (FORMAT/AD[0:1])/(FORMAT/DP)>=0.8' \
        -R "$BED_FILTER" \
        raw.vcf.gz -O v -o "${SAMPLE_NAME}.filtered.vcf"

    # Save raw and filtered results
    # gunzip -c raw.vcf.gz > "${SAMPLE_NAME}.raw.vcf"
    cp "${SAMPLE_NAME}.raw.vcf" "$OUT_DIR/"
    cp "${SAMPLE_NAME}.filtered.vcf" "$OUT_DIR/"

    echo "Results saved to: $OUT_DIR"
else
    echo "Error: DeepVariant execution failed."
    exit 1
fi