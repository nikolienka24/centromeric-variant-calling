#!/bin/bash
#PBS -N mashmap_simple_fasta
#PBS -l select=1:ncpus=16:mem=64gb:scratch_local=200gb
#PBS -l walltime=04:00:00
#PBS -j oe

# ==========================================
# 1. INPUT ARGUMENTS & USAGE EXAMPLE
# ==========================================
# Example Run:
# qsub mashmap.centromeres.sh "query_extracted.fasta" "hap1_extracted.fasta" "hap2_extracted.fasta" "/path/to/output"

if [ "$#" -lt 4 ]; then
    echo "Error: Missing arguments."
    echo "Usage: qsub $0 <query.fasta> <hap1.fasta> <hap2.fasta> <out_dir>"
    exit 1
fi

# Inputs (already extracted centromeres as FASTA files)
QUERY_FA=$1
HAP1_FA=$2
HAP2_FA=$3
OUT_DIR=$4

mkdir -p "$OUT_DIR"

# Construct a clean filename based on the query input name
QUERY_NAME=$(basename "$QUERY_FA" .fasta)
FINAL_OUT="$OUT_DIR/mashmap_best_${QUERY_NAME}.out"

# ==========================================
# 2. ENVIRONMENT SETUP
# ==========================================
# source /path/to/conda/etc/profile.d/conda.sh
# conda activate bioinf

# ==========================================
# 3. PREPARE TARGETS (SCRATCH)
# ==========================================
echo "Moving to scratch: $SCRATCHDIR"
cd "$SCRATCHDIR" || exit 1

# Combine Haplotypes into one reference file
# We add prefixes to the headers to know which is which in the results
sed "s/>/>HAP1_/" "$HAP1_FA" > targets.fasta
sed "s/>/>HAP2_/" "$HAP2_FA" >> targets.fasta

# ==========================================
# 4. RUN MASHMAP
# ==========================================
echo "Running MashMap alignment..."

# --noSplit: Treat the query as one unit
# -n 2: Get hits for both haplotypes
mashmap -r targets.fasta \
        -q "$QUERY_FA" \
        -t 16 \
        --noSplit \
        -n 2 \
        --pi 90 \
        -o mashmap_raw.out

# ==========================================
# 5. SELECT BEST HIT & SAVE
# ==========================================
if [ -s mashmap_raw.out ]; then
    # Sort by column 10 (Identity %) numerically, highest first
    sort -k10,10nr mashmap_raw.out | head -n 1 > best_match.out

    WINNER=$(awk '{print $6}' best_match.out)
    ID_VAL=$(awk '{print $10}' best_match.out)
    echo "Alignment complete. Winner: $WINNER ($ID_VAL% identity)"

    cp best_match.out "$FINAL_OUT"
else
    echo "Error: No alignment found between query and haplotypes."
    exit 1
fi

echo "Result saved to: $FINAL_OUT"