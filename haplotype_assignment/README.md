# Haplotype Assignment

This section covers haplotype assignment for centromeric regions using sequence-level alignment with MashMap. Two complementary strategies are used: direct alignment of **centromere sequences**, and alignment of **flanking regions** surrounding the centromere.

---

## Directory Structure

```
haplotype_assignment/
├── mashmap.centromeres.sh   # Assign haplotype by aligning centromere FASTA to H1/H2
└── mashmap.flanks.sh        # Assign haplotype by aligning left/right flanks to H1/H2
```

---

## 1. `mashmap.centromeres.sh`

Aligns an extracted centromere query sequence against both haplotypes (H1 and H2) simultaneously using MashMap. The haplotype producing the best-identity alignment is selected as the winner.

**Requirements:** `mashmap

**Usage (PBS/qsub):**
```bash
qsub mashmap.centromeres.sh \
    -v QRY="query.fasta",H1="hap1.fasta",H2="hap2.fasta",OUT="/path/to/output"
```

| Variable | Description |
|----------|-------------|
| `QRY` | Extracted centromere FASTA for the query sample |
| `H1` | Centromere FASTA for Haplotype 1 of the reference panel |
| `H2` | Centromere FASTA for Haplotype 2 of the reference panel |
| `OUT` | Output directory |

**Method:**
- H1 and H2 FASTAs are merged into a single target file with `HAP1_` / `HAP2_` header prefixes
- MashMap runs with `--noSplit -n 2 --pi 90` to retrieve hits against both haplotypes
- The single best hit (highest identity, column 10) is selected and saved

**Output:** `<OUT>/mashmap_best_<query_name>.out` — single-line best-hit MashMap record

---

## 2. `mashmap.flanks.sh`

Assigns a haplotype by aligning the **left and right flanking sequences** of the centromere against the corresponding flanks of H1 and H2. This is useful when the centromere itself is too repetitive for reliable mapping.

**Requirements:** `mashmap`

**Usage (PBS/qsub):**
```bash
qsub mashmap.flanks.sh \
    -v CHR="chrX",QRY="PAN027.paternal",H1="PAN028.h1",H2="PAN028.h2",IN="/path/flanks",OUT="/path/results"
```

| Variable | Description |
|----------|-------------|
| `CHR` | Chromosome name (e.g., `chrX`) — used to name the output subdirectory |
| `QRY` | Filename prefix for the query flank files |
| `H1` | Filename prefix for Haplotype 1 flank files |
| `H2` | Filename prefix for Haplotype 2 flank files |
| `IN` | Input directory containing all `*_left_flank.fasta` / `*_right_flank.fasta` files |
| `OUT` | Base output directory (`<OUT>/<CHR>/` will be created) |

**Expected input filenames** (inside `IN` directory):
```
<QRY>_left_flank.fasta    <QRY>_right_flank.fasta
<H1>_left_flank.fasta     <H1>_right_flank.fasta
<H2>_left_flank.fasta     <H2>_right_flank.fasta
```

**Method:**
- Left and right flanks are processed independently
- H1/H2 targets are merged with `H1_L_` / `H2_L_` and `H1_R_` / `H2_R_` header prefixes
- MashMap runs with `-n 2 --pi 95` for each flank
- The best hit per flank is selected by sorting on column 11 and taking the top result
- Both best hits are combined into a single summary file

**Outputs** (inside `<OUT>/<CHR>/`):

| File | Description |
|------|-------------|
| `<QRY>_best_hits_combined.out` | Best MashMap hit for left and right flank combined |
| `<QRY>_left.out` | Full MashMap output for the left flank |
| `<QRY>_right.out` | Full MashMap output for the right flank |

---

> **Note:** Edit the `ENVIRONMENT SETUP` section in each script to match your cluster's conda paths before submitting.