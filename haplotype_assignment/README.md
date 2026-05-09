# Haplotype Assignment

This section covers haplotype assignment for centromeric regions using sequence-level alignment with MashMap. Two complementary strategies are used: direct alignment of **centromeric sequences**, and alignment of **flanking regions** surrounding the centromere.

---

## Directory Structure

```
haplotype_assignment/
├── mashmap.centromeres.sh   # Assign haplotype by aligning centromere FASTA to H1/H2
├── mashmap.flanks.sh        # Assign haplotype by aligning left/right flanks to H1/H2
└── config.example.sh        # Template configuration file — copy to config.sh and fill in your paths
```

---

## 1. `mashmap.centromeres.sh`

Aligns an extracted centromere query sequence against both haplotypes (H1 and H2) simultaneously using MashMap. The haplotype producing the best-identity alignment is selected as the winner.

**Requirements:** `mashmap`

> **Note:** All input paths and sample identifiers are configured via `config.sh`. Copy `config.example.sh` to `config.sh` and fill in your paths before submitting.

**Usage (PBS/qsub):**
```bash
qsub mashmap.centromeres.sh
```

**Method:**
- H1 and H2 FASTAs are merged into a single target file with `HAP1_` / `HAP2_` header prefixes
- MashMap runs with `--noSplit -n 2 --pi 90` to retrieve hits against both haplotypes
- The single best hit (highest identity, column 10) is selected and saved

**Output:** `<OUT>/mashmap_best_<query_name>.out` — single-line best-hit MashMap record

---

## 2. `mashmap.flanks.sh`

Assigns a haplotype by aligning the **left and right flanking sequences** of the centromere against the corresponding flanks of H1 and H2. This is useful when the centromere itself is too repetitive for reliable mapping.

**Requirements:** `mashmap`

> **Note:** All input paths and sample identifiers are configured via `config.sh`. Copy `config.example.sh` to `config.sh` and fill in your paths before submitting.

**Usage (PBS/qsub):**
```bash
qsub mashmap.flanks.sh
```

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