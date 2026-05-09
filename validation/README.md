# Validation

This section contains the consensus variant calling logic used to cross-validate mutation calls across three independent bioinformatic tools: **Stretcher**, **Minimap2/Paftools**, and **Centrolign**. A variant is included in the final output only if it is independently supported by at least two tools.

---

## Directory Structure

```
validation/
└── intersect_variants.py   # Pairwise tool comparison and consensus variant calling
```

---

## `intersect_variants.py`

Loads variant calls from three tools, performs all pairwise comparisons using position proximity and sequence similarity, and assembles a final consensus set stratified by how many tools agree.

**Requirements:** `pandas`, `edlib`

**Usage:**
```bash
python intersect_variants.py \
    -c centrolign_variants.tsv \
    -s stretcher_variants.bedpe \
    -m minimap2_variants.vcf \
    -o /path/to/output/dir
```

| Argument | Description |
|----------|-------------|
| `-c` / `--centrolign` | Path to Centrolign TSV output file |
| `-s` / `--stretcher` | Path to Stretcher BEDPE output file |
| `-m` / `--minimap` | Path to Minimap2/Paftools VCF output file |
| `-o` / `--output_dir` | Directory where consensus TSV files will be saved |

---

## Input Files

The script expects one file per tool in its native output format:

| Tool | Format | Key Columns |
|------|--------|-------------|
| Centrolign | TSV | `Ref_Pos`, `Mut_Pos`, `Ref_Base`, `Alt_Base` |
| Stretcher | TSV (BEDPE-like) | `start1`, `start2`, `sequence1`, `sequence2` (0-based → converted to 1-based internally) |
| Minimap2/Paftools | VCF | Standard VCF columns (comment lines skipped) |

---

## Method

**1. Normalization**
All REF/ALT sequences are uppercased and stripped of common leading/trailing bases (left- and right-trimming) to put variants in a canonical form before comparison.

**2. Pairwise matching**
All three tool pairs are compared independently (Stretcher↔Minimap2, Stretcher↔Centrolign, Minimap2↔Centrolign). For each variant in tool A, candidates in tool B within a position window are considered:
- **±100 bp** window for short variants (≤50 bp)
- **±20,000 bp** window for structural variants (>50 bp)

Candidates are scored using a two-step similarity metric: length-ratio similarity followed by edit distance (Levenshtein via `edlib`, HW/infix mode). A pair is accepted if similarity ≥ 0.90.

**3. Consensus assembly**
Matched pairs are assembled into the final consensus in priority order:

| Priority | Condition | `tools_count` |
|----------|-----------|---------------|
| 1st | Triplet — same variant matched across all three tools | `3` |
| 2nd | Doublet — Stretcher + Minimap2 only | `2` |
| 3rd | Doublet — Stretcher + Centrolign only | `2` |
| 4th | Doublet — Minimap2 + Centrolign only | `2` |

Each variant index is marked as used after being assigned to a consensus entry to prevent double-counting.

---

## Outputs

Both files are saved to the directory specified by `--output_dir`:

| File | Description |
|------|-------------|
| `final_consensus_2_of_3.tsv` | All consensus variants supported by ≥2 tools |
| `final_consensus_3_of_3.tsv` | Subset of variants supported by all 3 tools |

**Output columns:**

| Column | Description |
|--------|-------------|
| `ref_pos_stretcher` | Reference position from Stretcher (`null` if not called) |
| `alt_pos_stretcher` | Alt position from Stretcher |
| `ref_pos_minimap` | Reference position from Minimap2 |
| `ref_pos_centrolign` | Reference position from Centrolign |
| `alt_pos_centrolign` | Alt position from Centrolign |
| `ref_seq` | Normalized REF sequence |
| `alt_seq` | Normalized ALT sequence |
| `tools_count` | Number of tools supporting this variant (2 or 3) |
| `sim_S_M` | Similarity score: Stretcher vs Minimap2 |
| `sim_S_C` | Similarity score: Stretcher vs Centrolign |
| `sim_M_C` | Similarity score: Minimap2 vs Centrolign |

Variants are sorted by the earliest reference position found across the three tools.