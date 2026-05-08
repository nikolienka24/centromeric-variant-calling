# Assembly — Variant Calling Pipeline

This folder contains three independent assembly-to-reference alignment pipelines used to call variants in centromeric regions. Each tool is applied to the same input sequences but uses a different alignment strategy. Their outputs are later cross-validated in the `validation/` step.

---

## Directory Structure

```
assembly/
├── centrolign/
│   ├── align_call.sh                  # Centrolign pairwise alignment + CIGAR parsing + filtering
│   ├── intersect_3gen.py              # De novo filter: subtract shared positions across generations
│   ├── parse.cigar_to_bed.all.py      # Convert full CIGAR output to BED format (all operations)
│   └── parse.cigar_to_tsv.only_mut.py # Convert CIGAR to TSV (mutations only, with offsets)
│
├── minimap2+paftools/
│   ├── align_call.sh                  # Minimap2 asm5 alignment + Paftools VCF calling + filtering
│   └── intersect_3gen.py              # De novo filter: subtract shared POS across generations
│
└── stretcher/
    ├── align_parse.sh                 # EMBOSS Stretcher alignment + straln parsing + filtering
    └── intersect_3gen.py              # De novo filter: subtract shared start1 positions
```

---

## `centrolign/`

### `align_call.sh`

Aligns two centromeric FASTA sequences using Centrolign and converts the CIGAR output to a variant TSV. Applies a three-stage BEDtools filtering pipeline to remove variants in problematic and gap regions.

**Requirements:** `centrolign`, `bedtools`, Python 3, 1 CPU, 64 GB RAM

> **Note:** All input paths (FASTAs, BED files, tool binaries, helper scripts) are **hardcoded**. Edit the `PATH CONFIGURATION` block at the top of the script before submitting.

**Usage (PBS/qsub):**
```bash
qsub align_call.sh
```
*(No `-v` arguments — sample IDs and paths are set directly in the script)*

**Pipeline steps:**
1. Extracts genomic offsets for both sequences from the centromere BED file (matched by FASTA filename against column 4)
2. Concatenates REF and MUT FASTAs and runs Centrolign to produce a CIGAR string
3. Calls `parse.cigar_to_tsv.only_mut.py` to convert the CIGAR into a variant TSV with absolute genomic coordinates
4. Applies three sequential BEDtools filters — reference problematic regions → mutation problematic regions → genomic gaps

**Outputs** (saved to `OUTPUT_FOLDER`):

| File | Description |
|------|-------------|
| `<CHR_ID>.cigar.txt` | Raw Centrolign CIGAR output |
| `<CHR_ID>.results.tsv` | Unfiltered variant TSV |
| `<CHR_ID>.results.filtered.tsv` | Final filtered variants |

---

### `parse.cigar_to_tsv.only_mut.py`

Parses a Centrolign CIGAR string and emits only positions where the reference and query differ (SNPs, insertions, deletions). Applies genomic offsets to report absolute coordinates.

**Requirements:** Python 3 (stdlib only)

**Usage:**
```bash
python parse.cigar_to_tsv.only_mut.py \
    -r ref.fasta -m mut.fasta -c cigar.txt \
    -o variants.tsv --chrom <name> --off_ref <int> --off_mut <int>
```

| Argument | Description |
|----------|-------------|
| `-r` / `--ref` | Reference FASTA |
| `-m` / `--mut` | Query (mutated) FASTA |
| `-c` / `--cigar` | Input CIGAR text file |
| `-o` / `--output` | Output TSV path |
| `--chrom` | Chromosome label for output records |
| `--off_ref` | Genomic offset for the reference sequence (default: 0) |
| `--off_mut` | Genomic offset for the query sequence (default: 0) |

**Output columns:** `Name`, `Type`, `Ref_Start`, `Ref_End`, `Mut_Pos`, `Ref_Base`, `Alt_Base`, `Length`

---

### `parse.cigar_to_bed.all.py`

Converts a full CIGAR string into a BED-like file recording every alignment operation (matches, deletions, insertions, soft-clips, hard-clips). Useful for inspecting the complete alignment structure rather than just mutations.

**Requirements:** Python 3 (stdlib only)

**Usage:**
```bash
python parse.cigar_to_bed.all.py \
    -i cigar.txt -o output.bed --chrom <name>
```

| Argument | Description |
|----------|-------------|
| `-i` / `--input` | Input CIGAR text file |
| `-o` / `--output` | Output BED file |
| `--chrom` | Chromosome name for output records |

**Output columns:** `#Chrom`, `Ref_Start`, `Ref_End`, `Mut_Start`, `Mut_End`, `Op`, `Length`, `Description`

---

### `intersect_3gen.py`

Identifies de novo mutations by subtracting positions found in a second-generation comparison from those in a first-generation comparison. Keeps only `Ref_Start` positions in the GP (generation 1) file that are absent from the GD (generation 2) file.

**Requirements:** `pandas`

**Usage:**
```bash
python intersect_3gen.py \
    -g1 gen1_vs_proband.tsv \
    -g2 gen2_vs_proband.tsv \
    -o de_novo_output.tsv
```

| Argument | Description |
|----------|-------------|
| `-g1` / `--input_gp` | Generation 1 variant TSV (Centrolign output) |
| `-g2` / `--input_gd` | Generation 2 variant TSV (Centrolign output) |
| `-o` / `--output` | Output TSV with unique (de novo) variants |

**Key column:** `Ref_Start`

---

## `minimap2+paftools/`

### `align_call.sh`

Aligns two centromeric assemblies using Minimap2 (`asm5` preset) and calls variants with Paftools. Applies coordinate offset correction to both the VCF `POS` field and the `QSTART` INFO tag, then filters on both reference- and query-side problematic/gap regions.

**Requirements:** `minimap2`, `paftools.js`, `bedtools`, 8 CPUs, 16 GB RAM

> **Note:** Input FASTA paths and BED annotation paths are **hardcoded**. Edit the `ENVIRONMENT SETUP` and `PREPARE INPUTS` sections before submitting.

**Usage (PBS/qsub):**
```bash
qsub align_call.sh \
    -v CHR_ID="chrX",REF_ID="PAN010",QRY_ID="PAN027",OUT="/path/to/results"
```

| Variable | Description |
|----------|-------------|
| `CHR_ID` | Chromosome identifier for output file naming |
| `REF_ID` | Reference sequence ID (generation 1) |
| `QRY_ID` | Query sequence ID (generation 2 / proband) |
| `OUT` | Output directory |

**Pipeline steps:**
1. FASTA headers are truncated to the first word to prevent a `paftools.js` TypeError
2. Genomic offsets for REF and QRY are looked up from `offsets.bed`; defaulting to 0 if not found
3. Minimap2 aligns with `-cx asm5 --cs` and Paftools calls variants into a raw VCF
4. An `awk` script corrects VCF `CHROM`, `POS` (+REF offset), and `QSTART` in INFO (+QRY offset)
5. Double-sided BEDtools filtering: reference-side (CHROM/POS) then query-side (extracted from `QSTART`)

**Outputs:**

| File | Description |
|------|-------------|
| `<CHR_ID>.paf` | Raw Minimap2 alignment |
| `<CHR_ID>_raw.vcf` | Unfiltered Paftools VCF |
| `<CHR_ID>.filtered.vcf` | Final filtered VCF |
| `alignment.log` | Minimap2 stderr log |

---

### `intersect_3gen.py`

Filters de novo variants by subtracting `POS` values found in a generation-3-vs-generation-2 comparison from a generation-1-vs-generation-2 comparison.

**Requirements:** `pandas`

**Usage:**
```bash
python intersect_3gen.py \
    -v1 vs_gen1.filtered.vcf \
    -v2 vs_gen3.filtered.vcf \
    -o /path/to/output_dir/
```

| Argument | Description |
|----------|-------------|
| `-v1` / `--vs_gen1` | Gen1-vs-Gen2 variant file (tab-delimited, must have `POS` column) |
| `-v2` / `--vs_gen2` | Gen3-vs-Gen2 variant file |
| `-o` / `--output_dir` | Output directory — result saved as `final_de_novo.tsv` |

**Key column:** `POS`

---

## `stretcher/`

### `align_parse.sh`

Performs global pairwise alignment using EMBOSS Stretcher, parses the alignment with `straln` to produce a variant BEDPE, and applies BEDtools filtering against problematic and gap regions.

**Requirements:** EMBOSS (`stretcher`, `straln`), `bedtools`, 1 CPU, 4 GB RAM

> **Note:** Input FASTA and BED paths are **hardcoded**. Edit the `PREPARE INPUTS` section before submitting. Gap open/extend penalties are set to 16/4.

**Usage (PBS/qsub):**
```bash
qsub align_parse.sh \
    -v CHR_ID="PAN027.chr3.maternal",S1="PAN027.chr3.maternal",S2="PAN010.chr3.haplotype2",OUT="/path/to/results"
```

| Variable | Description |
|----------|-------------|
| `CHR_ID` | Unique output file identifier |
| `S1` | Exact sequence name for sequence 1 (must match BED and FASTA filenames) |
| `S2` | Exact sequence name for sequence 2 |
| `OUT` | Output directory |

**Pipeline steps:**
1. Genomic offsets for S1 and S2 are extracted from `offsets.bed`
2. Stretcher performs global alignment (`markx0` format, gap open 16, extend 4)
3. `straln` parses the `.aln` file into a variant TSV using sequence names and offsets
4. BEDtools filters variants against combined problematic and gap BED files for both sequences

**Outputs:**

| File | Description |
|------|-------------|
| `<CHR_ID>.aln` | Raw Stretcher alignment file |
| `straln_raw_all/` | Directory with all raw `straln` TSV outputs |
| `<CHR_ID>.filtered.tsv` | Final filtered variant TSV |

---

### `intersect_3gen.py`

Identifies de novo variants by subtracting `start1` positions seen in a generation-2-vs-proband comparison from a generation-1-vs-proband comparison.

**Requirements:** `pandas`

**Usage:**
```bash
python intersect_3gen.py \
    -gp gen1.bedpe \
    -gd gen2.bedpe \
    -o /path/to/output_dir/
```

| Argument | Description |
|----------|-------------|
| `-gp` / `--input_gp` | Generation 1 BEDPE/TSV file |
| `-gd` / `--input_gd` | Generation 2 BEDPE/TSV file |
| `-o` / `--output_dir` | Output directory — result saved as `final_de_novo.tsv` |

**Key column:** `start1`

---

## Comparison of the Three Approaches

| | Centrolign | Minimap2 + Paftools | Stretcher |
|---|---|---|---|
| **Alignment type** | Centromere-aware whole-sequence | Assembly-to-reference (asm5) | Global pairwise (Needleman-Wunsch) |
| **Output format** | CIGAR → TSV | VCF | BEDPE via `straln` |
| **Offset handling** | Per-sequence BED lookup | Per-sequence BED lookup + awk VCF correction | Per-sequence BED lookup |
| **Coordinate key** | `Ref_Start` | `POS` | `start1` |
| **De novo filter** | `intersect_3gen.py` | `intersect_3gen.py` | `intersect_3gen.py` |