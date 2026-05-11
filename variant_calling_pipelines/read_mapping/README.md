# Read Mapping — Variant Calling Pipeline

This pipeline aligns Oxford Nanopore reads to a centromeric reference assembly and calls variants. It includes two variant calling approaches: DeepVariant for SNP/indel calling and Sniffles2 for structural variant calling.

---

## Directory Structure

```
read_mapping/
├── map_align.deepvariant.sh        # Minimap2 alignment + DeepVariant variant calling + filtering
├── config.example.deepvariant.sh   # Template config for DeepVariant pipeline
├── map_align.sniffles2.sh          # Minimap2 alignment + Sniffles2 SV calling + filtering
└── config.example.sniffles2.sh     # Template config for Sniffles2 pipeline
```

---

## `map_align.deepvariant.sh`

A two-step pipeline that first aligns ONT reads to a reference haplotype using Minimap2, then calls variants with DeepVariant (ONT R10.4 model) run inside a Singularity container. The resulting VCF is filtered for high-confidence PASS variants within centromeric regions and further filtered against problematic and gap regions using bedtools.

**Requirements:** `minimap2`, `samtools`, `bcftools`, `bedtools`, `singularity`, DeepVariant `.sif` image (v1.10.0)

> **Note:** All input paths, sample name, and tool binaries are configured via `config.sh`. Copy `config.example.deepvariant.sh` to `config.sh` and fill in your paths before submitting.

**Usage (PBS/qsub):**
```bash
qsub map_align.deepvariant.sh
```

---

## `map_align.sniffles2.sh`

A two-step pipeline that first aligns ONT reads to a reference haplotype using Minimap2 (or reuses a pre-existing BAM), then calls structural variants with Sniffles2. The resulting VCF is filtered for high-confidence PASS calls within centromeric regions and further filtered against problematic and gap regions using bedtools.

**Requirements:** `minimap2`, `samtools`, `bcftools`, `bedtools`, `sniffles`

> **Note:** All input paths, sample name, and tool binaries are configured via `config.sh`. Copy `config.example.sniffles2.sh` to `config.sh` and fill in your paths before submitting.

**Usage (PBS/qsub):**
```bash
qsub map_align.sniffles2.sh
```

---

## Pipeline Steps

### DeepVariant (`map_align.deepvariant.sh`)

**Step 1 — Offset extraction**
The genomic offset for the reference sequence is extracted from `BED_OFFSETS` by matching `REF_ID` against column 1. If `BED_OFFSETS` is not set or the file does not exist, the offset defaults to 0.

**Step 2 — Alignment (conditional)**
If a pre-existing BAM for the sample is found at `READS_BAM`, it is reused directly (re-indexed if the index is missing or outdated). Otherwise, reads are extracted from the input BAM with `samtools fastq` and aligned to the reference using Minimap2 (`map-ont` preset, 16 threads), then sorted and indexed.

**Step 3 — Variant calling (DeepVariant)**
DeepVariant runs inside a Singularity container using the `ONT_R104` model with 16 shards. The scratch directory is bind-mounted into the container at `/data`.

**Step 4 — Filtering**
The raw VCF is filtered in three stages:
- `bcftools view` applies quality filters and restricts calls to centromeric regions (`BED_REGIONS`):
  - `FILTER = PASS`
  - `QUAL ≥ 20`
  - `FORMAT/DP ≥ 20` (minimum read depth)
  - Allele frequency `FORMAT/AD[0:1] / FORMAT/DP ≥ 0.8`
- Genomic offset correction: if `BED_OFFSETS` is set and exists, the offset is added to the VCF `POS` field to convert to absolute coordinates; otherwise this step is skipped
- `bedtools intersect` removes variants in problematic regions (`BED_PROBLEMATIC`); skipped if not set or file not found
- `bedtools intersect` removes variants in genomic gaps (`BED_GAPS`); skipped if not set or file not found

---

### Sniffles2 (`map_align.sniffles2.sh`)

**Step 1 — Offset extraction**
The genomic offset for the reference sequence is extracted from `BED_OFFSETS` by matching `REF_ID` against column 1. If `BED_OFFSETS` is not set or the file does not exist, the offset defaults to 0.

**Step 2 — Alignment (conditional)**
If a pre-existing BAM for the sample is found at `READS_BAM`, it is reused directly (re-indexed if the index is missing or outdated). Otherwise, reads are extracted from the input BAM with `samtools fastq` and aligned to the reference using Minimap2 (`map-ont` preset, 16 threads), then sorted and indexed.

**Step 3 — Variant calling (Sniffles2)**
Sniffles2 calls structural variants directly from the BAM using the reference FASTA. Minimum SV length is 50 bp and minimum mapping quality is 20. Read names are included in the output (`--output-rnames`).

**Step 4 — Filtering**
The raw VCF is filtered in three stages:
- `bcftools view` applies quality filters and restricts calls to centromeric regions (`BED_REGIONS`):
  - `FILTER = PASS`
  - `QUAL ≥ 20`
  - `INFO/SUPPORT ≥ 10` (minimum number of supporting reads)
- Genomic offset correction: if `BED_OFFSETS` is set and exists, the offset is added to the VCF `POS` field to convert to absolute coordinates; otherwise this step is skipped
- `bedtools intersect` removes variants in problematic regions (`BED_PROBLEMATIC`); skipped if not set or file not found
- `bedtools intersect` removes variants in genomic gaps (`BED_GAPS`); skipped if not set or file not found

---

## Outputs

Results are saved to `<OUT>/<SAMPLE_NAME>/`:

### DeepVariant

| File | Description |
|------|-------------|
| `<SAMPLE_NAME>.raw.vcf.gz` | Unfiltered DeepVariant output |
| `<SAMPLE_NAME>.raw.vcf.gz.tbi` | Tabix index for the raw VCF |
| `<SAMPLE_NAME>.filtered.vcf` | High-confidence SNPs/indels passing all filters |

### Sniffles2

| File | Description |
|------|-------------|
| `<SAMPLE_NAME>.raw.vcf.gz` | Unfiltered Sniffles2 output |
| `<SAMPLE_NAME>.raw.vcf.gz.tbi` | Tabix index for the raw VCF |
| `<SAMPLE_NAME>.filtered.vcf` | High-confidence structural variants passing all filters |