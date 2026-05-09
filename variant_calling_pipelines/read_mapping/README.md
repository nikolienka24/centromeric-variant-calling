# Read Mapping — Variant Calling Pipeline

This pipeline aligns Oxford Nanopore reads to a centromeric reference assembly and calls variants using DeepVariant. It represents the standard read-level variant calling approach in the benchmark.

---

## Directory Structure

```
read_mapping/
├── align_call.sh        # Minimap2 alignment + DeepVariant variant calling + filtering
└── config.example.sh    # Template configuration file — copy to config.sh and fill in your paths
```

---

## `align_call.sh`

A two-step pipeline that first aligns ONT reads to a reference haplotype using Minimap2, then calls variants with DeepVariant (ONT R10.4 model) run inside a Singularity container. The resulting VCF is filtered for high-confidence PASS variants within centromeric regions and further filtered against problematic and gap regions using bedtools.

**Requirements:** `minimap2`, `samtools`, `bcftools`, `bedtools`, `singularity`, DeepVariant `.sif` image (v1.10.0)

> **Note:** All input paths, sample name, and tool binaries are configured via `config.sh`. Copy `config.example.sh` to `config.sh` and fill in your paths before submitting.

**Usage (PBS/qsub):**
```bash
qsub align_call.sh
```

---

## Pipeline Steps

**Step 1 — Offset extraction**
The genomic offset for the reference sequence is extracted from `BED_OFFSETS` by matching `REF_ID` against column 1.

**Step 2 — Alignment (conditional)**
If a pre-existing BAM for the sample is found in `BAM_DIR`, it is reused directly. Otherwise, reads are extracted from the input BAM with `samtools fastq` and aligned to the reference using Minimap2 (`map-ont` preset, 16 threads), then sorted and indexed.

**Step 3 — Variant calling (DeepVariant)**
DeepVariant runs inside a Singularity container using the `ONT_R104` model with 16 shards. The scratch directory is bind-mounted into the container at `/data`.

**Step 4 — Filtering**
The raw VCF is filtered in three stages:
- `bcftools view` applies quality filters and restricts calls to centromeric regions (`BED_REGIONS`):
  - `FILTER = PASS`
  - `QUAL ≥ 20`
  - `FORMAT/DP ≥ 20` (minimum read depth)
  - Allele frequency `FORMAT/AD[0:1] / FORMAT/DP ≥ 0.8`
- Genomic offset is added to the VCF `POS` field to convert to absolute coordinates
- `bedtools intersect` removes variants in problematic regions (`BED_PROBLEMATIC`)
- `bedtools intersect` removes variants in genomic gaps (`BED_GAPS`)

---

## Outputs

Results are saved to `<OUT>/<SAMPLE_NAME>/`:

| File | Description |
|------|-------------|
| `<SAMPLE_NAME>.raw.vcf.gz` | Unfiltered DeepVariant output |
| `<SAMPLE_NAME>.raw.vcf.gz.tbi` | Tabix index for the raw VCF |
| `<SAMPLE_NAME>.filtered.vcf` | High-confidence variants passing all filters |