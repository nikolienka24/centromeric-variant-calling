# Read Mapping — Variant Calling Pipeline

This pipeline aligns Oxford Nanopore reads to a centromeric reference assembly and calls variants using DeepVariant. It represents the standard read-level variant calling approach in the benchmark.

---

## Directory Structure

```
read_mapping/
└── map_call.sh    # Minimap2 alignment + DeepVariant variant calling
```

---

## `map_call.sh`

A two-step pipeline that first aligns ONT reads to a reference haplotype using Minimap2, then calls variants with DeepVariant (ONT R10.4 model) run inside a Singularity container. The resulting VCF is filtered for high-confidence PASS variants within centromeric BED regions.

**Requirements:** `minimap2`, `samtools`, `bcftools`, `singularity`, DeepVariant `.sif` image (v1.10.0)

**Usage (PBS/qsub):**
```bash
qsub map_call.sh -v SAMPLE="PAN027",CHR="chr22",PARENT="paternal",OUT="/path/to/results"
```

| Variable | Description |
|----------|-------------|
| `SAMPLE` | Sample ID (e.g., `PAN027`) |
| `CHR` | Chromosome to process (e.g., `chr22`) |
| `PARENT` | Haplotype origin (`paternal` or `maternal`) |
| `OUT` | Base output directory — results go into `<OUT>/<SAMPLE>.<CHR>.<PARENT>/` |

---

## Pipeline Steps

**Step 1 — Alignment (conditional)**  
If a pre-existing BAM for the sample is found at the expected path, it is reused directly. Otherwise, reads are extracted from the input BAM with `samtools fastq` and aligned to the reference using Minimap2 (`map-ont` preset, 16 threads), then sorted and indexed.

**Step 2 — Variant Calling (DeepVariant)**  
DeepVariant runs inside a Singularity container using the `ONT_R104` model with 16 shards. The scratch directory is bind-mounted into the container at `/data`.

**Step 3 — Filtering**  
The raw VCF is filtered with `bcftools view` applying the following criteria:
- `FILTER = PASS`
- `QUAL ≥ 20`
- `FORMAT/DP ≥ 20` (minimum read depth)
- Allele frequency `FORMAT/AD[0:1] / FORMAT/DP ≥ 0.8`
- Restricted to centromeric regions defined by the BED file

---

## Outputs

Results are saved to `<OUT>/<SAMPLE>.<CHR>.<PARENT>/`:

| File | Description |
|------|-------------|
| `<SAMPLE_NAME>.raw.vcf` | Unfiltered DeepVariant output |
| `<SAMPLE_NAME>.filtered.vcf` | High-confidence variants passing all filters |

---

## Notes

- Input file paths (reference FASTA, reads BAM, centromere BED, DeepVariant image) are **hardcoded** in the `ENVIRONMENT SETUP` and `PREPARE INPUTS` sections. Update these to match your storage layout before submitting.
- The DeepVariant Singularity image (`deepvariant_1.10.0.sif`) must be present at the path defined by `CONTAINER_IMG`.
- Edit the `ENVIRONMENT SETUP` section to match your cluster's conda configuration.