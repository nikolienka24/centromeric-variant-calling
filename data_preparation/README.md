# Data Preparation

This section covers the data preparation pipeline for benchmarking variant callers on centromeric regions. It consists of two stages: **reference extraction** and **read simulation** (with optional synthetic mutation injection).

---

## Directory Structure

```
data_preparation/
├── extract_centromeres.sh          # Extract centromeric regions from a reference genome
└── simulation/
    ├── generate_mutations_to_vcf.py  # Generate synthetic SNPs/indels in VCF format
    └── simulate_reads.nanosim.sh     # Simulate Oxford Nanopore reads with NanoSim
```

---

## 1. `extract_centromeres.sh`

Extracts centromeric sequences from a reference FASTA using coordinates defined in a BED file. Produces both a combined multi-region FASTA and individual per-region FASTA files.

**Requirements:** `samtools`, `bedtools`

**Usage (PBS/qsub):**
```bash
qsub extract_centromeres.sh -v REF="genome.fasta",BED="centromeres.bed",OUT="/path/to/output"
```

| Variable | Description |
|----------|-------------|
| `REF` | Path to the reference FASTA file |
| `BED` | BED file with centromere coordinates |
| `OUT` | Output directory |

**Outputs:**
- `<OUT>/reference.centromeres.fasta` — combined centromeric FASTA (+ `.fai` index)
- `<OUT>/per_chromosome/<chrom>_<start>_<end>.fasta` — individual region FASTAs (+ `.fai` indexes)

---

## 2. `simulation/generate_mutations_to_vcf.py`

Generates a synthetic VCF file containing random mutations (SNPs, 2-bp substitutions, short indels, long indels) placed within BED-defined regions. All mutations are validated against the reference to ensure they fall on canonical bases (A/C/G/T).

**Requirements:** `pysam`, Python 3

**Usage:**
```bash
python generate_mutations_to_vcf.py \
    -r reference.centromeres.fasta \
    -b centromeres.bed \
    -o synthetic_mutations.vcf \
    [--snps 40] [--sub2 20] [--short_indels 20] [--long_indels 20]
```

| Argument | Default | Description |
|----------|---------|-------------|
| `-r` / `--reference` | — | Reference FASTA (indexed) |
| `-b` / `--bed` | — | BED file defining target regions |
| `-o` / `--output` | — | Output VCF path |
| `--snps` | 40 | Number of SNPs |
| `--sub2` | 20 | Number of 2-bp substitutions |
| `--short_indels` | 20 | Short indels (3–50 bp) |
| `--long_indels` | 20 | Long structural indels (51–4000 bp, max 1 per chromosome) |

**Output:** VCF v4.2 file with `PASS`-filtered mutations and `SVTYPE`/`LEN` INFO fields.

---

## 3. `simulation/simulate_reads.nanosim.sh`

Simulates Oxford Nanopore reads over the (mutated) centromeric reference using NanoSim. Uses a pre-trained error/read-length model to produce realistic FASTQ output at a specified coverage depth.

**Requirements:** NanoSim (`simulator.py`)

**Usage (PBS/qsub):**
```bash
qsub simulate_reads.nanosim.sh \
    -v REF="reference.centromeres.fasta",MODEL="/path/to/nanosim_model",OUT="/path/to/output",COV=30
```

| Variable | Default | Description |
|----------|---------|-------------|
| `REF` | — | Reference FASTA to simulate from |
| `MODEL` | — | Path prefix of a trained NanoSim model |
| `OUT` | — | Output directory |
| `COV` | 60 | Target sequencing coverage |

**Output:** `<OUT>/simulated_reads_aligned_reads.fastq` (and associated NanoSim output files)

---

> **Note:** Environment setup (conda paths) inside the PBS scripts is cluster-specific. Edit the `ENVIRONMENT SETUP` section of each `.sh` script before submitting.