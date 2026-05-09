# Plots

This section contains all visualization and quality-control scripts used across the benchmarking pipeline. Scripts are organized into the main `plots/` directory for variant analysis plots, a `qc/` subdirectory for read quality control, and a `benchmarking/` subdirectory for tool performance metrics.

---

## Directory Structure

```
plots/
├── distribution_variants_types.py        # Stacked bar chart: variant concordance by type
├── variant_distribution_per_chr.py       # Variant distribution per chromosome
├── venn_diagrams.tool_cross_validation.py # Venn diagram: caller cross-validation & homopolymer analysis
├── qc/
│   └── nanoplot_qc.sh                    # NanoPlot QC report for simulated FASTQ reads
└── benchmarking/
    ├── joined_stats.py                   # Dashboard: combined runtime & memory across all tools
    └── separate_stats.py                 # Individual plots: runtime, memory, and efficiency scatter
```

---

## Main Scripts

### `distribution_variants_types.py`

Visualizes how many variant callers agree on each type of variant. Takes a consensus TSV (output of a multi-tool intersection step) and produces a stacked bar chart with 2-tool vs. 3-tool concordance per variant class.

**Requirements:** `pandas`, `matplotlib`

**Usage:**
```bash
python distribution_variants_types.py \
    -i consensus_variants.tsv \
    -o plots/variant_concordance.png
```

| Argument | Description |
|----------|-------------|
| `-i` / `--input` | Combined consensus TSV file with `ref_seq`, `alt_seq`, and `tools_count` columns |
| `-o` / `--output` | Output PNG path |

**Variant classes:** SNP, Substitution, Indel (≤50 bp), Indel (>50 bp)  
**Output:** Single stacked bar plot (`300 dpi`)

---

### `variant_distribution_per_chr.py`

Parses a VCF file and generates two plots: a per-chromosome stacked bar chart of variant types, and boxplots showing the length distribution of short and long indels.

**Requirements:** `pandas`, `matplotlib`, `seaborn`

**Usage:**
```bash
python variant_distribution_per_chr.py \
    -i variants.vcf \
    -o /path/to/output_dir/
```

| Argument | Description |
|----------|-------------|
| `-i` / `--input` | Input VCF file |
| `-o` / `--out_dir` | Output directory |

**Outputs:**
- `variant_distribution.png` — per-chromosome variant type counts
- `indel_lengths_combined.png` — side-by-side boxplots for short and long indels

---

### `venn_diagrams.tool_cross_validation.py`

Analyzes concordance across three variant callers (Stretcher, Minimap2, Centrolign) and visualizes the overlap in a 3-way Venn diagram. Each intersection region is annotated with both the regular variant count and the homopolymer variant count (`+N` format). Chr22 paternal/haplotype2 data is excluded automatically.

**Requirements:** `pandas`, `matplotlib`, `matplotlib_venn`

**Usage:**
```bash
python venn_diagrams.tool_cross_validation.py \
    -i 2gen.combined.2_of_3.tsv \
    -m master.2gen.tsv \
    -o /path/to/output/venn_concordance.png
```

| Argument | Description |
|----------|-------------|
| `-i` / `--intersect` | Intersection TSV file (variants called by ≥2 tools) |
| `-m` / `--master` | Master variant TSV file with per-tool calls and positions |
| `-o` / `--output` | Output PNG path |

**Output:** Venn diagram saved as PNG (`300 dpi`)

---

## `qc/` — Read Quality Control

### `qc/nanoplot_qc.sh`

Runs NanoPlot on a simulated or real FASTQ file to generate a full QC report including read length distribution, N50, and quality score plots.

**Requirements:** `NanoPlot`, 2 CPUs

> **Note:** Update `CONDA_BASE` and `ENV_PATH` in the `ENVIRONMENT SETUP` section of this script to match your cluster's conda configuration before submitting.

**Usage (PBS/qsub):**
```bash
qsub nanoplot_qc.sh -v IN="simulated_reads.fastq",OUT="/path/to/qc_results"
```

| Variable | Description |
|----------|-------------|
| `IN` | Input FASTQ file (long-read, Nanopore) |
| `OUT` | Output directory for NanoPlot reports and plots |

**Output:** Full NanoPlot HTML report and supporting figures in `<OUT>/`

---

## `benchmarking/` — Tool Performance Metrics

Both scripts process TSV performance logs (one file per tool) with `RealTime_sec` and `MaxMem_MB` columns. They cover the same set of tools with slightly different output formats.

**Supported tools / expected filenames:**

| Tool Label | Expected File |
|------------|---------------|
| `Centrolign_Align` | `centrolign_align.tsv` |
| `Centrolign_Pang` | `centrolign_pangenomes.tsv` |
| `Minimap2/Deepvariant` | `map+deepvariant.tsv` |
| `Minimap2/Paftools` | `minimap2+paftools.tsv` |
| `Stretcher` | `stretcher_align.tsv` |

> **Note:** `Centrolign_Pang` is excluded from the runtime bar chart in both scripts to maintain a readable log scale, but is included in memory and scatter plots.

---

### `benchmarking/joined_stats.py`

Generates a single three-panel dashboard image combining runtime, memory, and efficiency plots side by side. Also exports a summary statistics TSV.

**Requirements:** `pandas`, `matplotlib`, `seaborn`

**Usage:**
```bash
python joined_stats.py \
    -i /path/to/performance_logs/ \
    -o /path/to/output_dir/
```

| Argument | Description |
|----------|-------------|
| `-i` / `--input_dir` | Directory containing tool `.tsv` log files |
| `-o` / `--output_dir` | Output directory for plots and stats |

**Outputs:**
- `benchmarking_dashboard.png` — 3-panel figure (runtime, memory, efficiency)
- `benchmarking_summary_stats.tsv` — aggregated statistics table

---

### `benchmarking/separate_stats.py`

Generates the same three plots as `joined_stats.py` but saves each as a separate file. Also accepts an explicit path for the summary statistics export.

**Requirements:** `pandas`, `matplotlib`, `seaborn`

**Usage:**
```bash
python separate_stats.py \
    -i /path/to/performance_logs/ \
    -o /path/to/output_dir/ \
    -s /path/to/summary_stats.tsv
```

| Argument | Description |
|----------|-------------|
| `-i` / `--input_dir` | Directory containing tool `.tsv` log files |
| `-o` / `--output_dir` | Output directory for plots |
| `-s` / `--stats_out` | Path to save the summary TSV |

**Outputs:**
- `1_total_runtime.png` — total runtime per tool (log scale, excluding `Centrolign_Pang`)
- `2_peak_memory.png` — peak memory per tool (log scale)
- `3_efficiency_scatter.png` — time vs. memory scatter (log-log)
- `<stats_out>` — summary statistics TSV

---

## Notes

- All bar charts and scatter plots use the `magma` colormap for visual consistency across the pipeline.