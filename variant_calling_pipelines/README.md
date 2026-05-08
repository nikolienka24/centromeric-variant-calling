# Variant Calling Pipelines

This folder contains the three independent variant calling strategies evaluated in the benchmark. Each approach uses a different underlying methodology to detect mutations in centromeric regions, allowing for cross-validation of results in the `validation/` step.

---

## Directory Structure

```
variant_calling_pipelines/
├── read_mapping/               # Read-to-reference alignment based calling
├── pangenomes/                 # Pangenome graph-based alignment and calling
└── assembly/                   # Assembly-to-reference alignment based calling
    ├── stretcher/              # EMBOSS Stretcher global alignment
    ├── centrolign/             # Centrolign centromere-aware aligner
    └── minimap2+paftools/      # Minimap2 assembly alignment + Paftools calling
```

---

## Overview of Strategies

### `read_mapping/`
Simulated long reads are aligned directly to the centromeric reference using a read-level aligner. Variants are called from the resulting alignment. This approach captures sequencing-level noise and is representative of a standard nanopore variant calling workflow.  
→ See [`read_mapping/README.md`](read_mapping/README.md)

### `pangenomes/`
Centromeric sequences are incorporated into a pangenome graph structure using Centrolign in pangenome mode. Variants are identified by comparing paths through the graph rather than through linear alignment.  
→ See [`pangenomes/README.md`](pangenomes/README.md)

### `assembly/`
Query assemblies are aligned to the reference centromere at the sequence level using three different alignment tools. Each tool produces an independent set of variant calls, which are later compared in the `validation/` step.

| Tool | Approach                                                                     |
|------|------------------------------------------------------------------------------|
| **Stretcher** | Global pairwise alignment (EMBOSS), best suited for similar-length sequences |
| **Centrolign** | Centromere-aware pairwise-aligner, designed for repetitive regions           |
| **Minimap2 + Paftools** | Assembly-to-reference mapping followed by VCF variant extraction             |

→ See [`assembly/README.md`](assembly/README.md)

---

## Role in the Broader Pipeline

```
data_preparation/          → reference extraction & read simulation
haplotype_assignment/      → assign query to correct H1/H2 haplotype
        │
        ▼
variant_calling_pipelines/
    ├── read_mapping/
    ├── pangenomes/
    └── assembly/
            │
            ▼
        validation/        → cross-tool consensus calling
            │
            ▼
        plots/             → visualization & benchmarking
```

Each pipeline in this folder produces a set of variant calls in TSV or VCF format that feed directly into `validation/intersect_variants.py`.