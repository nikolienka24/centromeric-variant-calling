# Pangenomes — Variant Calling Pipeline

This pipeline constructs a multi-sequence pangenome graph from centromeric assemblies using Centrolign and extracts pairwise variant calls by comparing paths through the graph's variant matrix. It is specifically designed for the repetitive, structurally complex nature of centromeric regions.

---

## Directory Structure

```
pangenomes/
├── build_graph_call.sh      # Pangenome construction + variant matrix extraction
└── convert_var_matrix.py    # Pairwise variant comparison from the matrix
```

---

## `build_graph_call.sh`

Builds a centromere pangenome graph from four input sequences (one reference + three pedigree samples), extracts a variant matrix from the graph, and runs pairwise comparisons between the reference and each relative.

**Requirements:** `centrolign`, `make_var_mat`, `mashtree`, Python 3

> **Note:** All input paths, sample IDs, and tool binaries are configured via `config.sh`. Copy `config.example.sh` to `config.sh` and fill in your paths before submitting.

**Usage (PBS/qsub):**
```bash
qsub build_graph_call.sh
```

### Pipeline Steps

**Step 1 — FASTA preparation**
Each input FASTA is re-headered to a clean sample label and localized to scratch. All four sequences are concatenated into `all_sequences.fasta`.

**Step 2 — Guide tree**
A pre-computed Newick-format guide tree is provided via `GUIDE_TREE` in `config.sh` and copied to scratch. This replaces the previous `mashtree`-based approach, giving the user full control over the tree topology.
**Step 3 — Pangenome graph construction**
Centrolign uses the guide tree (`-T`) to progressively align all sequences and outputs a pangenome graph in GFA format. Intermediate subproblem files are stored in `subproblems_<CHR>/`.

**Step 4 — Variant matrix extraction**
`make_var_mat` processes the GFA graph with `--base --indels --mnvs` flags to produce a matrix TSV where rows are sequences and columns are variant positions.

**Step 5 — Pairwise comparison**
`convert_var_matrix.py` is called twice — once comparing the reference (`S2_ID`) against each relative (`S1_ID`, `S3_ID`) — to produce generation-specific variant TSVs.

### Outputs

| File | Description                                        |
|------|----------------------------------------------------|
| `<CHR>_pangenome.gfa` | Full pangenome graph in GFA format                 |
| `<CHR>_guide_tree.nwk` | Newick guide tree used for Centrolign              |
| `<CHR>_matrix.tsv` | Raw variant matrix (all sequences × all positions) |
| `<CHR>_variants_final.1gen.tsv` | Reference vs. generation 1 relative variants       |
| `<CHR>_variants_final.2gen.tsv` | Reference vs. generation 2 relative variants       |

---

## `convert_var_matrix.py`

Compares two named rows in the variant matrix TSV and writes out only the positions where they differ. Called internally by `build_graph_call.sh` but can also be run independently.

**Requirements:** `centrolign`, `make_var_mat`, Python 3

**Usage:**
```bash
python convert_var_matrix.py \
    -i matrix.tsv \
    -o variants_output.tsv \
    -r1 <reference_id> \
    -r2 <relative_id>
```

| Argument | Description |
|----------|-------------|
| `-i` / `--input` | Input variant matrix TSV (`make_var_mat` output) |
| `-o` / `--output` | Output variants TSV |
| `-r1` / `--row1` | Name of the first sequence row (e.g., reference) |
| `-r2` / `--row2` | Name of the second sequence row (e.g., relative) |

**Output format** (`variants_output.tsv`):

| Column | Description |
|--------|-------------|
| `COL_IDX` | Column index in the original matrix (variant position) |
| `<row1>` | Value at this position for sequence 1 |
| `<row2>` | Value at this position for sequence 2 |

Only positions where the two sequences differ are written. The total count of differences is printed to stdout.