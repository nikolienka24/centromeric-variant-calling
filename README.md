# Benchmarking Methods for Variant Detection in Human Centromeric Regions

This repository contains computational pipelines and custom tools designed for the systematic detection and validation of genetic variants within human centromeric regions.

## Project Overview
Centromeres are essential for cell division but remain one of the most difficult 
regions of the human genome to analyze due to their highly repetitive alpha-satellite 
DNA structure. This project contains pipelines and tools for systematical evaluation of three primary 
approaches for variant calling in these complex regions:

1.  **Mapping-based approach**: Utilizing Minimap2 and the DeepVariant deep-learning framework.
2.  **Assembly-based approach**: Utilizing EMBOSS Stretcher, Centrolign (pairwise mode), and the Minimap2 + Paftools.js workflows.
3.  **Pangenome-based approach**: Utilizing Centrolign in pangenome mode.

---

## Directory Structure
The repository is organized into subdirectories representing different stages of the analytical workflow. Each subdirectory contains its own detailed README with script descriptions.

* **`data_preparation/`**: Contains scripts for extracting centromeric sequence and simulating mutations.
* **`haplotype_assignment/`**: Scripts for identifying inherited parental haplotypes.
* **`plots/`**: Utility scripts for generating quality control plots and other visualizations.
* **`validation/`**: Scripts for cross-method consensus validation and automated validation scripts.
* **`variant_calling_pipelines/read_mapping/`**: Implementation of the mapping-based pipeline.
* **`variant_calling_pipelines/assembly/`**: Workflows for assembly-based variant detection. 
* **`variant_calling_pipelines/pangenomes/`**: Scripts for constructing pangenome graphs with Centrolign and extracting variants from pangenome matrices.
---

## Technical Infrastructure
All scripts are prepared to be executed on the **MetaCentrum** computational grid (CESNET).