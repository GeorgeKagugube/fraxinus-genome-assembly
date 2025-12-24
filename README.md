# F. mandshurica Genome Assembly & Polishing Pipeline

Snakemake workflow for assembling a **Fraxinus mandshurica** reference genome from ONT long reads and Illumina short reads.  
It runs QC, genome size estimation, assembly (**Shasta**), polishing (**Medaka** + **NextPolish×2**), and evaluation (**BUSCO**, 
**QUAST**, optional **Merqury**).

---

## Contents
- [Quick start](#quick-start)
- [Inputs & outputs](#inputs--outputs)
- [Run critical sections only](#run-critical-sections-only)
- [DAG (workflow graph)](#dag-workflow-graph)
- [Troubleshooting](#troubleshooting)
- [Environments/containers](#environmentscontainers)
- [Reproducing the DAG image](#reproducing-the-dag-image)

---

## Quick start

```bash
#  (recommended) create a clean snakemake env
mamba create -n snakemake snakemake=7 -c conda-forge -y
conda activate snakemake

# 1) run a dry-run to verify the DAG
snakemake -n --use-conda

# 2) execute with 32 cores (adjust as needed)
snakemake --use-conda --cores 32
```
---
##  Overview of Steps

1. **ONT Read QC** – NanoPlot summarises quality metrics.
2. **Genome Size Estimation** – k-mer analysis with Jellyfish.
3. **Assembly** – using **Shasta**.
4. **Initial QC** – BUSCO & QUAST on raw assembly.
5. **Polishing**:
   - **Medaka** (ONT-based polishing)
   - **NextPolish** (Illumina-based polishing, two rounds)
6. **Post-polishing QC** – BUSCO & QUAST after each round.
7. **Optional** – Merqury completeness & QV analysis.

---

## Input Files

- **ONT reads**: A single `.fastq.gz` or `.fq.gz` file.
- **Illumina reads**: Multiple paired-end `R1`/`R2` `.fastq.gz` files (one pair per sample).
- **Configurable paths** in the Snakefile or by environment variables.

---

## 🖥️ System Requirements

- **Snakemake** ≥ 7.x
- **mamba** or **conda**
- Optional: **Singularity** or **Apptainer** for containerised runs.

---

## 🚀 Usage

### 1. Clone the Repository
```bash
git clone https://github.com/<yourusername>/Fmandshurica-genome-pipeline.git
cd Fmandshurica-genome-pipeline

# Annotation Workflow (Structural + Functional)

This workflow annotates the polished genome with **BRAKER3** (RNA‑seq + protein homology), then adds functional layers via 
**InterProScan** and **eggNOG‑mapper**. It also supports optional **Liftoff** transfer from a close relative.
```
## Inputs
- Genome FASTA (soft‑masked is best).
- RNA‑seq paired‑end FASTQs (edit `RNASEQ_LIBS` in Snakefile).
- Close‑relative protein FASTA for homology hints (`PROTEINS_FA`).
- (Optional) Liftoff donor genome FASTA + GFF3.

## Key Outputs
- **Structural**: `03_braker/braker.gff3`, `03_braker/braker.proteins.faa`, `03_braker/braker.cds.fna`
- **Functional**:
  - InterPro: `04_function/interproscan.tsv`
  - eggNOG: `04_function/eggnog/eggnog.emapper.annotations`
  - **Merged**: `04_function/functional_merged.tsv`
  - Coverage plot: `04_function/functional_summary.png`
- **QC**: BUSCO on predicted proteome: `06_qc/busco_proteins/short_summary...proteins.txt`

## Common commands
```bash
# Dry run
snakemake -n -p --use-conda

# Full run
snakemake --use-conda --cores 32

# Structural only (masking + BRAKER3)
snakemake --use-conda --cores 32 03_braker/braker.gff3 03_braker/braker.proteins.faa

# Functional only (after BRAKER3)
snakemake --use-conda --cores 16 04_function/functional_merged.tsv 04_function/functional_summary.png
```

