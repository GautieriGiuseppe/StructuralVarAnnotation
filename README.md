# ImmuneVariantCalling

Snakemake-based pipeline for alignment and variant calling, designed for HPC execution.

## Overview

The workflow is divided in three modules:

- Alignment
- Alignment Quality Control
- Variant Calling

The varian calling module is further split into two independent workflows:

- Benchmarking (CHM13-centered)
- Annotation (GRCh38-centered)

The pipeline is implemented using Snakemake and designed for scalable execution via SLURM (sbatch).

Repository structure:

```
.
|__ align.smk                           # Alignment workflow
|__ alignqc.smk	                        # Alignment QC workflow
|__ config.yml                          # General config containing references and parameters
|__ alignqc_env.yml                     # Conda environment for QC
|__ run_snakemake_align.sh              # SBATCH script to run the alignment
|__ run_snakemake_alignqc.sh            # SBATCH script to run the alignemnt QC
|
|__ variant_calling/                    # Variant calling module
|  |--benchmark_master.smk              # CHM13 benchmarking workflow
|  |--annotation_master.smk             # GRCh38 annotation workflow
|  |--sample_variant_calling.smk        # Shared per-sample SV calling
|  |--cohort_merge.smk                  # CHM13 cohort merging
|  |--cohort_merge_grch38.smk           # GRCh38 cohort merging
|  |--needLR_grch38.smk                 # needLR annotation
|  |--crossref_confirmation_grch38smk   # CHM13 -> GRCh38 confirmation
|  |--run_snakemake_benchmarking.sh
|  |--run_snakemake_annotation.sh
|  |_envs/                              # Conda environments
|
|__ README.md
```


## Workflow Description

### Alignment

Reads are aligned againts both:

- GRCh38
- CHM13 (T2T)

This dual reference allows: 

- Improved variant discovery
- Cross-reference comparison
- Benchmarking across genome builds

### Variant Calling

Structural variants are detected using:

- Sniffles2
- CuteSV
- Delly lr

Per-sample calls are filtered and normalized, then merged at cohort level.

## Requirements

- Conda/Mamba
- Snakemake >= 7.32.4

These are the only requirements since the environments and tools used are built automatically when needed.

## Usage

Specify the input paths of the files through a samples tsv file before running the workflow.

samples tsv file example

```
sample_id	ubam	summary	fastq	batch_id
A_1	/path/to/A_1.ubam	/path/to/A_1_sequencing_summary.tsv.gz	/path/to/A_1.fastq.gz	A
A_2	/path/to/A_2.ubam	/path/to/A_2_sequencing_summary.tsv.gz	/path/to/A_2.fastq.gz	A
A_3	/path/tp/A_3.ubam	/path/to/A_3_sequencing_summary.tsv.gz	/path/to/A_3.fastq.gz	A
```
