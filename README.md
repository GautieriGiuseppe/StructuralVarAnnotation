# ImmuneVariantCalling

Snakemake-based pipeline for alignment and variant calling, designed for HPC execution.

## Overview

The workflow is divided in three modules:

- Alignment
- Alignment Quality Control
- Variant Calling

It makes usage of snakemake and sbatch scripts.

'''
.
|__ align.smk 			# Alignment workflow
|__ alignqc.smk			# Alignment QC workflow
|__ config.yml			# General config containing references and parameters
|__ alignqc_env.yml 		# Conda environment for QC
|__ run_snakemake_align.sh 	# SBATCH script to run the alignment
|__ run_snakemake_alignqc.sh 	# SBATCH script to run the alignemnt QC
|__ variant_calling/ 		# Variant calling module
|__ README.md
'''

## Workflow Description

The purpose of the workflow is to identify all the variants in the cohort present in the samples tsv file.

To achieve this, the alignment is performed using both Grch38 and chm13 (T2T) references. 

After the alignment and QC to compare the quality of it across references, the variant calling process is executed.

Here the workflow employes three different tools (Sniffles, Delly and CuteSV) to capture all the variants, merging with union all the resulting variants to collect them.

To achieve this, the liftover among references is done using Crossmap and the final merging by OctopuSV.

## Requirements

- Conda/Mamba
- Snakemake >= 7.32.4

These are the only requirements since the environments and tools used are built automatically when needed.

## Usage

Specify the input paths of the files through a samples tsv file before running the workflow.
