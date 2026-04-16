# Variant Calling Module

This module implements structural variant discovery, benchmarking, and annotation.

## Workflows

The two workflows are distinct but share some of the execution.

### Benchmarking (CHM13-centered)

In this workflow all the variants detected by tools on both reference are collected using merging and the Grch38 are lifted on chm13 in order to evaluate them
### Annotation (Grch38-centered)

In this workflow the variants are collected on Grch38 in order to annotate them using the needLR tool. 
As a layer of confirm the chm13 liftover is performed to validate on both references.

Run:
```bash
sbatch run_snakemake_annotation.sh
``` 
The needLR tool is already integrated in the workflow, but is also available as a docker container in the folder. 
