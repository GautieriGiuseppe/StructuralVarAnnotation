# needLR v3.5 Container (GRCh38 Annotation)

## Overview

This container provides a reproducible environment to run **needLR v3.5**
for structural variant annotation on **GRCh38 cohort VCFs**.

It is designed to be used:

- standalone via Singularity
- integrated in Snakemake pipelines

---

## Requirements

- Singularity / Apptainer
- Input files:
  - bgzipped VCF (`.vcf.gz`)
  - VCF index (`.tbi`)
  - GRCh38 reference genome (`.fa`)
  - sample names file (`.txt`)

---

## Container

needlr.sif

---

## Standard Usage

### Run needLR (cohort mode)

```bash
singularity exec needlr.sif \
bash -c "
cd /opt/needLR_v3.5_local && \
./needLR_v3.5_cohort.sh \
  -g /path/to/reference.fa \
  -v /path/to/input.vcf.gz \
  -a /path/to/samples.txt \
  -d 8
"
```

### Snakemake integration

```
rule needlr_grch38:
    input:
        vcf = config["output"] + "/cohort_results/GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz",
        tbi = config["output"] + "/cohort_results/GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz.tbi"
    output:
        directory(
            config["output"] + "/needLR_output/GRCh38_final_cohort_survivor_genotyped_matrix_needLR_v3.5_cohort"
        )
    singularity:
        "needlr.sif"
    threads: config["hc"]
    resources:
        mem_mb = config["hm"],
        time = config["ht"]
    shell:
        r"""
        mkdir -p {config["output"]}/needLR_output

        singularity exec {snakemake.singularity.image} \
        bash -c "
        cd /opt/needLR_v3.5_local && \
        ./needLR_v3.5_cohort.sh \
            -g {config["reference"]["grch38"]} \
            -v {input.vcf} \
            -a {config["output"]}/samples_names.txt \
            -d {threads}
        "
        """
```
