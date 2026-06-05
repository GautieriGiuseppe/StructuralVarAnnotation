import os

# ==============================================================================
# 1. PARSING CONFIG & METADATA
# ==============================================================================

samples_file = config["samples"]
OUTDIR = config["output"]

print(f"Loading samples from: {samples_file}")

VALID_PAIRS = []

with open(samples_file, "r") as f:
    header = f.readline().strip().split("\t")

    sample_idx = header.index("sample_id")
    batch_idx = header.index("batch_id")

    for line in f:
        if line.strip():
            fields = line.strip().split("\t")
            sample = fields[sample_idx]
            batch = fields[batch_idx]

            VALID_PAIRS.append((batch, sample))

print(f"Loaded {len(VALID_PAIRS)} sample/batch pairs.")


# ==============================================================================
# 2. RULE ALL
# ==============================================================================

rule all_alignqc:
    input:
        grch38_alfred=[
            f"{OUTDIR}/{batch}/{sample}/02.alignqc/grch38/{sample}.alfred.tsv.gz"
            for batch, sample in VALID_PAIRS
        ],
        grch38_mosdepth=[
            f"{OUTDIR}/{batch}/{sample}/02.alignqc/grch38/{sample}_grch38.mosdepth.global.dist.txt"
            for batch, sample in VALID_PAIRS
        ],
        chm13_alfred=[
            f"{OUTDIR}/{batch}/{sample}/02.alignqc/chm13/{sample}.alfred.tsv.gz"
            for batch, sample in VALID_PAIRS
        ],
        chm13_mosdepth=[
            f"{OUTDIR}/{batch}/{sample}/02.alignqc/chm13/{sample}_chm13.mosdepth.global.dist.txt"
            for batch, sample in VALID_PAIRS
        ]


# ==============================================================================
# 3. ANALYSIS RULES
# ==============================================================================

rule alfred_stats_grch38:
    input:
        reference=config["reference"]["grch38"],
        bam=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/grch38/{{sample}}.srt.bam"
    output:
        f"{OUTDIR}/{{batch}}/{{sample}}/02.alignqc/grch38/{{sample}}.alfred.tsv.gz"
    threads:
        config["lc"]
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["mm"],
        time=lambda wildcards, attempt: attempt * config["ht"]
    params:
        sample=lambda wc: wc.sample
    conda:
        "variant_calling/envs/snakemake.yml"
    shell:
        r"""
        mkdir -p $(dirname {output})

        alfred qc \
            -r {input.reference} \
            -a {params.sample} \
            -o {output} \
            {input.bam}
        """


rule mosdepth_stats_grch38:
    input:
        bam=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/grch38/{{sample}}.srt.bam"
    output:
        f"{OUTDIR}/{{batch}}/{{sample}}/02.alignqc/grch38/{{sample}}_grch38.mosdepth.global.dist.txt"
    threads:
        config["lc"]
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["mm"],
        time=lambda wildcards, attempt: attempt * config["mt"]
    params:
        prefix=f"{OUTDIR}/{{batch}}/{{sample}}/02.alignqc/grch38/{{sample}}_grch38"
    conda:
        "variant_calling/envs/snakemake.yml"
    shell:
        r"""
        mkdir -p $(dirname {output})

        mosdepth \
            {params.prefix} \
            {input.bam} \
            -n \
            -x \
            --by 1000
        """


rule alfred_stats_chm13:
    input:
        reference=config["reference"]["chm13"],
        bam=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/chm13/{{sample}}.srt.bam"
    output:
        f"{OUTDIR}/{{batch}}/{{sample}}/02.alignqc/chm13/{{sample}}.alfred.tsv.gz"
    threads:
        config["lc"]
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["mm"],
        time=lambda wildcards, attempt: attempt * config["ht"]
    params:
        sample=lambda wc: wc.sample
    conda:
        "variant_calling/envs/snakemake.yml"
    shell:
        r"""
        mkdir -p $(dirname {output})

        alfred qc \
            -r {input.reference} \
            -a {params.sample} \
            -o {output} \
            {input.bam}
        """


rule mosdepth_stats_chm13:
    input:
        bam=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/chm13/{{sample}}.srt.bam"
    output:
        f"{OUTDIR}/{{batch}}/{{sample}}/02.alignqc/chm13/{{sample}}_chm13.mosdepth.global.dist.txt"
    threads:
        config["lc"]
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["mm"],
        time=lambda wildcards, attempt: attempt * config["mt"]
    params:
        prefix=f"{OUTDIR}/{{batch}}/{{sample}}/02.alignqc/chm13/{{sample}}_chm13"
    conda:
        "variant_calling/envs/snakemake.yml"
    shell:
        r"""
        mkdir -p $(dirname {output})

        mosdepth \
            {params.prefix} \
            {input.bam} \
            -n \
            -x \
            --by 1000
        """