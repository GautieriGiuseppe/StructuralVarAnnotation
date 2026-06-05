import os

samples_file = config["samples"]
OUTDIR = config["output"]

print(f"Loading samples from: {samples_file}")

VALID_PAIRS = []
ubam_map = {}

with open(samples_file, "r") as f:
    header = f.readline().strip().split("\t")

    sample_idx = header.index("sample_id")
    ubam_idx = header.index("ubam")
    batch_idx = header.index("batch_id")

    for line in f:
        if line.strip():
            fields = line.strip().split("\t")

            sample = fields[sample_idx]
            ubam = fields[ubam_idx]
            batch = fields[batch_idx]

            VALID_PAIRS.append((batch, sample))
            ubam_map[(batch, sample)] = ubam

            print(f"  {batch}/{sample} -> {ubam}")

SAMPLES = sorted(set(sample for batch, sample in VALID_PAIRS))
BATCHES = sorted(set(batch for batch, sample in VALID_PAIRS))

print(f"Loaded {len(SAMPLES)} samples: {SAMPLES}")
print(f"Batches: {BATCHES}")


rule all_align:
    input:
        grch38=[
            f"{OUTDIR}/{batch}/{sample}/01.align/grch38/{sample}.srt.bam"
            for batch, sample in VALID_PAIRS
        ],
        chm13=[
            f"{OUTDIR}/{batch}/{sample}/01.align/chm13/{sample}.srt.bam"
            for batch, sample in VALID_PAIRS
        ]

rule minimap2_GRCh38:
    input:
        reference=config["reference"]["grch38"],
        sample=lambda wc: ubam_map[(wc.batch, wc.sample)]
    output:
        f"{OUTDIR}/{{batch}}/{{sample}}/01.align/grch38/{{sample}}.srt.bam"
    threads:
        config["mc"]
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["hm"],
        time=lambda wildcards, attempt: attempt * config["vht"]
    params:
        read_group=lambda wc: f"@RG\\tID:{wc.sample}\\tPL:ONT\\tSM:{wc.sample}",
        prefix=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/grch38/{{sample}}"
    conda:
        "variant_calling/envs/snakemake.yml"
    shell:
        r"""
        mkdir -p $(dirname {output})

        samtools fastq \
            -T MM,ML \
            {input.sample} | \
        minimap2 \
            -a \
            -x map-ont \
            -t {threads} \
            --MD \
            -y \
            --rmq=yes \
            --cs \
            -R '{params.read_group}' \
            {input.reference} \
            - | \
        samtools sort \
            -m 3G \
            -@ {threads} \
            --write-index \
            -o {output} \
            -T {params.prefix} \
            -
        """

rule minimap2_chm13:
    input:
        reference=config["reference"]["chm13"],
        sample=lambda wc: ubam_map[(wc.batch, wc.sample)]
    output:
        f"{OUTDIR}/{{batch}}/{{sample}}/01.align/chm13/{{sample}}.srt.bam"
    threads:
        config["mc"]
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["hm"],
        time=lambda wildcards, attempt: attempt * config["vht"]
    params:
        read_group=lambda wc: f"@RG\\tID:{wc.sample}\\tPL:ONT\\tSM:{wc.sample}",
        prefix=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/chm13/{{sample}}"
    conda:
        "variant_calling/envs/snakemake.yml"
    shell:
        r"""
        mkdir -p $(dirname {output})

        samtools fastq \
            -T MM,ML \
            {input.sample} | \
        minimap2 \
            -a \
            -x map-ont \
            -t {threads} \
            --MD \
            -y \
            --rmq=yes \
            --cs \
            -R '{params.read_group}' \
            {input.reference} \
            - | \
        samtools sort \
            -m 3G \
            -@ {threads} \
            --write-index \
            -o {output} \
            -T {params.prefix} \
            -
        """