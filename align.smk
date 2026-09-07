import os

samples_file = config["samples"]
OUTDIR = config["output"].rstrip("/")

print(f"Loading samples from: {samples_file}")

VALID_PAIRS = []
input_read_map = {}
input_read_type_map = {}


def is_missing_path(x):
    return x is None or str(x).strip() in {"", "NA", "NaN", "nan", "None", "."}


def detect_read_type(path):
    p = str(path).lower()

    # If multiple files are comma-separated, inspect the first one.
    first = p.split(",")[0].strip()

    if first.endswith((".bam", ".ubam", ".cram")):
        return "bam"

    if first.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
        return "fastq"

    return "fastq"


with open(samples_file, "r") as f:
    header = f.readline().strip().split("\t")
    sample_idx = header.index("sample_id")
    batch_idx = header.index("batch_id")

    ubam_idx = header.index("ubam") if "ubam" in header else None
    fastq_idx = header.index("fastq") if "fastq" in header else None

    for line in f:
        if line.strip():
            fields = line.strip().split("\t")
            sample = fields[sample_idx]
            batch = fields[batch_idx]

            ubam = fields[ubam_idx] if ubam_idx is not None else ""
            fastq = fields[fastq_idx] if fastq_idx is not None else ""

            VALID_PAIRS.append((batch, sample))

            if not is_missing_path(ubam):
                input_read_map[(batch, sample)] = [ubam]
                input_read_type_map[(batch, sample)] = detect_read_type(ubam)

            elif not is_missing_path(fastq):
                fastq_files = [x.strip() for x in fastq.split(",") if x.strip()]
                input_read_map[(batch, sample)] = fastq_files
                input_read_type_map[(batch, sample)] = "fastq"
            else:
                raise ValueError(f"No ubam or fastq provided for {batch}/{sample}")


SAMPLES = sorted(set(sample for batch, sample in VALID_PAIRS))
BATCHES = sorted(set(batch for batch, sample in VALID_PAIRS))

print(f"Loaded {len(SAMPLES)} samples: {SAMPLES}")
print(f"Batches: {BATCHES}")


def get_reads(wc):
    return input_read_map[(wc.batch, wc.sample)]


def get_read_type(wc):
    return input_read_type_map[(wc.batch, wc.sample)]


rule all_align:
    input:
        grch38_bam=[
            f"{OUTDIR}/{batch}/{sample}/01.align/grch38/{sample}.srt.bam"
            for batch, sample in VALID_PAIRS
        ],
        grch38_bai=[
            f"{OUTDIR}/{batch}/{sample}/01.align/grch38/{sample}.srt.bam.bai"
            for batch, sample in VALID_PAIRS
        ],
        chm13_bam=[
            f"{OUTDIR}/{batch}/{sample}/01.align/chm13/{sample}.srt.bam"
            for batch, sample in VALID_PAIRS
        ],
        chm13_bai=[
            f"{OUTDIR}/{batch}/{sample}/01.align/chm13/{sample}.srt.bam.bai"
            for batch, sample in VALID_PAIRS
        ]


rule minimap2_GRCh38:
    input:
        reference=config["reference"]["grch38"],
        reads=get_reads
    output:
        bam=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/grch38/{{sample}}.srt.bam",
        bai=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/grch38/{{sample}}.srt.bam.bai"
    threads:
        config["mc"]
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["vhm"],
        time=lambda wildcards, attempt: attempt * config["vht"]
    params:
        read_group=lambda wc: f"@RG\\tID:{wc.sample}\\tPL:ONT\\tSM:{wc.sample}",
        prefix=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/grch38/{{sample}}",
        read_type=get_read_type
    conda:
        "variant_calling/envs/snakemake.yml"
    shell:
        r"""
        mkdir -p $(dirname {output.bam})

        if [[ "{params.read_type}" == "bam" ]]; then
            samtools fastq \
                -T MM,ML \
                {input.reads} | \
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
                -m 1G \
                -@ {threads} \
                -o {output.bam} \
                -T {params.prefix}.${{SLURM_JOB_ID:-$$}} \
                -
        else
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
                {input.reads} | \
            samtools sort \
                -m 1G \
                -@ {threads} \
                -o {output.bam} \
                -T {params.prefix}.${{SLURM_JOB_ID:-$$}} \
                -
        fi

        samtools index {output.bam}
        test -s {output.bam}
        test -s {output.bai}
        """

rule minimap2_chm13:
    input:
        reference=config["reference"]["chm13"],
        reads=get_reads
    output:
        bam=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/chm13/{{sample}}.srt.bam",
        bai=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/chm13/{{sample}}.srt.bam.bai"
    threads:
        config["mc"]
    resources:
        mem_mb=lambda wildcards, attempt: attempt * config["vhm"],
        time=lambda wildcards, attempt: attempt * config["vht"]
    params:
        read_group=lambda wc: f"@RG\\tID:{wc.sample}\\tPL:ONT\\tSM:{wc.sample}",
        prefix=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/chm13/{{sample}}",
        read_type=get_read_type
    conda:
        "variant_calling/envs/snakemake.yml"
    shell:
        r"""
        mkdir -p $(dirname {output.bam})

        if [[ "{params.read_type}" == "bam" ]]; then
            samtools fastq \
                -T MM,ML \
                {input.reads} | \
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
                -m 1G \
                -@ {threads} \
                -o {output.bam} \
                -T {params.prefix}.${{SLURM_JOB_ID:-$$}} \
                -
        else
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
                {input.reads} | \
            samtools sort \
                -m 1G \
                -@ {threads} \
                -o {output.bam} \
                -T {params.prefix}.${{SLURM_JOB_ID:-$$}} \
                -
        fi

        samtools index {output.bam}
        test -s {output.bam}
        test -s {output.bai}
        """