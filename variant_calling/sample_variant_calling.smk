import os

OUTDIR = config["output"]

# ==============================================================================
#  VARIANT CALLING RULES
# ==============================================================================

# ------------------------------------------------------------------------------
# SNIFFLES
# ------------------------------------------------------------------------------

rule sniffles_grch38:
    input:
        sample=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/grch38/{{sample}}.srt.bam"
    output:
        vcf=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/{{sample}}_grch38_sniffles.vcf",
        snf=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/{{sample}}_grch38_sniffles.snf"
    conda:
        "envs/sniffles.yml"
    threads:
        config["mc"]
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        r"""
        mkdir -p $(dirname {output.vcf})

        sniffles \
            --input {input.sample} \
            --vcf {output.vcf} \
            --snf {output.snf} \
            --threads {threads}
        """


rule sniffles_chm13:
    input:
        sample=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/chm13/{{sample}}.srt.bam"
    output:
        vcf=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/{{sample}}_chm13_sniffles.vcf",
        snf=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/{{sample}}_chm13_sniffles.snf"
    conda:
        "envs/sniffles.yml"
    threads:
        config["mc"]
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        r"""
        mkdir -p $(dirname {output.vcf})

        sniffles \
            --input {input.sample} \
            --vcf {output.vcf} \
            --snf {output.snf} \
            --threads {threads}
        """


# ------------------------------------------------------------------------------
# DELLY
# ------------------------------------------------------------------------------

rule delly_grch38:
    input:
        reference=config["reference"]["grch38"],
        sample=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/grch38/{{sample}}.srt.bam"
    output:
        bcf=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/{{sample}}_grch38_delly.bcf"
    conda:
        "envs/delly.yml"
    threads:
        config["mc"]
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        r"""
        mkdir -p $(dirname {output.bcf})

        delly lr \
            -y ont \
            -g {input.reference} \
            {input.sample} \
            -o {output.bcf}
        """


rule delly_chm13:
    input:
        reference=config["reference"]["chm13"],
        sample=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/chm13/{{sample}}.srt.bam"
    output:
        bcf=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/{{sample}}_chm13_delly.bcf"
    conda:
        "envs/delly.yml"
    threads:
        config["lc"]
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        r"""
        mkdir -p $(dirname {output.bcf})

        delly lr \
            -y ont \
            -g {input.reference} \
            {input.sample} \
            -o {output.bcf}
        """


# ------------------------------------------------------------------------------
# CUTESV
# ------------------------------------------------------------------------------

rule cuteSV_grch38:
    input:
        bam=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/grch38/{{sample}}.srt.bam",
        ref=config["reference"]["grch38"]
    output:
        vcf=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/{{sample}}_grch38_cuteSV.vcf"
    conda:
        "envs/cutesv.yml"
    threads:
        config["hc"]
    resources:
        mem_mb=config["hm"],
        time=config["ht"]
    shell:
        r"""
        mkdir -p $(dirname {output.vcf})

        TEMP_DIR="$(dirname {output.vcf})/temp_cutesv_grch38_{wildcards.sample}"

        rm -rf "$TEMP_DIR"
        mkdir -p "$TEMP_DIR"

        cuteSV \
            --threads {threads} \
            --sample {wildcards.sample} \
            --genotype \
            {input.bam} \
            {input.ref} \
            {output.vcf} \
            "$TEMP_DIR"

        rm -rf "$TEMP_DIR"
        """


rule cuteSV_chm13:
    input:
        bam=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/chm13/{{sample}}.srt.bam",
        ref=config["reference"]["chm13"]
    output:
        vcf=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/{{sample}}_chm13_cuteSV.vcf"
    conda:
        "envs/cutesv.yml"
    threads:
        config["hc"]
    resources:
        mem_mb=config["vhm"],
        time=config["vht"]
    shell:
        r"""
        mkdir -p $(dirname {output.vcf})

        TEMP_DIR="$(dirname {output.vcf})/temp_cutesv_chm13_{wildcards.sample}"

        rm -rf "$TEMP_DIR"
        mkdir -p "$TEMP_DIR"

        cuteSV \
            --threads {threads} \
            --sample {wildcards.sample} \
            --genotype \
            {input.bam} \
            {input.ref} \
            {output.vcf} \
            "$TEMP_DIR"

        rm -rf "$TEMP_DIR"
        """


# ==============================================================================
#  BCFTOOLS / FORMAT CONVERSION
# ==============================================================================

rule bcf_to_vcf:
    input:
        "{path}/{sample}_{ref}_{tool}.bcf"
    output:
        "{path}/{sample}_{ref}_{tool}.vcf"
    wildcard_constraints:
        tool="delly"
    conda:
        "envs/snakemake.yml"
    threads:
        2
    resources:
        mem_mb=config["hm"],
        time=config["ht"]
    shell:
        r"""
        mkdir -p $(dirname {output})

        bcftools view \
            -O v \
            -o {output} \
            {input}
        """


# ==============================================================================
#  DELLY SYMBOLIC CONVERSION
# ==============================================================================

rule delly_to_symbolic:
    input:
        vcf=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/{{sample}}_{{ref}}_delly.vcf"
    output:
        vcf=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/{{sample}}_{{ref}}_delly_symbolic.vcf"
    conda:
        "envs/snakemake.yml"
    threads:
        1
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        r"""
        mkdir -p $(dirname {output.vcf})

        python3 variant_calling/delly_to_symbolic.py \
            {input.vcf} \
            {output.vcf}
        """