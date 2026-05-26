import os

OUTDIR = config["output"]
GRCH38_REF = config["reference"]["grch38"]

PAIRS = VALID_PAIRS
BATCHES = [b for b, s in PAIRS]
SAMPLES = [s for b, s in PAIRS]

# ------------------------------------------------------------------------------
rule all_grch38_force_genotype:
    input:
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz",
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz.tbi"

# ------------------------------------------------------------------------------
rule grch38_force_genotype:
    input:
        bam=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/grch38/{{sample}}.srt.bam",
        vcf=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz"
    output:
        f"{OUTDIR}/{{batch}}/{{sample}}/04.force_calling/{{sample}}_GRCh38_final_cohort_survivor_genotyped.vcf"
    conda:
        "variant_calling/envs/sniffles.yml"
    threads: config["hc"]
    resources:
        mem_mb=config["hm"],
        time=config["ht"]
    shell:
        """
        sniffles \
          --input {input.bam} \
          --reference {GRCH38_REF} \
          --genotype-vcf {input.vcf} \
          --vcf {output} \
          --sample {wildcards.sample} \
          --threads {threads}
        """

# ------------------------------------------------------------------------------
rule grch38_fix_genotype_header:
    input:
        f"{OUTDIR}/{{batch}}/{{sample}}/04.force_calling/{{sample}}_GRCh38_final_cohort_survivor_genotyped.vcf"
    output:
        f"{OUTDIR}/{{batch}}/{{sample}}/04.force_calling/{{sample}}_GRCh38_final_cohort_survivor_genotyped.fixed.vcf"
    conda:
        "variant_calling/envs/sniffles.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        """
        python variant_calling/fix_genotype_header.py {input} {output} {wildcards.sample}
        """

# ------------------------------------------------------------------------------
rule grch38_compress_genotyped:
    input:
        f"{OUTDIR}/{{batch}}/{{sample}}/04.force_calling/{{sample}}_GRCh38_final_cohort_survivor_genotyped.fixed.vcf"
    output:
        bcf=f"{OUTDIR}/{{batch}}/{{sample}}/04.force_calling/{{sample}}_GRCh38_final_cohort_survivor_genotyped.bcf",
        csi=f"{OUTDIR}/{{batch}}/{{sample}}/04.force_calling/{{sample}}_GRCh38_final_cohort_survivor_genotyped.bcf.csi"
    conda:
        "variant_calling/envs/truvari_env.yml"
    threads: 4
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        """
        bcftools view -e 'ALT="."' -O b -o {output.bcf} {input}
        bcftools index --threads {threads} {output.bcf}
        """

# ------------------------------------------------------------------------------
rule grch38_merge_genotyped_matrix:
    input:
        expand(
            f"{OUTDIR}/{{batch}}/{{sample}}/04.force_calling/{{sample}}_GRCh38_final_cohort_survivor_genotyped.bcf",
            zip,
            batch=BATCHES,
            sample=SAMPLES
        )
    output:
        vcf=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz",
        tbi=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz.tbi"
    conda:
        "variant_calling/envs/truvari_env.yml"
    threads: 4
    resources:
        mem_mb=config["hm"],
        time=config["ht"]
    shell:
        """
        bcftools merge -m none -O z -o {output.vcf} {input}
        tabix -p vcf {output.vcf}
        """