import os

OUTDIR = config["output"]

PAIRS = VALID_PAIRS

HAPLOTYPE_DIR = f"{OUTDIR}/cohort_results/haplotypes_grch38"

GRCH38_GENOTYPED_VCF = (
    f"{OUTDIR}/cohort_results/"
    "GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz"
)

GRCH38_GENOTYPED_TBI = (
    f"{OUTDIR}/cohort_results/"
    "GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz.tbi"
)

GRCH38_REF = config["reference_uncompressed"]["grch38"]


# ==============================================================================
# Main target for haplotyping
# ==============================================================================

rule all_haplotype_grch38:
    input:
        phased_vcfs=[
            f"{OUTDIR}/{batch}/{sample}/05.haplotyping/grch38/{sample}.whatshap.phased.vcf.gz"
            for batch, sample in PAIRS
        ],
        phased_tbis=[
            f"{OUTDIR}/{batch}/{sample}/05.haplotyping/grch38/{sample}.whatshap.phased.vcf.gz.tbi"
            for batch, sample in PAIRS
        ],
        haplotagged_bams=[
            f"{OUTDIR}/{batch}/{sample}/05.haplotyping/grch38/{sample}.haplotagged.bam"
            for batch, sample in PAIRS
        ],
        haplotagged_bais=[
            f"{OUTDIR}/{batch}/{sample}/05.haplotyping/grch38/{sample}.haplotagged.bam.bai"
            for batch, sample in PAIRS
        ],
        haplotag_lists=[
            f"{OUTDIR}/{batch}/{sample}/05.haplotyping/grch38/{sample}.haplotag_list.tsv"
            for batch, sample in PAIRS
        ]


# ==============================================================================
# Subset cohort genotyped VCF to one sample
# ==============================================================================

rule subset_sample_grch38_genotyped_vcf_for_whatshap:
    input:
        vcf=GRCH38_GENOTYPED_VCF,
        tbi=GRCH38_GENOTYPED_TBI
    output:
        vcf=f"{OUTDIR}/{{batch}}/{{sample}}/05.haplotyping/grch38/{{sample}}.unphased.input.vcf.gz",
        tbi=f"{OUTDIR}/{{batch}}/{{sample}}/05.haplotyping/grch38/{{sample}}.unphased.input.vcf.gz.tbi"
    conda:
        "variant_calling/envs/whatshap.yml"
    threads:
        1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    shell:
        r"""
        mkdir -p $(dirname {output.vcf})

        bcftools view \
            -s {wildcards.sample} \
            -O z \
            -o {output.vcf} \
            {input.vcf}

        tabix -f -p vcf {output.vcf}
        """


# ==============================================================================
# WhatsHap phase
# ==============================================================================

rule whatshap_phase_grch38:
    input:
        vcf=f"{OUTDIR}/{{batch}}/{{sample}}/05.haplotyping/grch38/{{sample}}.unphased.input.vcf.gz",
        tbi=f"{OUTDIR}/{{batch}}/{{sample}}/05.haplotyping/grch38/{{sample}}.unphased.input.vcf.gz.tbi",
        bam=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/grch38/{{sample}}.srt.bam",
        ref=GRCH38_REF
    output:
        vcf=f"{OUTDIR}/{{batch}}/{{sample}}/05.haplotyping/grch38/{{sample}}.whatshap.phased.vcf.gz",
        tbi=f"{OUTDIR}/{{batch}}/{{sample}}/05.haplotyping/grch38/{{sample}}.whatshap.phased.vcf.gz.tbi"
    conda:
        "variant_calling/envs/whatshap.yml"
    threads:
        config["lc"]
    resources:
        mem_mb=config["hm"],
        time=config["ht"]
    params:
        tmp_vcf=f"{OUTDIR}/{{batch}}/{{sample}}/05.haplotyping/grch38/{{sample}}.whatshap.phased.tmp.vcf"
    shell:
        r"""
        mkdir -p $(dirname {output.vcf})

        if [ ! -f {input.bam}.bai ] && [ ! -f {input.bam}.csi ]; then
            samtools index {input.bam}
        fi

        whatshap phase \
            --reference {input.ref} \
            --sample {wildcards.sample} \
            --output {params.tmp_vcf} \
            {input.vcf} \
            {input.bam}

        bgzip -f -c {params.tmp_vcf} > {output.vcf}
        tabix -f -p vcf {output.vcf}

        rm -f {params.tmp_vcf}
        """


# ==============================================================================
# WhatsHap haplotag BAM
# ==============================================================================

rule whatshap_haplotag_grch38:
    input:
        vcf=f"{OUTDIR}/{{batch}}/{{sample}}/05.haplotyping/grch38/{{sample}}.whatshap.phased.vcf.gz",
        tbi=f"{OUTDIR}/{{batch}}/{{sample}}/05.haplotyping/grch38/{{sample}}.whatshap.phased.vcf.gz.tbi",
        bam=f"{OUTDIR}/{{batch}}/{{sample}}/01.align/grch38/{{sample}}.srt.bam",
        ref=GRCH38_REF
    output:
        bam=f"{OUTDIR}/{{batch}}/{{sample}}/05.haplotyping/grch38/{{sample}}.haplotagged.bam",
        bai=f"{OUTDIR}/{{batch}}/{{sample}}/05.haplotyping/grch38/{{sample}}.haplotagged.bam.bai",
        haplotag_list=f"{OUTDIR}/{{batch}}/{{sample}}/05.haplotyping/grch38/{{sample}}.haplotag_list.tsv"
    conda:
        "variant_calling/envs/whatshap.yml"
    threads:
        config["lc"]
    resources:
        mem_mb=config["hm"],
        time=config["ht"]
    shell:
        r"""
        mkdir -p $(dirname {output.bam})

        if [ ! -f {input.bam}.bai ] && [ ! -f {input.bam}.csi ]; then
            samtools index {input.bam}
        fi

        whatshap haplotag \
            --reference {input.ref} \
            --output-haplotag-list {output.haplotag_list} \
            -o {output.bam} \
            {input.vcf} \
            {input.bam}

        samtools index {output.bam}
        """