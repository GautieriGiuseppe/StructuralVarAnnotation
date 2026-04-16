import os

OUTDIR = config["output"]

PAIRS = VALID_PAIRS
BATCHES = [b for b, s in PAIRS]
SAMPLES = [s for b, s in PAIRS]

# ------------------------------------------------------------------------------
rule all_grch38_cohort_merge:
    input:
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz",
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz.tbi",
        f"{OUTDIR}/cohort_results/GRCh38_cohort_support_table.tsv"

# ------------------------------------------------------------------------------
rule pass_sniffles_grch38:
    input:
        f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/{{sample}}_grch38_sniffles.vcf"
    output:
        f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/{{sample}}_grch38_sniffles_pass.vcf"
    conda:
        "envs/truvari_env.yml"
    threads: 2
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        "bcftools view -f PASS {input} -O v -o {output}"

rule pass_delly_grch38:
    input:
        f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/{{sample}}_grch38_delly.vcf"
    output:
        f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/{{sample}}_grch38_delly_pass.vcf"
    conda:
        "envs/truvari_env.yml"
    threads: 2
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        "bcftools view -f PASS {input} -O v -o {output}"

rule pass_cutesv_grch38:
    input:
        f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/{{sample}}_grch38_cuteSV.vcf"
    output:
        f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/{{sample}}_grch38_cuteSV_pass.vcf"
    conda:
        "envs/truvari_env.yml"
    threads: 2
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        "bcftools view -f PASS {input} -O v -o {output}"

# ------------------------------------------------------------------------------
rule grch38_sample_survivor_list:
    input:
        sniffles=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/{{sample}}_grch38_sniffles_pass.vcf",
        delly=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/{{sample}}_grch38_delly_pass.vcf",
        cutesv=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/{{sample}}_grch38_cuteSV_pass.vcf"
    output:
        f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/survivor/{{sample}}_grch38_vcf_list.txt"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(output[0], "w") as out:
            out.write(f"{input.sniffles}\n{input.delly}\n{input.cutesv}\n")

# ------------------------------------------------------------------------------
rule grch38_sample_survivor_merge:
    input:
        listfile=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/survivor/{{sample}}_grch38_vcf_list.txt"
    output:
        vcf=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/survivor/{{sample}}_grch38_survivor_merged.vcf",
        table=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/survivor/{{sample}}_grch38_support_table.tsv"
    conda:
        "envs/truvari_env.yml"
    threads: 4
    resources:
        mem_mb=config["hm"],
        time=config["ht"]
    shell:
        """
        SURVIVOR merge {input.listfile} 500 1 1 0 1 30 {output.vcf}
        bcftools query -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/SUPP\t%INFO/SUPP_VEC\n' {output.vcf} > {output.table}
        """

# ------------------------------------------------------------------------------
rule grch38_cohort_survivor_list:
    input:
        expand(
            f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/survivor/{{sample}}_grch38_survivor_merged.vcf",
            zip,
            batch=BATCHES,
            sample=SAMPLES
        )
    output:
        f"{OUTDIR}/cohort_results/GRCh38_cohort_survivor_vcf_list.txt"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(output[0], "w") as out:
            for vcf in input:
                out.write(f"{vcf}\n")

# ------------------------------------------------------------------------------
rule grch38_cohort_survivor_merge:
    input:
        listfile=f"{OUTDIR}/cohort_results/GRCh38_cohort_survivor_vcf_list.txt"
    output:
        vcf=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf",
        vcfgz=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz",
        tbi=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz.tbi",
        table=f"{OUTDIR}/cohort_results/GRCh38_cohort_support_table.tsv"
    conda:
        "envs/truvari_env.yml"
    threads: 4
    resources:
        mem_mb=config["hm"],
        time=config["ht"]
    shell:
        """
        SURVIVOR merge {input.listfile} 500 1 1 0 1 30 {output.vcf}
        bcftools sort {output.vcf} -O z -o {output.vcfgz}
        tabix -p vcf {output.vcfgz}
        bcftools query -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/SUPP\t%INFO/SUPP_VEC\n' {output.vcfgz} > {output.table}
        """