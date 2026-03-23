# ==============================================================================
#  PASS FILTERING
# ==============================================================================

rule pass_filter_sniffles_chm13:
    input:
        config["output"] + "/{batch}/{sample}/03.variant_calling/{sample}_chm13_sniffles.vcf"
    output:
        config["output"] + "/{batch}/{sample}/03.variant_calling/{sample}_chm13_sniffles_pass.vcf"
    conda:
        "envs/snakemake.yml"
    threads: 2
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        "bcftools view -f PASS {input} -O v -o {output}"


rule pass_filter_delly_chm13:
    input:
        config["output"] + "/{batch}/{sample}/03.variant_calling/{sample}_chm13_delly.vcf"
    output:
        config["output"] + "/{batch}/{sample}/03.variant_calling/{sample}_chm13_delly_pass.vcf"
    conda:
        "envs/snakemake.yml"
    threads: 2
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        "bcftools view -f PASS {input} -O v -o {output}"


rule pass_filter_cutesv_chm13:
    input:
        config["output"] + "/{batch}/{sample}/03.variant_calling/{sample}_chm13_cuteSV.vcf"
    output:
        config["output"] + "/{batch}/{sample}/03.variant_calling/{sample}_chm13_cuteSV_pass.vcf"
    conda:
        "envs/snakemake.yml"
    threads: 2
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        "bcftools view -f PASS {input} -O v -o {output}"


rule pass_filter_sniffles_gr38:
    input:
        config["output"] + "/{batch}/{sample}/03.variant_calling/liftover/sniffles_{sample}_grch38-to-chm13.vcf"
    output:
        config["output"] + "/{batch}/{sample}/03.variant_calling/liftover/sniffles_{sample}_grch38-to-chm13_pass.vcf"
    conda:
        "envs/snakemake.yml"
    threads: 2
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        "bcftools view -f PASS {input} -O v -o {output}"


rule pass_filter_delly_gr38:
    input:
        config["output"] + "/{batch}/{sample}/03.variant_calling/liftover/delly_{sample}_grch38-to-chm13.vcf"
    output:
        config["output"] + "/{batch}/{sample}/03.variant_calling/liftover/delly_{sample}_grch38-to-chm13_pass.vcf"
    conda:
        "envs/snakemake.yml"
    threads: 2
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        "bcftools view -f PASS {input} -O v -o {output}"


rule pass_filter_cutesv_gr38:
    input:
        config["output"] + "/{batch}/{sample}/03.variant_calling/liftover/cuteSV_{sample}_grch38-to-chm13.vcf"
    output:
        config["output"] + "/{batch}/{sample}/03.variant_calling/liftover/cuteSV_{sample}_grch38-to-chm13_pass.vcf"
    conda:
        "envs/snakemake.yml"
    threads: 2
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        "bcftools view -f PASS {input} -O v -o {output}"


# ==============================================================================
#  SURVIVOR PER-SAMPLE MERGE
# ==============================================================================

rule survivor_sample_merge:
    input:
        sniffles_chm13=config["output"] + "/{batch}/{sample}/03.variant_calling/{sample}_chm13_sniffles_pass.vcf",
        delly_chm13=config["output"] + "/{batch}/{sample}/03.variant_calling/{sample}_chm13_delly_pass.vcf",
        cutesv_chm13=config["output"] + "/{batch}/{sample}/03.variant_calling/{sample}_chm13_cuteSV_pass.vcf",
        sniffles_gr38=config["output"] + "/{batch}/{sample}/03.variant_calling/liftover/sniffles_{sample}_grch38-to-chm13_pass.vcf",
        delly_gr38=config["output"] + "/{batch}/{sample}/03.variant_calling/liftover/delly_{sample}_grch38-to-chm13_pass.vcf",
        cutesv_gr38=config["output"] + "/{batch}/{sample}/03.variant_calling/liftover/cuteSV_{sample}_grch38-to-chm13_pass.vcf"
    output:
        merged_vcf=config["output"] + "/{batch}/{sample}/03.variant_calling/survivor/{sample}_survivor_merged.vcf",
        support_table=config["output"] + "/{batch}/{sample}/03.variant_calling/survivor/{sample}_support_table.tsv",
        file_list=config["output"] + "/{batch}/{sample}/03.variant_calling/survivor/{sample}_vcf_list.txt"
    params:
        survivor=config.get("survivor", "SURVIVOR"),
        max_dist=500,
        min_support=1,
        use_type=1,
        use_strand=0,
        use_size=1,
        min_size=30
    conda:
        "envs/snakemake.yml"
    threads: 4
    resources:
        mem_mb=config["hm"],
        time=config["ht"]
    shell:
        r"""
        mkdir -p $(dirname {output.merged_vcf})

        printf "%s\n" \
            {input.sniffles_chm13} \
            {input.delly_chm13} \
            {input.cutesv_chm13} \
            {input.sniffles_gr38} \
            {input.delly_gr38} \
            {input.cutesv_gr38} > {output.file_list}

        {params.survivor} merge \
            {output.file_list} \
            {params.max_dist} \
            {params.min_support} \
            {params.use_type} \
            {params.use_strand} \
            {params.use_size} \
            {params.min_size} \
            {output.merged_vcf}

        bcftools query \
            -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/SUPP\t%INFO/SUPP_VEC\n' \
            {output.merged_vcf} > {output.support_table}
        """


# ==============================================================================
#  COHORT MERGE
# ==============================================================================

rule survivor_cohort_merge:
    input:
        sample_vcfs=[
            config["output"] + f"/{b}/{s}/03.variant_calling/survivor/{s}_survivor_merged.vcf"
            for b, s in VALID_PAIRS
        ]
    output:
        file_list=config["output"] + "/cohort_results/survivor_cohort_vcf_list.txt",
        cohort_vcf=config["output"] + "/cohort_results/CHM13_final_cohort_survivor.vcf",
        cohort_vcf_gz=config["output"] + "/cohort_results/CHM13_final_cohort_survivor.vcf.gz",
        support_table=config["output"] + "/cohort_results/CHM13_cohort_support_table.tsv"
    params:
        survivor=config.get("survivor", "SURVIVOR"),
        max_dist=500,
        min_support=1,
        use_type=1,
        use_strand=0,
        use_size=1,
        min_size=30
    conda:
        "envs/snakemake.yml"
    threads: 4
    resources:
        mem_mb=config["hm"],
        time=config["ht"]
    shell:
        r"""
        mkdir -p $(dirname {output.cohort_vcf})

        printf "%s\n" {input.sample_vcfs} > {output.file_list}

        {params.survivor} merge \
            {output.file_list} \
            {params.max_dist} \
            {params.min_support} \
            {params.use_type} \
            {params.use_strand} \
            {params.use_size} \
            {params.min_size} \
            {output.cohort_vcf}

        bcftools sort -O z -o {output.cohort_vcf_gz} {output.cohort_vcf}
        bcftools index -t {output.cohort_vcf_gz}

        bcftools query \
            -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/SUPP\t%INFO/SUPP_VEC\n' \
            {output.cohort_vcf} > {output.support_table}
        """


# ==============================================================================
#  TIER FILTERING
# ==============================================================================

rule filter_tier1:
    input:
        config["output"] + "/cohort_results/CHM13_final_cohort_survivor.vcf.gz"
    output:
        config["output"] + "/cohort_results/CHM13_tier1_supported.vcf"
    conda:
        "envs/snakemake.yml"
    threads: 2
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        r"""
        bcftools view -h {input} > {output}
        bcftools view -H {input} | \
        awk -F'\t' 'BEGIN{{OFS="\t"}} {{
            supp="";
            suppvec="";
            n=split($8, info, ";");
            for (i=1; i<=n; i++) {{
                if (info[i] ~ /^SUPP=/) {{
                    split(info[i], a, "=");
                    supp=a[2];
                }}
                if (info[i] ~ /^SUPP_VEC=/) {{
                    split(info[i], b, "=");
                    suppvec=b[2];
                }}
            }}
            if ((supp+0) > 1 && substr(suppvec,1,3) ~ /1/) print $0;
        }}' >> {output}
        """


rule filter_tier2:
    input:
        config["output"] + "/cohort_results/CHM13_final_cohort_survivor.vcf.gz"
    output:
        config["output"] + "/cohort_results/CHM13_genotype_input.vcf"
    conda:
        "envs/snakemake.yml"
    threads: 2
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        r"""
        bcftools view -h {input} > {output}
        bcftools view -H {input} | \
        awk -F'\t' 'BEGIN{{OFS="\t"}} {{
            supp="";
            n=split($8, info, ";");
            for (i=1; i<=n; i++) {{
                if (info[i] ~ /^SUPP=/) {{
                    split(info[i], a, "=");
                    supp=a[2];
                }}
            }}
            if ((supp+0) > 1) print $0;
        }}' >> {output}
        """


# ==============================================================================
#  GENOTYPING
# ==============================================================================

rule genotype:
    input:
        bam=config["output"] + "/{batch}/{sample}/01.align/chm13/{sample}.srt.bam",
        ref=config["reference"]["chm13"],
        sites=config["output"] + "/cohort_results/{cohort}.vcf"
    output:
        vcf=config["output"] + "/{batch}/{sample}/04.force_calling/{sample}_{cohort}_genotyped.vcf"
    conda:
        "envs/sniffles.yml"
    threads: config["hc"]
    resources:
        mem_mb=config["hm"],
        time=config["ht"]
    shell:
        r"""
        mkdir -p $(dirname {output.vcf})

        sniffles \
            --input {input.bam} \
            --reference {input.ref} \
            --genotype-vcf {input.sites} \
            --vcf {output.vcf} \
            --sample {wildcards.sample} \
            --threads {threads}
        """


# ==============================================================================
#  FIX BROKEN SNIFFLES GENOTYPE HEADER
# ==============================================================================

rule fix_genotype:
    input:
        config["output"] + "/{batch}/{sample}/04.force_calling/{sample}_{cohort}_genotyped.vcf"
    output:
        config["output"] + "/{batch}/{sample}/04.force_calling/{sample}_{cohort}_genotyped.fixed.vcf"
    conda:
        "envs/snakemake.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        r"""
        python3 fix_genotype_header.py \
            {input} \
            {output} \
            {wildcards.sample}
        """


# ==============================================================================
#  COMPRESS FIXED GENOTYPES
# ==============================================================================

rule compress:
    input:
        config["output"] + "/{batch}/{sample}/04.force_calling/{sample}_{cohort}_genotyped.fixed.vcf"
    output:
        bcf=config["output"] + "/{batch}/{sample}/04.force_calling/{sample}_{cohort}_genotyped.bcf",
        csi=config["output"] + "/{batch}/{sample}/04.force_calling/{sample}_{cohort}_genotyped.bcf.csi"
    conda:
        "envs/snakemake.yml"
    threads: 4
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        r"""
        bcftools view -e 'ALT="."' -O u {input} | \
        bcftools sort -O b -o {output.bcf}
        bcftools index --threads {threads} {output.bcf}
        """


# ==============================================================================
#  MERGE GENOTYPED COHORT
# ==============================================================================

rule merge_genotypes:
    input:
        bcfs=lambda wc: [
            config["output"] + f"/{b}/{s}/04.force_calling/{s}_{wc.cohort}_genotyped.bcf"
            for b, s in VALID_PAIRS
        ]
    output:
        vcf=config["output"] + "/cohort_results/{cohort}_genotyped_matrix.vcf.gz",
        tbi=config["output"] + "/cohort_results/{cohort}_genotyped_matrix.vcf.gz.tbi"
    conda:
        "envs/snakemake.yml"
    threads: 4
    resources:
        mem_mb=config["hm"],
        time=config["ht"]
    shell:
        r"""
        bcftools merge -m none -O z -o {output.vcf} {input.bcfs}
        bcftools index -t {output.vcf}
        """


# ==============================================================================
#  GENOTYPABILITY
# ==============================================================================

rule genotypability:
    input:
        config["output"] + "/{batch}/{sample}/04.force_calling/{sample}_{cohort}_genotyped.fixed.vcf"
    output:
        config["output"] + "/{batch}/{sample}/04.force_calling/{sample}_{cohort}_genotypability.txt"
    conda:
        "envs/snakemake.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        r"""
        bcftools view -e 'ALT="."' -H {input} | \
        awk 'BEGIN{{valid=0; total=0}} {{
            split($10, f, ":");
            if (f[1] ~ /^[0-9.]+\/[0-9.]+$/ && f[1] != "./.") valid++;
            total++;
        }} END {{
            if (total > 0) print valid/total;
            else print 0;
        }}' > {output}
        """