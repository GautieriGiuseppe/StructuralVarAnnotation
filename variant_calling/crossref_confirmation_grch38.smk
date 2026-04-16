import os

OUTDIR = config["output"]

# ------------------------------------------------------------------------------
# Convenience target
# ------------------------------------------------------------------------------
rule all_grch38_confirmation:
    input:
        f"{OUTDIR}/cohort_results/CHM13_final_cohort_survivor_to_GRCh38.vcf.gz",
        f"{OUTDIR}/cohort_results/CHM13_final_cohort_survivor_to_GRCh38.vcf.gz.tbi",
        f"{OUTDIR}/cohort_results/grch38_vs_chm13lifted_truvari/summary.txt",
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_confirmation.tsv"

# ------------------------------------------------------------------------------
# Lift CHM13 cohort VCF to GRCh38
# ------------------------------------------------------------------------------
rule crossmap_chm13_cohort_to_grch38:
    input:
        vcf=f"{OUTDIR}/cohort_results/CHM13_final_cohort_survivor.vcf.gz",
        tbi=f"{OUTDIR}/cohort_results/CHM13_final_cohort_survivor.vcf.gz.tbi",
        chain=config["chain_file"],
        ref=config["reference_uncompressed"]["grch38"]
    output:
        vcf=f"{OUTDIR}/cohort_results/CHM13_final_cohort_survivor_to_GRCh38.raw.vcf"
    conda:
        "envs/crossmap_truvari.yml"
    threads: 2
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        r"""
        mkdir -p $(dirname {output.vcf})

        CrossMap vcf \
            {input.chain} \
            {input.vcf} \
            {input.ref} \
            {output.vcf}
        """

# ------------------------------------------------------------------------------
# Sort + bgzip + index lifted VCF
# ------------------------------------------------------------------------------
rule sort_index_chm13_lifted_grch38:
    input:
        vcf=f"{OUTDIR}/cohort_results/CHM13_final_cohort_survivor_to_GRCh38.raw.vcf"
    output:
        vcf=f"{OUTDIR}/cohort_results/CHM13_final_cohort_survivor_to_GRCh38.vcf.gz",
        tbi=f"{OUTDIR}/cohort_results/CHM13_final_cohort_survivor_to_GRCh38.vcf.gz.tbi"
    conda:
        "envs/crossmap_truvari.yml"
    threads: 4
    resources:
        mem_mb=config["hm"],
        time=config["ht"]
    shell:
        r"""
        bcftools sort {input.vcf} -O z -o {output.vcf}
        tabix -f -p vcf {output.vcf}
        """

# ------------------------------------------------------------------------------
# Compare native GRCh38 vs lifted CHM13 using truvari
#
# base = native GRCh38 set
# comp = CHM13 lifted to GRCh38
# ------------------------------------------------------------------------------
rule truvari_match_grch38_vs_chm13lifted:
    input:
        base=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz",
        comp=f"{OUTDIR}/cohort_results/CHM13_final_cohort_survivor_to_GRCh38.vcf.gz",
        ref=config["reference"]["grch38"]
    output:
        summary=f"{OUTDIR}/cohort_results/grch38_vs_chm13lifted_truvari/summary.txt",
        tp_base=f"{OUTDIR}/cohort_results/grch38_vs_chm13lifted_truvari/tp-base.vcf.gz",
        tp_comp=f"{OUTDIR}/cohort_results/grch38_vs_chm13lifted_truvari/tp-comp.vcf.gz",
        fn=f"{OUTDIR}/cohort_results/grch38_vs_chm13lifted_truvari/fn.vcf.gz",
        fp=f"{OUTDIR}/cohort_results/grch38_vs_chm13lifted_truvari/fp.vcf.gz"
    conda:
        "envs/crossmap_truvari.yml"
    threads: 4
    resources:
        mem_mb=config["hm"],
        time=config["ht"]
    shell:
        r"""
        rm -rf {OUTDIR}/cohort_results/grch38_vs_chm13lifted_truvari

        truvari bench \
            -b {input.base} \
            -c {input.comp} \
            -f {input.ref} \
            -o {OUTDIR}/cohort_results/grch38_vs_chm13lifted_truvari \
            --passonly \
            --pctseq 0 \
            --pctsize 0 \
            --pctovl 0 \
            -r 500 \
            --multimatch \
        """

# ------------------------------------------------------------------------------
# Build confirmation table for native GRCh38 variants
#
# tp-base = GRCh38 variants matched by lifted CHM13
# fn      = GRCh38 variants not matched by lifted CHM13
# ------------------------------------------------------------------------------
rule build_grch38_confirmation_table:
    input:
        tp_base=f"{OUTDIR}/cohort_results/grch38_vs_chm13lifted_truvari/tp-base.vcf.gz",
        fn=f"{OUTDIR}/cohort_results/grch38_vs_chm13lifted_truvari/fn.vcf.gz"
    output:
        tsv=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_confirmation.tsv"
    conda:
        "envs/crossmap_truvari.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    shell:
        r"""
        (
            bcftools query -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\tconfirmed_by_CHM13\n' {input.tp_base}
            bcftools query -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\tGRCh38_only\n' {input.fn}
        ) | \
        awk 'BEGIN{{OFS="\t"; print "CHROM","POS","SVTYPE","SVLEN","CrossRef_support"}} {{print $1,$2,$3,$4,$5}}' \
        > {output.tsv}
        """