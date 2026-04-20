import os

OUTDIR = config["output"]
NEEDLR_SCRIPT = config["needlr"]["cohort_script"]
NEEDLR_OUTDIR = os.path.join(OUTDIR, config["needlr"].get("outdir", "needLR_output"))

# ------------------------------------------------------------------------------
# Prepare a clean, sorted, indexed VCF for needLR
# ------------------------------------------------------------------------------
rule grch38_needlr_prepare_input:
    input:
        vcf=os.path.join(
            OUTDIR,
            "cohort_results",
            "GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz"
        )
    output:
        vcf=os.path.join(
            NEEDLR_OUTDIR,
            "GRCh38_final_cohort_survivor_genotyped_matrix.needlr_input.vcf.gz"
        ),
        tbi=os.path.join(
            NEEDLR_OUTDIR,
            "GRCh38_final_cohort_survivor_genotyped_matrix.needlr_input.vcf.gz.tbi"
        )
    conda:
        "envs/truvari_env.yml"
    threads: 4
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        r"""
        mkdir -p {NEEDLR_OUTDIR}

        bcftools view -e 'ALT="."' {input.vcf} -O u | \
        bcftools sort -O z -o {output.vcf}

        tabix -f -p vcf {output.vcf}
        """

# ------------------------------------------------------------------------------
# Run needLR on the GRCh38 full cohort genotyped matrix
# ------------------------------------------------------------------------------

rule grch38_needlr:
    input:
        vcf=os.path.join(
            NEEDLR_OUTDIR,
            "GRCh38_final_cohort_survivor_genotyped_matrix.needlr_input.vcf.gz"
        ),
        tbi=os.path.join(
            NEEDLR_OUTDIR,
            "GRCh38_final_cohort_survivor_genotyped_matrix.needlr_input.vcf.gz.tbi"
        )
    output:
        outdir=directory(os.path.join(
            NEEDLR_OUTDIR,
            "GRCh38_final_cohort_survivor_genotyped_matrix_needLR_cohort"
        ))
    conda:
        "envs/needlr_env.yml"
    threads: config["hc"]
    resources:
        mem_mb=config["hm"],
        time=config["ht"]
    shell:
        r"""
        rm -rf {output.outdir}
        mkdir -p {output.outdir}

        needLR annotate -Q {input.vcf} -O {output.outdir} --all
        """