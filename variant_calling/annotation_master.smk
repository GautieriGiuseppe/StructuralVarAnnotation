import os

samples_file = config["samples"]
VALID_PAIRS = []

with open(samples_file, "r") as f:
    header = f.readline().strip().split("\t")
    sample_idx = header.index("sample_id")
    batch_idx = header.index("batch_id")

    for line in f:
        if line.strip():
            fields = line.strip().split("\t")
            VALID_PAIRS.append((fields[batch_idx], fields[sample_idx]))

OUTDIR = config["output"]
NEEDLR_OUTDIR = os.path.join(OUTDIR, config["needlr"].get("outdir", "needLR_output"))

include: "sample_variant_calling.smk"
include: "cohort_merge_grch38.smk"
include: "force_genotype_grch38.smk"
include: "needLR_grch38.smk"
include: "crossref_confirmation_grch38.smk"

rule all:
    input:
        # GRCh38 cohort
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz",
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz.tbi",
        f"{OUTDIR}/cohort_results/GRCh38_cohort_support_table.tsv",
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz",
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz.tbi",

        # needLR
        os.path.join(
            NEEDLR_OUTDIR,
            "GRCh38_final_cohort_survivor_genotyped_matrix_needLR_cohort"
        ),

        # cross-reference confirmation
        f"{OUTDIR}/cohort_results/CHM13_final_cohort_survivor_to_GRCh38.vcf.gz",
        f"{OUTDIR}/cohort_results/CHM13_final_cohort_survivor_to_GRCh38.vcf.gz.tbi",
        f"{OUTDIR}/cohort_results/grch38_vs_chm13lifted_truvari/summary.json",
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_confirmation.tsv"