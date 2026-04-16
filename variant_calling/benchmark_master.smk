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

# Paths
OUTDIR = config["output"]
NEEDLR_OUTDIR = os.path.join(OUTDIR, config["needlr"].get("outdir", "needLR_output"))

# Includes
include: "sample_variant_calling.smk"
include: "cohort_merge.smk"
include: "cohort_merge_grch38.smk"
include: "force_genotype_grch38.smk"
include: "needLR_grch38.smk"

# ==============================================================================
# CHM13 BENCHMARKING
# ==============================================================================
rule all:
    input:
        f"{OUTDIR}/cohort_results/CHM13_final_cohort_survivor.vcf.gz",
        expand(
            f"{OUTDIR}/cohort_results/{{cohort}}_genotyped_matrix.vcf.gz",
            cohort=[
                "CHM13_final_cohort_survivor",
                "CHM13_tier1_supported",
                "CHM13_genotype_input"
            ]
        ),
        expand(
            f"{OUTDIR}/{{batch}}/{{sample}}/04.force_calling/{{sample}}_{{cohort}}_genotypability.txt",
            batch=[b for b, s in VALID_PAIRS],
            sample=[s for b, s in VALID_PAIRS],
            cohort=[
                "CHM13_final_cohort_survivor",
                "CHM13_tier1_supported",
                "CHM13_genotype_input"
            ]
        )

# ==============================================================================
# GRCh38 PREP
# ==============================================================================
rule all_grch38:
    input:
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz",
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz.tbi",
        f"{OUTDIR}/cohort_results/GRCh38_cohort_support_table.tsv",
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz",
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz.tbi"

# ==============================================================================
# GRCh38 ANNOTATION (needLR)
# ==============================================================================
rule all_grch38_annotation:
    input:
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz",
        os.path.join(
            NEEDLR_OUTDIR,
            "GRCh38_final_cohort_survivor_genotyped_matrix_needLR_v3.5_cohort",
            "GRCh38_final_cohort_survivor_genotyped_matrix_needLR_v3.5_cohort_RESULTS.txt"
        )