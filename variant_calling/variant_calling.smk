import os

# ==============================================================================
# PARSING CONFIG & METADATA
# ==============================================================================
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

# ==============================================================================
# INCLUDING RULE FILES
# ==============================================================================
include: "sample_variant_calling.smk"
include: "cohort_merge.smk"

# ==============================================================================
# RULE ALL
# ==============================================================================
rule all:
    input:
        config["output"] + "/cohort_results/CHM13_final_cohort_survivor.vcf.gz",
        expand(
            config["output"] + "/cohort_results/{cohort}_genotyped_matrix.vcf.gz",
            cohort=[
                "CHM13_final_cohort_survivor",
                "CHM13_tier1_supported",
                "CHM13_genotype_input"
            ]
        ),
        expand(
            config["output"] + "/{batch}/{sample}/04.force_calling/{sample}_{cohort}_genotypability.txt",
            batch=[b for b, s in VALID_PAIRS],
            sample=[s for b, s in VALID_PAIRS],
            cohort=[
                "CHM13_final_cohort_survivor",
                "CHM13_tier1_supported",
                "CHM13_genotype_input"
            ]
        )