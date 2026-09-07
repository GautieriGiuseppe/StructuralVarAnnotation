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

OUTDIR = config["output"].rstrip("/")

NEEDLR_OUTDIR = os.path.join(
    OUTDIR,
    config.get("needlr", {}).get("outdir", "needLR_output")
)

include: "align.smk"
include: "alignqc.smk"
include: "variant_calling/sample_variant_calling.smk"
include: "variant_calling/cohort_merge_grch38.smk"
include: "variant_calling/force_genotype_grch38.smk"
include: "variant_calling/haplotype_grch38.smk"
include: "variant_calling/needLR_grch38.smk"
include: "variant_calling/needLR_trio_grch38.smk"
include: "variant_calling/crossref_confirmation_grch38.smk"
include: "variant_calling/qc_report_generator.smk"


rule all:
    input:
        rules.all_alignqc.input,

        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz",
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz.tbi",
        f"{OUTDIR}/cohort_results/GRCh38_cohort_support_table.tsv",

        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz",
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz.tbi",

        os.path.join(
            NEEDLR_OUTDIR,
            "GRCh38_final_cohort_survivor_genotyped_matrix_needLR_cohort"
        ),

        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_confirmation.tsv",
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_confirmation_summary.json",

        f"{OUTDIR}/cohort_results/read_qc/multiqc/multiqc_report.html",

        f"{OUTDIR}/cohort_results/qc_report/GRCh38_full_pipeline_QC_report.html",
        f"{OUTDIR}/cohort_results/qc_report/GRCh38_full_pipeline_QC_summary.tsv",


rule all_trio:
    input:
        rules.all_alignqc.input,

        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz",
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz.tbi",
        f"{OUTDIR}/cohort_results/GRCh38_cohort_support_table.tsv",

        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz",
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz.tbi",

        rules.all_needlr_trio_grch38.input,

        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_confirmation.tsv",
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_confirmation_summary.json",

        f"{OUTDIR}/cohort_results/read_qc/multiqc/multiqc_report.html",

        rules.build_full_grch38_trio_qc_report.output[0],
        rules.build_full_grch38_trio_qc_report.output[1]