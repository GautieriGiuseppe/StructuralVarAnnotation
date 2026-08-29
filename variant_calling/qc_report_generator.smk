import os

OUTDIR = config["output"].rstrip("/")

QC_REPORT_DIR = f"{OUTDIR}/cohort_results/qc_report"
READ_QC_DIR = f"{OUTDIR}/cohort_results/read_qc"

N_SAMPLES = len(VALID_PAIRS)

SAMPLE_IDS = sorted(set([sample for batch, sample in VALID_PAIRS]))
BATCH_IDS = sorted(set([batch for batch, sample in VALID_PAIRS]))

DEFAULT_BATCH_FOR_QC = BATCH_IDS[0]
SAMPLES_CSV_FOR_QC = ",".join(SAMPLE_IDS)

NEEDLR_OUTDIR = os.path.join(
    OUTDIR,
    config.get("needlr", {}).get("outdir", "needLR_output")
)


# ============================================================
# Trio metadata, parsed safely
# ============================================================

TRIO_FILE = config.get("needlr", {}).get("trio_file", None)
FAMILY_IDS = []

if TRIO_FILE is not None and os.path.exists(TRIO_FILE):
    with open(TRIO_FILE, "r") as f:
        header = f.readline().rstrip("\n").split("\t")

        if "family_id" not in header:
            raise ValueError(
                f"Trio file {TRIO_FILE} is missing required column: family_id"
            )

        fam_idx = header.index("family_id")

        for line in f:
            if line.strip():
                fields = line.rstrip("\n").split("\t")
                FAMILY_IDS.append(fields[fam_idx])


# ============================================================
# Read QC: NanoPlot + MultiQC only
#
# Supports:
#   1. sequencing_summary.tsv / sequencing_summary.tsv.gz
#   2. single FASTQ
#   3. comma-separated FASTQs
#
# samples.tsv examples:
#
# sample_id  ubam  summary                      fastq                         batch_id
# CHILD      NA    child_summary.tsv.gz          child.fastq.gz                 FAM001
# CHILD      NA    NA                            child_1.fq.gz,child_2.fq.gz    FAM001
# ============================================================

SAMPLE_FASTQ = {}
SAMPLE_SUMMARY = {}


def is_missing_path(x):
    return x is None or str(x).strip() in {"", "NA", "NaN", "nan", "None", "."}


with open(config["samples"], "r") as f:
    header = f.readline().rstrip("\n").split("\t")

    sample_idx = header.index("sample_id")

    fastq_idx = header.index("fastq") if "fastq" in header else None
    summary_idx = header.index("summary") if "summary" in header else None

    for line in f:
        if line.strip():
            fields = line.rstrip("\n").split("\t")
            sample = fields[sample_idx]

            fastq = fields[fastq_idx] if fastq_idx is not None else ""
            summary = fields[summary_idx] if summary_idx is not None else ""

            SAMPLE_FASTQ[sample] = fastq
            SAMPLE_SUMMARY[sample] = summary


def get_sample_read_qc_input(wc):
    """
    Prefer ONT sequencing summary if available.
    Fall back to FASTQ if summary is missing.
    FASTQ can be a single path or comma-separated paths.
    """

    summary = SAMPLE_SUMMARY.get(wc.sample, "")
    fastq = SAMPLE_FASTQ.get(wc.sample, "")

    if not is_missing_path(summary):
        return [x.strip() for x in summary.split(",") if x.strip()]

    if not is_missing_path(fastq):
        return [x.strip() for x in fastq.split(",") if x.strip()]

    raise ValueError(
        f"No summary or FASTQ file provided for sample {wc.sample}. "
        "NanoPlot requires either summary or fastq in samples.tsv."
    )


def get_sample_read_qc_mode(wc):
    """
    Return 'summary' if summary is available, otherwise 'fastq'.
    """

    summary = SAMPLE_SUMMARY.get(wc.sample, "")

    if not is_missing_path(summary):
        return "summary"

    return "fastq"


READ_QC_NANOPLOT_HTMLS = [
    f"{READ_QC_DIR}/nanoplot/{batch}/{sample}/NanoPlot-report.html"
    for batch, sample in VALID_PAIRS
]

READ_QC_NANOPLOT_STATS = [
    f"{READ_QC_DIR}/nanoplot/{batch}/{sample}/NanoStats.txt"
    for batch, sample in VALID_PAIRS
]


rule nanoplot_reads:
    input:
        reads=get_sample_read_qc_input
    output:
        html=f"{READ_QC_DIR}/nanoplot/{{batch}}/{{sample}}/NanoPlot-report.html",
        stats=f"{READ_QC_DIR}/nanoplot/{{batch}}/{{sample}}/NanoStats.txt"
    conda:
        "envs/read_qc.yml"
    threads: 4
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    params:
        outdir=f"{READ_QC_DIR}/nanoplot/{{batch}}/{{sample}}",
        mode=get_sample_read_qc_mode
    shell:
        r"""
        mkdir -p {params.outdir}

        if [[ "{params.mode}" == "summary" ]]; then
            tmp_summary="{params.outdir}/{wildcards.sample}.sequencing_summary.tsv"
            rm -f "$tmp_summary"

            python - <<'PY'
from pathlib import Path
import gzip

summary_files = "{input.reads}".split()
out = Path("{params.outdir}") / "{wildcards.sample}.sequencing_summary.tsv"

header = None
expected_ncols = None
n_written = 0

def open_maybe_gzip(file_path):
    if str(file_path).endswith(".gz"):
        return gzip.open(file_path, "rt")
    return open(file_path, "r")

with out.open("w") as fout:
    for file_path in summary_files:
        file_path = file_path.strip()

        if not file_path:
            continue

        with open_maybe_gzip(file_path) as fin:
            local_header = fin.readline().rstrip("\n")

            if not local_header:
                continue

            local_cols = local_header.split("\t")

            if header is None:
                header = local_header
                expected_ncols = len(local_cols)
                fout.write(header + "\n")
            else:
                if local_header != header:
                    raise ValueError(
                        "Header mismatch in %s. Expected the same sequencing_summary columns."
                        % file_path
                    )

            for line in fin:
                line = line.rstrip("\n")

                if not line:
                    continue

                if line == header:
                    continue

                fields = line.split("\t")

                if len(fields) != expected_ncols:
                    continue

                fout.write(line + "\n")
                n_written += 1

if header is None or n_written == 0:
    raise ValueError(
        "No valid sequencing_summary rows written to %s" % out
    )

print("Wrote merged sequencing summary: %s" % out)
print("Valid rows: %s" % n_written)
PY

            NanoPlot \
                --summary "$tmp_summary" \
                --threads {threads} \
                --outdir {params.outdir} \
                --prefix "" \
                --N50

        else
            NanoPlot \
                --fastq {input.reads} \
                --threads {threads} \
                --outdir {params.outdir} \
                --prefix "" \
                --N50
        fi

        test -s {output.html}
        test -s {output.stats}
        """


rule multiqc_read_qc:
    input:
        nanoplot=READ_QC_NANOPLOT_HTMLS,
        nanostats=READ_QC_NANOPLOT_STATS
    output:
        html=f"{READ_QC_DIR}/multiqc/multiqc_report.html"
    conda:
        "envs/read_qc.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    params:
        outdir=f"{READ_QC_DIR}/multiqc",
        search_dir=READ_QC_DIR
    shell:
        r"""
        mkdir -p {params.outdir}

        multiqc \
            {params.search_dir} \
            --outdir {params.outdir} \
            --filename multiqc_report.html \
            --force

        test -s {output.html}
        """


# ============================================================
# Tool/reference support UpSet plot
# ============================================================

rule plot_upset_qc:
    input:
        vcf=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz",
        tbi=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz.tbi"
    output:
        png=f"{OUTDIR}/cohort_results/tool_reference_upset_new_cohort/GRCh38_integrated_toolref_upset_mean_sample_frequency.png",
        pdf=f"{OUTDIR}/cohort_results/tool_reference_upset_new_cohort/GRCh38_integrated_toolref_upset_mean_sample_frequency.pdf",
        variant_table=f"{OUTDIR}/cohort_results/tool_reference_upset_new_cohort/toolref_variant_support_table.tsv",
        intersections=f"{OUTDIR}/cohort_results/tool_reference_upset_new_cohort/toolref_upset_intersections.tsv",
        set_sizes=f"{OUTDIR}/cohort_results/tool_reference_upset_new_cohort/toolref_upset_set_sizes.tsv"
    conda:
        "envs/plot_qc.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    params:
        outdir=f"{OUTDIR}/cohort_results/tool_reference_upset_new_cohort",
        prefix="GRCh38_integrated_toolref_upset_mean_sample_frequency",
        top_n=40,
        cmap="viridis",
        sample_support_max=N_SAMPLES
    shell:
        r"""
        mkdir -p {params.outdir}

        python plots/plot_upset.py \
            --vcf {input.vcf} \
            --out-dir {params.outdir} \
            --out-prefix {params.prefix} \
            --top-n {params.top_n} \
            --cmap {params.cmap} \
            --sample-support-min 1 \
            --sample-support-max {params.sample_support_max}
        """


# ============================================================
# needLR cohort annotation plots
#
# These are standard cohort-needLR plots.
# They are intentionally NOT required by the trio report.
# ============================================================

rule plot_needlr_qc:
    input:
        vcf=os.path.join(
            NEEDLR_OUTDIR,
            "GRCh38_final_cohort_survivor_genotyped_matrix_needLR_cohort",
            "GRCh38_final_cohort_survivor_genotyped_matrix.needlr_input_needLR_1kg_v4.0",
            "GRCh38_final_cohort_survivor_genotyped_matrix.needlr_input.needLR.4.0.vcf.gz"
        )
    output:
        png=f"{OUTDIR}/cohort_results/needlr_annotation_plots/needlr_annotation_burden_and_support.png",
        pdf=f"{OUTDIR}/cohort_results/needlr_annotation_plots/needlr_annotation_burden_and_support.pdf",
        variant_table=f"{OUTDIR}/cohort_results/needlr_annotation_plots/needlr_variant_annotation_table.tsv",
        burden=f"{OUTDIR}/cohort_results/needlr_annotation_plots/needlr_annotation_burden_summary.tsv",
        known_novel=f"{OUTDIR}/cohort_results/needlr_annotation_plots/needlr_known_vs_novel_annotation_burden.tsv",
        context=f"{OUTDIR}/cohort_results/needlr_annotation_plots/needlr_genomic_context_burden.tsv",
        carrier=f"{OUTDIR}/cohort_results/needlr_annotation_plots/needlr_carrier_count_distribution.tsv"
    conda:
        "envs/plot_qc.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    params:
        outdir=f"{OUTDIR}/cohort_results/needlr_annotation_plots",
        prefix="needlr_annotation_burden_and_support",
        max_carriers=N_SAMPLES
    shell:
        r"""
        mkdir -p {params.outdir}

        python plots/plot_NeedLR.py \
            --vcf {input.vcf} \
            --out-dir {params.outdir} \
            --out-prefix {params.prefix} \
            --title "needLR annotation burden and support" \
            --max-carriers {params.max_carriers}
        """


rule plot_allele_pop_frequency_qc:
    input:
        vcf=os.path.join(
            NEEDLR_OUTDIR,
            "GRCh38_final_cohort_survivor_genotyped_matrix_needLR_cohort",
            "GRCh38_final_cohort_survivor_genotyped_matrix.needlr_input_needLR_1kg_v4.0",
            "GRCh38_final_cohort_survivor_genotyped_matrix.needlr_input.needLR.4.0.vcf.gz"
        )
    output:
        png=f"{OUTDIR}/cohort_results/needlr_annotation_plots/needlr_control_population_frequency_summary.png",
        pdf=f"{OUTDIR}/cohort_results/needlr_annotation_plots/needlr_control_population_frequency_summary.pdf",
        variant_table=f"{OUTDIR}/cohort_results/needlr_annotation_plots/needlr_control_population_frequency_variant_table.tsv",
        mean_af=f"{OUTDIR}/cohort_results/needlr_annotation_plots/needlr_mean_control_af_among_present.tsv",
        presence=f"{OUTDIR}/cohort_results/needlr_annotation_plots/needlr_control_population_presence.tsv",
        summary=f"{OUTDIR}/cohort_results/needlr_annotation_plots/needlr_control_population_frequency_summary.tsv"
    conda:
        "envs/plot_qc.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    params:
        outdir=f"{OUTDIR}/cohort_results/needlr_annotation_plots",
        prefix="needlr_control_population_frequency_summary"
    shell:
        r"""
        mkdir -p {params.outdir}

        python plots/plot_allele_pop_frequency.py \
            --vcf {input.vcf} \
            --out-dir {params.outdir} \
            --out-prefix {params.prefix} \
            --title "needLR control population frequency summary"
        """


rule plot_needlr_carrier_dynamic_qc:
    input:
        vcf=os.path.join(
            NEEDLR_OUTDIR,
            "GRCh38_final_cohort_survivor_genotyped_matrix_needLR_cohort",
            "GRCh38_final_cohort_survivor_genotyped_matrix.needlr_input_needLR_1kg_v4.0",
            "GRCh38_final_cohort_survivor_genotyped_matrix.needlr_input.needLR.4.0.vcf.gz"
        )
    output:
        present_png=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_{N_SAMPLES}/needlr_popfreq_violin_carrier_counts_1_{N_SAMPLES}_present_only.png",
        present_pdf=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_{N_SAMPLES}/needlr_popfreq_violin_carrier_counts_1_{N_SAMPLES}_present_only.pdf",
        absent_png=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_{N_SAMPLES}/needlr_popfreq_violin_carrier_counts_1_{N_SAMPLES}_including_absent_variants.png",
        absent_pdf=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_{N_SAMPLES}/needlr_popfreq_violin_carrier_counts_1_{N_SAMPLES}_including_absent_variants.pdf",
        summary_png=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_{N_SAMPLES}/needlr_popfreq_carrier_counts_1_{N_SAMPLES}_summary_lines.png",
        summary_pdf=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_{N_SAMPLES}/needlr_popfreq_carrier_counts_1_{N_SAMPLES}_summary_lines.pdf",
        variant_tsv=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_{N_SAMPLES}/needlr_variants_carrier_counts_1_{N_SAMPLES}.tsv",
        long_tsv=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_{N_SAMPLES}/needlr_popfreq_long_carrier_counts_1_{N_SAMPLES}.tsv",
        summary_tsv=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_{N_SAMPLES}/needlr_popfreq_summary_carrier_counts_1_{N_SAMPLES}.tsv"
    conda:
        "envs/plot_qc.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    params:
        outdir=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_{N_SAMPLES}",
        out_prefix="needlr_popfreq_violin_carrier_counts",
        summary_prefix="needlr_popfreq_carrier_counts",
        max_carriers=N_SAMPLES
    shell:
        r"""
        mkdir -p {params.outdir}

        python plots/plot_needLR_popfreq_carriers.py \
            --vcf {input.vcf} \
            --out-dir {params.outdir} \
            --out-prefix {params.out_prefix} \
            --summary-prefix {params.summary_prefix} \
            --title-prefix "needLR control population frequency" \
            --min-carriers 1 \
            --max-carriers {params.max_carriers}
        """


# ============================================================
# GRCh38/CHM13 cross-reference confirmation plots
# ============================================================

rule plot_grch38_confirmation_qc:
    input:
        vcf=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz",
        tbi=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz.tbi"
    output:
        fig1_png=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_fig1_summary.png",
        fig1_pdf=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_fig1_summary.pdf",
        fig2_png=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_fig2_confirmation_patterns.png",
        fig2_pdf=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_fig2_confirmation_patterns.pdf",
        fig3_png=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_fig3_chromosome_confirmation.png",
        fig3_pdf=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_fig3_chromosome_confirmation.pdf",
        table=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_infofield_confirmation_table.tsv",
        metrics=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_infofield_summary_metrics.tsv"
    conda:
        "envs/plot_qc.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    params:
        outdir=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield",
        workflow_outdir=OUTDIR,
        batch=DEFAULT_BATCH_FOR_QC,
        samples=SAMPLES_CSV_FOR_QC,
        title_prefix="Cross-reference confirmation"
    shell:
        r"""
        mkdir -p {params.outdir}

        python plots/plot_grch38_confirmation.py \
            --integrated-vcf {input.vcf} \
            --out-dir {params.outdir} \
            --workflow-outdir {params.workflow_outdir} \
            --batch {params.batch} \
            --samples {params.samples} \
            --title-prefix "{params.title_prefix}"
        """


rule plot_svlen_breakpoint_qc:
    input:
        vcf=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz",
        tbi=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz.tbi",
        metadata=f"{OUTDIR}/cohort_results/GRCh38_toolref_48way_vcf_metadata.tsv"
    output:
        png=f"{OUTDIR}/cohort_results/crossref_confirmation_from_integrated_GRCh38/native_vs_lifted_support_distance/native_grch38_vs_lifted_chm13_support_distance.png",
        pdf=f"{OUTDIR}/cohort_results/crossref_confirmation_from_integrated_GRCh38/native_vs_lifted_support_distance/native_grch38_vs_lifted_chm13_support_distance.pdf",
        table=f"{OUTDIR}/cohort_results/crossref_confirmation_from_integrated_GRCh38/native_vs_lifted_support_distance/native_grch38_vs_lifted_chm13_support_distance_table.tsv",
        summary=f"{OUTDIR}/cohort_results/crossref_confirmation_from_integrated_GRCh38/native_vs_lifted_support_distance/native_grch38_vs_lifted_chm13_support_distance_summary.tsv"
    conda:
        "envs/plot_qc.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    params:
        outdir=f"{OUTDIR}/cohort_results/crossref_confirmation_from_integrated_GRCh38/native_vs_lifted_support_distance",
        prefix="native_grch38_vs_lifted_chm13_support_distance"
    shell:
        r"""
        mkdir -p {params.outdir}

        python plots/plot_svlen_breakpoint.py \
            --vcf {input.vcf} \
            --metadata {input.metadata} \
            --out-dir {params.outdir} \
            --out-prefix {params.prefix}
        """


# ============================================================
# Optional liftover-cost plots
#
# This rule is available, but it is not required by the trio
# report unless you explicitly add its outputs to the report input.
# ============================================================

rule plot_liftover_cost_qc:
    input:
        cohort_vcf=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz",
        metadata=f"{OUTDIR}/cohort_results/GRCh38_toolref_48way_vcf_metadata.tsv"
    output:
        counts=f"{OUTDIR}/cohort_results/liftover_cost_counts_chm13_to_grch38.tsv",
        sizes=f"{OUTDIR}/cohort_results/liftover_cost_variant_sizes_chm13_to_grch38.tsv",
        retention=f"{OUTDIR}/cohort_results/liftover_retention_chm13_to_grch38.tsv",
        retention_by_svtype=f"{OUTDIR}/cohort_results/liftover_retention_by_svtype_chm13_to_grch38.tsv",
        fig1_png=f"{OUTDIR}/cohort_results/figures_liftover_cost_chm13_to_grch38/fig1_liftover_cohort_total_counts_chm13_to_grch38.png",
        fig1_pdf=f"{OUTDIR}/cohort_results/figures_liftover_cost_chm13_to_grch38/fig1_liftover_cohort_total_counts_chm13_to_grch38.pdf",
        fig2_png=f"{OUTDIR}/cohort_results/figures_liftover_cost_chm13_to_grch38/fig2_liftover_cohort_svtype_composition_chm13_to_grch38.png",
        fig2_pdf=f"{OUTDIR}/cohort_results/figures_liftover_cost_chm13_to_grch38/fig2_liftover_cohort_svtype_composition_chm13_to_grch38.pdf",
        fig3_png=f"{OUTDIR}/cohort_results/figures_liftover_cost_chm13_to_grch38/fig3_liftover_retention_by_svtype_chm13_to_grch38.png",
        fig3_pdf=f"{OUTDIR}/cohort_results/figures_liftover_cost_chm13_to_grch38/fig3_liftover_retention_by_svtype_chm13_to_grch38.pdf",
        fig4_png=f"{OUTDIR}/cohort_results/figures_liftover_cost_chm13_to_grch38/fig4_liftover_raw_counts_per_tool_chm13_to_grch38.png",
        fig4_pdf=f"{OUTDIR}/cohort_results/figures_liftover_cost_chm13_to_grch38/fig4_liftover_raw_counts_per_tool_chm13_to_grch38.pdf",
        fig5_png=f"{OUTDIR}/cohort_results/figures_liftover_cost_chm13_to_grch38/fig5_liftover_svtype_breakdown_per_tool_chm13_to_grch38.png",
        fig5_pdf=f"{OUTDIR}/cohort_results/figures_liftover_cost_chm13_to_grch38/fig5_liftover_svtype_breakdown_per_tool_chm13_to_grch38.pdf",
        fig6_png=f"{OUTDIR}/cohort_results/figures_liftover_cost_chm13_to_grch38/fig6_liftover_size_distribution_chm13_to_grch38.png",
        fig6_pdf=f"{OUTDIR}/cohort_results/figures_liftover_cost_chm13_to_grch38/fig6_liftover_size_distribution_chm13_to_grch38.pdf"
    conda:
        "envs/plot_qc.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    params:
        samples=config["samples"],
        outdir=OUTDIR,
        figure_outdir=f"{OUTDIR}/cohort_results/figures_liftover_cost_chm13_to_grch38",
        title_suffix=f"{N_SAMPLES}-sample cohort"
    shell:
        r"""
        mkdir -p {params.figure_outdir}

        python plots/plot_liftover_cost_chm13_to_grch38.py \
            --samples {params.samples} \
            --outdir {params.outdir} \
            --out-dir {params.figure_outdir} \
            --title-suffix "{params.title_suffix}"
        """


# ============================================================
# Full cohort QC report
#
# This report includes standard cohort needLR plots.
# Use this for normal cohort mode, not trio-comparator mode.
# ============================================================

rule build_full_grch38_qc_report:
    input:
        cohort_vcf=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz",
        support_table=f"{OUTDIR}/cohort_results/GRCh38_cohort_support_table.tsv",
        genotyped_vcf=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz",
        confirmation_tsv=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_confirmation.tsv",
        confirmation_summary=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_confirmation_summary.json",

        read_qc_multiqc=f"{OUTDIR}/cohort_results/read_qc/multiqc/multiqc_report.html",

        upset_png=f"{OUTDIR}/cohort_results/tool_reference_upset_new_cohort/GRCh38_integrated_toolref_upset_mean_sample_frequency.png",
        upset_pdf=f"{OUTDIR}/cohort_results/tool_reference_upset_new_cohort/GRCh38_integrated_toolref_upset_mean_sample_frequency.pdf",

        crossref_fig1_png=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_fig1_summary.png",
        crossref_fig2_png=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_fig2_confirmation_patterns.png",
        crossref_fig3_png=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_fig3_chromosome_confirmation.png",
        crossref_table=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_infofield_confirmation_table.tsv",
        crossref_metrics=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_infofield_summary_metrics.tsv",

        svlen_breakpoint_png=f"{OUTDIR}/cohort_results/crossref_confirmation_from_integrated_GRCh38/native_vs_lifted_support_distance/native_grch38_vs_lifted_chm13_support_distance.png",
        svlen_breakpoint_pdf=f"{OUTDIR}/cohort_results/crossref_confirmation_from_integrated_GRCh38/native_vs_lifted_support_distance/native_grch38_vs_lifted_chm13_support_distance.pdf",
        svlen_breakpoint_summary=f"{OUTDIR}/cohort_results/crossref_confirmation_from_integrated_GRCh38/native_vs_lifted_support_distance/native_grch38_vs_lifted_chm13_support_distance_summary.tsv",

        needlr_png=f"{OUTDIR}/cohort_results/needlr_annotation_plots/needlr_annotation_burden_and_support.png",
        needlr_popfreq_png=f"{OUTDIR}/cohort_results/needlr_annotation_plots/needlr_control_population_frequency_summary.png",
        needlr_carrier_png=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_{N_SAMPLES}/needlr_popfreq_violin_carrier_counts_1_{N_SAMPLES}_present_only.png"
    output:
        html=f"{QC_REPORT_DIR}/GRCh38_full_pipeline_QC_report.html",
        summary=f"{QC_REPORT_DIR}/GRCh38_full_pipeline_QC_summary.tsv"
    conda:
        "envs/plot_qc.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    params:
        cohort_results_dir=f"{OUTDIR}/cohort_results",
        samples=config["samples"]
    shell:
        r"""
        mkdir -p $(dirname {output.html})

        python plots/build_full_qc_report.py \
            --mode cohort \
            --samples {params.samples} \
            --cohort-vcf {input.cohort_vcf} \
            --support-table {input.support_table} \
            --genotyped-vcf {input.genotyped_vcf} \
            --confirmation-tsv {input.confirmation_tsv} \
            --confirmation-summary {input.confirmation_summary} \
            --cohort-results-dir {params.cohort_results_dir} \
            --out-html {output.html} \
            --out-summary {output.summary}
        """


# ============================================================
# Full trio QC report
#
# This report excludes standard cohort needLR plots and instead
# requires needLR comparator outputs from the trio branch.
# ============================================================

rule build_full_grch38_trio_qc_report:
    input:
        # Main cohort outputs
        cohort_vcf=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz",
        support_table=f"{OUTDIR}/cohort_results/GRCh38_cohort_support_table.tsv",
        genotyped_vcf=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz",
        confirmation_tsv=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_confirmation.tsv",
        confirmation_summary=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_confirmation_summary.json",

        # Read QC
        read_qc_multiqc=f"{OUTDIR}/cohort_results/read_qc/multiqc/multiqc_report.html",

        # UpSet
        upset_png=f"{OUTDIR}/cohort_results/tool_reference_upset_new_cohort/GRCh38_integrated_toolref_upset_mean_sample_frequency.png",
        upset_pdf=f"{OUTDIR}/cohort_results/tool_reference_upset_new_cohort/GRCh38_integrated_toolref_upset_mean_sample_frequency.pdf",

        # GRCh38/CHM13 confirmation plots
        crossref_fig1_png=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_fig1_summary.png",
        crossref_fig2_png=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_fig2_confirmation_patterns.png",
        crossref_fig3_png=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_fig3_chromosome_confirmation.png",
        crossref_table=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_infofield_confirmation_table.tsv",
        crossref_metrics=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_infofield_summary_metrics.tsv",

        # Native/lifted breakpoint and SVLEN distance
        svlen_breakpoint_png=f"{OUTDIR}/cohort_results/crossref_confirmation_from_integrated_GRCh38/native_vs_lifted_support_distance/native_grch38_vs_lifted_chm13_support_distance.png",
        svlen_breakpoint_pdf=f"{OUTDIR}/cohort_results/crossref_confirmation_from_integrated_GRCh38/native_vs_lifted_support_distance/native_grch38_vs_lifted_chm13_support_distance.pdf",
        svlen_breakpoint_summary=f"{OUTDIR}/cohort_results/crossref_confirmation_from_integrated_GRCh38/native_vs_lifted_support_distance/native_grch38_vs_lifted_chm13_support_distance_summary.tsv",

        # Trio needLR comparator
        trio_done=[
            os.path.join(
                NEEDLR_OUTDIR,
                "trio",
                family_id,
                f"{family_id}.needLR_comparator.done"
            )
            for family_id in FAMILY_IDS
        ]
    output:
        html=f"{QC_REPORT_DIR}/GRCh38_full_pipeline_QC_report_trio.html",
        summary=f"{QC_REPORT_DIR}/GRCh38_full_pipeline_QC_summary_trio.tsv"
    conda:
        "envs/plot_qc.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    params:
        cohort_results_dir=f"{OUTDIR}/cohort_results",
        trio_dir=os.path.join(NEEDLR_OUTDIR, "trio"),
        samples=config["samples"]
    shell:
        r"""
        mkdir -p $(dirname {output.html})

        python plots/build_full_qc_report.py \
            --mode trio \
            --samples {params.samples} \
            --cohort-vcf {input.cohort_vcf} \
            --support-table {input.support_table} \
            --genotyped-vcf {input.genotyped_vcf} \
            --confirmation-tsv {input.confirmation_tsv} \
            --confirmation-summary {input.confirmation_summary} \
            --cohort-results-dir {params.cohort_results_dir} \
            --needlr-trio-dir {params.trio_dir} \
            --out-html {output.html} \
            --out-summary {output.summary}
        """