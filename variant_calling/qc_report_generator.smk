import os

OUTDIR = config["output"]

QC_REPORT_DIR = f"{OUTDIR}/cohort_results/qc_report"

READ_QC_DIR = f"{OUTDIR}/cohort_results/read_qc"

SAMPLE_FASTQ = {}

with open(config["samples"], "r") as f:
    header = f.readline().rstrip("\n").split("\t")
    sample_idx = header.index("sample_id")
    fastq_idx = header.index("fastq")

    for line in f:
        if line.strip():
            fields = line.rstrip("\n").split("\t")
            SAMPLE_FASTQ[fields[sample_idx]] = fields[fastq_idx]


def get_sample_fastq(wc):
    return SAMPLE_FASTQ[wc.sample]

SAMPLE_FASTQ = {}
SAMPLE_SUMMARY = {}

with open(config["samples"], "r") as f:
    header = f.readline().rstrip("\n").split("\t")
    sample_idx = header.index("sample_id")
    fastq_idx = header.index("fastq")
    summary_idx = header.index("summary")

    for line in f:
        if line.strip():
            fields = line.rstrip("\n").split("\t")
            SAMPLE_FASTQ[fields[sample_idx]] = fields[fastq_idx]
            SAMPLE_SUMMARY[fields[sample_idx]] = fields[summary_idx]


def get_sample_fastq(wc):
    return SAMPLE_FASTQ[wc.sample]


def get_sample_summary(wc):
    return SAMPLE_SUMMARY[wc.sample]

READ_QC_FASTQC_HTMLS = [
    f"{READ_QC_DIR}/fastqc/{batch}/{sample}/{sample}_fastqc.html"
    for batch, sample in VALID_PAIRS
]

READ_QC_NANOPLOT_HTMLS = [
    f"{READ_QC_DIR}/nanoplot/{batch}/{sample}/NanoPlot-report.html"
    for batch, sample in VALID_PAIRS
]

READ_QC_NANOPLOT_STATS = [
    f"{READ_QC_DIR}/nanoplot/{batch}/{sample}/NanoStats.txt"
    for batch, sample in VALID_PAIRS
]

rule fastqc_reads:
    input:
        fastq=get_sample_fastq
    output:
        html=f"{READ_QC_DIR}/fastqc/{{batch}}/{{sample}}/{{sample}}_fastqc.html",
        zip=f"{READ_QC_DIR}/fastqc/{{batch}}/{{sample}}/{{sample}}_fastqc.zip"
    conda:
        "envs/read_qc.yml"
    threads: 8
    resources:
        mem_mb=config["hm"],
        time=config["ht"]
    params:
        outdir=f"{READ_QC_DIR}/fastqc/{{batch}}/{{sample}}"
    shell:
        r"""
        mkdir -p {params.outdir}

        fastqc \
            --threads {threads} \
            --outdir {params.outdir} \
            {input.fastq}

        html=$(find {params.outdir} -maxdepth 1 -name "*_fastqc.html" | head -1)
        zipf=$(find {params.outdir} -maxdepth 1 -name "*_fastqc.zip" | head -1)

        if [ -z "$html" ]; then
            echo "ERROR: FastQC HTML output not found in {params.outdir}" >&2
            exit 1
        fi

        if [ -z "$zipf" ]; then
            echo "ERROR: FastQC ZIP output not found in {params.outdir}" >&2
            exit 1
        fi

        if [ "$html" != "{output.html}" ]; then
            cp "$html" {output.html}
        fi

        if [ "$zipf" != "{output.zip}" ]; then
            cp "$zipf" {output.zip}
        fi

        test -s {output.html}
        test -s {output.zip}
        """


rule nanoplot_reads:
    input:
        summary=get_sample_summary
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
        outdir=f"{READ_QC_DIR}/nanoplot/{{batch}}/{{sample}}"
    shell:
        r"""
        mkdir -p {params.outdir}

        NanoPlot \
            --summary {input.summary} \
            --threads {threads} \
            --outdir {params.outdir} \
            --prefix "" \
            --N50

        test -s {output.html}
        test -s {output.stats}
        """


rule multiqc_read_qc:
    input:
        fastqc=READ_QC_FASTQC_HTMLS,
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
        """

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
    #conda:
    #    "envs/plot_qc.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    params:
        outdir=f"{OUTDIR}/cohort_results/tool_reference_upset_new_cohort",
        prefix="GRCh38_integrated_toolref_upset_mean_sample_frequency",
        top_n=40,
        cmap="viridis"
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
            --sample-support-max 8
        """


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
    #conda:
    #    "envs/plot_qc.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    params:
        outdir=f"{OUTDIR}/cohort_results/needlr_annotation_plots",
        prefix="needlr_annotation_burden_and_support"
    shell:
        r"""
        mkdir -p {params.outdir}

        python plots/plot_NeedLR.py \
            --vcf {input.vcf} \
            --out-dir {params.outdir} \
            --out-prefix {params.prefix}
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
    #conda:
    #    "envs/plot_qc.yml"
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
            --out-prefix {params.prefix}
        """

rule plot_needlr_carrier_1_8_qc:
    input:
        vcf=os.path.join(
            NEEDLR_OUTDIR,
            "GRCh38_final_cohort_survivor_genotyped_matrix_needLR_cohort",
            "GRCh38_final_cohort_survivor_genotyped_matrix.needlr_input_needLR_1kg_v4.0",
            "GRCh38_final_cohort_survivor_genotyped_matrix.needlr_input.needLR.4.0.vcf.gz"
        )
    output:
        present_png=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_8/needlr_popfreq_violin_carrier_counts_1_8_present_only.png",
        present_pdf=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_8/needlr_popfreq_violin_carrier_counts_1_8_present_only.pdf",
        absent_png=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_8/needlr_popfreq_violin_carrier_counts_1_8_including_absent_variants.png",
        absent_pdf=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_8/needlr_popfreq_violin_carrier_counts_1_8_including_absent_variants.pdf",
        summary_png=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_8/needlr_popfreq_carrier_counts_1_8_summary_lines.png",
        summary_pdf=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_8/needlr_popfreq_carrier_counts_1_8_summary_lines.pdf",
        variant_tsv=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_8/needlr_variants_carrier_counts_1_8.tsv",
        long_tsv=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_8/needlr_popfreq_long_carrier_counts_1_8.tsv",
        summary_tsv=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_8/needlr_popfreq_summary_carrier_counts_1_8.tsv"
    #conda:
    #    "envs/plot_qc.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    params:
        outdir=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_8",
        prefix="needlr_popfreq_violin_carrier_counts_1_8",
        summary_prefix="needlr_popfreq_carrier_counts_1_8"
    shell:
        r"""
        mkdir -p {params.outdir}

        python plots/plot_needLR_popfreq_carriers_1-8.py \
            --vcf {input.vcf} \
            --out-dir {params.outdir} \
            --out-prefix {params.prefix} \
            --summary-prefix {params.summary_prefix}
        """

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
    #conda:
    #    "envs/plot_qc.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    params:
        outdir=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield",
        workflow_outdir=OUTDIR,
        batch=lambda wc: VALID_PAIRS[0][0],
        samples=lambda wc: ",".join(sorted(set([s for b, s in VALID_PAIRS]))),
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
    #conda:
    #    "envs/plot_qc.yml"
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


rule build_full_grch38_qc_report:
    input:
        # Read QC
        read_qc_multiqc=f"{OUTDIR}/cohort_results/read_qc/multiqc/multiqc_report.html",

        # UpSet
        upset_png=f"{OUTDIR}/cohort_results/tool_reference_upset_new_cohort/GRCh38_integrated_toolref_upset_mean_sample_frequency.png",
        upset_pdf=f"{OUTDIR}/cohort_results/tool_reference_upset_new_cohort/GRCh38_integrated_toolref_upset_mean_sample_frequency.pdf",

        # needLR annotation burden
        needlr_png=f"{OUTDIR}/cohort_results/needlr_annotation_plots/needlr_annotation_burden_and_support.png",
        needlr_pdf=f"{OUTDIR}/cohort_results/needlr_annotation_plots/needlr_annotation_burden_and_support.pdf",

        # needLR population frequency
        needlr_popfreq_png=f"{OUTDIR}/cohort_results/needlr_annotation_plots/needlr_control_population_frequency_summary.png",
        needlr_popfreq_pdf=f"{OUTDIR}/cohort_results/needlr_annotation_plots/needlr_control_population_frequency_summary.pdf",

        # needLR carrier-count 1-8 plots
        carrier_1_8_present_png=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_8/needlr_popfreq_violin_carrier_counts_1_8_present_only.png",
        carrier_1_8_absent_png=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_8/needlr_popfreq_violin_carrier_counts_1_8_including_absent_variants.png",
        carrier_1_8_summary_png=f"{OUTDIR}/cohort_results/needlr_population_frequency_carriers_1_8/needlr_popfreq_carrier_counts_1_8_summary_lines.png",

        # GRCh38 confirmation plots
        crossref_fig1_png=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_fig1_summary.png",
        crossref_fig2_png=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_fig2_confirmation_patterns.png",
        crossref_fig3_png=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_fig3_chromosome_confirmation.png",
        crossref_table=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_infofield_confirmation_table.tsv",
        crossref_metrics=f"{OUTDIR}/cohort_results/crossref_confirmation_infofield/crossref_infofield_summary_metrics.tsv",

        # Native/lifted breakpoint and SVLEN distance
        svlen_breakpoint_png=f"{OUTDIR}/cohort_results/crossref_confirmation_from_integrated_GRCh38/native_vs_lifted_support_distance/native_grch38_vs_lifted_chm13_support_distance.png",
        svlen_breakpoint_pdf=f"{OUTDIR}/cohort_results/crossref_confirmation_from_integrated_GRCh38/native_vs_lifted_support_distance/native_grch38_vs_lifted_chm13_support_distance.pdf",
        svlen_breakpoint_summary=f"{OUTDIR}/cohort_results/crossref_confirmation_from_integrated_GRCh38/native_vs_lifted_support_distance/native_grch38_vs_lifted_chm13_support_distance_summary.tsv",

        # Main cohort outputs
        cohort_vcf=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz",
        support_table=f"{OUTDIR}/cohort_results/GRCh38_cohort_support_table.tsv",
        genotyped_vcf=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz",
        confirmation_tsv=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_confirmation.tsv",
        confirmation_summary=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_confirmation_summary.json"
    output:
        html=f"{QC_REPORT_DIR}/GRCh38_full_pipeline_QC_report.html",
        summary=f"{QC_REPORT_DIR}/GRCh38_full_pipeline_QC_summary.tsv"
    params:
        outdir=QC_REPORT_DIR,
        samples=config["samples"],
        config_file="config.yml"
    conda:
        "envs/crossmap_truvari.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    shell:
        r"""
        mkdir -p {params.outdir}

        python plots/build_full_qc_report.py \
            --samples {params.samples} \
            --cohort-vcf {input.cohort_vcf} \
            --support-table {input.support_table} \
            --genotyped-vcf {input.genotyped_vcf} \
            --confirmation-tsv {input.confirmation_tsv} \
            --confirmation-summary {input.confirmation_summary} \
            --cohort-results-dir {OUTDIR}/cohort_results \
            --out-html {output.html} \
            --out-summary {output.summary}
        """