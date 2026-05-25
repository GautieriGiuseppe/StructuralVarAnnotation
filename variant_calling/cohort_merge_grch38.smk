import os

OUTDIR = config["output"]

PAIRS = VALID_PAIRS
BATCHES = [b for b, s in PAIRS]
SAMPLES = [s for b, s in PAIRS]

# ==============================================================================
# CHAIN
# ==============================================================================

CHM13_TO_GRCH38_CHAIN = config["chain_file"]

# ==============================================================================
# TOOL / REFERENCE SET ORDERS
# ==============================================================================

# GRCh38-space 48-way merge:
# 0 = Sniffles native GRCh38
# 1 = Delly native GRCh38
# 2 = CuteSV native GRCh38
# 3 = Sniffles CHM13 lifted to GRCh38
# 4 = Delly CHM13 lifted to GRCh38
# 5 = CuteSV CHM13 lifted to GRCh38
TOOLREF_SET_ORDER_GRCH38 = [
    ("Sniffles_GR38", "sniffles", "native_grch38"),
    ("Delly_GR38", "delly", "native_grch38"),
    ("CuteSV_GR38", "cuteSV", "native_grch38"),
    ("Sniffles_CHM13_to_GR38", "sniffles", "lifted_chm13_to_grch38"),
    ("Delly_CHM13_to_GR38", "delly", "lifted_chm13_to_grch38"),
    ("CuteSV_CHM13_to_GR38", "cuteSV", "lifted_chm13_to_grch38"),
]


# ==============================================================================
# MAIN TARGETS PRODUCED BY THIS FILE
# ==============================================================================

rule all_grch38_cohort_merge:
    input:
        # GRCh38 tool/reference-supported cohort
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz",
        f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz.tbi",
        f"{OUTDIR}/cohort_results/GRCh38_cohort_support_table.tsv",


# ==============================================================================
# NATIVE GRCh38 PASS FILTERING
# ==============================================================================

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


# ==============================================================================
# PER-TOOL CHM13 -> GRCh38 LIFTOVER
# ==============================================================================

rule crossmap_tool_chm13_to_grch38:
    input:
        vcf=lambda wc: (
            f"{OUTDIR}/{wc.batch}/{wc.sample}/03.variant_calling/{wc.sample}_chm13_delly_symbolic.vcf"
            if wc.tool == "delly"
            else f"{OUTDIR}/{wc.batch}/{wc.sample}/03.variant_calling/{wc.sample}_chm13_{wc.tool}.vcf"
        ),
        chain=CHM13_TO_GRCH38_CHAIN,
        ref=config["reference_uncompressed"]["grch38"]
    output:
        raw=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/liftover/{{tool}}_{{sample}}_chm13-to-grch38.raw.vcf"
    conda:
        "envs/crossmap_truvari.yml"
    threads: 2
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        r"""
        mkdir -p $(dirname {output.raw})

        CrossMap vcf \
            {input.chain} \
            {input.vcf} \
            {input.ref} \
            {output.raw}
        """


rule clean_sort_tool_chm13_to_grch38:
    input:
        raw=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/liftover/{{tool}}_{{sample}}_chm13-to-grch38.raw.vcf"
    output:
        vcf=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/liftover/{{tool}}_{{sample}}_chm13-to-grch38.toolref.vcf"
    conda:
        "envs/crossmap_truvari.yml"
    threads: 2
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        r"""
        tmp=$(dirname {output.vcf})/{wildcards.tool}_{wildcards.sample}_chm13-to-grch38.toolref.cleaned.vcf

        awk 'BEGIN{{FS=OFS="\t"}}
             /^#/ {{print; next}}
             ($2 ~ /^[0-9]+$/) && ($2 > 0) {{print}}
        ' {input.raw} > $tmp

        bcftools sort $tmp -O v -o {output.vcf}

        rm -f $tmp
        """


rule pass_filter_tool_chm13_lifted_to_grch38:
    input:
        vcf=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/liftover/{{tool}}_{{sample}}_chm13-to-grch38.toolref.vcf"
    output:
        vcf=f"{OUTDIR}/{{batch}}/{{sample}}/03.variant_calling/liftover/{{tool}}_{{sample}}_chm13-to-grch38.toolref_pass.vcf"
    conda:
        "envs/crossmap_truvari.yml"
    threads: 2
    resources:
        mem_mb=config["mm"],
        time=config["ht"]
    shell:
        "bcftools view -f PASS {input.vcf} -O v -o {output.vcf}"


# ==============================================================================
# GRCh38 48-WAY INPUT LIST
# ==============================================================================

def grch38_toolref_vcf_path(batch, sample, tool, source):
    if source == "native_grch38":
        return f"{OUTDIR}/{batch}/{sample}/03.variant_calling/{sample}_grch38_{tool}_pass.vcf"

    if source == "lifted_chm13_to_grch38":
        return f"{OUTDIR}/{batch}/{sample}/03.variant_calling/liftover/{tool}_{sample}_chm13-to-grch38.toolref_pass.vcf"

    raise ValueError(f"Unknown GRCh38 source: {source}")


def grch38_toolref_48way_paths():
    paths = []

    for set_label, tool, source in TOOLREF_SET_ORDER_GRCH38:
        for batch, sample in VALID_PAIRS:
            paths.append(grch38_toolref_vcf_path(batch, sample, tool, source))

    return paths


ALL_GRCH38_TOOLREF_48WAY_VCFS = grch38_toolref_48way_paths()


rule grch38_toolref_48way_list:
    input:
        vcfs=ALL_GRCH38_TOOLREF_48WAY_VCFS
    output:
        listfile=f"{OUTDIR}/cohort_results/GRCh38_toolref_48way_vcf_list.txt",
        metadata=f"{OUTDIR}/cohort_results/GRCh38_toolref_48way_vcf_metadata.tsv"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    run:
        os.makedirs(os.path.dirname(output.listfile), exist_ok=True)

        with open(output.listfile, "w") as out:
            for vcf in input.vcfs:
                out.write(f"{vcf}\n")

        with open(output.metadata, "w") as out:
            out.write("vector_index\tset_index\tset_label\ttool\tsource\tbatch\tsample\tpath\n")

            vector_index = 0
            for set_index, (set_label, tool, source) in enumerate(TOOLREF_SET_ORDER_GRCH38):
                for batch, sample in VALID_PAIRS:
                    path = grch38_toolref_vcf_path(batch, sample, tool, source)
                    out.write(
                        f"{vector_index}\t{set_index}\t{set_label}\t{tool}\t"
                        f"{source}\t{batch}\t{sample}\t{path}\n"
                    )
                    vector_index += 1


# ==============================================================================
# GRCh38 DIRECT 48-WAY SURVIVOR MERGE
# ==============================================================================

rule grch38_toolref_48way_merge:
    input:
        listfile=f"{OUTDIR}/cohort_results/GRCh38_toolref_48way_vcf_list.txt"
    output:
        vcf=f"{OUTDIR}/cohort_results/GRCh38_toolref_48way_merged.vcf"
    log:
        f"{OUTDIR}/cohort_results/logs/GRCh38_toolref_48way_merge.log"
    params:
        survivor="SURVIVOR",
        max_dist=500,
        min_support=1,
        use_type=1,
        use_strand=0,
        use_size=1,
        min_size=30
    conda:
        "envs/crossmap_truvari.yml"
    threads: 4
    resources:
        mem_mb=config["hm"],
        time=config["ht"]
    shell:
        r"""
        mkdir -p $(dirname {log})

        {params.survivor} merge \
            {input.listfile} \
            {params.max_dist} \
            {params.min_support} \
            {params.use_type} \
            {params.use_strand} \
            {params.use_size} \
            {params.min_size} \
            {output.vcf} > {log} 2>&1
        """


# ==============================================================================
# ANNOTATE GRCh38 WITH TOOL/REFERENCE SUPPORT
# ==============================================================================

rule annotate_grch38_toolref_support:
    input:
        vcf=f"{OUTDIR}/cohort_results/GRCh38_toolref_48way_merged.vcf",
        metadata=f"{OUTDIR}/cohort_results/GRCh38_toolref_48way_vcf_metadata.tsv"
    output:
        vcf=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf"
    conda:
        "envs/crossmap_truvari.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    shell:
        r"""
        python add_toolref_support_info.py \
            --vcf {input.vcf} \
            --metadata {input.metadata} \
            --out {output.vcf} \
            --native-source native_grch38 \
            --lifted-source lifted_chm13_to_grch38 \
            --native-info-prefix NATIVE_GRCH38 \
            --lifted-info-prefix LIFTED_CHM13_GRCH38 \
            --any-native-flag ANY_NATIVE_GRCH38 \
            --any-lifted-flag ANY_LIFTED_CHM13_GRCH38
        """


# ==============================================================================
# SORT, INDEX, AND WRITE GRCh38 SUPPORT TABLE
# ==============================================================================

rule grch38_cohort_survivor_merge:
    input:
        vcf=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf"
    output:
        vcfgz=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz",
        tbi=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz.tbi",
        table=f"{OUTDIR}/cohort_results/GRCh38_cohort_support_table.tsv"
    conda:
        "envs/crossmap_truvari.yml"
    threads: 4
    resources:
        mem_mb=config["hm"],
        time=config["ht"]
    shell:
        r"""
        bcftools sort {input.vcf} -O z -o {output.vcfgz}
        tabix -f -p vcf {output.vcfgz}

        bcftools query \
            -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/SUPP\t%INFO/SUPP_VEC\t%INFO/TOOLREF_SUPP\t%INFO/TOOLREF_SUPP_VEC\t%INFO/SAMPLE_SUPP\t%INFO/SAMPLE_SUPP_VEC\t%INFO/NATIVE_GRCH38_SUPP\t%INFO/LIFTED_CHM13_GRCH38_SUPP\n' \
            {output.vcfgz} > {output.table}
        """