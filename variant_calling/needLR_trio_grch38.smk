import os

OUTDIR = config["output"]

NEEDLR_OUTDIR = os.path.join(
    OUTDIR,
    config.get("needlr", {}).get("outdir", "needLR_output")
)

TRIO_FILE = config.get("needlr", {}).get("trio_file", None)

TRIOS = []

if TRIO_FILE is not None:
    if os.path.exists(TRIO_FILE):
        with open(TRIO_FILE, "r") as f:
            header = f.readline().rstrip("\n").split("\t")

            required = [
                "family_id",
                "proband",
                "mother",
                "father",
                "proband_batch",
                "mother_batch",
                "father_batch",
            ]

            missing = [c for c in required if c not in header]
            if missing:
                raise ValueError(
                    f"Trio file {TRIO_FILE} is missing columns: {missing}. "
                    f"Available columns: {header}"
                )

            fam_idx = header.index("family_id")
            proband_idx = header.index("proband")
            mother_idx = header.index("mother")
            father_idx = header.index("father")
            proband_batch_idx = header.index("proband_batch")
            mother_batch_idx = header.index("mother_batch")
            father_batch_idx = header.index("father_batch")

            for line in f:
                if line.strip():
                    fields = line.rstrip("\n").split("\t")
                    TRIOS.append({
                        "family_id": fields[fam_idx],
                        "proband": fields[proband_idx],
                        "mother": fields[mother_idx],
                        "father": fields[father_idx],
                        "proband_batch": fields[proband_batch_idx],
                        "mother_batch": fields[mother_batch_idx],
                        "father_batch": fields[father_batch_idx],
                    })
    else:
        print(
            f"WARNING: needLR trio file configured but not found: {TRIO_FILE}. "
            "Trio rules will have no family targets unless --trio is run with a valid trio_file."
        )


FAMILY_IDS = [t["family_id"] for t in TRIOS]
TRIO_BY_FAMILY = {t["family_id"]: t for t in TRIOS}

GENOTYPED_VCF = (
    f"{OUTDIR}/cohort_results/"
    "GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz"
)

GENOTYPED_TBI = (
    f"{OUTDIR}/cohort_results/"
    "GRCh38_final_cohort_survivor_genotyped_matrix.vcf.gz.tbi"
)


rule all_needlr_trio_grch38:
    input:
        GENOTYPED_VCF,
        GENOTYPED_TBI,
        [
            os.path.join(
                NEEDLR_OUTDIR,
                "trio",
                family_id,
                f"{family_id}.needLR_comparator.done"
            )
            for family_id in FAMILY_IDS
        ]


rule subset_trio_sample_vcf_grch38:
    input:
        vcf=GENOTYPED_VCF,
        tbi=GENOTYPED_TBI
    output:
        vcf=os.path.join(
            NEEDLR_OUTDIR,
            "trio",
            "{family_id}",
            "vcfs",
            "{sample}.vcf.gz"
        ),
        tbi=os.path.join(
            NEEDLR_OUTDIR,
            "trio",
            "{family_id}",
            "vcfs",
            "{sample}.vcf.gz.tbi"
        )
    conda:
        "envs/needlr_env.yml"
    threads:
        1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    shell:
        r"""
        mkdir -p $(dirname {output.vcf})

        bcftools view \
            -s {wildcards.sample} \
            -O z \
            -o {output.vcf} \
            {input.vcf}

        tabix -f -p vcf {output.vcf}
        """


def trio_vcfs(wildcards):
    trio = TRIO_BY_FAMILY[wildcards.family_id]

    return {
        "proband": os.path.join(
            NEEDLR_OUTDIR,
            "trio",
            wildcards.family_id,
            "vcfs",
            f"{trio['proband']}.vcf.gz"
        ),
        "mother": os.path.join(
            NEEDLR_OUTDIR,
            "trio",
            wildcards.family_id,
            "vcfs",
            f"{trio['mother']}.vcf.gz"
        ),
        "father": os.path.join(
            NEEDLR_OUTDIR,
            "trio",
            wildcards.family_id,
            "vcfs",
            f"{trio['father']}.vcf.gz"
        ),
    }


rule needlr_trio_comparator_grch38:
    input:
        unpack(trio_vcfs)
    output:
        done=os.path.join(
            NEEDLR_OUTDIR,
            "trio",
            "{family_id}",
            "{family_id}.needLR_comparator.done"
        )
    conda:
        "envs/needlr_env.yml"
    threads:
        config["mc"]
    resources:
        mem_mb=config["hm"],
        time=config["ht"]
    params:
        outdir=lambda wc: os.path.join(
            NEEDLR_OUTDIR,
            "trio",
            wc.family_id,
            "needLR_comparator"
        )
    shell:
        r"""
        mkdir -p {params.outdir}

        needLR comparator \
            -O {params.outdir} \
            -T {threads} \
            -P {input.mother},{input.father} \
            {input.proband}

        touch {output.done}
        """