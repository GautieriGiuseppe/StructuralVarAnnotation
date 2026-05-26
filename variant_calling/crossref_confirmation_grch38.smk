import os

OUTDIR = config["output"]

# ==============================================================================
# Build confirmation table directly from the integrated GRCh38 cohort INFO fields
# ==============================================================================

rule build_grch38_confirmation_from_info:
    input:
        vcf=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz",
        tbi=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor.vcf.gz.tbi"
    output:
        tsv=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_confirmation.tsv",
        summary=f"{OUTDIR}/cohort_results/GRCh38_final_cohort_survivor_confirmation_summary.json"
    conda:
        "variant_calling/envs/crossmap_truvari.yml"
    threads: 1
    resources:
        mem_mb=config["mm"],
        time=config["mt"]
    shell:
        r"""
        python - <<'PY'
import gzip
import json

vcf_path = "{input.vcf}"
tsv_path = "{output.tsv}"
summary_path = "{output.summary}"

def parse_info(info_str):
    info = {{}}
    if info_str in ("", "."):
        return info
    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = True
    return info

def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

counts = {{
    "confirmed_by_CHM13": 0,
    "GRCh38_only": 0,
    "CHM13_lifted_only": 0,
    "unsupported_unknown": 0,
}}

total = 0
native_records = 0
lifted_records = 0

with open_text(vcf_path) as inp, open(tsv_path, "w") as out:
    out.write(
        "CHROM\\tPOS\\tID\\tSVTYPE\\tSVLEN\\t"
        "TOOLREF_SUPP\\tTOOLREF_SUPP_VEC\\t"
        "SAMPLE_SUPP\\tSAMPLE_SUPP_VEC\\t"
        "NATIVE_GRCH38_SUPP\\tLIFTED_CHM13_GRCH38_SUPP\\t"
        "CrossRef_support\\n"
    )

    for line in inp:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\\n").split("\\t")
        if len(fields) < 8:
            continue

        chrom = fields[0]
        pos = fields[1]
        var_id = fields[2]
        info = parse_info(fields[7])

        svtype = info.get("SVTYPE", ".")
        svlen = info.get("SVLEN", ".")
        toolref_supp = info.get("TOOLREF_SUPP", "0")
        toolref_vec = info.get("TOOLREF_SUPP_VEC", ".")
        sample_supp = info.get("SAMPLE_SUPP", "0")
        sample_vec = info.get("SAMPLE_SUPP_VEC", ".")
        native = int(info.get("NATIVE_GRCH38_SUPP", "0"))
        lifted = int(info.get("LIFTED_CHM13_GRCH38_SUPP", "0"))

        if native > 0:
            native_records += 1
        if lifted > 0:
            lifted_records += 1

        if native > 0 and lifted > 0:
            label = "confirmed_by_CHM13"
        elif native > 0 and lifted == 0:
            label = "GRCh38_only"
        elif native == 0 and lifted > 0:
            label = "CHM13_lifted_only"
        else:
            label = "unsupported_unknown"

        counts[label] += 1
        total += 1

        out.write(
            f"{{chrom}}\\t{{pos}}\\t{{var_id}}\\t{{svtype}}\\t{{svlen}}\\t"
            f"{{toolref_supp}}\\t{{toolref_vec}}\\t"
            f"{{sample_supp}}\\t{{sample_vec}}\\t"
            f"{{native}}\\t{{lifted}}\\t{{label}}\\n"
        )

native_denominator = counts["confirmed_by_CHM13"] + counts["GRCh38_only"]
if native_denominator > 0:
    confirmation_rate_among_native_grch38 = counts["confirmed_by_CHM13"] / native_denominator
else:
    confirmation_rate_among_native_grch38 = None

summary = {{
    "total_records": total,
    "confirmed_by_CHM13": counts["confirmed_by_CHM13"],
    "GRCh38_only": counts["GRCh38_only"],
    "CHM13_lifted_only": counts["CHM13_lifted_only"],
    "unsupported_unknown": counts["unsupported_unknown"],
    "native_grch38_records": native_records,
    "lifted_chm13_grch38_records": lifted_records,
    "confirmation_rate_among_native_grch38": confirmation_rate_among_native_grch38,
    "definition": {{
        "confirmed_by_CHM13": "NATIVE_GRCH38_SUPP > 0 and LIFTED_CHM13_GRCH38_SUPP > 0",
        "GRCh38_only": "NATIVE_GRCH38_SUPP > 0 and LIFTED_CHM13_GRCH38_SUPP == 0",
        "CHM13_lifted_only": "NATIVE_GRCH38_SUPP == 0 and LIFTED_CHM13_GRCH38_SUPP > 0"
    }}
}}

with open(summary_path, "w") as out:
    json.dump(summary, out, indent=4)

print(json.dumps(summary, indent=4))
PY
        """