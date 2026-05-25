
"""
needLR annotation burden and support plots
------------------------------------------

Input:
  Final needLR-annotated GRCh38 genotyped cohort VCF.

Output:
  - needlr_annotation_burden_and_support.png
  - needlr_annotation_burden_and_support.pdf
  - needlr_annotation_burden_summary.tsv
  - needlr_carrier_count_distribution.tsv

This script uses pure Python parsing of the VCF.
No bcftools required.
"""

import os
import gzip
import re
from collections import Counter

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter


# =============================================================================
# CONFIG
# =============================================================================

OUTDIR = "/group/dominguez/shared_notebooks/Immune_variation/mapping"

NEEDLR_ANNOTATED_VCF = (
    OUTDIR
    + "/needLR_output/"
    + "GRCh38_final_cohort_survivor_genotyped_matrix_needLR_cohort/"
    + "GRCh38_final_cohort_survivor_genotyped_matrix.needlr_input_needLR_1kg_v4.0/"
    + "GRCh38_final_cohort_survivor_genotyped_matrix.needlr_input.needLR.4.0.vcf.gz"
)

OUT_DIR = OUTDIR + "/cohort_results/needlr_annotation_plots"
OUT_PREFIX = "needlr_annotation_burden_and_support"

CANONICAL_SVTYPES = {"DEL", "INS", "DUP", "INV"}

COLORS = {
    "Non Genic": "#BDBDBD",
    "Genic": "#5B83B1",
    "Exonic": "#FF8C2A",
    "OMIM": "#63A95B",
    "Context": "#5E55C8",
    "Carrier": "#2FA882",
}


# =============================================================================
# HELPERS
# =============================================================================

def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_info(info_str):
    info = {}
    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = True
    return info


def is_present(value):
    if value is None:
        return False
    if value is True:
        return True

    s = str(value).strip()
    if s in {"", ".", "NA", "NaN", "nan", "NAN", "None", "none"}:
        return False

    return True


def as_float(value, default=np.nan):
    try:
        if not is_present(value):
            return default
        return float(str(value).split(",")[0])
    except Exception:
        return default


def as_int(value, default=np.nan):
    try:
        if not is_present(value):
            return default
        return int(float(str(value).split(",")[0]))
    except Exception:
        return default


def any_info_present(info, keys):
    return any(is_present(info.get(k)) for k in keys)


def any_flag_true(info, keys):
    for k in keys:
        if k in info:
            v = info[k]
            if v is True:
                return True
            if str(v).strip() not in {"0", "False", "false", ".", "", "NA"}:
                return True
    return False


def fmt_int(x, _):
    return f"{int(x):,}"


def style_ax(ax, grid_axis="y"):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#CCCCCC")
    ax.spines["bottom"].set_color("#CCCCCC")
    ax.tick_params(labelsize=10, colors="#444444")
    ax.grid(axis=grid_axis, color="#EAEAEA", linewidth=0.8)
    ax.set_axisbelow(True)


def add_bar_labels(ax, bars, values, suffix="", fontsize=9):
    ymax = max([b.get_height() for b in bars] + [1])
    ax.set_ylim(0, ymax * 1.18)

    for b, v in zip(bars, values):
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() + ymax * 0.015,
            f"{v:.1f}{suffix}" if isinstance(v, float) else f"{int(v):,}{suffix}",
            ha="center",
            va="bottom",
            fontsize=fontsize,
            clip_on=False,
        )


def normalize_svtype(raw):
    if raw is None:
        return "OTHER"

    s = str(raw).upper()

    if "DEL" in s:
        return "DEL"
    if "INS" in s:
        return "INS"
    if "DUP" in s:
        return "DUP"
    if "INV" in s:
        return "INV"

    return s


def infer_control_af(info):
    """
    Try several possible needLR control frequency fields.
    Returns ALL/global control allele/population frequency when present.
    """
    candidate_keys = [
        "Pop_Freq_ALL",
        "Allele_Freq_ALL",
        "AF_ALL",
        "AF_1KGP_ALL",
        "gnomAD_AF",
        "gnomad_AF",
    ]

    for k in candidate_keys:
        val = as_float(info.get(k))
        if pd.notna(val):
            return val

    return np.nan


def infer_carrier_count(info):
    """
    Prefer cohort population count.
    Fall back to cohort allele count if needed.
    """
    candidate_keys = [
        "Pop_Count_Cohort",
        "Cohort_Pop_Count",
        "N_CARRIERS",
        "Carrier_Count",
        "AC_Cohort",
        "Allele_Count_Cohort",
        "Cohort_Allele_Count",
    ]

    for k in candidate_keys:
        val = as_int(info.get(k))
        if pd.notna(val):
            return int(val)

    return np.nan


# =============================================================================
# LOAD VCF
# =============================================================================

def load_needlr_vcf(path):
    rows = []

    if not os.path.exists(path):
        raise FileNotFoundError(path)

    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, vid, ref, alt, qual, filt, info_str = cols[:8]
            info = parse_info(info_str)

            svtype = normalize_svtype(info.get("SVTYPE"))

            if svtype not in CANONICAL_SVTYPES:
                continue

            svlen = abs(as_int(info.get("SVLEN"), default=0))

            carrier_count = infer_carrier_count(info)
            control_af = infer_control_af(info)

            genes = info.get("Genes")
            cds = info.get("CDS")
            omim = info.get("OMIM_phenotype", info.get("OMIM"))

            is_genic = is_present(genes)
            is_exonic = is_present(cds)
            is_omim = is_present(omim)

            is_repeat = any_info_present(info, ["Repeat", "Repeats", "RepeatMasker"])
            is_segdup = any_info_present(info, ["Segdup", "SegDup", "Segmental_duplication"])
            is_tre = any_info_present(info, ["TRE", "TandemRepeat", "Tandem_Repeats", "STR"])
            is_telomeric = any_info_present(info, ["Telomeric", "Telomere"])
            is_centromeric = any_info_present(info, ["Centromeric", "Centromere"])
            is_hiconf = any_info_present(info, ["HiConf", "HighConfidence", "High_confidence"])

            rows.append(
                {
                    "CHROM": chrom,
                    "POS": int(pos),
                    "ID": vid,
                    "SVTYPE": svtype,
                    "SVLEN": svlen,
                    "FILTER": filt,

                    "carrier_count": carrier_count,
                    "control_af": control_af,

                    "is_known": pd.notna(control_af) and control_af > 0,
                    "is_novel": pd.isna(control_af) or control_af == 0,

                    "is_singleton": pd.notna(carrier_count) and carrier_count == 1,
                    "is_rare": pd.notna(control_af) and control_af < 0.01,

                    "is_genic": is_genic,
                    "is_non_genic": not is_genic,
                    "is_exonic": is_exonic,
                    "is_omim": is_omim,

                    "is_repeat": is_repeat,
                    "is_segdup": is_segdup,
                    "is_tre": is_tre,
                    "is_telomeric": is_telomeric,
                    "is_centromeric": is_centromeric,
                    "is_hiconf": is_hiconf,
                    "has_genes_field": is_genic,
                    "has_cds_field": is_exonic,
                }
            )

    df = pd.DataFrame(rows)

    if df.empty:
        raise RuntimeError("No canonical SV records loaded from needLR VCF.")

    return df


# =============================================================================
# SUMMARY TABLES
# =============================================================================

def pct(mask, denom):
    if denom == 0:
        return np.nan
    return 100.0 * mask.sum() / denom


def burden_for_subset(df, label):
    n = len(df)

    return {
        "Group": label,
        "N": n,
        "Non Genic": pct(df["is_non_genic"], n),
        "Genic": pct(df["is_genic"], n),
        "Exonic": pct(df["is_exonic"], n),
        "OMIM": pct(df["is_omim"], n),
    }


def build_burden_tables(df):
    rare_df = df[df["is_rare"]].copy()
    novel_df = df[df["is_novel"]].copy()
    singleton_df = df[df["is_singleton"]].copy()
    known_df = df[df["is_known"]].copy()

    burden = pd.DataFrame(
        [
            burden_for_subset(df, "All"),
            burden_for_subset(rare_df, "Rare"),
            burden_for_subset(novel_df, "Novel"),
            burden_for_subset(singleton_df, "Singleton"),
        ]
    )

    known_novel = pd.DataFrame(
        [
            burden_for_subset(known_df, "Known"),
            burden_for_subset(novel_df, "Novel"),
        ]
    )

    context_items = [
        ("Repeat", "is_repeat"),
        ("Segdup", "is_segdup"),
        ("TRE", "is_tre"),
        ("Telomeric", "is_telomeric"),
        ("Centromeric", "is_centromeric"),
        ("HiConf", "is_hiconf"),
        ("Genes", "has_genes_field"),
        ("CDS", "has_cds_field"),
    ]

    context = []
    n = len(df)
    for label, col in context_items:
        context.append(
            {
                "Context": label,
                "Percent": pct(df[col], n),
                "Count": int(df[col].sum()),
                "N": n,
            }
        )

    context = pd.DataFrame(context)

    carrier = (
        df[df["carrier_count"].notna()]
        .assign(carrier_count=lambda x: x["carrier_count"].astype(int))
        .groupby("carrier_count")
        .size()
        .reset_index(name="count")
        .sort_values("carrier_count")
    )

    return burden, known_novel, context, carrier


# =============================================================================
# PLOT
# =============================================================================

def plot_needlr_burden(burden, known_novel, context, carrier, out_prefix):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.subplots_adjust(
        left=0.08,
        right=0.98,
        bottom=0.08,
        top=0.88,
        wspace=0.32,
        hspace=0.38,
    )

    categories = ["Non Genic", "Genic", "Exonic", "OMIM"]
    cat_colors = [COLORS[c] for c in categories]

    # -------------------------------------------------------------------------
    # A. Annotation burden
    # -------------------------------------------------------------------------
    ax = axes[0, 0]
    style_ax(ax)

    x = np.arange(len(burden))
    width = 0.18

    for i, cat in enumerate(categories):
        vals = burden[cat].values
        ax.bar(
            x + (i - 1.5) * width,
            vals,
            width=width,
            color=cat_colors[i],
            edgecolor="none",
            alpha=0.9,
            label=cat,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(burden["Group"])
    ax.set_ylabel("% of SVs")
    ax.set_title("A  Annotation burden", loc="left", fontweight="bold", fontsize=14)
    ax.set_ylim(0, max(75, np.nanmax(burden[categories].values) * 1.15))

    # -------------------------------------------------------------------------
    # B. Known vs novel annotation burden
    # -------------------------------------------------------------------------
    ax = axes[0, 1]
    style_ax(ax)

    x = np.arange(len(known_novel))

    for i, cat in enumerate(categories):
        vals = known_novel[cat].values
        ax.bar(
            x + (i - 1.5) * width,
            vals,
            width=width,
            color=cat_colors[i],
            edgecolor="none",
            alpha=0.9,
            label=cat,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(known_novel["Group"])
    ax.set_ylabel("% of SVs")
    ax.set_title("B  Known vs novel annotation burden", loc="left", fontweight="bold", fontsize=14)
    ax.set_ylim(0, max(75, np.nanmax(known_novel[categories].values) * 1.15))

    # -------------------------------------------------------------------------
    # C. Genomic context burden
    # -------------------------------------------------------------------------
    ax = axes[1, 0]
    style_ax(ax, grid_axis="x")

    context_plot = context.copy()
    context_plot["Context"] = pd.Categorical(
        context_plot["Context"],
        categories=[
            "Repeat",
            "Segdup",
            "TRE",
            "Telomeric",
            "Centromeric",
            "HiConf",
            "Genes",
            "CDS",
        ],
        ordered=True,
    )
    context_plot = context_plot.sort_values("Context")

    y = np.arange(len(context_plot))
    bars = ax.barh(
        y,
        context_plot["Percent"],
        color=COLORS["Context"],
        edgecolor="none",
        alpha=0.9,
    )

    ax.set_yticks(y)
    ax.set_yticklabels(context_plot["Context"])
    ax.invert_yaxis()
    ax.set_xlabel("% of SVs")
    ax.set_title("C  Genomic context burden", loc="left", fontweight="bold", fontsize=14)
    ax.set_xlim(0, max(90, context_plot["Percent"].max() * 1.15))

    # -------------------------------------------------------------------------
    # D. Carrier count distribution
    # -------------------------------------------------------------------------
    ax = axes[1, 1]
    style_ax(ax)

    carrier_plot = carrier.copy()
    carrier_plot = carrier_plot[
        (carrier_plot["carrier_count"] >= 1) &
        (carrier_plot["carrier_count"] <= 8)
    ]

    x = np.arange(len(carrier_plot))
    bars = ax.bar(
        x,
        carrier_plot["count"],
        color=COLORS["Carrier"],
        edgecolor="none",
        alpha=0.95,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(carrier_plot["carrier_count"].astype(str))
    ax.set_xlabel("Number of carriers")
    ax.set_ylabel("Number of SVs")
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_int))
    ax.set_title("D  Carrier count distribution", loc="left", fontweight="bold", fontsize=14)

    # -------------------------------------------------------------------------
    # Shared legend and title
    # -------------------------------------------------------------------------
    handles = [
        mpatches.Patch(color=COLORS["Non Genic"], label="Non Genic"),
        mpatches.Patch(color=COLORS["Genic"], label="Genic"),
        mpatches.Patch(color=COLORS["Exonic"], label="Exonic"),
        mpatches.Patch(color=COLORS["OMIM"], label="OMIM"),
    ]

    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 0.94),
        fontsize=11,
    )

    fig.suptitle(
        "needLR annotation burden and support",
        fontsize=18,
        fontweight="bold",
        y=0.985,
    )

    fig.savefig(out_prefix + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(out_prefix + ".pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("Loading needLR VCF:")
    print(NEEDLR_ANNOTATED_VCF)

    df = load_needlr_vcf(NEEDLR_ANNOTATED_VCF)

    print()
    print("Loaded canonical SVs:", f"{len(df):,}")
    print(df["SVTYPE"].value_counts().to_string())

    burden, known_novel, context, carrier = build_burden_tables(df)

    out_prefix = os.path.join(OUT_DIR, OUT_PREFIX)

    df.to_csv(
        os.path.join(OUT_DIR, "needlr_variant_annotation_table.tsv"),
        sep="\t",
        index=False,
    )

    burden.to_csv(
        os.path.join(OUT_DIR, "needlr_annotation_burden_summary.tsv"),
        sep="\t",
        index=False,
    )

    known_novel.to_csv(
        os.path.join(OUT_DIR, "needlr_known_vs_novel_annotation_burden.tsv"),
        sep="\t",
        index=False,
    )

    context.to_csv(
        os.path.join(OUT_DIR, "needlr_genomic_context_burden.tsv"),
        sep="\t",
        index=False,
    )

    carrier.to_csv(
        os.path.join(OUT_DIR, "needlr_carrier_count_distribution.tsv"),
        sep="\t",
        index=False,
    )

    plot_needlr_burden(burden, known_novel, context, carrier, out_prefix)

    print()
    print("Saved:")
    print(out_prefix + ".png")
    print(out_prefix + ".pdf")
    print()
    print("Annotation burden:")
    print(burden.to_string(index=False))
    print()
    print("Known vs novel:")
    print(known_novel.to_string(index=False))
    print()
    print("Genomic context:")
    print(context.to_string(index=False))
    print()
    print("Carrier count distribution:")
    print(carrier.to_string(index=False))


if __name__ == "__main__":
    main()