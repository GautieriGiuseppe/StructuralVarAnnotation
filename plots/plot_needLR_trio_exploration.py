#!/usr/bin/env python3

"""
needLR trio comparator exploration plots
---------------------------------------

Input:
  needLR comparator VCF from trio mode.

Expected INFO fields include:
  Inheritance, SVTYPE, SVLEN, End_Pos, Query_ID,
  Genotype, Alt_Reads, Ref_Reads, Total_Reads,
  Maternal_GT, Maternal_Alt_Reads, Maternal_Ref_Reads, Maternal_Total_Reads,
  Paternal_GT, Paternal_Alt_Reads, Paternal_Ref_Reads, Paternal_Total_Reads,
  Pop_Freq_ALL, Allele_Freq_ALL, Allele_Freq_ALL_Control,
  Genes, CDS, OMIM_phenotype, GENCC_phenotype, HPO_terms,
  Repeat, Segdup, TRE, Centromeric, Pericentromeric, Telomeric, Gap, HiConf.

Output:
  - needlr_trio_exploration_summary.png
  - needlr_trio_exploration_summary.pdf
  - needlr_trio_population_and_priority.png
  - needlr_trio_population_and_priority.pdf
  - needlr_trio_variant_table.tsv
  - needlr_trio_inheritance_counts.tsv
  - needlr_trio_inheritance_svtype_counts.tsv
  - needlr_trio_annotation_summary.tsv
  - needlr_trio_high_priority_candidates.tsv

This script uses pure Python parsing of the VCF.
No bcftools required.
"""

import os
import gzip
import argparse

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter


# =============================================================================
# DEFAULTS
# =============================================================================

DEFAULT_OUT_PREFIX = "needlr_trio_exploration"

CANONICAL_SVTYPES = ["DEL", "INS", "DUP", "INV"]

INHERITANCE_ORDER = [
    "de_novo",
    "maternal",
    "paternal",
    "inherited",
    "unknown",
]

INHERITANCE_COLORS = {
    "de_novo": "#D95F02",
    "maternal": "#1B9E77",
    "paternal": "#7570B3",
    "inherited": "#4C78A8",
    "unknown": "#BDBDBD",
}

SVTYPE_COLORS = {
    "DEL": "#4C78A8",
    "INS": "#F58518",
    "DUP": "#54A24B",
    "INV": "#B279A2",
}

ANNOT_COLORS = {
    "Genic": "#5B83B1",
    "CDS": "#FF8C2A",
    "OMIM": "#63A95B",
    "GENCC": "#5E55C8",
    "HPO": "#2FA882",
    "Rare": "#D95F02",
    "Absent controls": "#A51C30",
}

CONTEXT_COLORS = {
    "HiConf": "#2FA882",
    "Repeat": "#8C8C8C",
    "Segdup": "#B279A2",
    "TRE": "#F58518",
    "Centromeric": "#D62728",
    "Telomeric": "#9467BD",
    "Gap": "#A51C30",
}


plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


# =============================================================================
# ARGPARSE
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot needLR trio comparator exploration summaries."
    )

    parser.add_argument(
        "--vcf",
        required=True,
        help="Input needLR trio comparator VCF or VCF.gz."
    )

    parser.add_argument(
        "--out-dir",
        required=True,
        help="Output directory."
    )

    parser.add_argument(
        "--out-prefix",
        default=DEFAULT_OUT_PREFIX,
        help="Output prefix without extension."
    )

    parser.add_argument(
        "--title-prefix",
        default="needLR trio comparator",
        help="Title prefix used in figures."
    )

    parser.add_argument(
        "--rare-af-threshold",
        type=float,
        default=0.001,
        help="Population allele/frequency threshold used to define rare SVs."
    )

    parser.add_argument(
        "--min-alt-reads",
        type=int,
        default=3,
        help="Minimum proband alternate supporting reads for support-aware candidate counts."
    )

    return parser.parse_args()


# =============================================================================
# HELPERS
# =============================================================================

def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_info(info_str):
    info = {}

    if info_str in {"", "."}:
        return info

    for item in info_str.split(";"):
        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True

    return info


def is_missing(x):
    return x is None or str(x).strip() in {
        "",
        ".",
        "NA",
        "NaN",
        "nan",
        "NAN",
        "None",
        "none",
    }


def is_present(x):
    return not is_missing(x)


def as_int(x):
    try:
        if is_missing(x):
            return np.nan
        return int(float(str(x).split(",")[0]))
    except Exception:
        return np.nan


def as_float(x):
    try:
        if is_missing(x):
            return np.nan
        return float(str(x).split(",")[0])
    except Exception:
        return np.nan


def first_float(info, keys):
    for key in keys:
        value = as_float(info.get(key))
        if pd.notna(value):
            return value
    return np.nan


def first_int(info, keys):
    for key in keys:
        value = as_int(info.get(key))
        if pd.notna(value):
            return value
    return np.nan


def normalize_svtype(x):
    s = str(x).upper()

    if "DEL" in s:
        return "DEL"
    if "INS" in s:
        return "INS"
    if "DUP" in s:
        return "DUP"
    if "INV" in s:
        return "INV"

    return "OTHER"


def normalize_inheritance(x):
    if is_missing(x):
        return "unknown"

    s = str(x).strip()

    if s in {"de_novo", "denovo", "de-novo"}:
        return "de_novo"

    if s in {"maternal", "mother"}:
        return "maternal"

    if s in {"paternal", "father"}:
        return "paternal"

    if s in {"inherited", "parental"}:
        return "inherited"

    return s


def flag_present(info, keys):
    for key in keys:
        value = info.get(key)

        if value is True:
            return True

        if is_missing(value):
            continue

        s = str(value).strip().lower()

        if s in {"true", "yes", "y", "1"}:
            return True

        if s in {"false", "no", "n", "0", "."}:
            continue

        return True

    return False


def alt_fraction(alt_reads, total_reads):
    if pd.isna(alt_reads) or pd.isna(total_reads) or total_reads == 0:
        return np.nan
    return float(alt_reads) / float(total_reads)


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


# =============================================================================
# LOAD VCF
# =============================================================================

def load_trio_vcf(vcf_path, rare_af_threshold):
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(vcf_path)

    rows = []

    with open_text(vcf_path) as fh:
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

            inheritance = normalize_inheritance(info.get("Inheritance"))

            proband_alt = as_int(info.get("Alt_Reads"))
            proband_ref = as_int(info.get("Ref_Reads"))
            proband_total = as_int(info.get("Total_Reads"))

            maternal_alt = as_int(info.get("Maternal_Alt_Reads"))
            maternal_ref = as_int(info.get("Maternal_Ref_Reads"))
            maternal_total = as_int(info.get("Maternal_Total_Reads"))

            paternal_alt = as_int(info.get("Paternal_Alt_Reads"))
            paternal_ref = as_int(info.get("Paternal_Ref_Reads"))
            paternal_total = as_int(info.get("Paternal_Total_Reads"))

            pop_freq_all = first_float(
                info,
                [
                    "Pop_Freq_ALL",
                    "Allele_Freq_ALL",
                    "Allele_Freq_ALL_Control",
                    "AF_ALL",
                    "AF_1KGP_ALL",
                ],
            )

            pop_count_all = first_int(
                info,
                [
                    "Pop_Count_ALL",
                    "Allele_Count_ALL",
                    "Pop_Count_Cohort",
                    "Allele_Count_Cohort",
                ],
            )

            svlen = as_int(info.get("SVLEN"))
            if pd.notna(svlen):
                svlen = abs(int(svlen))

            genes = info.get("Genes", ".")
            cds = info.get("CDS", ".")
            omim = info.get("OMIM_phenotype", ".")
            gencc = info.get("GENCC_phenotype", ".")
            hpo = info.get("HPO_terms", ".")

            row = {
                "CHROM": chrom,
                "POS": int(pos),
                "ID": vid,
                "REF": ref,
                "ALT": alt,
                "FILTER": filt,

                "SVTYPE": svtype,
                "SVLEN": svlen,
                "End_Pos": as_int(info.get("End_Pos")),

                "Query_ID": info.get("Query_ID", "."),
                "Inheritance": inheritance,
                "Control_support": info.get("Control_support", "."),
                "Genotype": info.get("Genotype", "."),

                "Alt_Reads": proband_alt,
                "Ref_Reads": proband_ref,
                "Total_Reads": proband_total,
                "Alt_Fraction": alt_fraction(proband_alt, proband_total),

                "Maternal_GT": info.get("Maternal_GT", "."),
                "Maternal_Alt_Reads": maternal_alt,
                "Maternal_Ref_Reads": maternal_ref,
                "Maternal_Total_Reads": maternal_total,
                "Maternal_Alt_Fraction": alt_fraction(maternal_alt, maternal_total),

                "Paternal_GT": info.get("Paternal_GT", "."),
                "Paternal_Alt_Reads": paternal_alt,
                "Paternal_Ref_Reads": paternal_ref,
                "Paternal_Total_Reads": paternal_total,
                "Paternal_Alt_Fraction": alt_fraction(paternal_alt, paternal_total),

                "Pop_Count_ALL": pop_count_all,
                "Pop_Freq_ALL": pop_freq_all,

                "Genes": genes,
                "CDS": cds,
                "OMIM_phenotype": omim,
                "OMIM_MOI": info.get("OMIM_MOI", "."),
                "GENCC_phenotype": gencc,
                "GENCC_support": info.get("GENCC_support", "."),
                "GENCC_MOI": info.get("GENCC_MOI", "."),
                "HPO_terms": hpo,
                "pLI": info.get("pLI", "."),

                "is_genic": is_present(genes),
                "is_coding": is_present(cds),
                "is_omim": is_present(omim),
                "is_gencc": is_present(gencc),
                "is_hpo": is_present(hpo),

                "HiConf": flag_present(info, ["HiConf", "HighConfidence", "High_confidence"]),
                "Repeat": flag_present(info, ["Repeat", "Repeats", "RepeatMasker"]),
                "Segdup": flag_present(info, ["Segdup", "SegDup", "Segmental_duplication"]),
                "TRE": flag_present(info, ["TRE", "TandemRepeat", "Tandem_Repeats", "STR"]),
                "Centromeric": flag_present(info, ["Centromeric", "Centromere"]),
                "Pericentromeric": flag_present(info, ["Pericentromeric"]),
                "Telomeric": flag_present(info, ["Telomeric", "Telomere"]),
                "Gap": flag_present(info, ["Gap"]),
                "Homopolymer": flag_present(info, ["Homopolymer"]),
            }

            row["is_rare"] = pd.isna(pop_freq_all) or pop_freq_all <= rare_af_threshold
            row["is_absent_controls"] = pd.isna(pop_freq_all) or pop_freq_all == 0

            rows.append(row)

    df = pd.DataFrame(rows)

    if df.empty:
        raise RuntimeError("No canonical trio SV records loaded from VCF.")

    return df


# =============================================================================
# SUMMARIES
# =============================================================================

def add_priority_score(df, min_alt_reads):
    d = df.copy()

    d["priority_score"] = 0

    d.loc[d["Inheritance"].eq("de_novo"), "priority_score"] += 5
    d.loc[d["Inheritance"].isin(["maternal", "paternal", "inherited"]), "priority_score"] += 1

    d.loc[d["is_absent_controls"], "priority_score"] += 3
    d.loc[d["is_rare"], "priority_score"] += 2

    d.loc[d["is_genic"], "priority_score"] += 1
    d.loc[d["is_coding"], "priority_score"] += 3
    d.loc[d["is_omim"], "priority_score"] += 3
    d.loc[d["is_gencc"], "priority_score"] += 3
    d.loc[d["is_hpo"], "priority_score"] += 2

    d.loc[d["Alt_Reads"].fillna(0) >= min_alt_reads, "priority_score"] += 1
    d.loc[d["Alt_Reads"].fillna(0) >= 5, "priority_score"] += 1
    d.loc[d["Total_Reads"].fillna(0) >= 10, "priority_score"] += 1

    d.loc[d["HiConf"], "priority_score"] += 1

    d.loc[d["Gap"], "priority_score"] -= 2
    d.loc[d["Centromeric"], "priority_score"] -= 1
    d.loc[d["Pericentromeric"], "priority_score"] -= 1

    return d


def build_summary_tables(df):
    inheritance_counts = (
        df.groupby("Inheritance")
        .size()
        .reset_index(name="n_variants")
    )

    inheritance_counts["Inheritance"] = pd.Categorical(
        inheritance_counts["Inheritance"],
        categories=[x for x in INHERITANCE_ORDER if x in set(inheritance_counts["Inheritance"])],
        ordered=True,
    )

    inheritance_counts = inheritance_counts.sort_values("Inheritance")

    inheritance_svtype = (
        df.groupby(["Inheritance", "SVTYPE"])
        .size()
        .reset_index(name="n_variants")
    )

    annotation_rows = []

    for inheritance, sub in df.groupby("Inheritance"):
        n = len(sub)

        annotation_rows.append({
            "Inheritance": inheritance,
            "n_variants": n,
            "pct_genic": 100 * sub["is_genic"].mean(),
            "pct_coding": 100 * sub["is_coding"].mean(),
            "pct_omim": 100 * sub["is_omim"].mean(),
            "pct_gencc": 100 * sub["is_gencc"].mean(),
            "pct_hpo": 100 * sub["is_hpo"].mean(),
            "pct_rare": 100 * sub["is_rare"].mean(),
            "pct_absent_controls": 100 * sub["is_absent_controls"].mean(),
            "median_svlen": sub["SVLEN"].median(),
            "median_alt_reads": sub["Alt_Reads"].median(),
            "median_total_reads": sub["Total_Reads"].median(),
            "median_pop_freq_all": sub["Pop_Freq_ALL"].median(),
        })

    annotation_summary = pd.DataFrame(annotation_rows)

    context_items = [
        ("HiConf", "HiConf"),
        ("Repeat", "Repeat"),
        ("Segdup", "Segdup"),
        ("TRE", "TRE"),
        ("Centromeric", "Centromeric"),
        ("Telomeric", "Telomeric"),
        ("Gap", "Gap"),
    ]

    context_rows = []

    for label, col in context_items:
        context_rows.append({
            "Context": label,
            "n_variants": int(df[col].sum()),
            "pct_variants": 100 * df[col].mean(),
        })

    context_summary = pd.DataFrame(context_rows)

    return inheritance_counts, inheritance_svtype, annotation_summary, context_summary


def get_ordered_inheritance_values(df):
    values = list(df["Inheritance"].dropna().unique())
    ordered = [x for x in INHERITANCE_ORDER if x in values]
    extra = sorted([x for x in values if x not in ordered])
    return ordered + extra


# =============================================================================
# PLOT 1: MAIN TRIO EXPLORATION
# =============================================================================

def plot_main_trio_summary(
    df,
    inheritance_counts,
    inheritance_svtype,
    annotation_summary,
    out_prefix,
    title_prefix,
):
    inheritance_order = get_ordered_inheritance_values(df)

    fig, axes = plt.subplots(2, 2, figsize=(15, 10))

    fig.subplots_adjust(
        left=0.08,
        right=0.98,
        bottom=0.08,
        top=0.88,
        wspace=0.30,
        hspace=0.36,
    )

    # -------------------------------------------------------------------------
    # A. Inheritance counts
    # -------------------------------------------------------------------------
    ax = axes[0, 0]
    style_ax(ax)

    counts = inheritance_counts.copy()
    counts["Inheritance"] = counts["Inheritance"].astype(str)
    counts = counts.set_index("Inheritance").reindex(inheritance_order).reset_index()
    counts["n_variants"] = counts["n_variants"].fillna(0)

    colors = [
        INHERITANCE_COLORS.get(x, "#BDBDBD")
        for x in counts["Inheritance"]
    ]

    x = np.arange(len(counts))

    ax.bar(
        x,
        counts["n_variants"],
        color=colors,
        edgecolor="none",
        alpha=0.95,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(counts["Inheritance"], rotation=25, ha="right")
    ax.set_ylabel("Number of SVs")
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_int))
    ax.set_title("A  Inheritance class counts", loc="left", fontweight="bold", fontsize=14)

    ymax = max(counts["n_variants"].max() * 1.15, 1)
    ax.set_ylim(0, ymax)

    for xi, val in zip(x, counts["n_variants"]):
        ax.text(
            xi,
            val + ymax * 0.015,
            f"{int(val):,}",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    # -------------------------------------------------------------------------
    # B. SVTYPE composition by inheritance
    # -------------------------------------------------------------------------
    ax = axes[0, 1]
    style_ax(ax)

    pivot = (
        inheritance_svtype
        .pivot(index="Inheritance", columns="SVTYPE", values="n_variants")
        .fillna(0)
    )

    pivot = pivot.reindex(inheritance_order).fillna(0)

    for svtype in CANONICAL_SVTYPES:
        if svtype not in pivot.columns:
            pivot[svtype] = 0

    pivot = pivot[CANONICAL_SVTYPES]

    bottom = np.zeros(len(pivot))

    x = np.arange(len(pivot))

    for svtype in CANONICAL_SVTYPES:
        vals = pivot[svtype].values

        ax.bar(
            x,
            vals,
            bottom=bottom,
            color=SVTYPE_COLORS[svtype],
            edgecolor="none",
            alpha=0.95,
            label=svtype,
        )

        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(pivot.index.astype(str), rotation=25, ha="right")
    ax.set_ylabel("Number of SVs")
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_int))
    ax.set_title("B  SVTYPE composition by inheritance", loc="left", fontweight="bold", fontsize=14)
    ax.legend(frameon=False, ncol=4)

    # -------------------------------------------------------------------------
    # C. Proband and parental alt fractions
    # -------------------------------------------------------------------------
    ax = axes[1, 0]
    style_ax(ax)

    frac_rows = []

    for _, row in df.iterrows():
        frac_rows.append({
            "Inheritance": row["Inheritance"],
            "sample_role": "proband",
            "alt_fraction": row["Alt_Fraction"],
        })
        frac_rows.append({
            "Inheritance": row["Inheritance"],
            "sample_role": "mother",
            "alt_fraction": row["Maternal_Alt_Fraction"],
        })
        frac_rows.append({
            "Inheritance": row["Inheritance"],
            "sample_role": "father",
            "alt_fraction": row["Paternal_Alt_Fraction"],
        })

    frac_df = pd.DataFrame(frac_rows)
    frac_df = frac_df[frac_df["alt_fraction"].notna()].copy()

    role_order = ["proband", "mother", "father"]
    role_colors = {
        "proband": "#D95F02",
        "mother": "#1B9E77",
        "father": "#7570B3",
    }

    positions = []
    labels = []
    data = []
    colors = []

    pos = 0

    for inh in inheritance_order:
        for role in role_order:
            vals = frac_df[
                (frac_df["Inheritance"] == inh)
                & (frac_df["sample_role"] == role)
            ]["alt_fraction"].values

            if len(vals) == 0:
                vals = [np.nan]

            positions.append(pos)
            labels.append(role)
            data.append(vals)
            colors.append(role_colors[role])
            pos += 1

        pos += 1

    bp = ax.boxplot(
        data,
        positions=positions,
        widths=0.65,
        patch_artist=True,
        showfliers=False,
    )

    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
        patch.set_edgecolor("none")

    for median in bp["medians"]:
        median.set_color("black")
        median.set_linewidth(1.2)

    group_centers = []
    group_labels = []

    pos = 0
    for inh in inheritance_order:
        group_centers.append(pos + 1)
        group_labels.append(inh)
        pos += 4

    ax.set_xticks(group_centers)
    ax.set_xticklabels(group_labels, rotation=25, ha="right")
    ax.set_ylabel("Alt-read fraction")
    ax.set_ylim(0, 1.05)
    ax.set_title("C  Proband and parental alt-read fractions", loc="left", fontweight="bold", fontsize=14)

    handles = [
        mpatches.Patch(color=role_colors["proband"], label="proband"),
        mpatches.Patch(color=role_colors["mother"], label="mother"),
        mpatches.Patch(color=role_colors["father"], label="father"),
    ]

    ax.legend(handles=handles, frameon=False, ncol=3, loc="upper right")

    # -------------------------------------------------------------------------
    # D. Annotation burden by inheritance
    # -------------------------------------------------------------------------
    ax = axes[1, 1]
    style_ax(ax)

    ann = annotation_summary.copy()
    ann["Inheritance"] = ann["Inheritance"].astype(str)
    ann = ann.set_index("Inheritance").reindex(inheritance_order).reset_index()

    categories = [
        ("Genic", "pct_genic"),
        ("CDS", "pct_coding"),
        ("OMIM", "pct_omim"),
        ("GENCC", "pct_gencc"),
        ("Rare", "pct_rare"),
        ("Absent controls", "pct_absent_controls"),
    ]

    x = np.arange(len(ann))
    width = 0.12

    for i, (label, col) in enumerate(categories):
        vals = ann[col].fillna(0).values

        ax.bar(
            x + (i - 2.5) * width,
            vals,
            width=width,
            color=ANNOT_COLORS[label],
            edgecolor="none",
            alpha=0.92,
            label=label,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(ann["Inheritance"], rotation=25, ha="right")
    ax.set_ylabel("% of SVs")
    ax.set_ylim(0, 105)
    ax.set_title("D  Annotation burden by inheritance", loc="left", fontweight="bold", fontsize=14)
    ax.legend(frameon=False, ncol=3, fontsize=9)

    fig.suptitle(
        f"{title_prefix}: trio SV inheritance and annotation summary",
        fontsize=18,
        fontweight="bold",
        y=0.985,
    )

    fig.savefig(out_prefix + "_summary.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_prefix + "_summary.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# PLOT 2: FREQUENCY, SIZE, PRIORITY, CONTEXT
# =============================================================================

def plot_population_priority_context(
    df,
    context_summary,
    out_prefix,
    title_prefix,
):
    inheritance_order = get_ordered_inheritance_values(df)

    fig, axes = plt.subplots(2, 2, figsize=(15, 10))

    fig.subplots_adjust(
        left=0.08,
        right=0.98,
        bottom=0.08,
        top=0.88,
        wspace=0.30,
        hspace=0.36,
    )

    # -------------------------------------------------------------------------
    # A. Population frequency
    # -------------------------------------------------------------------------
    ax = axes[0, 0]
    style_ax(ax)

    data = []
    labels = []
    colors = []

    for inh in inheritance_order:
        vals = df[df["Inheritance"] == inh]["Pop_Freq_ALL"].copy()
        vals = vals.fillna(0)
        vals = vals.clip(lower=0)

        data.append(vals.values * 100)
        labels.append(inh)
        colors.append(INHERITANCE_COLORS.get(inh, "#BDBDBD"))

    bp = ax.boxplot(
        data,
        labels=labels,
        patch_artist=True,
        showfliers=False,
    )

    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
        patch.set_edgecolor("none")

    for median in bp["medians"]:
        median.set_color("black")
        median.set_linewidth(1.2)

    ax.set_ylabel("Control population frequency (%)")
    ax.set_xticklabels(labels, rotation=25, ha="right")
    ax.set_title("A  Control population frequency", loc="left", fontweight="bold", fontsize=14)

    # -------------------------------------------------------------------------
    # B. SVLEN distribution
    # -------------------------------------------------------------------------
    ax = axes[0, 1]
    style_ax(ax)

    data = []
    labels = []
    colors = []

    for inh in inheritance_order:
        vals = df[
            (df["Inheritance"] == inh)
            & (df["SVLEN"].notna())
            & (df["SVLEN"] > 0)
        ]["SVLEN"].copy()

        vals = np.log10(vals.values)

        data.append(vals)
        labels.append(inh)
        colors.append(INHERITANCE_COLORS.get(inh, "#BDBDBD"))

    bp = ax.boxplot(
        data,
        labels=labels,
        patch_artist=True,
        showfliers=False,
    )

    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
        patch.set_edgecolor("none")

    for median in bp["medians"]:
        median.set_color("black")
        median.set_linewidth(1.2)

    ax.set_ylabel("log10(SVLEN bp)")
    ax.set_xticklabels(labels, rotation=25, ha="right")
    ax.set_title("B  SV size distribution", loc="left", fontweight="bold", fontsize=14)

    # -------------------------------------------------------------------------
    # C. Priority score distribution
    # -------------------------------------------------------------------------
    ax = axes[1, 0]
    style_ax(ax)

    score_counts = (
        df.groupby(["Inheritance", "priority_score"])
        .size()
        .reset_index(name="n_variants")
    )

    for inh in inheritance_order:
        sub = score_counts[score_counts["Inheritance"] == inh].copy()

        if sub.empty:
            continue

        ax.plot(
            sub["priority_score"],
            sub["n_variants"],
            marker="o",
            linewidth=2.0,
            color=INHERITANCE_COLORS.get(inh, "#BDBDBD"),
            label=inh,
        )

    ax.set_xlabel("Candidate priority score")
    ax.set_ylabel("Number of SVs")
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_int))
    ax.set_title("C  Candidate priority score distribution", loc="left", fontweight="bold", fontsize=14)
    ax.legend(frameon=False)

    # -------------------------------------------------------------------------
    # D. Genomic context burden
    # -------------------------------------------------------------------------
    ax = axes[1, 1]
    style_ax(ax, grid_axis="x")

    context = context_summary.copy()
    context = context.sort_values("pct_variants", ascending=True)

    y = np.arange(len(context))

    colors = [
        CONTEXT_COLORS.get(x, "#8C8C8C")
        for x in context["Context"]
    ]

    ax.barh(
        y,
        context["pct_variants"],
        color=colors,
        edgecolor="none",
        alpha=0.92,
    )

    ax.set_yticks(y)
    ax.set_yticklabels(context["Context"])
    ax.set_xlabel("% of SVs")
    ax.set_title("D  Genomic context flags", loc="left", fontweight="bold", fontsize=14)

    for yi, val in zip(y, context["pct_variants"]):
        ax.text(
            val + 1,
            yi,
            f"{val:.1f}%",
            va="center",
            fontsize=9,
        )

    xmax = max(context["pct_variants"].max() * 1.20, 10)
    ax.set_xlim(0, xmax)

    fig.suptitle(
        f"{title_prefix}: population frequency, SV size and candidate context",
        fontsize=18,
        fontweight="bold",
        y=0.985,
    )

    fig.savefig(out_prefix + "_population_priority.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_prefix + "_population_priority.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# WRITE TABLES
# =============================================================================

def write_outputs(
    df,
    inheritance_counts,
    inheritance_svtype,
    annotation_summary,
    context_summary,
    out_dir,
    out_prefix,
):
    variant_table = os.path.join(out_dir, f"{out_prefix}_variant_table.tsv")
    inheritance_counts_tsv = os.path.join(out_dir, f"{out_prefix}_inheritance_counts.tsv")
    inheritance_svtype_tsv = os.path.join(out_dir, f"{out_prefix}_inheritance_svtype_counts.tsv")
    annotation_summary_tsv = os.path.join(out_dir, f"{out_prefix}_annotation_summary.tsv")
    context_summary_tsv = os.path.join(out_dir, f"{out_prefix}_context_summary.tsv")
    high_priority_tsv = os.path.join(out_dir, f"{out_prefix}_high_priority_candidates.tsv")
    rare_denovo_tsv = os.path.join(out_dir, f"{out_prefix}_rare_de_novo_candidates.tsv")

    df.to_csv(variant_table, sep="\t", index=False)
    inheritance_counts.to_csv(inheritance_counts_tsv, sep="\t", index=False)
    inheritance_svtype.to_csv(inheritance_svtype_tsv, sep="\t", index=False)
    annotation_summary.to_csv(annotation_summary_tsv, sep="\t", index=False)
    context_summary.to_csv(context_summary_tsv, sep="\t", index=False)

    high_priority = df[
        (df["priority_score"] >= 8)
        & (
            df["Inheritance"].eq("de_novo")
            | df["is_coding"]
            | df["is_omim"]
            | df["is_gencc"]
            | df["is_hpo"]
        )
    ].copy()

    high_priority = high_priority.sort_values(
        [
            "priority_score",
            "Inheritance",
            "is_coding",
            "is_omim",
            "is_gencc",
            "Alt_Reads",
        ],
        ascending=[False, True, False, False, False, False],
    )

    rare_denovo = df[
        df["Inheritance"].eq("de_novo")
        & df["is_rare"]
        & (
            df["is_genic"]
            | df["is_coding"]
            | df["is_omim"]
            | df["is_gencc"]
            | df["is_hpo"]
        )
    ].copy()

    rare_denovo = rare_denovo.sort_values(
        [
            "priority_score",
            "is_coding",
            "is_omim",
            "is_gencc",
            "Alt_Reads",
        ],
        ascending=[False, False, False, False, False],
    )

    high_priority.to_csv(high_priority_tsv, sep="\t", index=False)
    rare_denovo.to_csv(rare_denovo_tsv, sep="\t", index=False)

    return {
        "variant_table": variant_table,
        "inheritance_counts": inheritance_counts_tsv,
        "inheritance_svtype": inheritance_svtype_tsv,
        "annotation_summary": annotation_summary_tsv,
        "context_summary": context_summary_tsv,
        "high_priority": high_priority_tsv,
        "rare_denovo": rare_denovo_tsv,
        "n_high_priority": len(high_priority),
        "n_rare_denovo": len(rare_denovo),
    }


# =============================================================================
# MAIN
# =============================================================================

def main():
    args = parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print("Loading needLR trio comparator VCF:")
    print(args.vcf)

    df = load_trio_vcf(
        vcf_path=args.vcf,
        rare_af_threshold=args.rare_af_threshold,
    )

    df = add_priority_score(
        df=df,
        min_alt_reads=args.min_alt_reads,
    )

    inheritance_counts, inheritance_svtype, annotation_summary, context_summary = build_summary_tables(df)

    out_prefix_path = os.path.join(args.out_dir, args.out_prefix)

    plot_main_trio_summary(
        df=df,
        inheritance_counts=inheritance_counts,
        inheritance_svtype=inheritance_svtype,
        annotation_summary=annotation_summary,
        out_prefix=out_prefix_path,
        title_prefix=args.title_prefix,
    )

    plot_population_priority_context(
        df=df,
        context_summary=context_summary,
        out_prefix=out_prefix_path,
        title_prefix=args.title_prefix,
    )

    outputs = write_outputs(
        df=df,
        inheritance_counts=inheritance_counts,
        inheritance_svtype=inheritance_svtype,
        annotation_summary=annotation_summary,
        context_summary=context_summary,
        out_dir=args.out_dir,
        out_prefix=args.out_prefix,
    )

    print()
    print("Loaded canonical trio SVs:", f"{len(df):,}")
    print()
    print("Inheritance counts:")
    print(inheritance_counts.to_string(index=False))
    print()
    print("Annotation summary:")
    print(annotation_summary.to_string(index=False))
    print()
    print("High-priority candidates:", outputs["n_high_priority"])
    print("Rare annotated de novo candidates:", outputs["n_rare_denovo"])
    print()
    print("Saved:")
    print(out_prefix_path + "_summary.png")
    print(out_prefix_path + "_summary.pdf")
    print(out_prefix_path + "_population_priority.png")
    print(out_prefix_path + "_population_priority.pdf")
    print(outputs["variant_table"])
    print(outputs["inheritance_counts"])
    print(outputs["inheritance_svtype"])
    print(outputs["annotation_summary"])
    print(outputs["context_summary"])
    print(outputs["high_priority"])
    print(outputs["rare_denovo"])


if __name__ == "__main__":
    main()