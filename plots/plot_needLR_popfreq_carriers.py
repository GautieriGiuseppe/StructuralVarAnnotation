#!/usr/bin/env python3

import os
import gzip
import argparse
import math

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


# =============================================================================
# DEFAULTS
# =============================================================================

DEFAULT_BASE_DIR = "/group/dominguez/shared_notebooks/Immune_variation/mapping"

DEFAULT_NEEDLR_VCF = (
    DEFAULT_BASE_DIR
    + "/needLR_output/"
    + "GRCh38_final_cohort_survivor_genotyped_matrix_needLR_cohort/"
    + "GRCh38_final_cohort_survivor_genotyped_matrix.needlr_input_needLR_1kg_v4.0/"
    + "GRCh38_final_cohort_survivor_genotyped_matrix.needlr_input.needLR.4.0.vcf.gz"
)

DEFAULT_OUT_DIR = (
    DEFAULT_BASE_DIR
    + "/cohort_results/"
    + "needlr_population_frequency_carriers_dynamic"
)

DEFAULT_OUT_PREFIX = "needlr_popfreq_violin_carrier_counts"
DEFAULT_SUMMARY_PREFIX = "needlr_popfreq_carrier_counts"

CANONICAL_SVTYPES = {"DEL", "INS", "DUP", "INV"}
ANCESTRIES = ["AFR", "AMR", "EAS", "EUR", "SAS", "ALL"]

ANCESTRY_COLORS = {
    "AFR": "#9ecae1",
    "AMR": "#9ecae1",
    "EAS": "#9ecae1",
    "EUR": "#08519c",
    "SAS": "#9ecae1",
    "ALL": "#9ecae1",
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
        description=(
            "Plot needLR control population frequency distributions "
            "stratified by cohort carrier count."
        )
    )

    parser.add_argument(
        "--vcf",
        default=DEFAULT_NEEDLR_VCF,
        help="Input needLR-annotated VCF.gz."
    )

    parser.add_argument(
        "--out-dir",
        default=DEFAULT_OUT_DIR,
        help="Output directory."
    )

    parser.add_argument(
        "--out-prefix",
        default=DEFAULT_OUT_PREFIX,
        help="Output prefix for violin plots without extension."
    )

    parser.add_argument(
        "--summary-prefix",
        default=DEFAULT_SUMMARY_PREFIX,
        help="Output prefix for summary line plots without extension."
    )

    parser.add_argument(
        "--title-prefix",
        default="needLR control population frequency",
        help="Title prefix used in figures."
    )

    parser.add_argument(
        "--max-carriers",
        type=int,
        default=8,
        help=(
            "Maximum cohort carrier count to include. "
            "Use 3 for a trio, 8 for the original 8-sample cohort, etc."
        )
    )

    parser.add_argument(
        "--min-carriers",
        type=int,
        default=1,
        help="Minimum cohort carrier count to include."
    )

    return parser.parse_args()


# =============================================================================
# HELPERS
# =============================================================================

def open_vcf(path):
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
    }


def as_float(x):
    try:
        if is_missing(x):
            return np.nan
        return float(str(x).split(",")[0])
    except Exception:
        return np.nan


def as_int(x):
    try:
        if is_missing(x):
            return np.nan
        return int(float(str(x).split(",")[0]))
    except Exception:
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

    return s


def get_first_float(info, keys):
    for key in keys:
        val = as_float(info.get(key))
        if pd.notna(val):
            return val
    return np.nan


def get_first_int(info, keys):
    for key in keys:
        val = as_int(info.get(key))
        if pd.notna(val):
            return val
    return np.nan


def get_pop_count(info, ancestry):
    keys = [
        f"Pop_Count_{ancestry}",
        f"Population_Count_{ancestry}",
        f"Control_Pop_Count_{ancestry}",
        f"Allele_Count_{ancestry}",
        f"Control_Allele_Count_{ancestry}",
        f"AC_{ancestry}",
        f"AC_1KGP_{ancestry}",
        f"1KGP_AC_{ancestry}",
    ]

    return get_first_int(info, keys)


def get_pop_freq(info, ancestry):
    keys = [
        f"Pop_Freq_{ancestry}",
        f"Population_Freq_{ancestry}",
        f"Control_Pop_Freq_{ancestry}",
        f"Allele_Freq_{ancestry}",
        f"Control_Allele_Freq_{ancestry}",
        f"AF_{ancestry}",
        f"AF_1KGP_{ancestry}",
        f"1KGP_AF_{ancestry}",
    ]

    return get_first_float(info, keys)


def get_cohort_pop_count(info):
    keys = [
        "Pop_Count_Cohort",
        "Cohort_Pop_Count",
        "Population_Count_Cohort",
        "Allele_Count_Cohort",
        "Cohort_Allele_Count",
    ]

    return get_first_int(info, keys)


def get_cohort_pop_freq(info):
    keys = [
        "Pop_Freq_Cohort",
        "Cohort_Pop_Freq",
        "Population_Freq_Cohort",
        "Allele_Freq_Cohort",
        "Cohort_Allele_Freq",
    ]

    return get_first_float(info, keys)


def write_empty_outputs(out_dir, out_prefix, summary_prefix, carrier_range_label, reason):
    """
    Write placeholder outputs so Snakemake does not fail when there are no records.
    """

    os.makedirs(out_dir, exist_ok=True)

    variant_tsv = os.path.join(
        out_dir,
        f"needlr_variants_carrier_counts_{carrier_range_label}.tsv"
    )

    long_tsv = os.path.join(
        out_dir,
        f"needlr_popfreq_long_carrier_counts_{carrier_range_label}.tsv"
    )

    summary_tsv = os.path.join(
        out_dir,
        f"needlr_popfreq_summary_carrier_counts_{carrier_range_label}.tsv"
    )

    pd.DataFrame({"reason": [reason]}).to_csv(
        variant_tsv,
        sep="\t",
        index=False,
    )

    pd.DataFrame({"reason": [reason]}).to_csv(
        long_tsv,
        sep="\t",
        index=False,
    )

    pd.DataFrame({"metric": ["status"], "value": [reason]}).to_csv(
        summary_tsv,
        sep="\t",
        index=False,
    )

    figure_paths = [
        os.path.join(out_dir, f"{out_prefix}_present_only.png"),
        os.path.join(out_dir, f"{out_prefix}_present_only.pdf"),
        os.path.join(out_dir, f"{out_prefix}_including_absent_variants.png"),
        os.path.join(out_dir, f"{out_prefix}_including_absent_variants.pdf"),
        os.path.join(out_dir, f"{summary_prefix}_summary_lines.png"),
        os.path.join(out_dir, f"{summary_prefix}_summary_lines.pdf"),
    ]

    for path in figure_paths:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.text(
            0.5,
            0.5,
            reason,
            ha="center",
            va="center",
            wrap=True,
            fontsize=12,
        )
        ax.axis("off")
        fig.tight_layout()
        fig.savefig(path, dpi=300, bbox_inches="tight")
        plt.close(fig)

    print("No data available. Placeholder outputs written:")
    print(variant_tsv)
    print(long_tsv)
    print(summary_tsv)
    for path in figure_paths:
        print(path)


# =============================================================================
# LOAD NEEDLR VCF
# =============================================================================

def load_needlr_carriers(vcf_path, carrier_counts):
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(vcf_path)

    rows = []

    with open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, vid, ref, alt, qual, filt, info_str = cols[:8]
            info = parse_info(info_str)

            svtype = normalize_svtype(info.get("SVTYPE", "OTHER"))
            if svtype not in CANONICAL_SVTYPES:
                continue

            cohort_pop_count = get_cohort_pop_count(info)
            if pd.isna(cohort_pop_count):
                continue

            cohort_pop_count = int(cohort_pop_count)

            if cohort_pop_count not in carrier_counts:
                continue

            svlen = as_int(info.get("SVLEN"))
            if pd.notna(svlen):
                svlen = abs(int(svlen))

            row = {
                "CHROM": chrom,
                "POS": int(pos),
                "ID": vid,
                "SVTYPE": svtype,
                "SVLEN": svlen,
                "Cohort_Pop_Count": cohort_pop_count,
                "Cohort_Pop_Freq": get_cohort_pop_freq(info),
            }

            for anc in ANCESTRIES:
                row[f"Pop_Count_{anc}"] = get_pop_count(info, anc)
                row[f"Pop_Freq_{anc}"] = get_pop_freq(info, anc)

            rows.append(row)

    df = pd.DataFrame(rows)

    if df.empty:
        return df

    return df


# =============================================================================
# LONG TABLE AND SUMMARY
# =============================================================================

def make_long_popfreq_table(df):
    rows = []

    for _, r in df.iterrows():
        for anc in ANCESTRIES:
            pop_count = r.get(f"Pop_Count_{anc}", np.nan)
            pop_freq = r.get(f"Pop_Freq_{anc}", np.nan)

            present = False
            if pd.notna(pop_count):
                present = pop_count > 0
            elif pd.notna(pop_freq):
                present = pop_freq > 0

            rows.append({
                "CHROM": r["CHROM"],
                "POS": r["POS"],
                "ID": r["ID"],
                "SVTYPE": r["SVTYPE"],
                "SVLEN": r["SVLEN"],
                "Cohort_Pop_Count": r["Cohort_Pop_Count"],
                "ancestry": anc,
                "Pop_Count": pop_count,
                "Pop_Freq": pop_freq,
                "present_in_ancestry": present,
            })

    return pd.DataFrame(rows)


def summarize_pop_freq(df, carrier_counts):
    summary_rows = []

    for carrier_count in carrier_counts:
        d = df[df["Cohort_Pop_Count"] == carrier_count].copy()

        for anc in ANCESTRIES:
            count_col = f"Pop_Count_{anc}"
            freq_col = f"Pop_Freq_{anc}"

            if len(d) == 0:
                summary_rows.append({
                    "Cohort_Pop_Count": carrier_count,
                    "ancestry": anc,
                    "n_variants": 0,
                    "n_present": 0,
                    "pct_present": np.nan,
                    "mean_pop_freq_all_variants": np.nan,
                    "mean_pop_freq_present_only": np.nan,
                    "median_pop_freq_present_only": np.nan,
                    "max_pop_freq": np.nan,
                })
                continue

            present = pd.Series(False, index=d.index)

            if count_col in d.columns:
                present = d[count_col].fillna(0) > 0
            elif freq_col in d.columns:
                present = d[freq_col].fillna(0) > 0

            vals_all = d[freq_col].fillna(0)
            vals_present = d.loc[present, freq_col].dropna()

            summary_rows.append({
                "Cohort_Pop_Count": carrier_count,
                "ancestry": anc,
                "n_variants": len(d),
                "n_present": int(present.sum()),
                "pct_present": 100 * present.mean(),
                "mean_pop_freq_all_variants": vals_all.mean(),
                "mean_pop_freq_present_only": vals_present.mean() if len(vals_present) else np.nan,
                "median_pop_freq_present_only": vals_present.median() if len(vals_present) else np.nan,
                "max_pop_freq": vals_all.max(),
            })

    return pd.DataFrame(summary_rows)


# =============================================================================
# PLOTS
# =============================================================================

def plot_popfreq_violin(
    long_df,
    out_prefix,
    title_prefix,
    carrier_counts,
    present_only=True,
):
    d = long_df.copy()

    if present_only:
        d = d[d["present_in_ancestry"]].copy()
        title_suffix = "among present variants"
        filename_suffix = "present_only"
    else:
        d["Pop_Freq"] = d["Pop_Freq"].fillna(0)
        title_suffix = "including absent variants"
        filename_suffix = "including_absent_variants"

    d = d[d["Pop_Freq"].notna()].copy()

    if d.empty:
        out_png = f"{out_prefix}_{filename_suffix}.png"
        out_pdf = f"{out_prefix}_{filename_suffix}.pdf"

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.text(
            0.5,
            0.5,
            f"No population-frequency data available {title_suffix}",
            ha="center",
            va="center",
            wrap=True,
        )
        ax.axis("off")
        fig.tight_layout()
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
        plt.close(fig)

        print("Saved:", out_png)
        print("Saved:", out_pdf)
        return

    d["Pop_Freq_percent"] = d["Pop_Freq"] * 100

    d["ancestry"] = pd.Categorical(
        d["ancestry"],
        categories=ANCESTRIES,
        ordered=True,
    )

    d["Cohort_Pop_Count"] = pd.Categorical(
        d["Cohort_Pop_Count"],
        categories=carrier_counts,
        ordered=True,
    )

    n = len(carrier_counts)
    ncols = min(4, n)
    nrows = int(math.ceil(n / ncols))

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.5 * ncols, 4.0 * nrows),
        sharey=True,
    )

    axes = np.atleast_1d(axes).flatten()

    for ax, carrier_count in zip(axes, carrier_counts):
        sub = d[d["Cohort_Pop_Count"] == carrier_count].copy()

        if sub.empty:
            ax.set_visible(False)
            continue

        sns.violinplot(
            data=sub,
            x="ancestry",
            y="Pop_Freq_percent",
            order=ANCESTRIES,
            inner=None,
            cut=0,
            linewidth=0.8,
            palette=ANCESTRY_COLORS,
            ax=ax,
        )

        ax.set_title(
            f"Cohort carrier count = {carrier_count}",
            fontsize=11,
            fontweight="bold",
        )

        ax.set_xlabel("")
        ax.set_ylabel("Pop_Freq (%)")

        ax.grid(axis="y", color="#EAEAEA", linewidth=0.8)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        n_variants = (
            sub[["CHROM", "POS", "SVTYPE", "SVLEN"]]
            .drop_duplicates()
            .shape[0]
        )

        ax.text(
            0.02,
            0.96,
            f"n={n_variants:,}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8,
        )

    for ax in axes[len(carrier_counts):]:
        ax.set_visible(False)

    carrier_label = f"{min(carrier_counts)}-{max(carrier_counts)}"

    fig.suptitle(
        f"{title_prefix} distributions {title_suffix}\n"
        f"stratified by cohort carrier count {carrier_label} and 1KGP ancestry",
        fontsize=16,
        fontweight="bold",
        y=0.99,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.92])

    out_png = f"{out_prefix}_{filename_suffix}.png"
    out_pdf = f"{out_prefix}_{filename_suffix}.pdf"

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("Saved:", out_png)
    print("Saved:", out_pdf)


def plot_summary_lines(
    summary_df,
    out_prefix,
    title_prefix,
    carrier_counts,
):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]

    for anc in ANCESTRIES:
        d = summary_df[
            summary_df["ancestry"] == anc
        ].sort_values("Cohort_Pop_Count")

        color = ANCESTRY_COLORS[anc]
        lw = 2.8 if anc == "EUR" else 1.8

        ax.plot(
            d["Cohort_Pop_Count"],
            d["mean_pop_freq_present_only"] * 100,
            marker="o",
            linewidth=lw,
            color=color,
            label=anc,
        )

    ax.set_xlabel("Cohort carrier count")
    ax.set_ylabel("Mean Pop_Freq among present SVs (%)")
    ax.set_title("A  Mean control population frequency", loc="left", fontweight="bold")
    ax.set_xticks(carrier_counts)
    ax.grid(axis="y", color="#EAEAEA")
    ax.legend(frameon=False, ncol=3)

    ax = axes[1]

    for anc in ANCESTRIES:
        d = summary_df[
            summary_df["ancestry"] == anc
        ].sort_values("Cohort_Pop_Count")

        color = ANCESTRY_COLORS[anc]
        lw = 2.8 if anc == "EUR" else 1.8

        ax.plot(
            d["Cohort_Pop_Count"],
            d["pct_present"],
            marker="o",
            linewidth=lw,
            color=color,
            label=anc,
        )

    ax.set_xlabel("Cohort carrier count")
    ax.set_ylabel("% of SVs present")
    ax.set_title("B  Presence across control populations", loc="left", fontweight="bold")
    ax.set_xticks(carrier_counts)
    ax.grid(axis="y", color="#EAEAEA")
    ax.legend(frameon=False, ncol=3)

    carrier_label = f"{min(carrier_counts)}-{max(carrier_counts)}"

    fig.suptitle(
        f"{title_prefix} context for cohort carrier counts {carrier_label}",
        fontsize=16,
        fontweight="bold",
        y=0.98,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.92])

    out_png = f"{out_prefix}_summary_lines.png"
    out_pdf = f"{out_prefix}_summary_lines.pdf"

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("Saved:", out_png)
    print("Saved:", out_pdf)


# =============================================================================
# MAIN
# =============================================================================

def main():
    args = parse_args()

    if args.min_carriers < 1:
        raise ValueError("--min-carriers must be >= 1")

    if args.max_carriers < args.min_carriers:
        raise ValueError("--max-carriers must be >= --min-carriers")

    carrier_counts = list(range(args.min_carriers, args.max_carriers + 1))
    carrier_range_label = f"{args.min_carriers}_{args.max_carriers}"

    os.makedirs(args.out_dir, exist_ok=True)

    print("Loading:")
    print(args.vcf)
    print("Carrier counts:", carrier_counts)

    out_prefix = f"{args.out_prefix}_{carrier_range_label}"
    summary_prefix = f"{args.summary_prefix}_{carrier_range_label}"

    violin_prefix = os.path.join(
        args.out_dir,
        out_prefix,
    )

    summary_plot_prefix = os.path.join(
        args.out_dir,
        summary_prefix,
    )

    try:
        df = load_needlr_carriers(
            vcf_path=args.vcf,
            carrier_counts=carrier_counts,
        )
    except Exception:
        raise

    if df.empty:
        reason = (
            f"No canonical variants found with cohort carrier count "
            f"{args.min_carriers}-{args.max_carriers}. "
            "Check the needLR VCF path and INFO field names."
        )

        write_empty_outputs(
            out_dir=args.out_dir,
            out_prefix=out_prefix,
            summary_prefix=summary_prefix,
            carrier_range_label=carrier_range_label,
            reason=reason,
        )
        return

    long_df = make_long_popfreq_table(df)
    summary_df = summarize_pop_freq(
        df=df,
        carrier_counts=carrier_counts,
    )

    variant_tsv = os.path.join(
        args.out_dir,
        f"needlr_variants_carrier_counts_{carrier_range_label}.tsv"
    )

    long_tsv = os.path.join(
        args.out_dir,
        f"needlr_popfreq_long_carrier_counts_{carrier_range_label}.tsv"
    )

    summary_tsv = os.path.join(
        args.out_dir,
        f"needlr_popfreq_summary_carrier_counts_{carrier_range_label}.tsv"
    )

    df.to_csv(variant_tsv, sep="\t", index=False)
    long_df.to_csv(long_tsv, sep="\t", index=False)
    summary_df.to_csv(summary_tsv, sep="\t", index=False)

    plot_popfreq_violin(
        long_df=long_df,
        out_prefix=violin_prefix,
        title_prefix=args.title_prefix,
        carrier_counts=carrier_counts,
        present_only=True,
    )

    plot_popfreq_violin(
        long_df=long_df,
        out_prefix=violin_prefix,
        title_prefix=args.title_prefix,
        carrier_counts=carrier_counts,
        present_only=False,
    )

    plot_summary_lines(
        summary_df=summary_df,
        out_prefix=summary_plot_prefix,
        title_prefix=args.title_prefix,
        carrier_counts=carrier_counts,
    )

    print()
    print("Saved tables:")
    print(variant_tsv)
    print(long_tsv)
    print(summary_tsv)

    print()
    print("Variant counts by cohort carrier count:")
    print(df["Cohort_Pop_Count"].value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()