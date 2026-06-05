#!/usr/bin/env python3

"""
Cross-reference confirmation plots from integrated GRCh38 cohort.

This version uses the INFO fields already present in:
    GRCh38_final_cohort_survivor.vcf.gz

Confirmation logic:
    confirmed_by_CHM13 = ANY_LIFTED_CHM13_GRCH38
                         or LIFTED_CHM13_GRCH38_SUPP > 0

    GRCh38_only        = no lifted CHM13 support

Canonical SVs:
    DEL, INS, DUP, INV

Outputs:
    crossref_infofield_confirmation_table.tsv
    crossref_infofield_summary_metrics.tsv
    crossref_fig1_summary.png/pdf
    crossref_fig2_confirmation_patterns.png/pdf
    crossref_fig3_chromosome_confirmation.png/pdf
"""

import os
import gzip
import math
import argparse

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


# ==============================================================================
# DEFAULTS
# ==============================================================================

DEFAULT_WORKFLOW_OUTDIR = "/group/dominguez/shared_notebooks/Immune_variation/mapping"
DEFAULT_BASE_DIR = os.path.join(DEFAULT_WORKFLOW_OUTDIR, "cohort_results")

DEFAULT_GRCH38_INTEGRATED = os.path.join(
    DEFAULT_BASE_DIR,
    "GRCh38_final_cohort_survivor.vcf.gz"
)

DEFAULT_OUT_DIR = os.path.join(
    DEFAULT_BASE_DIR,
    "crossref_confirmation_infofield"
)

DEFAULT_BATCH = "NP057"

DEFAULT_SAMPLES = [
    "CCH_095", "CCH_096", "CCH_097", "CCH_098",
    "CCH_099", "CCH_100", "CCH_101", "CCH_102",
]

CANONICAL_SVTYPES = {"DEL", "INS", "DUP", "INV"}
SVTYPE_ORDER = ["DEL", "INS", "DUP", "INV"]

SVTYPE_COLORS = {
    "DEL": "#2980B9",
    "INS": "#27AE60",
    "DUP": "#8E44AD",
    "INV": "#E67E22",
}

STATUS_COLORS = {
    "confirmed_by_CHM13": "#1D9E75",
    "GRCh38_only": "#B8B8B8",
    "native_total": "#7F7F7F",
    "lifted_total": "#4C78A8",
    "not_lifted": "#D85A30",
    "lifted_only": "#4C78A8",
}

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
})


# ==============================================================================
# ARGPARSE
# ==============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate GRCh38 cross-reference confirmation plots using "
            "integrated cohort INFO fields."
        )
    )

    parser.add_argument(
        "--integrated-vcf",
        default=DEFAULT_GRCH38_INTEGRATED,
        help="Input GRCh38 integrated cohort VCF.gz."
    )

    parser.add_argument(
        "--out-dir",
        default=DEFAULT_OUT_DIR,
        help="Output directory for plots and summary tables."
    )

    parser.add_argument(
        "--workflow-outdir",
        default=DEFAULT_WORKFLOW_OUTDIR,
        help=(
            "Main workflow output directory containing batch/sample folders. "
            "Used only for the CHM13 liftover summary panel."
        )
    )

    parser.add_argument(
        "--batch",
        default=DEFAULT_BATCH,
        help="Batch ID used for per-sample paths, e.g. NP057."
    )

    parser.add_argument(
        "--samples",
        default=",".join(DEFAULT_SAMPLES),
        help="Comma-separated sample IDs."
    )

    parser.add_argument(
        "--title-prefix",
        default="Cross-reference confirmation",
        help="Title prefix used in figures."
    )

    parser.add_argument(
        "--skip-liftover-panel",
        action="store_true",
        help=(
            "Skip counting per-tool CHM13 native/lifted files. "
            "If set, the liftover panel is filled with NA metrics."
        )
    )

    return parser.parse_args()


# ==============================================================================
# HELPERS
# ==============================================================================

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
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = True

    return info


def safe_int(x, default=0):
    try:
        if x in [None, ".", "", True]:
            return default
        return int(float(str(x).split(",")[0]))
    except Exception:
        return default


def safe_float(x, default=np.nan):
    try:
        if x in [None, ".", "", True]:
            return default
        return float(str(x).split(",")[0])
    except Exception:
        return default


def normalize_svtype(raw):
    raw = str(raw).upper()

    if "DEL" in raw:
        return "DEL"
    if "INS" in raw:
        return "INS"
    if "DUP" in raw:
        return "DUP"
    if "INV" in raw:
        return "INV"

    return raw


def abs_svlen(raw):
    try:
        if raw in [None, ".", "", True]:
            return 0
        return abs(int(float(str(raw).split(",")[0])))
    except Exception:
        return 0


def has_info_flag(info, flag):
    return flag in info and info[flag] is True


def is_canonical_record(svtype, alt, info):
    svtype = normalize_svtype(svtype)
    alt = str(alt)

    if svtype not in CANONICAL_SVTYPES:
        return False

    if alt == ".":
        return False

    if "BND" in alt or "TRA" in alt:
        return False

    if svtype in {"BND", "TRA"}:
        return False

    end = safe_int(info.get("END"), default=None)
    pos = safe_int(info.get("_POS"), default=None)

    if end is not None and pos is not None:
        if end < pos:
            return False

    return True


def style_ax(ax, grid_axis="y"):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#CCCCCC")
    ax.spines["bottom"].set_color("#CCCCCC")
    ax.tick_params(labelsize=9, colors="#444444")
    ax.grid(axis=grid_axis, color="#EAEAEA", linewidth=0.8)
    ax.set_axisbelow(True)


def fmt_int(x, _):
    return f"{int(x):,}"


def save_figure(fig, prefix, out_dir):
    png = os.path.join(out_dir, prefix + ".png")
    pdf = os.path.join(out_dir, prefix + ".pdf")

    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("Saved:", png)
    print("Saved:", pdf)


def add_bar_labels(ax, bars, values, fontsize=9):
    valid_values = [
        float(v)
        for v in values
        if pd.notna(v)
    ]

    ymax = max([b.get_height() for b in bars] + valid_values + [1])
    ax.set_ylim(0, ymax * 1.18)

    for b, v in zip(bars, values):
        if pd.isna(v):
            label = "NA"
        else:
            label = f"{int(v):,}"

        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() + ymax * 0.015,
            label,
            ha="center",
            va="bottom",
            fontsize=fontsize,
            clip_on=False
        )


def add_percent_labels(ax, bars, values, counts=None, fontsize=8):
    valid_values = [
        float(v)
        for v in values
        if pd.notna(v)
    ]

    ymax = max([b.get_height() for b in bars] + valid_values + [1])
    ax.set_ylim(0, min(115, ymax * 1.25))

    for i, (b, v) in enumerate(zip(bars, values)):
        if pd.isna(v):
            label = "NA"
        elif counts is None:
            label = f"{v:.1f}%"
        else:
            label = f"{v:.1f}%\n(n={int(counts[i]):,})"

        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() + ymax * 0.025,
            label,
            ha="center",
            va="bottom",
            fontsize=fontsize,
            clip_on=False
        )


# ==============================================================================
# LOAD VCF TABLES
# ==============================================================================

def load_integrated_grch38(vcf_path):
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
            info["_POS"] = pos

            svtype = normalize_svtype(info.get("SVTYPE", "OTHER"))
            svlen = abs_svlen(info.get("SVLEN"))

            if filt != "PASS":
                continue

            if not is_canonical_record(svtype, alt, info):
                continue

            native_supp = safe_int(info.get("NATIVE_GRCH38_SUPP"), default=0)
            lifted_supp = safe_int(info.get("LIFTED_CHM13_GRCH38_SUPP"), default=0)

            any_native = has_info_flag(info, "ANY_NATIVE_GRCH38") or native_supp > 0
            any_lifted = has_info_flag(info, "ANY_LIFTED_CHM13_GRCH38") or lifted_supp > 0

            confirmed = any_lifted

            rows.append({
                "CHROM": chrom,
                "POS": int(pos),
                "ID": vid,
                "REF": ref,
                "ALT": alt,
                "FILTER": filt,
                "SVTYPE": svtype,
                "SVLEN": svlen,

                "SUPP": safe_int(info.get("SUPP")),
                "SUPP_VEC": info.get("SUPP_VEC", "."),

                "TOOLREF_SUPP": safe_int(info.get("TOOLREF_SUPP")),
                "TOOLREF_SUPP_VEC": info.get("TOOLREF_SUPP_VEC", "."),

                "SAMPLE_SUPP": safe_int(info.get("SAMPLE_SUPP")),
                "SAMPLE_SUPP_VEC": info.get("SAMPLE_SUPP_VEC", "."),

                "NATIVE_GRCH38_SUPP": native_supp,
                "LIFTED_CHM13_GRCH38_SUPP": lifted_supp,

                "ANY_NATIVE_GRCH38": any_native,
                "ANY_LIFTED_CHM13_GRCH38": any_lifted,

                "confirmed": confirmed,
                "CrossRef_support": (
                    "confirmed_by_CHM13" if confirmed else "GRCh38_only"
                ),
            })

    df = pd.DataFrame(rows)

    if df.empty:
        raise RuntimeError(f"No canonical PASS variants loaded from {vcf_path}")

    return df


def count_vcf_records(vcf_path, canonical_only=False):
    if not os.path.exists(vcf_path):
        return np.nan

    n = 0

    with open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")

            if len(cols) < 8:
                continue

            chrom, pos, vid, ref, alt, qual, filt, info_str = cols[:8]
            info = parse_info(info_str)
            info["_POS"] = pos

            if filt != "PASS":
                continue

            svtype = normalize_svtype(info.get("SVTYPE", "OTHER"))

            if canonical_only:
                if not is_canonical_record(svtype, alt, info):
                    continue

            n += 1

    return n


# ==============================================================================
# TABLE + METRICS
# ==============================================================================

def add_bins(df):
    df = df.copy()

    size_bins = [0, 50, 100, 500, 1000, 5000, 10000, 50000, np.inf]
    size_labels = [
        "1-50",
        "51-100",
        "101-500",
        "501-1k",
        "1-5k",
        "5-10k",
        "10-50k",
        ">50k",
    ]

    df["size_bin"] = pd.cut(
        df["SVLEN"],
        bins=size_bins,
        labels=size_labels,
        right=True
    )

    def carrier_bin(x):
        if pd.isna(x):
            return np.nan
        if x <= 1:
            return "1"
        if x == 2:
            return "2"
        if 3 <= x <= 5:
            return "3-5"
        return ">5"

    df["carrier_bin"] = df["SAMPLE_SUPP"].apply(carrier_bin)
    df["carrier_bin"] = pd.Categorical(
        df["carrier_bin"],
        categories=["1", "2", "3-5", ">5"],
        ordered=True
    )

    return df


def count_existing_per_tool_chm13_liftover_inputs(workflow_outdir, batch, samples):
    tools = {
        "sniffles": "sniffles",
        "delly": "delly",
        "cuteSV": "cuteSV",
    }

    native_total = 0
    lifted_total = 0

    native_by_tool = {}
    lifted_by_tool = {}

    for tool, tool_file in tools.items():
        native_by_tool[tool] = 0
        lifted_by_tool[tool] = 0

        for sample in samples:
            vcf_root = os.path.join(
                workflow_outdir,
                batch,
                sample,
                "03.variant_calling"
            )

            if tool == "delly":
                native_vcf = os.path.join(
                    vcf_root,
                    f"{sample}_chm13_delly_symbolic.vcf"
                )
            else:
                native_vcf = os.path.join(
                    vcf_root,
                    f"{sample}_chm13_{tool_file}.vcf"
                )

            lifted_vcf = os.path.join(
                vcf_root,
                "liftover",
                f"{tool_file}_{sample}_chm13-to-grch38.toolref_pass.vcf"
            )

            n_native = count_vcf_records(native_vcf, canonical_only=True)
            n_lifted = count_vcf_records(lifted_vcf, canonical_only=True)

            if pd.notna(n_native):
                native_total += int(n_native)
                native_by_tool[tool] += int(n_native)

            if pd.notna(n_lifted):
                lifted_total += int(n_lifted)
                lifted_by_tool[tool] += int(n_lifted)

    return native_total, lifted_total, native_by_tool, lifted_by_tool


def build_metrics(df, workflow_outdir, batch, samples, skip_liftover_panel=False):
    if skip_liftover_panel:
        chm13_native_raw = np.nan
        chm13_lifted_raw = np.nan
        chm13_not_lifted_raw = np.nan
        pct_liftover_success_raw = np.nan
        native_by_tool = {
            "sniffles": np.nan,
            "delly": np.nan,
            "cuteSV": np.nan,
        }
        lifted_by_tool = {
            "sniffles": np.nan,
            "delly": np.nan,
            "cuteSV": np.nan,
        }
    else:
        chm13_native_raw, chm13_lifted_raw, native_by_tool, lifted_by_tool = (
            count_existing_per_tool_chm13_liftover_inputs(
                workflow_outdir=workflow_outdir,
                batch=batch,
                samples=samples,
            )
        )

        chm13_not_lifted_raw = chm13_native_raw - chm13_lifted_raw

        pct_liftover_success_raw = (
            100 * chm13_lifted_raw / chm13_native_raw
            if chm13_native_raw > 0
            else np.nan
        )

    lifted_integrated = df["ANY_LIFTED_CHM13_GRCH38"]
    native_integrated = df["ANY_NATIVE_GRCH38"]

    metrics = {
        "n_grch38_canonical": len(df),
        "n_confirmed_by_chm13": int(df["confirmed"].sum()),
        "n_grch38_only": int((~df["confirmed"]).sum()),
        "pct_confirmed_by_chm13": 100 * df["confirmed"].mean(),

        "n_chm13_native_raw": chm13_native_raw,
        "n_chm13_lifted_raw": chm13_lifted_raw,
        "n_chm13_not_lifted_raw": chm13_not_lifted_raw,
        "pct_chm13_liftover_success_raw": pct_liftover_success_raw,

        "n_integrated_with_lifted_chm13_support": int(lifted_integrated.sum()),
        "n_integrated_shared_native_and_lifted": int((lifted_integrated & native_integrated).sum()),
        "n_integrated_lifted_only": int((lifted_integrated & ~native_integrated).sum()),
        "n_integrated_native_only": int((native_integrated & ~lifted_integrated).sum()),
    }

    for tool in native_by_tool:
        metrics[f"chm13_native_{tool}"] = native_by_tool[tool]
        metrics[f"chm13_lifted_{tool}"] = lifted_by_tool[tool]

        if pd.notna(native_by_tool[tool]) and pd.notna(lifted_by_tool[tool]):
            metrics[f"chm13_lost_{tool}"] = native_by_tool[tool] - lifted_by_tool[tool]
        else:
            metrics[f"chm13_lost_{tool}"] = np.nan

    return metrics


# ==============================================================================
# FIGURE 1: SUMMARY
# ==============================================================================

def plot_summary(df, metrics, out_dir, title_prefix):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # --------------------------------------------------------------------------
    # A. GRCh38 confirmation
    # --------------------------------------------------------------------------
    ax = axes[0]
    style_ax(ax)

    labels = [
        "GRCh38\ncanonical",
        "Confirmed",
        "GRCh38\nonly",
    ]

    values = [
        metrics["n_grch38_canonical"],
        metrics["n_confirmed_by_chm13"],
        metrics["n_grch38_only"],
    ]

    colors = [
        STATUS_COLORS["native_total"],
        STATUS_COLORS["confirmed_by_CHM13"],
        STATUS_COLORS["GRCh38_only"],
    ]

    bars = ax.bar(labels, values, color=colors, edgecolor="none", alpha=0.9)
    ax.set_ylabel("Number of SVs")
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_int))
    ax.set_title(
        "A  GRCh38 confirmation\ncanonical SVs only",
        loc="left",
        fontweight="bold"
    )
    add_bar_labels(ax, bars, values)

    pct = metrics["pct_confirmed_by_chm13"]
    ax.text(
        0.5,
        0.93,
        f"{pct:.1f}% confirmed",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=11,
        fontweight="bold",
        color=STATUS_COLORS["confirmed_by_CHM13"]
    )

    # --------------------------------------------------------------------------
    # B. CHM13 liftover outcome from existing per-tool CrossMap outputs
    # --------------------------------------------------------------------------
    ax = axes[1]
    style_ax(ax)

    labels = [
        "Native CHM13\nper-tool SVs",
        "Lifted to\nGRCh38",
        "Not lifted",
    ]

    values = [
        metrics["n_chm13_native_raw"],
        metrics["n_chm13_lifted_raw"],
        metrics["n_chm13_not_lifted_raw"],
    ]

    plot_values = [
        0 if pd.isna(v) else v
        for v in values
    ]

    colors = [
        STATUS_COLORS["native_total"],
        STATUS_COLORS["lifted_total"],
        STATUS_COLORS["not_lifted"],
    ]

    bars = ax.bar(labels, plot_values, color=colors, edgecolor="none", alpha=0.9)
    ax.set_ylabel("Number of SVs")
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_int))
    ax.set_title(
        "B  CHM13-to-GRCh38 liftover outcome\nexisting per-tool CrossMap outputs",
        loc="left",
        fontweight="bold"
    )
    add_bar_labels(ax, bars, values)

    if pd.notna(metrics["pct_chm13_liftover_success_raw"]):
        pct_lost = 100 - metrics["pct_chm13_liftover_success_raw"]
        text = f"{pct_lost:.1f}% not lifted"
    else:
        text = "Liftover counts not available"

    ax.text(
        0.5,
        0.93,
        text,
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=11,
        fontweight="bold",
        color=STATUS_COLORS["not_lifted"]
    )

    fig.suptitle(
        f"{title_prefix} and liftover summary",
        fontsize=16,
        fontweight="bold",
        y=0.98
    )

    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save_figure(fig, "crossref_fig1_summary", out_dir)


# ==============================================================================
# FIGURE 2: CONFIRMATION PATTERNS
# ==============================================================================

def plot_confirmation_patterns(df, out_dir, title_prefix):
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10))

    # --------------------------------------------------------------------------
    # A. SV type
    # --------------------------------------------------------------------------
    ax = axes[0, 0]
    style_ax(ax)

    rate = (
        df.groupby("SVTYPE")["confirmed"]
        .mean()
        .reindex(SVTYPE_ORDER)
        .dropna() * 100
    )

    counts = (
        df.groupby("SVTYPE")
        .size()
        .reindex(rate.index)
    )

    bars = ax.bar(
        rate.index,
        rate.values,
        color=[SVTYPE_COLORS[t] for t in rate.index],
        edgecolor="none",
        alpha=0.9
    )

    ax.set_ylabel("Confirmed GRCh38 SVs (%)")
    ax.set_xlabel("SV type")
    ax.set_title(
        "A  Confirmation rate by SV type",
        loc="left",
        fontweight="bold"
    )
    add_percent_labels(ax, bars, rate.values, counts.values)

    # --------------------------------------------------------------------------
    # B. Size
    # --------------------------------------------------------------------------
    ax = axes[0, 1]
    style_ax(ax)

    d = df[df["size_bin"].notna()].copy()

    rate = d.groupby("size_bin", observed=False)["confirmed"].mean() * 100
    counts = d.groupby("size_bin", observed=False).size()

    keep = counts[counts > 0].index
    rate = rate.loc[keep]
    counts = counts.loc[keep]

    x = np.arange(len(rate))

    ax.plot(
        x,
        rate.values,
        marker="o",
        linewidth=2.2,
        color="#534AB7"
    )

    ax.set_xticks(x)
    ax.set_xticklabels(rate.index.astype(str), rotation=30)
    ax.set_ylabel("Confirmed GRCh38 SVs (%)")
    ax.set_xlabel("SV length bin (bp)")
    ax.set_title(
        "B  Confirmation rate by size",
        loc="left",
        fontweight="bold"
    )

    ymax = max(rate.values) if len(rate) else 1
    ax.set_ylim(0, min(115, ymax * 1.25))

    for xi, yi, n in zip(x, rate.values, counts.values):
        ax.text(
            xi,
            yi + ymax * 0.035,
            f"n={int(n):,}",
            ha="center",
            va="bottom",
            fontsize=8
        )

    # --------------------------------------------------------------------------
    # C. Counts by type and status
    # --------------------------------------------------------------------------
    ax = axes[1, 0]
    style_ax(ax)

    tab = (
        df.groupby(["SVTYPE", "CrossRef_support"])
        .size()
        .unstack(fill_value=0)
        .reindex(SVTYPE_ORDER)
        .fillna(0)
    )

    x = np.arange(len(tab.index))
    width = 0.38

    confirmed_vals = (
        tab["confirmed_by_CHM13"].values
        if "confirmed_by_CHM13" in tab.columns
        else np.zeros(len(tab))
    )

    only_vals = (
        tab["GRCh38_only"].values
        if "GRCh38_only" in tab.columns
        else np.zeros(len(tab))
    )

    ax.bar(
        x - width / 2,
        confirmed_vals,
        width=width,
        color=STATUS_COLORS["confirmed_by_CHM13"],
        edgecolor="none",
        alpha=0.9,
        label="Confirmed"
    )

    ax.bar(
        x + width / 2,
        only_vals,
        width=width,
        color=STATUS_COLORS["GRCh38_only"],
        edgecolor="none",
        alpha=0.9,
        label="GRCh38 only"
    )

    ax.set_xticks(x)
    ax.set_xticklabels(tab.index)
    ax.set_ylabel("Number of SVs")
    ax.set_xlabel("SV type")
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_int))
    ax.set_title(
        "C  Counts by SV type",
        loc="left",
        fontweight="bold"
    )
    ax.legend(frameon=False)

    # --------------------------------------------------------------------------
    # D. Carrier count
    # --------------------------------------------------------------------------
    ax = axes[1, 1]
    style_ax(ax)

    d = df[df["carrier_bin"].notna()].copy()

    rate = d.groupby("carrier_bin", observed=False)["confirmed"].mean() * 100
    counts = d.groupby("carrier_bin", observed=False).size()

    keep = counts[counts > 0].index
    rate = rate.loc[keep]
    counts = counts.loc[keep]

    x = np.arange(len(rate))

    ax.plot(
        x,
        rate.values,
        marker="o",
        linewidth=2.2,
        color=STATUS_COLORS["confirmed_by_CHM13"]
    )

    ax.set_xticks(x)
    ax.set_xticklabels(rate.index.astype(str))
    ax.set_ylabel("Confirmed GRCh38 SVs (%)")
    ax.set_xlabel("Sample carrier count bin")
    ax.set_title(
        "D  Confirmation rate by carrier count",
        loc="left",
        fontweight="bold"
    )

    ymax = max(rate.values) if len(rate) else 1
    ax.set_ylim(0, min(115, ymax * 1.25))

    for xi, yi, n in zip(x, rate.values, counts.values):
        ax.text(
            xi,
            yi + ymax * 0.035,
            f"n={int(n):,}",
            ha="center",
            va="bottom",
            fontsize=8
        )

    fig.suptitle(
        f"{title_prefix} patterns from integrated GRCh38 cohort",
        fontsize=16,
        fontweight="bold",
        y=0.98
    )

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save_figure(fig, "crossref_fig2_confirmation_patterns", out_dir)


# ==============================================================================
# FIGURE 3: CHROMOSOMES
# ==============================================================================

def plot_chromosome_confirmation(df, out_dir, title_prefix):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    d = df[df["CHROM"].isin(chrom_order)].copy()

    # --------------------------------------------------------------------------
    # A. Counts
    # --------------------------------------------------------------------------
    ax = axes[0]
    style_ax(ax)

    tab = (
        d.groupby(["CHROM", "CrossRef_support"])
        .size()
        .unstack(fill_value=0)
        .reindex(chrom_order)
        .fillna(0)
    )

    tab = tab.loc[tab.sum(axis=1) > 0]

    x = np.arange(len(tab.index))
    bottom = np.zeros(len(tab.index))

    for status in ["GRCh38_only", "confirmed_by_CHM13"]:
        if status not in tab.columns:
            continue

        vals = tab[status].values

        ax.bar(
            x,
            vals,
            bottom=bottom,
            color=STATUS_COLORS[status],
            edgecolor="none",
            alpha=0.9,
            label=status.replace("_", " ")
        )

        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("chr", "") for c in tab.index], rotation=45)
    ax.set_ylabel("Number of GRCh38 SVs")
    ax.set_xlabel("Chromosome")
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_int))
    ax.set_title(
        "A  Confirmed vs GRCh38-only counts",
        loc="left",
        fontweight="bold"
    )
    ax.legend(frameon=False)

    # --------------------------------------------------------------------------
    # B. Rate
    # --------------------------------------------------------------------------
    ax = axes[1]
    style_ax(ax)

    rate = (
        d.groupby("CHROM")["confirmed"]
        .mean()
        .reindex(chrom_order)
        .dropna() * 100
    )

    x = np.arange(len(rate.index))

    ax.bar(
        x,
        rate.values,
        color="#534AB7",
        edgecolor="none",
        alpha=0.9
    )

    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("chr", "") for c in rate.index], rotation=45)
    ax.set_ylabel("Confirmed GRCh38 SVs (%)")
    ax.set_xlabel("Chromosome")
    ax.set_title(
        "B  Chromosome-level confirmation rate",
        loc="left",
        fontweight="bold"
    )

    fig.suptitle(
        f"Chromosomal distribution of {title_prefix.lower()}",
        fontsize=16,
        fontweight="bold",
        y=0.98
    )

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save_figure(fig, "crossref_fig3_chromosome_confirmation", out_dir)


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    args = parse_args()

    samples = [
        s.strip()
        for s in args.samples.split(",")
        if s.strip()
    ]

    os.makedirs(args.out_dir, exist_ok=True)

    if not os.path.exists(args.integrated_vcf):
        raise FileNotFoundError(args.integrated_vcf)

    print("Loading integrated GRCh38 cohort:")
    print(args.integrated_vcf)

    df = load_integrated_grch38(args.integrated_vcf)
    df = add_bins(df)

    metrics = build_metrics(
        df=df,
        workflow_outdir=args.workflow_outdir,
        batch=args.batch,
        samples=samples,
        skip_liftover_panel=args.skip_liftover_panel,
    )

    out_table = os.path.join(
        args.out_dir,
        "crossref_infofield_confirmation_table.tsv"
    )

    df.to_csv(out_table, sep="\t", index=False)

    out_metrics = os.path.join(
        args.out_dir,
        "crossref_infofield_summary_metrics.tsv"
    )

    pd.DataFrame({
        "metric": list(metrics.keys()),
        "value": list(metrics.values())
    }).to_csv(out_metrics, sep="\t", index=False)

    print("\nSummary:")
    for k, v in metrics.items():
        if isinstance(v, float) and not math.isnan(v):
            print(f"{k}: {v:.2f}")
        else:
            print(f"{k}: {v}")

    print("\nSaved:")
    print(out_table)
    print(out_metrics)

    plot_summary(
        df=df,
        metrics=metrics,
        out_dir=args.out_dir,
        title_prefix=args.title_prefix,
    )

    plot_confirmation_patterns(
        df=df,
        out_dir=args.out_dir,
        title_prefix=args.title_prefix,
    )

    plot_chromosome_confirmation(
        df=df,
        out_dir=args.out_dir,
        title_prefix=args.title_prefix,
    )

    print("\nDone.")
    print("Output directory:")
    print(args.out_dir)


if __name__ == "__main__":
    main()