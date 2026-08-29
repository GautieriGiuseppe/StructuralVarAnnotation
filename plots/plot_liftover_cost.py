#!/usr/bin/env python3

"""
Liftover cost analysis for the integrated GRCh38 workflow.

Comparison:
  1. Native GRCh38 per-tool calls
  2. Native CHM13 per-tool calls
  3. Lifted CHM13 to GRCh38 per-tool calls

Retention:
  lifted CHM13 to GRCh38 / native CHM13

This version is pipeline-ready:
  - reads sample_id and batch_id from samples.tsv
  - supports any number of samples, including trio mode
  - writes outputs under the selected workflow outdir
  - does not require bcftools
"""

import os
import gzip
import argparse
import warnings

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter


warnings.filterwarnings("ignore", category=FutureWarning)


# ==============================================================================
# CONSTANTS
# ==============================================================================

TOOLS = ["Sniffles", "Delly", "CuteSV"]

TOOL_FILE = {
    "Sniffles": "sniffles",
    "Delly": "delly",
    "CuteSV": "cuteSV",
}

SV_TYPES = ["DEL", "INS", "INV", "DUP", "SV"]

TYPE_COLORS = {
    "DEL": "#2980B9",
    "INS": "#27AE60",
    "INV": "#E67E22",
    "DUP": "#8E44AD",
    "SV": "#95A5A6",
}

CALLSET_LABELS = [
    "Native\nGRCh38",
    "Native\nCHM13",
    "Lifted\nCHM13→GRCh38",
]

CALLSET_COLORS = {
    "Native\nGRCh38": "#2C5F8A",
    "Native\nCHM13": "#C0392B",
    "Lifted\nCHM13→GRCh38": "#E8735A",
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
            "Analyze liftover cost from native CHM13 to lifted GRCh38 callsets "
            "using per-sample per-tool VCFs."
        )
    )

    parser.add_argument(
        "--samples",
        required=True,
        help=(
            "Samples TSV containing at least sample_id and batch_id columns. "
            "This is the same samples.tsv used by the workflow."
        ),
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help=(
            "Main workflow output directory, for example "
            "/scratch/giuseppe.gautieri/trio_output."
        ),
    )

    parser.add_argument(
        "--out-dir",
        default=None,
        help=(
            "Optional output directory for figures. If not provided, figures are "
            "written to <outdir>/cohort_results/figures_liftover_cost_chm13_to_grch38."
        ),
    )

    parser.add_argument(
        "--title-suffix",
        default=None,
        help=(
            "Optional text appended to figure titles, for example 'trio mode'."
        ),
    )

    return parser.parse_args()


# ==============================================================================
# HELPERS
# ==============================================================================

def fmt_int(x, _):
    return f"{int(x):,}"


def save(fig, out_dir, name):
    os.makedirs(out_dir, exist_ok=True)

    pdf = os.path.join(out_dir, f"{name}.pdf")
    png = os.path.join(out_dir, f"{name}.png")

    fig.savefig(pdf, bbox_inches="tight", dpi=300)
    fig.savefig(png, bbox_inches="tight", dpi=300)
    plt.close(fig)

    print(f"Saved: {pdf}")
    print(f"Saved: {png}")


def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_info_field(info_str):
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


def normalize_svtype(raw_type):
    raw_type = str(raw_type).upper()

    if "DEL" in raw_type:
        return "DEL"
    if "INS" in raw_type:
        return "INS"
    if "INV" in raw_type:
        return "INV"
    if "DUP" in raw_type:
        return "DUP"

    return "SV"


def parse_svlen(raw_len):
    if raw_len in [None, ".", "", True]:
        return np.nan

    try:
        first = str(raw_len).split(",")[0]
        return abs(int(float(first)))
    except Exception:
        return np.nan


def load_sample_pairs(samples_tsv):
    if not os.path.exists(samples_tsv):
        raise FileNotFoundError(f"Samples TSV not found: {samples_tsv}")

    pairs = []

    with open(samples_tsv, "r") as f:
        header = f.readline().rstrip("\n").split("\t")

        required = {"sample_id", "batch_id"}
        missing = required - set(header)

        if missing:
            raise ValueError(
                f"Samples TSV is missing required columns: {missing}\n"
                f"Available columns: {header}"
            )

        sample_idx = header.index("sample_id")
        batch_idx = header.index("batch_id")

        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            sample = fields[sample_idx].strip()
            batch = fields[batch_idx].strip()

            if sample and batch:
                pairs.append((batch, sample))

    if not pairs:
        raise RuntimeError(f"No sample/batch pairs loaded from {samples_tsv}")

    return pairs


def vcf_paths(outdir, batch, sample, tool):
    tool_file = TOOL_FILE[tool]
    root = os.path.join(outdir, batch, sample, "03.variant_calling")

    if tool == "Delly":
        chm13_native = os.path.join(
            root,
            f"{sample}_chm13_delly_symbolic.vcf"
        )
    else:
        chm13_native = os.path.join(
            root,
            f"{sample}_chm13_{tool_file}.vcf"
        )

    return {
        "Native\nGRCh38": os.path.join(
            root,
            f"{sample}_grch38_{tool_file}.vcf"
        ),
        "Native\nCHM13": chm13_native,
        "Lifted\nCHM13→GRCh38": os.path.join(
            root,
            "liftover",
            f"{tool_file}_{sample}_chm13-to-grch38.toolref_pass.vcf"
        ),
    }


def parse_vcf_counts(vcf_path):
    """
    Return one row per PASS variant with SVTYPE and absolute SVLEN.
    """

    if not os.path.exists(vcf_path):
        print(f"[MISSING] {vcf_path}")
        return pd.DataFrame(columns=["SVTYPE", "SVLEN"])

    rows = []

    try:
        with open_vcf(vcf_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue

                cols = line.rstrip("\n").split("\t")

                if len(cols) < 8:
                    continue

                filt = cols[6]
                info_str = cols[7]

                if filt != "PASS":
                    continue

                info = parse_info_field(info_str)

                raw_type = info.get("SVTYPE", ".")
                raw_len = info.get("SVLEN", ".")

                svtype = normalize_svtype(raw_type)
                svlen = parse_svlen(raw_len)

                rows.append({
                    "SVTYPE": svtype,
                    "SVLEN": svlen,
                })

    except Exception as e:
        print(f"[WARNING] failed to parse {vcf_path}: {e}")
        return pd.DataFrame(columns=["SVTYPE", "SVLEN"])

    if not rows:
        return pd.DataFrame(columns=["SVTYPE", "SVLEN"])

    return pd.DataFrame(rows)


def title_with_suffix(title, suffix):
    if suffix:
        return f"{title}\n({suffix})"
    return title


# ==============================================================================
# LOAD DATA
# ==============================================================================

def load_liftover_data(outdir, sample_pairs):
    print("Loading callsets...")

    summary_records = []
    size_records = []

    for batch, sample in sample_pairs:
        for tool in TOOLS:
            paths = vcf_paths(
                outdir=outdir,
                batch=batch,
                sample=sample,
                tool=tool,
            )

            for callset, path in paths.items():
                df = parse_vcf_counts(path)

                if df.empty:
                    print(f"[SKIP] {batch}/{sample} {tool} {callset.replace(chr(10), ' ')}")
                    continue

                print(
                    f"{batch}/{sample} {tool} "
                    f"{callset.replace(chr(10), ' ')}: {len(df):,}"
                )

                type_counts = df["SVTYPE"].value_counts().to_dict()

                row = {
                    "Batch": batch,
                    "Sample": sample,
                    "Tool": tool,
                    "Callset": callset,
                    "Total": int(len(df)),
                }

                for svtype in SV_TYPES:
                    row[f"n_{svtype}"] = int(type_counts.get(svtype, 0))

                summary_records.append(row)

                tmp = df.copy()
                tmp["Batch"] = batch
                tmp["Sample"] = sample
                tmp["Tool"] = tool
                tmp["Callset"] = callset
                size_records.append(tmp)

    master_df = pd.DataFrame(summary_records)

    if master_df.empty:
        raise RuntimeError(
            "No VCF records were loaded. Check file paths and whether VCFs contain PASS variants."
        )

    size_df = (
        pd.concat(size_records, ignore_index=True)
        if size_records
        else pd.DataFrame()
    )

    return master_df, size_df


# ==============================================================================
# TABLES
# ==============================================================================

def build_retention_tables(master_df):
    retention_rows = []

    for sample in sorted(master_df["Sample"].unique()):
        for tool in TOOLS:
            sub = master_df[
                (master_df["Sample"] == sample) &
                (master_df["Tool"] == tool)
            ]

            native = sub[sub["Callset"] == "Native\nCHM13"]["Total"].values
            lifted = sub[sub["Callset"] == "Lifted\nCHM13→GRCh38"]["Total"].values

            if len(native) and len(lifted) and native[0] > 0:
                retention_rows.append({
                    "Sample": sample,
                    "Tool": tool,
                    "Native_CHM13": int(native[0]),
                    "Lifted_CHM13_to_GRCh38": int(lifted[0]),
                    "Lost": int(native[0] - lifted[0]),
                    "Retention": lifted[0] / native[0],
                })

    retention_df = pd.DataFrame(retention_rows)

    type_retention_rows = []

    for sample in sorted(master_df["Sample"].unique()):
        for tool in TOOLS:
            sub = master_df[
                (master_df["Sample"] == sample) &
                (master_df["Tool"] == tool)
            ]

            native = sub[sub["Callset"] == "Native\nCHM13"]
            lifted = sub[sub["Callset"] == "Lifted\nCHM13→GRCh38"]

            if native.empty or lifted.empty:
                continue

            native = native.iloc[0]
            lifted = lifted.iloc[0]

            for svtype in SV_TYPES:
                n_native = native[f"n_{svtype}"]
                n_lifted = lifted[f"n_{svtype}"]

                if n_native > 0:
                    type_retention_rows.append({
                        "Sample": sample,
                        "Tool": tool,
                        "SVTYPE": svtype,
                        "Native_CHM13": int(n_native),
                        "Lifted_CHM13_to_GRCh38": int(n_lifted),
                        "Lost": int(n_native - n_lifted),
                        "Retention": n_lifted / n_native,
                    })

    type_retention_df = pd.DataFrame(type_retention_rows)

    return retention_df, type_retention_df


def write_tables(master_df, size_df, retention_df, type_retention_df, cohort_out):
    os.makedirs(cohort_out, exist_ok=True)

    paths = {
        "counts": os.path.join(
            cohort_out,
            "liftover_cost_counts_chm13_to_grch38.tsv",
        ),
        "sizes": os.path.join(
            cohort_out,
            "liftover_cost_variant_sizes_chm13_to_grch38.tsv",
        ),
        "retention": os.path.join(
            cohort_out,
            "liftover_retention_chm13_to_grch38.tsv",
        ),
        "retention_by_svtype": os.path.join(
            cohort_out,
            "liftover_retention_by_svtype_chm13_to_grch38.tsv",
        ),
    }

    master_df.to_csv(paths["counts"], sep="\t", index=False)
    size_df.to_csv(paths["sizes"], sep="\t", index=False)
    retention_df.to_csv(paths["retention"], sep="\t", index=False)
    type_retention_df.to_csv(paths["retention_by_svtype"], sep="\t", index=False)

    print("\nSaved raw tables:")
    for path in paths.values():
        print(path)


# ==============================================================================
# PLOTS
# ==============================================================================

def plot_fig1_cohort_total_counts(master_df, out_dir, title_suffix):
    cohort_counts = (
        master_df
        .groupby("Callset")["Total"]
        .sum()
        .reindex(CALLSET_LABELS)
        .fillna(0)
    )

    fig, ax = plt.subplots(figsize=(8.5, 6))

    x = np.arange(len(CALLSET_LABELS))
    colors = [CALLSET_COLORS[c] for c in CALLSET_LABELS]

    bars = ax.bar(
        x,
        cohort_counts.values,
        color=colors,
        edgecolor="white",
        linewidth=0.8,
        alpha=0.92,
    )

    ymax = max(cohort_counts.max(), 1)

    for bar, val in zip(bars, cohort_counts.values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + ymax * 0.015,
            f"{int(val):,}",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=10,
        )

    native_chm13 = cohort_counts["Native\nCHM13"]
    lifted = cohort_counts["Lifted\nCHM13→GRCh38"]

    if native_chm13 > 0:
        pct_lost = (native_chm13 - lifted) / native_chm13 * 100

        ax.annotate(
            "",
            xy=(2, lifted),
            xytext=(1, native_chm13),
            arrowprops=dict(arrowstyle="->", color="#E74C3C", lw=1.6),
        )

        ax.text(
            1.5,
            max(native_chm13, lifted) * 1.06,
            f"-{pct_lost:.1f}%",
            ha="center",
            color="#E74C3C",
            fontweight="bold",
            fontsize=10,
        )

    ax.set_title(
        title_with_suffix(
            "Variant count: native GRCh38 vs native CHM13 vs lifted CHM13→GRCh38",
            title_suffix,
        ),
        fontweight="bold",
        fontsize=13,
    )

    ax.set_ylabel("Variant count")
    ax.set_xticks(x)
    ax.set_xticklabels(CALLSET_LABELS)
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_int))
    ax.grid(axis="y", alpha=0.3)

    handles = [
        mpatches.Patch(color=CALLSET_COLORS[c], label=c.replace("\n", " "))
        for c in CALLSET_LABELS
    ]

    ax.legend(
        handles=handles,
        frameon=False,
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.5, -0.16),
    )

    save(fig, out_dir, "fig1_liftover_cohort_total_counts_chm13_to_grch38")


def plot_fig2_cohort_svtype_composition(master_df, out_dir, title_suffix):
    composition = (
        master_df
        .groupby("Callset")[[f"n_{t}" for t in SV_TYPES]]
        .sum()
        .reindex(CALLSET_LABELS)
        .fillna(0)
    )

    fig, ax = plt.subplots(figsize=(8.5, 6))

    x = np.arange(len(CALLSET_LABELS))
    bottom = np.zeros(len(CALLSET_LABELS))
    max_total = max(composition.sum(axis=1).max(), 1)

    for svtype in SV_TYPES:
        vals = composition[f"n_{svtype}"].values

        ax.bar(
            x,
            vals,
            bottom=bottom,
            color=TYPE_COLORS[svtype],
            edgecolor="white",
            linewidth=0.5,
            label=svtype,
            alpha=0.92,
        )

        for i, val in enumerate(vals):
            if val > max_total * 0.04:
                ax.text(
                    i,
                    bottom[i] + val / 2,
                    f"{int(val):,}",
                    ha="center",
                    va="center",
                    color="white",
                    fontweight="bold",
                    fontsize=9,
                )

        bottom += vals

    ax.set_title(
        title_with_suffix(
            "SV type composition: native GRCh38 vs native CHM13 vs lifted CHM13→GRCh38",
            title_suffix,
        ),
        fontweight="bold",
        fontsize=13,
    )

    ax.set_ylabel("Variant count")
    ax.set_xticks(x)
    ax.set_xticklabels(CALLSET_LABELS)
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_int))
    ax.grid(axis="y", alpha=0.3)

    ax.legend(
        title="SV type",
        frameon=False,
        ncol=5,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.18),
    )

    save(fig, out_dir, "fig2_liftover_cohort_svtype_composition_chm13_to_grch38")


def plot_fig3_retention_by_svtype(type_retention_df, out_dir, title_suffix):
    if type_retention_df.empty:
        print("[SKIP] fig3 retention by SV type: no retention data.")
        return

    retention_by_type = (
        type_retention_df
        .groupby("SVTYPE")
        .agg(
            Native_CHM13=("Native_CHM13", "sum"),
            Lifted_CHM13_to_GRCh38=("Lifted_CHM13_to_GRCh38", "sum"),
        )
        .reset_index()
    )

    retention_by_type = retention_by_type[
        retention_by_type["SVTYPE"].isin(SV_TYPES)
    ].copy()

    if retention_by_type.empty:
        print("[SKIP] fig3 retention by SV type: no SV types available.")
        return

    retention_by_type["Retention"] = (
        retention_by_type["Lifted_CHM13_to_GRCh38"] /
        retention_by_type["Native_CHM13"]
    )

    retention_by_type["SVTYPE"] = pd.Categorical(
        retention_by_type["SVTYPE"],
        categories=SV_TYPES,
        ordered=True,
    )

    retention_by_type = retention_by_type.sort_values("SVTYPE")

    fig, ax = plt.subplots(figsize=(8.5, 6))

    x_type = np.arange(len(retention_by_type))
    vals = retention_by_type["Retention"].values * 100
    types = retention_by_type["SVTYPE"].astype(str).values
    colors = [TYPE_COLORS[t] for t in types]

    bars = ax.bar(
        x_type,
        vals,
        color=colors,
        edgecolor="white",
        linewidth=0.8,
        alpha=0.92,
    )

    ymax = max(vals.max(), 1)

    for bar, val in zip(bars, vals):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + ymax * 0.02,
            f"{val:.1f}%",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=10,
        )

    ax.axhline(100, color="grey", linestyle="--", linewidth=1.0, alpha=0.6)
    ax.set_ylim(0, max(110, ymax * 1.15))
    ax.set_xticks(x_type)
    ax.set_xticklabels(types)
    ax.set_ylabel("Retention rate (%)")

    ax.set_title(
        title_with_suffix(
            "Liftover retention rate by SV type\n"
            "(% of native CHM13 variants surviving liftover to GRCh38)",
            title_suffix,
        ),
        fontweight="bold",
        fontsize=13,
    )

    ax.grid(axis="y", alpha=0.3)

    save(fig, out_dir, "fig3_liftover_retention_by_svtype_chm13_to_grch38")


def plot_fig4_raw_counts_per_tool(master_df, out_dir, n_samples, title_suffix):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 5), sharey=False)

    for ax, tool in zip(axes, TOOLS):
        sub = master_df[master_df["Tool"] == tool]

        means = (
            sub.groupby("Callset")["Total"]
            .mean()
            .reindex(CALLSET_LABELS)
            .fillna(0)
        )

        stds = (
            sub.groupby("Callset")["Total"]
            .std()
            .reindex(CALLSET_LABELS)
            .fillna(0)
        )

        x_tool = np.arange(len(CALLSET_LABELS))

        bars = ax.bar(
            x_tool,
            means.values,
            yerr=stds.values,
            capsize=3,
            color=[CALLSET_COLORS[c] for c in CALLSET_LABELS],
            edgecolor="white",
            linewidth=0.6,
            alpha=0.92,
            error_kw={"linewidth": 0.8},
        )

        ymax = max(means.max(), 1)

        for i, (bar, mean_val) in enumerate(zip(bars, means.values)):
            yerr = stds.values[i]

            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + yerr + ymax * 0.02,
                f"{int(mean_val):,}",
                ha="center",
                fontsize=8,
                fontweight="bold",
            )

        native_chm13_mean = means["Native\nCHM13"]
        lifted_mean = means["Lifted\nCHM13→GRCh38"]

        if native_chm13_mean > 0:
            pct_lost = (
                (native_chm13_mean - lifted_mean) /
                native_chm13_mean * 100
            )

            ax.annotate(
                "",
                xy=(2, lifted_mean),
                xytext=(1, native_chm13_mean),
                arrowprops=dict(arrowstyle="->", color="#E74C3C", lw=1.2),
            )

            ax.text(
                1.5,
                max(native_chm13_mean, lifted_mean) * 1.08,
                f"-{pct_lost:.1f}%\nlost",
                ha="center",
                fontsize=8,
                color="#E74C3C",
                fontweight="bold",
            )

        ax.set_title(tool, fontweight="bold")
        ax.set_xticks(x_tool)
        ax.set_xticklabels(CALLSET_LABELS)
        ax.set_ylabel("Mean variant count ± SD" if ax == axes[0] else "")
        ax.yaxis.set_major_formatter(FuncFormatter(fmt_int))
        ax.grid(axis="y", alpha=0.3)

    fig.suptitle(
        title_with_suffix(
            f"Per-tool variant count: native GRCh38 vs native CHM13 vs lifted CHM13→GRCh38\n"
            f"(mean ± SD across {n_samples} samples)",
            title_suffix,
        ),
        fontweight="bold",
        fontsize=13,
    )

    handles = [
        mpatches.Patch(color=CALLSET_COLORS[c], label=c.replace("\n", " "))
        for c in CALLSET_LABELS
    ]

    fig.legend(
        handles=handles,
        frameon=False,
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.5, -0.04),
    )

    plt.tight_layout()

    save(fig, out_dir, "fig4_liftover_raw_counts_per_tool_chm13_to_grch38")


def plot_fig5_svtype_breakdown_per_tool(master_df, out_dir, n_samples, title_suffix):
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 5), sharey=False)

    for ax, tool in zip(axes, TOOLS):
        sub = master_df[master_df["Tool"] == tool]

        type_means = {
            svtype: (
                sub.groupby("Callset")[f"n_{svtype}"]
                .mean()
                .reindex(CALLSET_LABELS)
                .fillna(0)
            )
            for svtype in SV_TYPES
        }

        x_tool = np.arange(len(CALLSET_LABELS))
        bottom = np.zeros(len(CALLSET_LABELS))

        for svtype in SV_TYPES:
            vals = type_means[svtype].values

            if vals.sum() == 0:
                continue

            ax.bar(
                x_tool,
                vals,
                bottom=bottom,
                color=TYPE_COLORS[svtype],
                edgecolor="white",
                linewidth=0.4,
                label=svtype,
                alpha=0.92,
            )

            bottom += vals

        ax.set_title(tool, fontweight="bold")
        ax.set_xticks(x_tool)
        ax.set_xticklabels(CALLSET_LABELS)
        ax.set_ylabel("Mean variant count" if ax == axes[0] else "")
        ax.yaxis.set_major_formatter(FuncFormatter(fmt_int))
        ax.grid(axis="y", alpha=0.3)

    handles = [mpatches.Patch(color=TYPE_COLORS[t], label=t) for t in SV_TYPES]

    fig.legend(
        handles=handles,
        title="SV type",
        frameon=False,
        ncol=5,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.04),
    )

    fig.suptitle(
        title_with_suffix(
            f"Per-tool SV type composition: native GRCh38 vs native CHM13 vs lifted CHM13→GRCh38\n"
            f"(mean across {n_samples} samples)",
            title_suffix,
        ),
        fontweight="bold",
        fontsize=13,
    )

    plt.tight_layout()

    save(fig, out_dir, "fig5_liftover_svtype_breakdown_per_tool_chm13_to_grch38")


def plot_fig6_size_distribution(size_df, out_dir, title_suffix):
    if size_df.empty:
        print("[SKIP] fig6 size distribution: size_df is empty.")
        return

    size_plot = size_df[
        (size_df["SVLEN"] > 30) &
        (size_df["SVLEN"] < 1e6) &
        (size_df["SVTYPE"].isin(["DEL", "INS"]))
    ].copy()

    if size_plot.empty:
        print("[SKIP] fig6 size distribution: no DEL/INS size records available.")
        return

    size_plot["log10_SVLEN"] = np.log10(size_plot["SVLEN"])

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 5), sharey=True)

    for ax, tool in zip(axes, TOOLS):
        sub_tool = size_plot[size_plot["Tool"] == tool].copy()

        data = [
            sub_tool[sub_tool["Callset"] == callset]["log10_SVLEN"].dropna().values
            for callset in CALLSET_LABELS
        ]

        if all(len(d) == 0 for d in data):
            ax.set_title(tool, fontweight="bold")
            ax.text(
                0.5,
                0.5,
                "No size data",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.axis("off")
            continue

        parts = ax.violinplot(
            data,
            positions=np.arange(len(CALLSET_LABELS)),
            showmeans=False,
            showmedians=True,
            showextrema=False,
        )

        for body, callset in zip(parts["bodies"], CALLSET_LABELS):
            body.set_facecolor(CALLSET_COLORS[callset])
            body.set_edgecolor("black")
            body.set_alpha(0.65)

        if "cmedians" in parts:
            parts["cmedians"].set_color("black")
            parts["cmedians"].set_linewidth(1.0)

        ax.set_title(tool, fontweight="bold")
        ax.set_xticks(np.arange(len(CALLSET_LABELS)))
        ax.set_xticklabels(CALLSET_LABELS)
        ax.set_ylabel("log10(SVLEN bp)" if ax == axes[0] else "")
        ax.grid(axis="y", alpha=0.3)

        for log_size, label in [
            (np.log10(50), "50 bp"),
            (np.log10(1000), "1 kb"),
            (np.log10(10000), "10 kb"),
        ]:
            ax.axhline(
                log_size,
                color="grey",
                linestyle=":",
                linewidth=0.7,
                alpha=0.6,
            )

            if ax == axes[0]:
                ax.text(
                    -0.45,
                    log_size,
                    label,
                    fontsize=7,
                    va="center",
                    color="grey",
                )

    fig.suptitle(
        title_with_suffix(
            "SV size distribution: native GRCh38 vs native CHM13 vs lifted CHM13→GRCh38\n"
            "(DEL and INS only, log scale)",
            title_suffix,
        ),
        fontweight="bold",
        fontsize=13,
    )

    plt.tight_layout()

    save(fig, out_dir, "fig6_liftover_size_distribution_chm13_to_grch38")


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    args = parse_args()

    outdir = args.outdir.rstrip("/")
    cohort_out = os.path.join(outdir, "cohort_results")

    out_dir = (
        args.out_dir
        if args.out_dir is not None
        else os.path.join(
            cohort_out,
            "figures_liftover_cost_chm13_to_grch38",
        )
    )

    sample_pairs = load_sample_pairs(args.samples)
    n_samples = len(sample_pairs)

    print("Samples TSV:", args.samples)
    print("Workflow outdir:", outdir)
    print("Cohort output:", cohort_out)
    print("Figure output:", out_dir)
    print(f"Loaded {n_samples} sample/batch pairs:")
    for batch, sample in sample_pairs:
        print(f"  {batch}\t{sample}")

    master_df, size_df = load_liftover_data(
        outdir=outdir,
        sample_pairs=sample_pairs,
    )

    retention_df, type_retention_df = build_retention_tables(master_df)

    write_tables(
        master_df=master_df,
        size_df=size_df,
        retention_df=retention_df,
        type_retention_df=type_retention_df,
        cohort_out=cohort_out,
    )

    print("\nMean retention by tool:")
    if not retention_df.empty:
        print(
            retention_df
            .groupby("Tool")[["Retention", "Lost"]]
            .mean()
            .round(3)
            .to_string()
        )
    else:
        print("No retention data available.")

    print("\nMean retention by SV type:")
    if not type_retention_df.empty:
        print(
            type_retention_df
            .groupby("SVTYPE")[["Retention", "Lost"]]
            .mean()
            .round(3)
            .to_string()
        )
    else:
        print("No SVTYPE retention data available.")

    plot_fig1_cohort_total_counts(
        master_df=master_df,
        out_dir=out_dir,
        title_suffix=args.title_suffix,
    )

    plot_fig2_cohort_svtype_composition(
        master_df=master_df,
        out_dir=out_dir,
        title_suffix=args.title_suffix,
    )

    plot_fig3_retention_by_svtype(
        type_retention_df=type_retention_df,
        out_dir=out_dir,
        title_suffix=args.title_suffix,
    )

    plot_fig4_raw_counts_per_tool(
        master_df=master_df,
        out_dir=out_dir,
        n_samples=n_samples,
        title_suffix=args.title_suffix,
    )

    plot_fig5_svtype_breakdown_per_tool(
        master_df=master_df,
        out_dir=out_dir,
        n_samples=n_samples,
        title_suffix=args.title_suffix,
    )

    plot_fig6_size_distribution(
        size_df=size_df,
        out_dir=out_dir,
        title_suffix=args.title_suffix,
    )

    print("\n" + "=" * 70)
    print("DONE")
    print("Figures saved to:", out_dir)
    print("Tables saved to:", cohort_out)
    print("=" * 70)


if __name__ == "__main__":
    main()