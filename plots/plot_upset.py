import os
import gzip
import argparse

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors


# =============================================================================
# DEFAULTS
# =============================================================================

DEFAULT_VCF = (
    "/group/dominguez/shared_notebooks/Immune_variation/mapping/cohort_results/"
    "GRCh38_final_cohort_survivor.vcf.gz"
)

DEFAULT_OUT_DIR = (
    "/group/dominguez/shared_notebooks/Immune_variation/mapping/cohort_results/"
    "tool_reference_upset_new_cohort"
)

DEFAULT_OUT_PREFIX = "GRCh38_integrated_toolref_upset_mean_sample_frequency"

DEFAULT_TOP_N = 40

# Actual TOOLREF_SUPP_VEC order from cohort_merge_grch38.smk:
# 0 Sniffles native GRCh38
# 1 Delly native GRCh38
# 2 CuteSV native GRCh38
# 3 Sniffles CHM13 lifted to GRCh38
# 4 Delly CHM13 lifted to GRCh38
# 5 CuteSV CHM13 lifted to GRCh38
VECTOR_ORDER = [
    "Sniffles\nGR38",
    "Delly\nGR38",
    "CuteSV\nGR38",
    "Sniffles\nCHM13→GR38",
    "Delly\nCHM13→GR38",
    "CuteSV\nCHM13→GR38",
]

# Display order in the plot, interleaved by tool
DISPLAY_ORDER = [0, 3, 1, 4, 2, 5]
SET_LABELS = [VECTOR_ORDER[i] for i in DISPLAY_ORDER]
N_SETS = len(SET_LABELS)

SAMPLE_SUPPORT_MIN = 1
SAMPLE_SUPPORT_MAX = 8


# =============================================================================
# ARGPARSE
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate an UpSet-style plot for the integrated GRCh38 "
            "tool/reference support vectors from GRCh38_final_cohort_survivor.vcf.gz."
        )
    )

    parser.add_argument(
        "--vcf",
        default=DEFAULT_VCF,
        help="Input integrated GRCh38 cohort VCF.gz containing TOOLREF_SUPP_VEC."
    )

    parser.add_argument(
        "--out-dir",
        default=DEFAULT_OUT_DIR,
        help="Output directory for plots and TSV summaries."
    )

    parser.add_argument(
        "--out-prefix",
        default=DEFAULT_OUT_PREFIX,
        help="Output file prefix without extension."
    )

    parser.add_argument(
        "--top-n",
        type=int,
        default=DEFAULT_TOP_N,
        help="Number of largest intersections to plot."
    )

    parser.add_argument(
        "--cmap",
        default="viridis",
        help="Matplotlib colormap used for mean sample support."
    )

    parser.add_argument(
        "--sample-support-min",
        type=float,
        default=SAMPLE_SUPPORT_MIN,
        help="Minimum value for color scale."
    )

    parser.add_argument(
        "--sample-support-max",
        type=float,
        default=SAMPLE_SUPPORT_MAX,
        help="Maximum value for color scale."
    )

    parser.add_argument(
        "--title",
        default=(
            "Tool Intersection UpSet — Bar Color = Mean Sample Frequency\n"
            "Integrated GRCh38 cohort: native GRCh38 + CHM13 calls lifted to GRCh38"
        ),
        help="Plot title."
    )

    return parser.parse_args()


# =============================================================================
# HELPERS
# =============================================================================

def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


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


def safe_int(x):
    try:
        if x in [None, ".", ""]:
            return np.nan
        return int(float(x))
    except Exception:
        return np.nan


def color_for_mean_sample_supp(mean_sample_supp, cmap_name, vmin, vmax):
    if pd.isna(mean_sample_supp):
        return "#BDBDBD"

    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.get_cmap(cmap_name)

    return cmap(norm(float(mean_sample_supp)))


def reorder_vec(vec):
    return "".join(vec[i] for i in DISPLAY_ORDER)


# =============================================================================
# LOAD
# =============================================================================

def load_toolref_patterns(vcf):
    rows = []

    with open_vcf(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, vid, ref, alt, qual, filt, info_str = cols[:8]
            info = parse_info(info_str)

            vec = info.get("TOOLREF_SUPP_VEC")
            if vec is None:
                continue

            vec = str(vec).strip()
            if len(vec) != 6:
                continue

            rows.append({
                "CHROM": chrom,
                "POS": int(pos),
                "ID": vid,
                "SVTYPE": info.get("SVTYPE", "."),
                "SVLEN": info.get("SVLEN", "."),
                "TOOLREF_SUPP_VEC_ORIGINAL": vec,
                "TOOLREF_SUPP_VEC": reorder_vec(vec),
                "TOOLREF_SUPP": safe_int(info.get("TOOLREF_SUPP")),
                "SAMPLE_SUPP": safe_int(info.get("SAMPLE_SUPP")),
                "NATIVE_GRCH38_SUPP": safe_int(info.get("NATIVE_GRCH38_SUPP")),
                "LIFTED_CHM13_GRCH38_SUPP": safe_int(info.get("LIFTED_CHM13_GRCH38_SUPP")),
            })

    return pd.DataFrame(rows)


# =============================================================================
# SUMMARIZE
# =============================================================================

def summarize_intersections(df, top_n, cmap_name, vmin, vmax):
    grouped = (
        df.groupby("TOOLREF_SUPP_VEC")
        .agg(
            intersection_size=("TOOLREF_SUPP_VEC", "size"),
            mean_sample_supp=("SAMPLE_SUPP", "mean"),
            median_sample_supp=("SAMPLE_SUPP", "median"),
            mean_toolref_supp=("TOOLREF_SUPP", "mean"),
            mean_native_grch38_supp=("NATIVE_GRCH38_SUPP", "mean"),
            mean_lifted_chm13_supp=("LIFTED_CHM13_GRCH38_SUPP", "mean"),
        )
        .reset_index()
    )

    grouped = grouped.sort_values(
        ["intersection_size", "mean_sample_supp"],
        ascending=[False, False]
    ).head(top_n).reset_index(drop=True)

    grouped["bar_color"] = grouped["mean_sample_supp"].apply(
        lambda x: color_for_mean_sample_supp(x, cmap_name, vmin, vmax)
    )

    set_sizes = []

    for i, label in enumerate(SET_LABELS):
        set_sizes.append({
            "set": label,
            "size": int(df["TOOLREF_SUPP_VEC"].str[i].eq("1").sum())
        })

    return grouped, pd.DataFrame(set_sizes)


# =============================================================================
# PLOT
# =============================================================================

def plot_upset(
    intersections,
    set_sizes,
    out_prefix,
    cmap_name,
    vmin,
    vmax,
    title,
):
    n = len(intersections)
    x = np.arange(n)

    fig = plt.figure(figsize=(24, 10))

    gs = fig.add_gridspec(
        nrows=2,
        ncols=2,
        width_ratios=[1.35, 8.65],
        height_ratios=[3.4, 2.3],
        wspace=0.12,
        hspace=0.08,
    )

    ax_empty = fig.add_subplot(gs[0, 0])
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_sets = fig.add_subplot(gs[1, 0])
    ax_matrix = fig.add_subplot(gs[1, 1])

    ax_empty.axis("off")

    bars = ax_bar.bar(
        x,
        intersections["intersection_size"],
        color=intersections["bar_color"],
        width=0.75,
        edgecolor="none",
        align="center"
    )

    ymax = intersections["intersection_size"].max() * 1.18
    ax_bar.set_ylim(0, ymax)

    left_pad = 0.85
    right_pad = 0.55
    ax_bar.set_xlim(-left_pad, n - right_pad)

    for i, bar in enumerate(bars):
        val = intersections.loc[i, "intersection_size"]
        ax_bar.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + ymax * 0.01,
            f"{int(val)}",
            ha="center",
            va="bottom",
            fontsize=8,
            clip_on=False
        )

    ax_bar.set_ylabel("Intersection size", fontsize=11)
    ax_bar.set_xticks([])
    ax_bar.grid(axis="y", color="#D0D0D0", linewidth=0.8)
    ax_bar.set_axisbelow(True)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)

    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    sm = cm.ScalarMappable(
        cmap=cm.get_cmap(cmap_name),
        norm=norm
    )
    sm.set_array([])

    cbar = fig.colorbar(
        sm,
        ax=ax_bar,
        fraction=0.025,
        pad=0.015
    )

    cbar.set_label("Mean sample support", fontsize=10)

    if vmin == 1 and vmax == 8:
        cbar.set_ticks(range(1, 9))

    y = np.arange(N_SETS)

    set_sizes = set_sizes.copy()
    set_sizes["set"] = pd.Categorical(
        set_sizes["set"],
        categories=SET_LABELS,
        ordered=True
    )
    set_sizes = set_sizes.sort_values("set")

    ax_sets.barh(y, set_sizes["size"], color="black", height=0.62)

    ax_sets.set_yticks(y)
    ax_sets.set_yticklabels([""] * N_SETS)
    ax_sets.invert_xaxis()
    ax_sets.invert_yaxis()
    ax_sets.grid(axis="x", color="#D0D0D0", linewidth=0.8)
    ax_sets.set_axisbelow(True)
    ax_sets.spines["top"].set_visible(False)
    ax_sets.spines["right"].set_visible(False)
    ax_sets.spines["left"].set_visible(False)

    xmax_sets = set_sizes["size"].max()
    ax_sets.set_xlim(xmax_sets * 1.15, 0)

    for yi, val in zip(y, set_sizes["size"]):
        ax_sets.text(
            val + xmax_sets * 0.03,
            yi,
            f"{int(val)}",
            ha="right",
            va="center",
            fontsize=9
        )

    ax_matrix.set_xlim(-left_pad, n - right_pad)
    ax_matrix.set_ylim(-0.5, N_SETS - 0.5)
    ax_matrix.invert_yaxis()

    for yi in range(N_SETS):
        if yi % 2 == 1:
            ax_matrix.axhspan(
                yi - 0.5,
                yi + 0.5,
                color="#F2F2F2",
                zorder=0
            )

    for i in range(n):
        ax_matrix.scatter(
            [i] * N_SETS,
            np.arange(N_SETS),
            s=95,
            color="#D0D0D0",
            zorder=1
        )

    for i, vec in enumerate(intersections["TOOLREF_SUPP_VEC"]):
        active = [j for j, v in enumerate(vec) if v == "1"]

        if len(active) > 1:
            ax_matrix.plot(
                [i, i],
                [min(active), max(active)],
                color="black",
                linewidth=1.6,
                zorder=2
            )

        if active:
            ax_matrix.scatter(
                [i] * len(active),
                active,
                s=120,
                color="black",
                zorder=3
            )

    ax_matrix.set_yticks(np.arange(N_SETS))
    ax_matrix.set_yticklabels(SET_LABELS, fontsize=10)
    ax_matrix.set_xticks([])
    ax_matrix.spines["top"].set_visible(False)
    ax_matrix.spines["right"].set_visible(False)
    ax_matrix.spines["bottom"].set_visible(False)
    ax_matrix.spines["left"].set_visible(False)

    fig.suptitle(
        title,
        fontsize=15,
        fontweight="bold",
        y=0.98
    )

    fig.savefig(out_prefix + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(out_prefix + ".pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================

def main():
    args = parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    df = load_toolref_patterns(args.vcf)

    if df.empty:
        raise RuntimeError(
            "No variants with TOOLREF_SUPP_VEC found. "
            "Check that the input VCF is the annotated integrated GRCh38 cohort."
        )

    intersections, set_sizes = summarize_intersections(
        df=df,
        top_n=args.top_n,
        cmap_name=args.cmap,
        vmin=args.sample_support_min,
        vmax=args.sample_support_max,
    )

    df.to_csv(
        os.path.join(args.out_dir, "toolref_variant_support_table.tsv"),
        sep="\t",
        index=False
    )

    intersections.to_csv(
        os.path.join(args.out_dir, "toolref_upset_intersections.tsv"),
        sep="\t",
        index=False
    )

    set_sizes.to_csv(
        os.path.join(args.out_dir, "toolref_upset_set_sizes.tsv"),
        sep="\t",
        index=False
    )

    out_prefix = os.path.join(args.out_dir, args.out_prefix)

    plot_upset(
        intersections=intersections,
        set_sizes=set_sizes,
        out_prefix=out_prefix,
        cmap_name=args.cmap,
        vmin=args.sample_support_min,
        vmax=args.sample_support_max,
        title=args.title,
    )

    print("Saved:")
    print(out_prefix + ".png")
    print(out_prefix + ".pdf")
    print()
    print("Total variants with TOOLREF_SUPP_VEC:", len(df))
    print()
    print("Top intersections:")
    print(intersections.head(15).to_string(index=False))
    print()
    print("Set sizes:")
    print(set_sizes.to_string(index=False))


if __name__ == "__main__":
    main()