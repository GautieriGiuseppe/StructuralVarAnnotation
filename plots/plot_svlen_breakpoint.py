import argparse
import os
import re
import gzip
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


# =============================================================================
# CONFIG
# =============================================================================

BASE_DIR = "/group/dominguez/shared_notebooks/Immune_variation/mapping/cohort_results"

VCF = os.path.join(
    BASE_DIR,
    "GRCh38_final_cohort_survivor.vcf.gz"
)

METADATA = os.path.join(
    BASE_DIR,
    "GRCh38_toolref_48way_vcf_metadata.tsv"
)

OUT_DIR = os.path.join(
    BASE_DIR,
    "crossref_confirmation_from_integrated_GRCh38",
    "native_vs_lifted_support_distance"
)

OUT_PREFIX = "native_grch38_vs_lifted_chm13_support_distance"

CANONICAL_SVTYPES = {"DEL", "INS", "DUP", "INV"}


# =============================================================================
# STYLE
# =============================================================================

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12,
    "axes.titlesize": 17,
    "axes.labelsize": 13,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.9,
})


def fmt_int(x, _):
    return f"{int(x):,}"


def style_ax(ax):
    ax.grid(axis="both", color="#E6E6E6", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines["left"].set_color("#CCCCCC")
    ax.spines["bottom"].set_color("#CCCCCC")
    ax.tick_params(colors="#444444")
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_int))


# =============================================================================
# VCF HELPERS
# =============================================================================

def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_info(info_str):
    info = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = True
    return info


def safe_int(x, default=None):
    try:
        if x in [None, "", "."]:
            return default
        return int(float(str(x).split(",")[0]))
    except Exception:
        return default


def parse_format_sample(fmt_keys, sample_value):
    vals = sample_value.split(":")
    d = {}

    for i, key in enumerate(fmt_keys):
        d[key] = vals[i] if i < len(vals) else "."

    return d


def is_supported_sample(sample_dict):
    gt = sample_dict.get("GT", ".")
    co = sample_dict.get("CO", ".")
    ln = sample_dict.get("LN", ".")

    if gt not in [".", "./.", ".|."]:
        return True

    if co not in [".", "NA", "NAN", "NaN", "nan", ""]:
        return True

    if ln not in [".", "0", "NA", "NAN", "NaN", "nan", ""]:
        return True

    return False


def parse_ln_values(ln_value):
    if ln_value in [None, "", ".", "NA", "NAN", "NaN", "nan"]:
        return []

    out = []
    for x in str(ln_value).split(","):
        try:
            val = abs(int(float(x)))
            if val > 0:
                out.append(val)
        except Exception:
            continue

    return out


def parse_co_positions(co_value, fallback_pos=None):
    """
    CO examples:
      chr1_10279-chr1_10279
      chr1_10801-chr1_10801,chr1_10801-chr1_10801

    Returns first coordinate from each CO token.
    """
    positions = []

    if co_value in [None, "", ".", "NA", "NAN", "NaN", "nan"]:
        if fallback_pos is not None:
            return [fallback_pos]
        return []

    for token in str(co_value).split(","):
        token = token.strip()

        m = re.search(r"[^_]+_(\d+)-", token)
        if m:
            positions.append(int(m.group(1)))
            continue

        m = re.search(r":(\d+)", token)
        if m:
            positions.append(int(m.group(1)))
            continue

    if not positions and fallback_pos is not None:
        positions = [fallback_pos]

    return positions


def min_abs_pairwise_diff(a, b):
    if not a or not b:
        return np.nan

    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]

    if len(a) == 0 or len(b) == 0:
        return np.nan

    a.sort()
    b.sort()

    i = 0
    j = 0
    best = np.inf

    while i < len(a) and j < len(b):
        diff = abs(a[i] - b[j])
        if diff < best:
            best = diff

        if a[i] < b[j]:
            i += 1
        else:
            j += 1

    return float(best)


# =============================================================================
# LOAD METADATA
# =============================================================================

def load_metadata(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Metadata file not found: {path}")

    md = pd.read_csv(path, sep="\t")

    required = {"vector_index", "set_index", "set_label", "tool", "source", "batch", "sample"}
    missing = required - set(md.columns)

    if missing:
        raise ValueError(
            f"Metadata file is missing required columns: {missing}\n"
            f"Available columns: {list(md.columns)}"
        )

    md["vector_index"] = md["vector_index"].astype(int)
    md = md.sort_values("vector_index").reset_index(drop=True)

    return md


# =============================================================================
# LOAD DISTANCES FROM INTEGRATED VCF
# =============================================================================

def load_native_lifted_distances(vcf_path, metadata_path):
    md = load_metadata(metadata_path)

    rows = []

    if not os.path.exists(vcf_path):
        raise FileNotFoundError(f"VCF not found: {vcf_path}")

    native_col_indices = md.index[md["source"].eq("native_grch38")].tolist()
    lifted_col_indices = md.index[md["source"].eq("lifted_chm13_to_grch38")].tolist()

    if not native_col_indices:
        raise RuntimeError("No native_grch38 entries found in metadata.")

    if not lifted_col_indices:
        raise RuntimeError("No lifted_chm13_to_grch38 entries found in metadata.")

    n_total = 0
    n_both = 0
    n_used = 0

    with open_text(vcf_path) as fh:
        sample_names = None

        for line in fh:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_names = header[9:]

                if len(sample_names) != len(md):
                    print(
                        "[WARNING] Number of VCF sample columns does not match metadata rows.\n"
                        f"          VCF sample columns: {len(sample_names)}\n"
                        f"          Metadata rows:     {len(md)}\n"
                        "          I will still use metadata row order against FORMAT columns."
                    )
                continue

            if line.startswith("#"):
                continue

            n_total += 1

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue

            chrom = cols[0]
            pos = safe_int(cols[1])
            vid = cols[2]
            info = parse_info(cols[7])
            fmt_keys = cols[8].split(":")
            sample_values = cols[9:]

            svtype = str(info.get("SVTYPE", "OTHER")).upper()
            if svtype not in CANONICAL_SVTYPES:
                continue

            native_supp = safe_int(info.get("NATIVE_GRCH38_SUPP"), 0)
            lifted_supp = safe_int(info.get("LIFTED_CHM13_GRCH38_SUPP"), 0)

            if native_supp <= 0 or lifted_supp <= 0:
                continue

            n_both += 1

            native_positions = []
            lifted_positions = []
            native_lengths = []
            lifted_lengths = []

            for idx in native_col_indices:
                if idx >= len(sample_values):
                    continue

                sd = parse_format_sample(fmt_keys, sample_values[idx])

                if not is_supported_sample(sd):
                    continue

                native_positions.extend(parse_co_positions(sd.get("CO"), fallback_pos=pos))
                native_lengths.extend(parse_ln_values(sd.get("LN")))

            for idx in lifted_col_indices:
                if idx >= len(sample_values):
                    continue

                sd = parse_format_sample(fmt_keys, sample_values[idx])

                if not is_supported_sample(sd):
                    continue

                lifted_positions.extend(parse_co_positions(sd.get("CO"), fallback_pos=pos))
                lifted_lengths.extend(parse_ln_values(sd.get("LN")))

            breakpoint_diff = min_abs_pairwise_diff(native_positions, lifted_positions)
            svlen_diff = min_abs_pairwise_diff(native_lengths, lifted_lengths)

            if np.isnan(breakpoint_diff) and np.isnan(svlen_diff):
                continue

            n_used += 1

            rows.append({
                "CHROM": chrom,
                "POS": pos,
                "ID": vid,
                "SVTYPE": svtype,
                "NATIVE_GRCH38_SUPP": native_supp,
                "LIFTED_CHM13_GRCH38_SUPP": lifted_supp,
                "n_native_positions": len(native_positions),
                "n_lifted_positions": len(lifted_positions),
                "n_native_lengths": len(native_lengths),
                "n_lifted_lengths": len(lifted_lengths),
                "breakpoint_diff": breakpoint_diff,
                "svlen_diff": svlen_diff,
            })

    df = pd.DataFrame(rows)

    print("Parsed VCF records:", n_total)
    print("Canonical records with native + lifted support:", n_both)
    print("Records with usable distance metrics:", n_used)

    if df.empty:
        raise RuntimeError("No usable native/lifted support distance rows were produced.")

    return df


# =============================================================================
# PLOT
# =============================================================================

def add_median_p90(ax, values, x_offset_frac=0.01):
    values = pd.Series(values).dropna().astype(float)

    if values.empty:
        return np.nan, np.nan

    median = values.median()
    p90 = values.quantile(0.90)

    ax.axvline(
        median,
        color="#222222",
        linestyle="--",
        linewidth=1.7,
        zorder=5
    )

    ax.axvline(
        p90,
        color="#D85A30",
        linestyle="--",
        linewidth=1.7,
        zorder=5
    )

    ymax = ax.get_ylim()[1]
    xspan = ax.get_xlim()[1] - ax.get_xlim()[0]

    ax.text(
        median + xspan * x_offset_frac,
        ymax * 0.92,
        f"median={median:.0f}",
        ha="left",
        va="top",
        fontsize=10,
        color="#222222"
    )

    ax.text(
        p90 + xspan * x_offset_frac,
        ymax * 0.78,
        f"p90={p90:.0f}",
        ha="left",
        va="top",
        fontsize=10,
        color="#222222"
    )

    return median, p90


def plot_distributions(df, out_dir, out_prefix):
    os.makedirs(out_dir, exist_ok=True)

    pos = df["breakpoint_diff"].dropna().astype(float)
    svlen = df["svlen_diff"].dropna().astype(float)

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.3))

    # -------------------------------------------------------------------------
    # A. Breakpoint distance
    # -------------------------------------------------------------------------
    ax = axes[0]
    style_ax(ax)

    pos_plot = pos[pos <= pos.quantile(0.995)]
    pos_max = max(100, np.ceil(pos_plot.max() / 10) * 10)

    bins = np.arange(0, pos_max + 5, 5)

    ax.hist(
        pos_plot,
        bins=bins,
        color="#5B55B5",
        edgecolor="none",
        alpha=0.9
    )

    ax.set_xlim(-2, pos_max)
    ax.set_xlabel("Absolute position difference (bp)")
    ax.set_ylabel("Matched native/lifted SV clusters")
    ax.set_title("A  Breakpoint distance", loc="left", fontweight="bold")

    pos_median, pos_p90 = add_median_p90(ax, pos)

    # -------------------------------------------------------------------------
    # B. SVLEN difference
    # -------------------------------------------------------------------------
    ax = axes[1]
    style_ax(ax)

    svlen_plot = svlen[svlen <= svlen.quantile(0.995)]
    svlen_max = max(250, np.ceil(svlen_plot.max() / 50) * 50)

    bins = np.arange(0, svlen_max + 10, 10)

    ax.hist(
        svlen_plot,
        bins=bins,
        color="#4C78A8",
        edgecolor="none",
        alpha=0.9
    )

    ax.set_xlim(-5, svlen_max)
    ax.set_xlabel("Absolute SVLEN difference (bp)")
    ax.set_ylabel("Matched native/lifted SV clusters")
    ax.set_title("B  Length difference", loc="left", fontweight="bold")

    svlen_median, svlen_p90 = add_median_p90(ax, svlen)

    fig.suptitle(
        "Native GRCh38 vs lifted CHM13 support distance within integrated GRCh38 cohort",
        fontsize=16,
        fontweight="bold",
        y=1.02
    )

    fig.tight_layout()

    out_png = os.path.join(out_dir, out_prefix + ".png")
    out_pdf = os.path.join(out_dir, out_prefix + ".pdf")

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)

    summary = pd.DataFrame([
        {
            "metric": "breakpoint_diff",
            "n": int(pos.notna().sum()),
            "median": pos_median,
            "p90": pos_p90,
            "mean": pos.mean(),
            "max": pos.max(),
        },
        {
            "metric": "svlen_diff",
            "n": int(svlen.notna().sum()),
            "median": svlen_median,
            "p90": svlen_p90,
            "mean": svlen.mean(),
            "max": svlen.max(),
        },
    ])

    out_table = os.path.join(out_dir, out_prefix + "_table.tsv")
    out_summary = os.path.join(out_dir, out_prefix + "_summary.tsv")

    df.to_csv(out_table, sep="\t", index=False)
    summary.to_csv(out_summary, sep="\t", index=False)

    print()
    print("Saved:")
    print(out_png)
    print(out_pdf)
    print(out_table)
    print(out_summary)
    print()
    print(summary.to_string(index=False))


# =============================================================================
# MAIN
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot native GRCh38 vs lifted CHM13 breakpoint/SVLEN support distances."
    )

    parser.add_argument(
        "--vcf",
        default=VCF,
        help="Integrated GRCh38 cohort VCF.gz."
    )

    parser.add_argument(
        "--metadata",
        default=METADATA,
        help="Tool/reference metadata TSV."
    )

    parser.add_argument(
        "--out-dir",
        default=OUT_DIR,
        help="Output directory."
    )

    parser.add_argument(
        "--out-prefix",
        default=OUT_PREFIX,
        help="Output file prefix."
    )

    return parser.parse_args()


def main():
    args = parse_args()

    df = load_native_lifted_distances(
        vcf_path=args.vcf,
        metadata_path=args.metadata
    )

    plot_distributions(
        df=df,
        out_dir=args.out_dir,
        out_prefix=args.out_prefix
    )


if __name__ == "__main__":
    main()