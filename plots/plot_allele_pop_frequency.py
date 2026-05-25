
"""
needLR control population frequency plots
-----------------------------------------

Plots:
  A. Mean control AF among present variants
  B. Presence across control populations

Input:
  Final needLR-annotated GRCh38 genotyped cohort VCF.

Outputs:
  - needlr_control_population_frequency_summary.png
  - needlr_control_population_frequency_summary.pdf
  - needlr_control_population_frequency_summary.tsv

Definitions:
  - "Mean control AF among present variants":
      mean Pop_Freq / Allele_Freq for variants where that ancestry frequency > 0.
  - "Presence across control populations":
      percentage of cohort SVs with frequency > 0 in each 1KGP ancestry.
"""

import os
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

OUTDIR = "/group/dominguez/shared_notebooks/Immune_variation/mapping"

NEEDLR_ANNOTATED_VCF = (
    OUTDIR
    + "/needLR_output/"
    + "GRCh38_final_cohort_survivor_genotyped_matrix_needLR_cohort/"
    + "GRCh38_final_cohort_survivor_genotyped_matrix.needlr_input_needLR_1kg_v4.0/"
    + "GRCh38_final_cohort_survivor_genotyped_matrix.needlr_input.needLR.4.0.vcf.gz"
)

OUT_DIR = OUTDIR + "/cohort_results/needlr_annotation_plots"
OUT_PREFIX = "needlr_control_population_frequency_summary"

CANONICAL_SVTYPES = {"DEL", "INS", "DUP", "INV"}

ANCESTRIES = ["AFR", "AMR", "EAS", "EUR", "SAS", "ALL"]
ANCESTRIES_FOR_PRESENCE = ["AFR", "AMR", "EAS", "EUR", "SAS"]

# same blue for most ancestries, darker blue for EUR
BAR_COLORS_AF = {
    "AFR": "#5B83B1",
    "AMR": "#5B83B1",
    "EAS": "#5B83B1",
    "EUR": "#2C5F8A",
    "SAS": "#5B83B1",
    "ALL": "#5B83B1",
}

BAR_COLORS_PRESENCE = {
    "AFR": "#2FA882",
    "AMR": "#2FA882",
    "EAS": "#2FA882",
    "EUR": "#1D7F68",
    "SAS": "#2FA882",
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


def style_ax(ax, grid_axis="y"):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#CCCCCC")
    ax.spines["bottom"].set_color("#CCCCCC")
    ax.tick_params(labelsize=11, colors="#444444")
    ax.grid(axis=grid_axis, color="#EAEAEA", linewidth=0.8)
    ax.set_axisbelow(True)


def fmt_percent(x, _):
    return f"{x:.0f}"


def find_first_float(info, candidate_keys):
    for key in candidate_keys:
        val = as_float(info.get(key))
        if pd.notna(val):
            return val
    return np.nan


def infer_freq(info, ancestry):
    """
    Flexible parser for possible needLR frequency field names.

    It tries Pop_Freq first, then Allele_Freq, then abbreviated AF-style fields.
    """
    candidate_keys = [
        f"Pop_Freq_{ancestry}",
        f"Population_Freq_{ancestry}",
        f"Control_Pop_Freq_{ancestry}",
        f"Allele_Freq_{ancestry}",
        f"Control_Allele_Freq_{ancestry}",
        f"AF_{ancestry}",
        f"AF_1KGP_{ancestry}",
        f"1KGP_AF_{ancestry}",
        f"Freq_{ancestry}",
    ]

    return find_first_float(info, candidate_keys)


def infer_count(info, ancestry):
    """
    Flexible parser for possible needLR count field names.
    Counts are not strictly required, but they help classify presence if AF is missing.
    """
    candidate_keys = [
        f"Pop_Count_{ancestry}",
        f"Population_Count_{ancestry}",
        f"Control_Pop_Count_{ancestry}",
        f"Allele_Count_{ancestry}",
        f"Control_Allele_Count_{ancestry}",
        f"AC_{ancestry}",
        f"AC_1KGP_{ancestry}",
        f"1KGP_AC_{ancestry}",
        f"Count_{ancestry}",
    ]

    for key in candidate_keys:
        val = as_int(info.get(key))
        if pd.notna(val):
            return int(val)

    return np.nan


# =============================================================================
# LOAD VCF
# =============================================================================

def load_population_frequency_table(vcf_path):
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

            row = {
                "CHROM": chrom,
                "POS": int(pos),
                "ID": vid,
                "SVTYPE": svtype,
            }

            for ancestry in ANCESTRIES:
                freq = infer_freq(info, ancestry)
                count = infer_count(info, ancestry)

                row[f"{ancestry}_freq"] = freq
                row[f"{ancestry}_count"] = count

                if pd.notna(freq):
                    row[f"{ancestry}_present"] = freq > 0
                elif pd.notna(count):
                    row[f"{ancestry}_present"] = count > 0
                else:
                    row[f"{ancestry}_present"] = False

            rows.append(row)

    df = pd.DataFrame(rows)

    if df.empty:
        raise RuntimeError("No canonical SVs loaded from needLR VCF.")

    return df


# =============================================================================
# SUMMARIZE
# =============================================================================

def summarize_population_frequency(df):
    total_n = len(df)

    mean_af_rows = []
    for ancestry in ANCESTRIES:
        freq_col = f"{ancestry}_freq"
        present_col = f"{ancestry}_present"

        d = df[df[present_col] & df[freq_col].notna()].copy()

        mean_freq = d[freq_col].mean() if len(d) else np.nan

        mean_af_rows.append(
            {
                "ancestry": ancestry,
                "n_total_variants": total_n,
                "n_present_with_frequency": len(d),
                "mean_frequency_fraction": mean_freq,
                "mean_frequency_percent": mean_freq * 100 if pd.notna(mean_freq) else np.nan,
            }
        )

    mean_af = pd.DataFrame(mean_af_rows)

    presence_rows = []
    for ancestry in ANCESTRIES_FOR_PRESENCE:
        present_col = f"{ancestry}_present"
        n_present = int(df[present_col].sum())
        pct_present = 100.0 * n_present / total_n if total_n else np.nan

        presence_rows.append(
            {
                "ancestry": ancestry,
                "n_total_variants": total_n,
                "n_present": n_present,
                "percent_present": pct_present,
            }
        )

    presence = pd.DataFrame(presence_rows)

    return mean_af, presence


# =============================================================================
# PLOT
# =============================================================================

def plot_population_frequency(mean_af, presence, out_prefix):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.2))
    fig.subplots_adjust(
        left=0.07,
        right=0.98,
        bottom=0.14,
        top=0.84,
        wspace=0.35,
    )

    # -------------------------------------------------------------------------
    # A. Mean control AF among present variants
    # -------------------------------------------------------------------------
    ax = axes[0]
    style_ax(ax)

    x = np.arange(len(mean_af))
    values = mean_af["mean_frequency_percent"].values
    labels = mean_af["ancestry"].tolist()
    colors = [BAR_COLORS_AF[a] for a in labels]

    bars = ax.bar(
        x,
        values,
        color=colors,
        edgecolor="none",
        alpha=0.95,
    )

    ymax = np.nanmax(values) if len(values) else 1
    ax.set_ylim(0, ymax * 1.18)

    for b, v in zip(bars, values):
        if pd.notna(v):
            ax.text(
                b.get_x() + b.get_width() / 2,
                b.get_height() + ymax * 0.02,
                f"{v:.2f}",
                ha="center",
                va="bottom",
                fontsize=11,
            )

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Mean control AF among present SVs (%)", fontsize=12)
    ax.set_title(
        "A  Mean control AF among present variants",
        loc="left",
        fontweight="bold",
        fontsize=15,
    )
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_percent))

    # -------------------------------------------------------------------------
    # B. Presence across control populations
    # -------------------------------------------------------------------------
    ax = axes[1]
    style_ax(ax)

    x = np.arange(len(presence))
    values = presence["percent_present"].values
    labels = presence["ancestry"].tolist()
    colors = [BAR_COLORS_PRESENCE[a] for a in labels]

    bars = ax.bar(
        x,
        values,
        color=colors,
        edgecolor="none",
        alpha=0.95,
    )

    ymax = np.nanmax(values) if len(values) else 1
    ax.set_ylim(0, min(100, ymax * 1.18))

    for b, v in zip(bars, values):
        if pd.notna(v):
            ax.text(
                b.get_x() + b.get_width() / 2,
                b.get_height() + ymax * 0.02,
                f"{v:.1f}",
                ha="center",
                va="bottom",
                fontsize=11,
            )

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("% of SVs present", fontsize=12)
    ax.set_title(
        "B  Presence across control populations",
        loc="left",
        fontweight="bold",
        fontsize=15,
    )
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_percent))

    fig.suptitle(
        "needLR control population frequency summary",
        fontsize=18,
        fontweight="bold",
        y=0.98,
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

    df = load_population_frequency_table(NEEDLR_ANNOTATED_VCF)

    print()
    print("Loaded canonical SVs:", f"{len(df):,}")
    print(df["SVTYPE"].value_counts().to_string())

    mean_af, presence = summarize_population_frequency(df)

    summary = mean_af.merge(
        presence,
        on=["ancestry", "n_total_variants"],
        how="outer",
    )

    out_prefix = os.path.join(OUT_DIR, OUT_PREFIX)

    df.to_csv(
        os.path.join(OUT_DIR, "needlr_control_population_frequency_variant_table.tsv"),
        sep="\t",
        index=False,
    )

    mean_af.to_csv(
        os.path.join(OUT_DIR, "needlr_mean_control_af_among_present.tsv"),
        sep="\t",
        index=False,
    )

    presence.to_csv(
        os.path.join(OUT_DIR, "needlr_control_population_presence.tsv"),
        sep="\t",
        index=False,
    )

    summary.to_csv(
        os.path.join(OUT_DIR, "needlr_control_population_frequency_summary.tsv"),
        sep="\t",
        index=False,
    )

    plot_population_frequency(mean_af, presence, out_prefix)

    print()
    print("Saved:")
    print(out_prefix + ".png")
    print(out_prefix + ".pdf")
    print()
    print("Mean control AF among present variants:")
    print(mean_af.to_string(index=False))
    print()
    print("Presence across control populations:")
    print(presence.to_string(index=False))


if __name__ == "__main__":
    main()