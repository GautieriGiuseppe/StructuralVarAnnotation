#!/usr/bin/env python3
"""
Plot alignment QC summaries from ALFRED and mosdepth outputs.

Designed for StructuralVarAnnotation output layout:

  <outdir>/<batch>/<sample>/01.align/<ref>/<sample>.srt.bam
  <outdir>/<batch>/<sample>/02.alignqc/<ref>/<sample>.alfred.tsv.gz
  <outdir>/<batch>/<sample>/02.alignqc/<ref>/<sample>_<ref>.mosdepth.summary.txt
  <outdir>/<batch>/<sample>/02.alignqc/<ref>/<sample>_<ref>.mosdepth.global.dist.txt

Outputs:
  <prefix>_alfred_mapping_summary.tsv
  <prefix>_mosdepth_summary.tsv
  <prefix>_mosdepth_contig_coverage.tsv
  <prefix>_mosdepth_global_distribution.tsv
  <prefix>_alfred_mapping_unique.png/pdf
  <prefix>_mosdepth_coverage_summary.png/pdf
  <prefix>_mosdepth_contig_coverage.png/pdf
  <prefix>_mosdepth_global_distribution.png/pdf
"""

import argparse
import gzip
import os
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


DEFAULT_REFS = ["grch38", "chm13"]


# =============================================================================
# CLI
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot alignment QC summaries from ALFRED and mosdepth outputs."
    )

    parser.add_argument(
        "--samples",
        required=True,
        help="Samples TSV with sample_id and batch_id columns.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Workflow output directory containing <batch>/<sample> folders.",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        help="Directory where figures and summary tables will be written.",
    )
    parser.add_argument(
        "--out-prefix",
        default="alignment_qc",
        help="Prefix for output files.",
    )
    parser.add_argument(
        "--refs",
        default=",".join(DEFAULT_REFS),
        help="Comma-separated references to summarize. Default: grch38,chm13.",
    )
    parser.add_argument(
        "--mapq-threshold",
        type=int,
        default=20,
        help="MAPQ threshold used for unique/confident mapping from BAM.",
    )
    parser.add_argument(
        "--skip-unique-mapping",
        "--skip-bam-unique",
        dest="skip_unique_mapping",
        action="store_true",
        help=(
            "Skip samtools BAM counting for unique MAPQ-based mapping rate. "
            "Recommended for pipeline reports because it avoids repeated full BAM scans."
        ),
    )
    parser.add_argument(
        "--title-prefix",
        default="Alignment QC",
        help="Title prefix for plots.",
    )
    parser.add_argument(
        "--example-sample",
        default=None,
        help=(
            "Sample ID to preserve for sample-specific example plots. "
            "If not provided and --trio-file is not provided, the first sample "
            "in the samples TSV is used."
        ),
    )
    parser.add_argument(
        "--trio-file",
        default=None,
        help=(
            "Optional trio TSV with a proband column. If provided and "
            "--example-sample is not set, the first proband is used as the "
            "preserved example sample."
        ),
    )
    parser.add_argument(
        "--contig-top-n",
        type=int,
        default=30,
        help="Maximum number of contigs shown per reference in the contig coverage plot.",
    )
    parser.add_argument(
        "--max-coverage-depth",
        type=float,
        default=120,
        help="Maximum coverage depth shown in the mosdepth global distribution plot.",
    )

    return parser.parse_args()


# =============================================================================
# UTILITIES
# =============================================================================

def read_samples(samples_tsv):
    df = pd.read_csv(samples_tsv, sep="\t")
    required = {"sample_id", "batch_id"}
    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"Samples file is missing required columns: {missing}. "
            f"Available columns: {list(df.columns)}"
        )

    return df[["batch_id", "sample_id"]].drop_duplicates().copy()


def get_example_sample_from_trio_file(trio_file):
    if trio_file is None or str(trio_file).strip() == "":
        return None

    trio_file = str(trio_file)

    if not os.path.exists(trio_file):
        print(f"Trio file not found: {trio_file}")
        return None

    trio_df = pd.read_csv(trio_file, sep="\t")

    if "proband" not in trio_df.columns:
        raise ValueError(
            f"Trio file must contain a 'proband' column. "
            f"Available columns: {list(trio_df.columns)}"
        )

    if trio_df.empty:
        print(f"Trio file is empty: {trio_file}")
        return None

    proband = str(trio_df.iloc[0]["proband"]).strip()

    if proband == "" or proband.lower() == "nan":
        return None

    return proband


def safe_float(x):
    try:
        if x is None or str(x) in {"", ".", "NA", "NaN", "nan"}:
            return np.nan
        return float(x)
    except Exception:
        return np.nan


def run_count(cmd):
    result = subprocess.run(
        cmd,
        check=True,
        text=True,
        capture_output=True,
    )

    value = result.stdout.strip()
    if value == "":
        return 0

    return int(value)


def natural_contig_key(contig):
    s = str(contig)
    s0 = s.replace("chr", "")

    special = {
        "X": 23,
        "Y": 24,
        "M": 25,
        "MT": 25,
    }

    if s0 in special:
        return (0, special[s0], s)

    if s0.isdigit():
        return (0, int(s0), s)

    match = re.match(r"^([A-Za-z_]+)(\d+)$", s0)
    if match:
        return (1, match.group(1), int(match.group(2)), s)

    return (2, s)


# =============================================================================
# ALFRED PARSING
# =============================================================================

def read_alfred_me(path):
    """
    Read only the ALFRED ME alignment summary section.
    """
    header = None
    values = None

    with gzip.open(path, "rt") as handle:
        for line in handle:
            line = line.rstrip("\n")

            if not line.startswith("ME\t"):
                continue

            fields = line.split("\t")

            if len(fields) < 3:
                continue

            if fields[1] == "Sample":
                header = fields
            elif header is not None:
                values = fields
                break

    if header is None or values is None:
        raise ValueError(f"Could not find ME summary section in {path}")

    n = min(len(header), len(values))
    return dict(zip(header[:n], values[:n]))


def get_optional_float(metrics, key):
    if key not in metrics:
        return np.nan

    return safe_float(metrics[key])


def get_required_float(metrics, key):
    if key not in metrics:
        raise KeyError(
            f"Metric {key} not found. Available metrics: {list(metrics.keys())}"
        )

    return safe_float(metrics[key])


def get_alignment_error_rate(metrics):
    if "ErrorRate" in metrics:
        return safe_float(metrics["ErrorRate"]), "ErrorRate"

    for key in metrics:
        if "ErrorRate" in key:
            return safe_float(metrics[key]), key

    mismatch = get_optional_float(metrics, "MismatchRate")
    deletion = get_optional_float(metrics, "DeletionRate")
    insertion = get_optional_float(metrics, "InsertionRate")

    if pd.notna(mismatch) and pd.notna(deletion) and pd.notna(insertion):
        return mismatch + deletion + insertion, "MismatchRate + DeletionRate + InsertionRate"

    return np.nan, "NA"


def find_bam(workflow_outdir, batch, sample, ref):
    sample_dir = Path(workflow_outdir) / str(batch) / str(sample)

    preferred = sample_dir / "01.align" / ref / f"{sample}.srt.bam"
    if preferred.exists():
        return preferred

    candidates = []
    candidates.extend((sample_dir / "01.align" / ref).glob("*.bam"))
    candidates.extend((sample_dir / "01.align" / ref).glob("*.cram"))
    candidates.extend((sample_dir / "02.alignqc" / ref).glob("*.bam"))
    candidates.extend((sample_dir / "02.alignqc" / ref).glob("*.cram"))
    candidates.extend(sample_dir.glob(f"**/{ref}/*.bam"))
    candidates.extend(sample_dir.glob(f"**/{ref}/*.cram"))
    candidates.extend(sample_dir.glob(f"**/*{ref}*.bam"))
    candidates.extend(sample_dir.glob(f"**/*{ref}*.cram"))

    clean = []
    seen = set()

    for path in candidates:
        if str(path).endswith((".bai", ".crai", ".csi")):
            continue
        if path not in seen:
            clean.append(path)
            seen.add(path)

    if not clean:
        return None

    return clean[0]


def get_unique_mapping_stats(bam_path, mapq_threshold=20):
    total_primary = run_count([
        "samtools", "view", "-c", "-F", "2304", str(bam_path)
    ])

    primary_mapped = run_count([
        "samtools", "view", "-c", "-F", "2308", str(bam_path)
    ])

    unique_mapped = run_count([
        "samtools", "view", "-c", "-F", "2308", "-q", str(mapq_threshold), str(bam_path)
    ])

    if total_primary == 0:
        primary_mapping_percent = np.nan
        unique_mapping_percent = np.nan
    else:
        primary_mapping_percent = primary_mapped / total_primary * 100
        unique_mapping_percent = unique_mapped / total_primary * 100

    return {
        "total_primary_reads_bam": total_primary,
        "primary_mapped_reads_bam": primary_mapped,
        "unique_mapped_reads_bam": unique_mapped,
        "primary_mapping_percent_from_bam": primary_mapping_percent,
        "unique_mapping_percent": unique_mapping_percent,
    }


def collect_alfred_metrics(samples_df, workflow_outdir, refs, mapq_threshold, skip_unique_mapping):
    rows = []

    for _, sample_row in samples_df.iterrows():
        batch = str(sample_row["batch_id"])
        sample = str(sample_row["sample_id"])

        for ref in refs:
            alfred_path = (
                Path(workflow_outdir)
                / batch
                / sample
                / "02.alignqc"
                / ref
                / f"{sample}.alfred.tsv.gz"
            )

            if not alfred_path.exists():
                print(f"Missing ALFRED file: {alfred_path}")
                continue

            try:
                metrics = read_alfred_me(alfred_path)
            except Exception as exc:
                print(f"Could not parse ALFRED file {alfred_path}: {exc}")
                continue

            mapped_reads = get_required_float(metrics, "#Mapped")
            mapped_fraction = get_required_float(metrics, "MappedFraction")
            unmapped_reads = get_required_float(metrics, "#Unmapped")
            unmapped_fraction = get_required_float(metrics, "UnmappedFraction")
            total_reads_alfred = mapped_reads + unmapped_reads
            mapping_percent_alfred = mapped_fraction * 100

            alignment_error_rate, error_source = get_alignment_error_rate(metrics)

            bam_path = find_bam(workflow_outdir, batch, sample, ref)

            unique_stats = {
                "total_primary_reads_bam": np.nan,
                "primary_mapped_reads_bam": np.nan,
                "unique_mapped_reads_bam": np.nan,
                "primary_mapping_percent_from_bam": np.nan,
                "unique_mapping_percent": np.nan,
            }

            if bam_path is not None and not skip_unique_mapping:
                try:
                    unique_stats = get_unique_mapping_stats(
                        bam_path=bam_path,
                        mapq_threshold=mapq_threshold,
                    )
                except Exception as exc:
                    print(f"Could not compute unique mapping from BAM {bam_path}: {exc}")

            row = {
                "batch_id": batch,
                "sample": sample,
                "reference": ref,
                "alfred_path": str(alfred_path),
                "bam_path": str(bam_path) if bam_path is not None else "NA",
                "mapped_reads_alfred": mapped_reads,
                "unmapped_reads_alfred": unmapped_reads,
                "total_reads_alfred": total_reads_alfred,
                "mapped_fraction_alfred": mapped_fraction,
                "mapping_percent_alfred": mapping_percent_alfred,
                "unmapped_fraction_alfred": unmapped_fraction,
                "alignment_error_rate": alignment_error_rate,
                "alignment_error_percent": (
                    alignment_error_rate * 100
                    if pd.notna(alignment_error_rate)
                    else np.nan
                ),
                "error_source": error_source,
                "secondary_alignments": get_optional_float(metrics, "#SecondaryAlignments"),
                "secondary_fraction": get_optional_float(metrics, "SecondaryAlignmentFraction"),
                "supplementary_alignments": get_optional_float(metrics, "#SupplementaryAlignments"),
                "supplementary_fraction": get_optional_float(metrics, "SupplementaryAlignmentFraction"),
                "median_mapq": get_optional_float(metrics, "MedianMAPQ"),
                "mapq_threshold_for_unique": mapq_threshold,
            }

            row.update(unique_stats)
            rows.append(row)

            print(f"{sample} {ref}")
            print(f"  ALFRED overall mapping: {mapping_percent_alfred:.4f}%")
            if pd.notna(row["unique_mapping_percent"]):
                print(f"  BAM unique MAPQ>={mapq_threshold}: {row['unique_mapping_percent']:.4f}%")
            print(f"  alignment error rate: {alignment_error_rate}")
            print(f"  ALFRED: {alfred_path}")
            print(f"  BAM: {bam_path}")
            print()

    return pd.DataFrame(rows)


# =============================================================================
# MOSDEPTH PARSING
# =============================================================================

def read_mosdepth_summary(summary_path):
    df = pd.read_csv(summary_path, sep="\t")

    for col in ["length", "bases", "mean", "min", "max"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def collect_mosdepth_summary(samples_df, workflow_outdir, refs):
    rows = []

    for _, sample_row in samples_df.iterrows():
        batch = str(sample_row["batch_id"])
        sample = str(sample_row["sample_id"])

        for ref in refs:
            summary_path = (
                Path(workflow_outdir)
                / batch
                / sample
                / "02.alignqc"
                / ref
                / f"{sample}_{ref}.mosdepth.summary.txt"
            )

            if not summary_path.exists():
                print(f"Missing mosdepth summary: {summary_path}")
                continue

            try:
                df = read_mosdepth_summary(summary_path)
            except Exception as exc:
                print(f"Could not read mosdepth summary {summary_path}: {exc}")
                continue

            if df.empty:
                continue

            total_rows = (
                df[df["chrom"].astype(str).isin(["total", "genome"])]
                if "chrom" in df.columns
                else pd.DataFrame()
            )

            if total_rows.empty:
                total_row = df.iloc[-1]
            else:
                total_row = total_rows.iloc[-1]

            row = total_row.to_dict()
            row.update({
                "batch_id": batch,
                "sample": sample,
                "reference": ref,
                "mosdepth_summary_path": str(summary_path),
            })

            rows.append(row)

    if not rows:
        return pd.DataFrame()

    return pd.DataFrame(rows)


def collect_mosdepth_contig_coverage(samples_df, workflow_outdir, refs):
    rows = []

    for _, sample_row in samples_df.iterrows():
        batch = str(sample_row["batch_id"])
        sample = str(sample_row["sample_id"])

        for ref in refs:
            summary_path = (
                Path(workflow_outdir)
                / batch
                / sample
                / "02.alignqc"
                / ref
                / f"{sample}_{ref}.mosdepth.summary.txt"
            )

            if not summary_path.exists():
                continue

            try:
                df = read_mosdepth_summary(summary_path)
            except Exception as exc:
                print(f"Could not read mosdepth contig coverage {summary_path}: {exc}")
                continue

            if df.empty or "chrom" not in df.columns:
                continue

            d = df.copy()
            d["chrom"] = d["chrom"].astype(str)
            d = d[~d["chrom"].isin(["total", "genome"])].copy()

            if d.empty:
                continue

            d["batch_id"] = batch
            d["sample"] = sample
            d["reference"] = ref
            d["mosdepth_summary_path"] = str(summary_path)

            rows.append(d)

    if not rows:
        return pd.DataFrame()

    out = pd.concat(rows, ignore_index=True)

    for col in ["length", "bases", "mean", "min", "max"]:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")

    return out


def collect_mosdepth_global_distribution(samples_df, workflow_outdir, refs):
    rows = []

    for _, sample_row in samples_df.iterrows():
        batch = str(sample_row["batch_id"])
        sample = str(sample_row["sample_id"])

        for ref in refs:
            global_path = (
                Path(workflow_outdir)
                / batch
                / sample
                / "02.alignqc"
                / ref
                / f"{sample}_{ref}.mosdepth.global.dist.txt"
            )

            if not global_path.exists():
                print(f"Missing mosdepth global distribution: {global_path}")
                continue

            try:
                df = pd.read_csv(
                    global_path,
                    sep="\t",
                    header=None,
                    names=["chrom", "coverage", "fraction"],
                )
            except Exception as exc:
                print(f"Could not read mosdepth global distribution {global_path}: {exc}")
                continue

            if df.empty:
                continue

            genome = df[df["chrom"].astype(str).isin(["total", "genome"])]

            if genome.empty:
                genome = df.copy()

            genome["batch_id"] = batch
            genome["sample"] = sample
            genome["reference"] = ref
            genome["mosdepth_global_dist_path"] = str(global_path)
            genome["coverage"] = pd.to_numeric(genome["coverage"], errors="coerce")
            genome["fraction"] = pd.to_numeric(genome["fraction"], errors="coerce")

            rows.append(genome)

    if not rows:
        return pd.DataFrame()

    return pd.concat(rows, ignore_index=True)


# =============================================================================
# PLOTTING
# =============================================================================

def save_fig(fig, path_prefix):
    png = f"{path_prefix}.png"
    pdf = f"{path_prefix}.pdf"

    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, dpi=300, bbox_inches="tight")

    plt.close(fig)

    print(f"Saved: {png}")
    print(f"Saved: {pdf}")


def save_placeholder(path_prefix, title, message):
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.axis("off")

    ax.text(
        0.5,
        0.56,
        title,
        ha="center",
        va="center",
        fontsize=16,
        fontweight="bold",
    )
    ax.text(
        0.5,
        0.40,
        message,
        ha="center",
        va="center",
        fontsize=11,
        wrap=True,
    )

    save_fig(fig, path_prefix)


def annotate_sample_labels(ax, df, column, ref_order):
    if df is None or df.empty or column not in df.columns:
        return

    d = df[df[column].notna()].copy()
    if d.empty:
        return

    samples = sorted(d["sample"].astype(str).unique())
    if not samples:
        return

    offsets = {}
    if len(samples) == 1:
        offsets[samples[0]] = 0.0
    else:
        values = np.linspace(-0.18, 0.18, len(samples))
        offsets = dict(zip(samples, values))

    ref_to_x = {ref: i for i, ref in enumerate(ref_order)}

    for _, row in d.iterrows():
        ref = str(row["reference"])
        sample = str(row["sample"])

        if ref not in ref_to_x:
            continue

        x = ref_to_x[ref] + offsets.get(sample, 0.0)
        y = row[column]

        ax.text(
            x + 0.03,
            y,
            sample,
            fontsize=7,
            alpha=0.8,
            va="center",
        )


def plot_alfred_mapping(df, out_prefix, title_prefix, mapq_threshold, example_sample=None):
    if df is None or df.empty:
        print("No ALFRED data available for plotting.")
        save_placeholder(
            f"{out_prefix}_alfred_mapping_unique",
            f"{title_prefix}: ALFRED mapping metrics",
            "No ALFRED ME records were available for plotting.",
        )
        return

    sns.set_theme(style="whitegrid")

    ref_order = [ref for ref in DEFAULT_REFS if ref in set(df["reference"].astype(str))]
    extra_refs = [
        ref for ref in sorted(df["reference"].astype(str).unique())
        if ref not in ref_order
    ]
    ref_order.extend(extra_refs)

    fig, axes = plt.subplots(1, 3, figsize=(19, 4.8))

    panels = [
        (
            "alignment_error_percent",
            "Alignment Error Rate",
            "Error rate (%)",
        ),
        (
            "mapping_percent_alfred",
            "Overall Mapping Rate",
            "Mapped reads (%)",
        ),
        (
            "unique_mapping_percent",
            f"Unique Mapping Rate, MAPQ >= {mapq_threshold}",
            "Unique/confident mapped reads (%)",
        ),
    ]

    for ax, (column, title, ylabel) in zip(axes, panels):
        if column not in df.columns or df[column].dropna().empty:
            ax.text(
                0.5,
                0.5,
                f"No data for {column}",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_axis_off()
            continue

        sns.boxplot(
            data=df,
            x="reference",
            y=column,
            order=ref_order,
            ax=ax,
            showfliers=False,
        )

        sns.stripplot(
            data=df,
            x="reference",
            y=column,
            order=ref_order,
            color="black",
            alpha=0.55,
            size=5,
            jitter=False,
            ax=ax,
        )

        annotate_sample_labels(ax, df, column, ref_order)

        ax.set_title(title, fontweight="bold")
        ax.set_xlabel("Reference genome")
        ax.set_ylabel(ylabel)
        ax.set_ylim(0, 100)
        ax.set_yticks([0, 20, 40, 60, 80, 100])

    if example_sample:
        subtitle = f"sample labels shown; preserved example sample: {example_sample}"
    else:
        subtitle = "sample labels shown"

    fig.suptitle(
        f"{title_prefix}: ALFRED mapping metrics ({subtitle})",
        fontsize=15,
        fontweight="bold",
        y=1.05,
    )

    fig.tight_layout()
    save_fig(fig, f"{out_prefix}_alfred_mapping_unique")


def plot_mosdepth_summary(df, out_prefix, title_prefix):
    if df is None or df.empty or "mean" not in df.columns:
        print("No mosdepth summary data available for plotting.")
        save_placeholder(
            f"{out_prefix}_mosdepth_coverage_summary",
            f"{title_prefix}: mosdepth coverage summary",
            "No mosdepth summary records were available for plotting.",
        )
        return

    sns.set_theme(style="whitegrid")

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

    sns.boxplot(
        data=df,
        x="reference",
        y="mean",
        ax=axes[0],
        showfliers=False,
    )
    sns.stripplot(
        data=df,
        x="reference",
        y="mean",
        color="black",
        alpha=0.55,
        size=5,
        ax=axes[0],
    )
    axes[0].set_title("Mean coverage", fontweight="bold")
    axes[0].set_xlabel("Reference genome")
    axes[0].set_ylabel("Mean coverage")

    if "bases" in df.columns:
        plot_df = df.copy()
        plot_df["bases_gb"] = plot_df["bases"] / 1e9

        sns.boxplot(
            data=plot_df,
            x="reference",
            y="bases_gb",
            ax=axes[1],
            showfliers=False,
        )
        sns.stripplot(
            data=plot_df,
            x="reference",
            y="bases_gb",
            color="black",
            alpha=0.55,
            size=5,
            ax=axes[1],
        )
        axes[1].set_title("Aligned bases", fontweight="bold")
        axes[1].set_xlabel("Reference genome")
        axes[1].set_ylabel("Aligned bases (Gb)")
    else:
        axes[1].text(0.5, 0.5, "No bases column", ha="center", va="center")
        axes[1].set_axis_off()

    fig.suptitle(
        f"{title_prefix}: mosdepth coverage summary",
        fontsize=16,
        fontweight="bold",
        y=1.05,
    )

    fig.tight_layout()
    save_fig(fig, f"{out_prefix}_mosdepth_coverage_summary")


def select_contigs_for_plot(df, example_sample, contig_top_n):
    d = df.copy()

    if example_sample is not None and "sample" in d.columns:
        example_df = d[d["sample"].astype(str) == str(example_sample)].copy()
        if not example_df.empty:
            d = example_df
            sample_label = str(example_sample)
        else:
            sample_label = "all samples"
    else:
        sample_label = "all samples"

    d = d[d["mean"].notna()].copy()

    if d.empty:
        return d, sample_label

    selected = []

    for ref, sub in d.groupby("reference", sort=False):
        s = sub.copy()

        if "length" in s.columns and s["length"].notna().any():
            s = s.sort_values("length", ascending=False).head(contig_top_n)
        else:
            s = s.head(contig_top_n)

        selected.append(s)

    out = pd.concat(selected, ignore_index=True)

    contigs = list(dict.fromkeys(out["chrom"].astype(str).tolist()))
    contigs = sorted(contigs, key=natural_contig_key)

    out["chrom"] = pd.Categorical(
        out["chrom"].astype(str),
        categories=contigs,
        ordered=True,
    )

    out = out.sort_values(["reference", "chrom"])

    return out, sample_label


def plot_mosdepth_contig_coverage(df, out_prefix, title_prefix, example_sample=None, contig_top_n=30):
    if df is None or df.empty or "mean" not in df.columns or "chrom" not in df.columns:
        print("No mosdepth contig coverage data available for plotting.")
        save_placeholder(
            f"{out_prefix}_mosdepth_contig_coverage",
            f"{title_prefix}: mosdepth coverage by contig",
            "No contig-level mosdepth summary records were available for plotting.",
        )
        return

    d, sample_label = select_contigs_for_plot(df, example_sample, contig_top_n)

    if d.empty:
        save_placeholder(
            f"{out_prefix}_mosdepth_contig_coverage",
            f"{title_prefix}: mosdepth coverage by contig",
            "No contig-level mosdepth summary records remained after filtering.",
        )
        return

    sns.set_theme(style="whitegrid")

    n_contigs = d["chrom"].nunique()
    width = max(12, min(26, n_contigs * 0.45))

    fig, ax = plt.subplots(figsize=(width, 5.2))

    # Use a line plot, consistent with the global coverage distribution panel.
    # Markers are kept so individual contigs remain readable.
    sns.lineplot(
        data=d,
        x="chrom",
        y="mean",
        hue="reference",
        marker="o",
        linewidth=1.8,
        markersize=4.5,
        ax=ax,
    )

    ax.set_title(
        f"Coverage by contig: {sample_label}",
        fontweight="bold",
    )
    ax.set_xlabel("Contig")
    ax.set_ylabel("Mean coverage")
    ax.tick_params(axis="x", rotation=75)

    fig.suptitle(
        f"{title_prefix}: mosdepth coverage by contig",
        fontsize=16,
        fontweight="bold",
        y=1.03,
    )

    fig.tight_layout()
    save_fig(fig, f"{out_prefix}_mosdepth_contig_coverage")


def plot_mosdepth_global_distribution(df, out_prefix, title_prefix, example_sample=None, max_coverage_depth=120):
    if df is None or df.empty:
        print("No mosdepth global distribution data available for plotting.")
        save_placeholder(
            f"{out_prefix}_mosdepth_global_distribution",
            f"{title_prefix}: mosdepth global distribution",
            "No mosdepth global distribution records were available for plotting.",
        )
        return

    d = df.copy()

    if example_sample is not None and "sample" in d.columns:
        example_df = d[d["sample"].astype(str) == str(example_sample)].copy()

        if not example_df.empty:
            d = example_df
            sample_label = str(example_sample)
        else:
            sample_label = "all samples"
    else:
        sample_label = "all samples"

    d = d[d["coverage"].notna() & d["fraction"].notna()]
    d = d[d["coverage"] <= max_coverage_depth]

    if d.empty:
        print("No mosdepth global distribution rows after filtering.")
        save_placeholder(
            f"{out_prefix}_mosdepth_global_distribution",
            f"{title_prefix}: mosdepth global distribution",
            "No mosdepth global distribution rows remained after filtering.",
        )
        return

    sns.set_theme(style="whitegrid")

    fig, ax = plt.subplots(figsize=(10, 5))

    sns.lineplot(
        data=d,
        x="coverage",
        y="fraction",
        hue="reference",
        units="sample",
        estimator=None,
        alpha=0.75,
        linewidth=1.6,
        ax=ax,
    )

    ax.set_title(
        f"Genome-wide coverage distribution: {sample_label}",
        fontweight="bold",
    )
    ax.set_xlabel("Coverage depth")
    ax.set_ylabel("Fraction of genome")
    ax.set_ylim(0, 1)
    ax.set_xlim(0, min(max_coverage_depth, d["coverage"].max()))

    fig.suptitle(
        f"{title_prefix}: mosdepth global distribution",
        fontsize=16,
        fontweight="bold",
        y=1.03,
    )

    fig.tight_layout()
    save_fig(fig, f"{out_prefix}_mosdepth_global_distribution")


# =============================================================================
# MAIN
# =============================================================================

def main():
    args = parse_args()

    workflow_outdir = args.outdir.rstrip("/")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    refs = [x.strip() for x in args.refs.split(",") if x.strip()]

    samples_df = read_samples(args.samples)

    if args.example_sample is not None:
        example_sample = str(args.example_sample)
        example_source = "--example-sample"
    else:
        example_sample = get_example_sample_from_trio_file(args.trio_file)
        if example_sample is not None:
            example_source = "--trio-file proband"
        else:
            example_sample = str(samples_df.iloc[0]["sample_id"])
            example_source = "first sample in samples TSV"

    print(
        "Using preserved example sample for Alignment QC: "
        f"{example_sample} ({example_source})"
    )

    if args.skip_unique_mapping:
        print(
            "Skipping BAM-derived unique mapping calculation "
            "(--skip-unique-mapping / --skip-bam-unique)."
        )

    alfred_df = collect_alfred_metrics(
        samples_df=samples_df,
        workflow_outdir=workflow_outdir,
        refs=refs,
        mapq_threshold=args.mapq_threshold,
        skip_unique_mapping=args.skip_unique_mapping,
    )

    mosdepth_summary_df = collect_mosdepth_summary(
        samples_df=samples_df,
        workflow_outdir=workflow_outdir,
        refs=refs,
    )

    mosdepth_contig_df = collect_mosdepth_contig_coverage(
        samples_df=samples_df,
        workflow_outdir=workflow_outdir,
        refs=refs,
    )

    mosdepth_global_df = collect_mosdepth_global_distribution(
        samples_df=samples_df,
        workflow_outdir=workflow_outdir,
        refs=refs,
    )

    alfred_tsv = out_dir / f"{args.out_prefix}_alfred_mapping_summary.tsv"
    mosdepth_summary_tsv = out_dir / f"{args.out_prefix}_mosdepth_summary.tsv"
    mosdepth_contig_tsv = out_dir / f"{args.out_prefix}_mosdepth_contig_coverage.tsv"
    mosdepth_global_tsv = out_dir / f"{args.out_prefix}_mosdepth_global_distribution.tsv"

    alfred_df.to_csv(alfred_tsv, sep="\t", index=False)
    mosdepth_summary_df.to_csv(mosdepth_summary_tsv, sep="\t", index=False)
    mosdepth_contig_df.to_csv(mosdepth_contig_tsv, sep="\t", index=False)
    mosdepth_global_df.to_csv(mosdepth_global_tsv, sep="\t", index=False)

    print(f"Saved: {alfred_tsv}")
    print(f"Saved: {mosdepth_summary_tsv}")
    print(f"Saved: {mosdepth_contig_tsv}")
    print(f"Saved: {mosdepth_global_tsv}")

    out_prefix = out_dir / args.out_prefix

    plot_alfred_mapping(
        df=alfred_df,
        out_prefix=str(out_prefix),
        title_prefix=args.title_prefix,
        mapq_threshold=args.mapq_threshold,
        example_sample=example_sample,
    )

    plot_mosdepth_summary(
        df=mosdepth_summary_df,
        out_prefix=str(out_prefix),
        title_prefix=args.title_prefix,
    )

    plot_mosdepth_contig_coverage(
        df=mosdepth_contig_df,
        out_prefix=str(out_prefix),
        title_prefix=args.title_prefix,
        example_sample=example_sample,
        contig_top_n=args.contig_top_n,
    )

    plot_mosdepth_global_distribution(
        df=mosdepth_global_df,
        out_prefix=str(out_prefix),
        title_prefix=args.title_prefix,
        example_sample=example_sample,
        max_coverage_depth=args.max_coverage_depth,
    )

    print()
    print("Alignment QC plotting complete.")
    print(f"ALFRED rows: {len(alfred_df)}")
    print(f"mosdepth summary rows: {len(mosdepth_summary_df)}")
    print(f"mosdepth contig rows: {len(mosdepth_contig_df)}")
    print(f"mosdepth global distribution rows: {len(mosdepth_global_df)}")


if __name__ == "__main__":
    main()