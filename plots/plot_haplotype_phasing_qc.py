#!/usr/bin/env python3

import argparse
import os
import glob
import subprocess
from collections import Counter, defaultdict

import pandas as pd
import matplotlib.pyplot as plt


def run_cmd(cmd):
    return subprocess.check_output(cmd, shell=True, text=True)


def count_gt(vcf):
    out = run_cmd(f"bcftools query -f '[%GT\\n]' {vcf}")
    counts = Counter(line.strip() for line in out.splitlines() if line.strip())
    return counts


def count_phase_sets(vcf):
    try:
        out = run_cmd(f"bcftools query -f '[%PS\\n]' {vcf} 2>/dev/null")
    except subprocess.CalledProcessError:
        return 0

    ps = {
        line.strip()
        for line in out.splitlines()
        if line.strip() and line.strip() != "."
    }
    return len(ps)


def count_haplotagged_reads(bam):
    hp1 = 0
    hp2 = 0
    tagged = 0

    # Count HP:i tags directly from SAM records.
    p = subprocess.Popen(
        f"samtools view {bam}",
        shell=True,
        stdout=subprocess.PIPE,
        text=True,
    )

    for line in p.stdout:
        if "HP:i:1" in line:
            hp1 += 1
            tagged += 1
        elif "HP:i:2" in line:
            hp2 += 1
            tagged += 1
        elif "HP:i:" in line:
            tagged += 1

    p.wait()
    if p.returncode != 0:
        raise RuntimeError(f"samtools view failed for {bam}")

    return hp1, hp2, tagged


def collect_sample_metrics(haplo_dir):
    rows = []

    vcfs = sorted(
        glob.glob(
            os.path.join(
                haplo_dir,
                "*",
                "*",
                "05.haplotyping",
                "grch38",
                "*.whatshap.phased.vcf.gz",
            )
        )
    )

    # Also support running directly on the output root where files may be deeper.
    if not vcfs:
        vcfs = sorted(
            glob.glob(
                os.path.join(
                    haplo_dir,
                    "**",
                    "05.haplotyping",
                    "grch38",
                    "*.whatshap.phased.vcf.gz",
                ),
                recursive=True,
            )
        )

    if not vcfs:
        raise FileNotFoundError(
            f"No WhatsHap phased VCFs found below: {haplo_dir}"
        )

    for vcf in vcfs:
        sample = os.path.basename(vcf).replace(".whatshap.phased.vcf.gz", "")
        sample_dir = os.path.dirname(vcf)
        bam = os.path.join(sample_dir, f"{sample}.haplotagged.bam")

        gt_counts = count_gt(vcf)

        phased_het = gt_counts.get("0|1", 0) + gt_counts.get("1|0", 0)
        unphased_het = gt_counts.get("0/1", 0) + gt_counts.get("1/0", 0)
        hom_ref = gt_counts.get("0/0", 0)
        hom_alt = gt_counts.get("1/1", 0)
        missing = gt_counts.get("./.", 0)

        total_genotypes = sum(gt_counts.values())
        total_het = phased_het + unphased_het

        pct_all_phased = (
            100.0 * phased_het / total_genotypes
            if total_genotypes > 0
            else 0.0
        )

        pct_het_phased = (
            100.0 * phased_het / total_het
            if total_het > 0
            else 0.0
        )

        phase_sets = count_phase_sets(vcf)

        hp1 = hp2 = tagged = 0
        if os.path.exists(bam):
            hp1, hp2, tagged = count_haplotagged_reads(bam)

        rows.append(
            {
                "sample": sample,
                "total_genotypes": total_genotypes,
                "hom_ref_0_0": hom_ref,
                "hom_alt_1_1": hom_alt,
                "missing": missing,
                "unphased_het": unphased_het,
                "phased_het": phased_het,
                "pct_all_genotypes_phased": pct_all_phased,
                "pct_heterozygous_phased": pct_het_phased,
                "phase_sets": phase_sets,
                "hp1_reads": hp1,
                "hp2_reads": hp2,
                "haplotagged_reads": tagged,
                "vcf": vcf,
                "haplotagged_bam": bam if os.path.exists(bam) else "",
            }
        )

    return pd.DataFrame(rows).sort_values("sample")


def save_table_png(df, out_png):
    show_cols = [
        "sample",
        "total_genotypes",
        "phased_het",
        "unphased_het",
        "pct_all_genotypes_phased",
        "pct_heterozygous_phased",
        "phase_sets",
        "haplotagged_reads",
    ]

    table_df = df[show_cols].copy()
    table_df["pct_all_genotypes_phased"] = table_df[
        "pct_all_genotypes_phased"
    ].map(lambda x: f"{x:.2f}")
    table_df["pct_heterozygous_phased"] = table_df[
        "pct_heterozygous_phased"
    ].map(lambda x: f"{x:.2f}")

    fig_h = 0.45 * len(table_df) + 1.2
    fig, ax = plt.subplots(figsize=(14, fig_h))
    ax.axis("off")

    table = ax.table(
        cellText=table_df.values,
        colLabels=[
            "Sample",
            "Total GT",
            "Phased het",
            "Unphased het",
            "% all phased",
            "% het phased",
            "Phase sets",
            "Haplotagged reads",
        ],
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.3)

    ax.set_title(
        "Haplotype phasing summary",
        fontsize=14,
        fontweight="bold",
        pad=12,
    )

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def plot_all_vs_het_rate(df, out_png):
    x = range(len(df))

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(
        [i - 0.2 for i in x],
        df["pct_all_genotypes_phased"],
        width=0.4,
        label="All genotypes",
    )
    ax.bar(
        [i + 0.2 for i in x],
        df["pct_heterozygous_phased"],
        width=0.4,
        label="Heterozygous genotypes",
    )

    ax.set_xticks(list(x))
    ax.set_xticklabels(df["sample"], rotation=45, ha="right")
    ax.set_ylabel("Phased genotypes (%)")
    ax.set_title("WhatsHap phasing rate")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def plot_phase_sets(df, out_png):
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(df["sample"], df["phase_sets"])
    ax.set_ylabel("Number of phase sets")
    ax.set_title("Phase sets per sample")
    ax.set_xticklabels(df["sample"], rotation=45, ha="right")
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def plot_haplotagged_reads(df, out_png):
    x = range(len(df))

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(x, df["hp1_reads"], label="HP1 reads")
    ax.bar(x, df["hp2_reads"], bottom=df["hp1_reads"], label="HP2 reads")

    ax.set_xticks(list(x))
    ax.set_xticklabels(df["sample"], rotation=45, ha="right")
    ax.set_ylabel("Haplotagged reads")
    ax.set_title("Haplotagged read counts")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def plot_genotype_composition(df, out_png):
    cols = [
        "hom_ref_0_0",
        "hom_alt_1_1",
        "unphased_het",
        "phased_het",
        "missing",
    ]

    labels = [
        "0/0",
        "1/1",
        "0/1 unphased",
        "0|1 or 1|0 phased",
        "./.",
    ]

    x = range(len(df))
    bottom = [0] * len(df)

    fig, ax = plt.subplots(figsize=(11, 5))

    for col, label in zip(cols, labels):
        vals = df[col].values
        ax.bar(x, vals, bottom=bottom, label=label)
        bottom = [b + v for b, v in zip(bottom, vals)]

    ax.set_xticks(list(x))
    ax.set_xticklabels(df["sample"], rotation=45, ha="right")
    ax.set_ylabel("Genotype records")
    ax.set_title("Genotype composition after WhatsHap phasing")
    ax.legend(frameon=False, ncol=3)
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Summarize WhatsHap SV phasing and haplotagging outputs."
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Pipeline output root, e.g. /scratch/.../8-batch_output",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory where tables and plots are written.",
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    df = collect_sample_metrics(args.outdir)

    summary_tsv = os.path.join(args.output_dir, "haplotype_phasing_summary.tsv")
    df.to_csv(summary_tsv, sep="\t", index=False)

    save_table_png(
        df,
        os.path.join(args.output_dir, "haplotype_phasing_summary_table.png"),
    )

    plot_all_vs_het_rate(
        df,
        os.path.join(args.output_dir, "haplotype_phasing_rate.png"),
    )

    plot_phase_sets(
        df,
        os.path.join(args.output_dir, "haplotype_phase_sets.png"),
    )

    plot_haplotagged_reads(
        df,
        os.path.join(args.output_dir, "haplotagged_reads.png"),
    )

    plot_genotype_composition(
        df,
        os.path.join(args.output_dir, "haplotype_genotype_composition.png"),
    )

    print(f"Wrote {summary_tsv}")


if __name__ == "__main__":
    main()