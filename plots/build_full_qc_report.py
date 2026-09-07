#!/usr/bin/env python3
"""
Build full GRCh38 pipeline QC report.

Designed for the current no-Truvari workflow.

Modes:
  cohort:
    Standard cohort needLR annotation and population-frequency summaries.

  trio:
    Trio/family needLR comparator annotation. The upstream SV calling,
    cohort merge, force genotyping, confirmation, read QC, and alignment QC
    are the same, but the needLR section summarizes per-family comparator
    outputs instead of cohort needLR annotation plots.

Main chapters:
  1. Run overview
  2. Read QC
  3. Alignment QC
  4. Cohort SV construction
  5. Tool/reference support UpSet
  6. needLR annotation / trio comparator annotation
  7. needLR trio exploration, in trio mode
  8. GRCh38/CHM13 cross-reference confirmation
  9. Native/lifted breakpoint and SVLEN distance
  10. Output files
"""

import os
import gzip
import json
import glob
import base64
import argparse
from datetime import datetime

import pandas as pd


# =============================================================================
# ARGPARSE
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Build full GRCh38 pipeline QC HTML report."
    )

    parser.add_argument("--samples", required=True)
    parser.add_argument("--cohort-vcf", required=True)
    parser.add_argument("--support-table", required=True)
    parser.add_argument("--genotyped-vcf", required=True)
    parser.add_argument("--confirmation-tsv", required=True)
    parser.add_argument("--confirmation-summary", required=True)
    parser.add_argument("--cohort-results-dir", required=True)
    parser.add_argument("--out-html", required=True)
    parser.add_argument("--out-summary", required=True)

    parser.add_argument(
        "--mode",
        choices=["cohort", "trio"],
        default="cohort",
        help="Report mode: cohort needLR annotation or trio needLR comparator."
    )

    parser.add_argument(
        "--needlr-trio-dir",
        default=None,
        help="Directory containing per-family needLR comparator outputs. Used in --mode trio."
    )

    parser.add_argument(
        "--max-table-rows",
        type=int,
        default=12,
        help="Maximum number of rows shown for preview tables."
    )

    return parser.parse_args()


# =============================================================================
# BASIC HELPERS
# =============================================================================

def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def file_exists(path):
    return path is not None and os.path.exists(path)


def relpath(path, start):
    try:
        return os.path.relpath(path, start)
    except Exception:
        return path


def html_escape(x):
    if x is None:
        return ""

    x = str(x)

    return (
        x.replace("&", "&amp;")
         .replace("<", "&lt;")
         .replace(">", "&gt;")
         .replace('"', "&quot;")
    )


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


def safe_int(x, default=0):
    try:
        if x in [None, ".", "", True]:
            return default
        return int(float(str(x).split(",")[0]))
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
    if "TRA" in raw or "BND" in raw:
        return "TRA"

    return raw


def image_to_base64(path):
    if not file_exists(path):
        return None

    with open(path, "rb") as fh:
        encoded = base64.b64encode(fh.read()).decode("utf-8")

    return encoded


def missing_box(title, message):
    return f"""
    <div class="missing">
      <strong>{html_escape(title)}</strong><br>
      {html_escape(message)}
    </div>
    """


def png_figure_html(path, title, caption=None, cohort_results_dir=None):
    if not file_exists(path):
        return missing_box(title, f"Missing figure: {path}")

    encoded = image_to_base64(path)

    if cohort_results_dir:
        path_label = relpath(path, cohort_results_dir)
    else:
        path_label = path

    caption_html = ""
    if caption:
        caption_html = f"<p class='caption'>{html_escape(caption)}</p>"

    return f"""
    <div class="figure-block">
      <h4>{html_escape(title)}</h4>
      <img src="data:image/png;base64,{encoded}" alt="{html_escape(title)}">
      {caption_html}
      <p class="path">Source: {html_escape(path_label)}</p>
    </div>
    """


def link_button(path, label, report_dir):
    if not file_exists(path):
        return missing_box(label, f"Missing file: {path}")

    href = relpath(path, report_dir)

    return f"""
    <p>
      <a class="button-link" href="{html_escape(href)}">
        {html_escape(label)}
      </a>
    </p>
    """


def read_table(path, max_rows=None):
    if not file_exists(path):
        return None

    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        try:
            df = pd.read_csv(path)
        except Exception:
            return None

    if max_rows is not None:
        return df.head(max_rows)

    return df


def scroll_table_html(table_html):
    return f"""
    <div class="table-scroll">
      {table_html}
    </div>
    """


def dataframe_to_html(df, max_rows=12, index=False):
    if df is None or df.empty:
        return "<p class='note'>No table available.</p>"

    display = df.head(max_rows).copy()

    table = display.to_html(
        index=index,
        border=0,
        classes="data-table",
        escape=True
    )

    return scroll_table_html(table)


def compact_table_html(df, max_rows=12, drop_columns=None):
    if df is None or df.empty:
        return "<p class='note'>No table available.</p>"

    d = df.copy()

    if drop_columns:
        d = d.drop(
            columns=[c for c in drop_columns if c in d.columns],
            errors="ignore"
        )

    if len(d) > max_rows:
        d = d.head(max_rows)

    table = d.to_html(
        index=False,
        classes="data-table",
        border=0,
        escape=True
    )

    return scroll_table_html(table)


def collapsed_table_html(title, df, max_rows=20, drop_columns=None):
    if df is None or df.empty:
        return ""

    d = df.copy()

    if drop_columns:
        d = d.drop(
            columns=[c for c in d.columns if c in drop_columns],
            errors="ignore"
        )

    if len(d) > max_rows:
        note = f"<p class='note'>Showing first {max_rows} rows of {len(d):,}.</p>"
        d = d.head(max_rows)
    else:
        note = ""

    table = d.to_html(
        index=False,
        classes="data-table",
        border=0,
        escape=True
    )

    return f"""
    <details>
      <summary>{html_escape(title)}</summary>
      {note}
      {scroll_table_html(table)}
    </details>
    """


def metric_cards(metrics):
    html = ["<div class='metric-grid'>"]

    for key, value in metrics:
        html.append(
            f"""
            <div class="metric-card">
              <div class="metric-value">{html_escape(value)}</div>
              <div class="metric-label">{html_escape(key)}</div>
            </div>
            """
        )

    html.append("</div>")
    return "\n".join(html)


def format_int(x):
    try:
        return f"{int(x):,}"
    except Exception:
        return "NA"


# =============================================================================
# SAMPLE AND VCF SUMMARIES
# =============================================================================

def load_samples(samples_path):
    if not file_exists(samples_path):
        raise FileNotFoundError(samples_path)

    df = pd.read_csv(samples_path, sep="\t")

    required = {"sample_id", "batch_id"}
    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"Sample sheet is missing required columns: {missing}. "
            f"Available columns: {list(df.columns)}"
        )

    return df


def summarize_vcf(vcf_path):
    if not file_exists(vcf_path):
        return {
            "total_records": 0,
            "pass_records": 0,
            "canonical_pass_records": 0,
            "svtype_counts": {},
            "native_supported": 0,
            "lifted_supported": 0,
            "confirmed_native_and_lifted": 0,
            "lifted_only": 0,
            "native_only": 0,
        }

    canonical = {"DEL", "INS", "DUP", "INV"}

    total = 0
    pass_records = 0
    canonical_pass = 0
    svtype_counts = {}

    native_supported = 0
    lifted_supported = 0
    confirmed_native_and_lifted = 0
    lifted_only = 0
    native_only = 0

    with open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            total += 1

            filt = fields[6]
            info = parse_info(fields[7])
            svtype = normalize_svtype(info.get("SVTYPE", "OTHER"))

            if filt == "PASS":
                pass_records += 1

            if filt == "PASS" and svtype in canonical:
                canonical_pass += 1
                svtype_counts[svtype] = svtype_counts.get(svtype, 0) + 1

            native = safe_int(info.get("NATIVE_GRCH38_SUPP"), default=0)
            lifted = safe_int(info.get("LIFTED_CHM13_GRCH38_SUPP"), default=0)

            any_native = ("ANY_NATIVE_GRCH38" in info) or native > 0
            any_lifted = ("ANY_LIFTED_CHM13_GRCH38" in info) or lifted > 0

            if any_native:
                native_supported += 1
            if any_lifted:
                lifted_supported += 1
            if any_native and any_lifted:
                confirmed_native_and_lifted += 1
            elif any_native and not any_lifted:
                native_only += 1
            elif any_lifted and not any_native:
                lifted_only += 1

    return {
        "total_records": total,
        "pass_records": pass_records,
        "canonical_pass_records": canonical_pass,
        "svtype_counts": svtype_counts,
        "native_supported": native_supported,
        "lifted_supported": lifted_supported,
        "confirmed_native_and_lifted": confirmed_native_and_lifted,
        "lifted_only": lifted_only,
        "native_only": native_only,
    }


def read_confirmation_summary(path):
    if not file_exists(path):
        return {}

    try:
        with open(path, "r") as fh:
            return json.load(fh)
    except Exception:
        return {}


def read_metrics_tsv(path):
    df = read_table(path)

    if df is None or df.empty:
        return {}

    if {"metric", "value"}.issubset(df.columns):
        return dict(zip(df["metric"], df["value"]))

    return {}


# =============================================================================
# FILE DISCOVERY
# =============================================================================

def discover_files(base_dir, patterns):
    hits = []

    for pattern in patterns:
        hits.extend(glob.glob(os.path.join(base_dir, pattern), recursive=True))

    hits = sorted(set(hits))
    return hits


def classify_sample_paths(samples_df, workflow_outdir):
    rows = []

    for _, row in samples_df.iterrows():
        sample = str(row["sample_id"])
        batch = str(row["batch_id"])

        sample_dir = os.path.join(workflow_outdir, batch, sample)

        rows.append({
            "batch_id": batch,
            "sample_id": sample,
            "sample_dir": sample_dir,
            "exists": os.path.exists(sample_dir),
        })

    return pd.DataFrame(rows)


def read_qc_patterns():
    return [
        "**/*multiqc*.html",
        "**/*multiqc*.txt",
        "**/*multiqc*.tsv",
        "**/*NanoPlot-report.html",
        "**/*NanoPlot*/*.html",
        "**/*NanoPlot*/*.png",
        "**/*nanoplot*/*.html",
        "**/*nanoplot*/*.png",
        "**/*NanoStats*.txt",
        "**/*nanostats*.txt",
        "**/*LengthvsQuality*.png",
        "**/*Yield_By_Length*.png",
        "**/*HistogramReadlength*.png",
        "**/*WeightedHistogramReadlength*.png",
        "**/*Non_weightedHistogramReadlength*.png",
        "**/*sequencing_summary*.txt",
        "**/*read*q*.png",
        "**/*read*q*.tsv",
        "**/*read*q*.txt",
        "**/*read*quality*.png",
        "**/*read*quality*.tsv",
        "**/*read*quality*.txt",
        "**/*read_qc*.png",
        "**/*read_qc*.tsv",
        "**/*read_qc*.txt",
    ]


def build_read_qc_discovery(samples_df, workflow_outdir, cohort_results_dir):
    patterns = read_qc_patterns()
    read_qc_root = os.path.join(cohort_results_dir, "read_qc")

    rows = []

    for _, row in samples_df.iterrows():
        batch = str(row["batch_id"])
        sample = str(row["sample_id"])

        sample_dir = os.path.join(workflow_outdir, batch, sample)

        hits = []

        if os.path.exists(sample_dir):
            hits.extend(discover_files(sample_dir, patterns))

        if os.path.exists(read_qc_root):
            sample_read_qc_dirs = [
                os.path.join(read_qc_root, "nanoplot", batch, sample),
            ]

            for d in sample_read_qc_dirs:
                if os.path.exists(d):
                    hits.extend(discover_files(d, patterns))

        hits = sorted(set(hits))

        rows.append({
            "batch_id": batch,
            "sample_id": sample,
            "sample_dir_exists": os.path.exists(sample_dir),
            "n_read_qc_files_detected": len(hits),
            "example_files": "; ".join(relpath(x, workflow_outdir) for x in hits[:5]),
        })

    return pd.DataFrame(rows)


def build_global_read_qc_discovery(cohort_results_dir):
    read_qc_dir = os.path.join(cohort_results_dir, "read_qc")
    patterns = read_qc_patterns() + ["**/*.png", "**/*.html", "**/*.txt", "**/*.tsv"]

    hits = discover_files(read_qc_dir, patterns) if os.path.exists(read_qc_dir) else []
    hits = sorted(set(hits))

    rows = []

    for path in hits:
        rows.append({
            "file": relpath(path, cohort_results_dir),
            "exists": os.path.exists(path),
        })

    return pd.DataFrame(rows)


def alignment_qc_patterns():
    return [
        # Current pipeline layout: <outdir>/<batch>/<sample>/02.alignqc/<ref>/...
        "**/02.alignqc/**/*.alfred.tsv.gz",
        "**/02.alignqc/**/*alfred*.tsv.gz",
        "**/02.alignqc/**/*alfred*.tsv",
        "**/02.alignqc/**/*alfred*.txt",
        "**/02.alignqc/**/*alfred*.png",

        "**/02.alignqc/**/*.mosdepth.global.dist.txt",
        "**/02.alignqc/**/*.mosdepth.region.dist.txt",
        "**/02.alignqc/**/*.mosdepth.summary.txt",
        "**/02.alignqc/**/*mosdepth*.txt",
        "**/02.alignqc/**/*.regions.bed.gz",
        "**/02.alignqc/**/*.regions.bed.gz.csi",

        # Optional samtools-style QC outputs
        "**/02.alignqc/**/*flagstat*",
        "**/02.alignqc/**/*idxstats*",
        "**/02.alignqc/**/*samtools*stats*.txt",
        "**/02.alignqc/**/*.bam.stats",
        "**/02.alignqc/**/*.bam.flagstat",
        "**/02.alignqc/**/*.bam.idxstats",

        # Backward-compatible / generic patterns
        "**/*alfred*.png",
        "**/*alfred*.tsv.gz",
        "**/*alfred*.tsv",
        "**/*alfred*.txt",
        "**/*flagstat*",
        "**/*idxstats*",
        "**/*.bam.stats",
        "**/*.bam.flagstat",
        "**/*.bam.idxstats",
        "**/*.bam.mosdepth.summary.txt",
        "**/*mosdepth*summary*.txt",
        "**/*mosdepth*.global.dist.txt",
        "**/*mosdepth*.region.dist.txt",
        "**/*mosdepth*.per-base.bed.gz",
        "**/*mosdepth*.thresholds.bed.gz",
        "**/*.regions.bed.gz",
        "**/*.regions.bed.gz.csi",
        "**/*coverage*.txt",
        "**/*coverage*.tsv",
        "**/*coverage*.png",
        "**/*depth*.txt",
        "**/*depth*.tsv",
        "**/*depth*.png",
        "**/*alignment*q*.txt",
        "**/*alignment*q*.tsv",
        "**/*alignment*q*.png",
        "**/*align*q*.txt",
        "**/*align*q*.tsv",
        "**/*align*q*.png",
        "**/*multiqc*.html",
        "**/*multiqc*.txt",
        "**/*multiqc*.tsv",
    ]

def discover_sample_alignment_qc_files(sample, batch, workflow_outdir, cohort_results_dir):
    patterns = alignment_qc_patterns()

    candidate_roots = [
        os.path.join(workflow_outdir, batch, sample),
        os.path.join(workflow_outdir, "alignment_qc", batch, sample),
        os.path.join(workflow_outdir, "alignment_qc", sample),
        os.path.join(cohort_results_dir, "alignment_qc", batch, sample),
        os.path.join(cohort_results_dir, "alignment_qc", sample),
        os.path.join(cohort_results_dir, "alignment_qc"),
        os.path.join(cohort_results_dir, "align_qc"),
        os.path.join(cohort_results_dir, "alignment"),
    ]

    hits = []

    for root in candidate_roots:
        if os.path.exists(root):
            hits.extend(discover_files(root, patterns))

    # Fallback for layouts where alignment QC files are placed outside a
    # per-sample directory, but still contain the sample ID in the filename/path.
    # This keeps the report robust across old and new pipeline layouts.
    if not hits:
        sample_patterns = [
            f"**/*{sample}*flagstat*",
            f"**/*{sample}*idxstats*",
            f"**/*{sample}*mosdepth*",
            f"**/*{sample}*coverage*",
            f"**/*{sample}*depth*",
            f"**/*{sample}*alfred*",
            f"**/*{sample}*alignment*q*",
            f"**/*{sample}*align*q*",
        ]

        fallback_roots = [
            workflow_outdir,
            cohort_results_dir,
            os.path.join(cohort_results_dir, "alignment_qc"),
            os.path.join(cohort_results_dir, "align_qc"),
        ]

        for root in fallback_roots:
            if os.path.exists(root):
                hits.extend(discover_files(root, sample_patterns))

    hits = sorted(set(x for x in hits if os.path.isfile(x)))
    return hits


def build_alignment_qc_discovery(samples_df, workflow_outdir, cohort_results_dir):
    rows = []

    for _, row in samples_df.iterrows():
        batch = str(row["batch_id"])
        sample = str(row["sample_id"])
        sample_dir = os.path.join(workflow_outdir, batch, sample)

        hits = discover_sample_alignment_qc_files(
            sample=sample,
            batch=batch,
            workflow_outdir=workflow_outdir,
            cohort_results_dir=cohort_results_dir,
        )

        rows.append({
            "batch_id": batch,
            "sample_id": sample,
            "sample_dir": sample_dir,
            "sample_dir_exists": os.path.exists(sample_dir),
            "n_alignment_qc_files_detected": len(hits),
            "example_files": "; ".join(relpath(x, workflow_outdir) for x in hits[:8]),
        })

    return pd.DataFrame(rows)





def discover_mosdepth_summary_tables(workflow_outdir):
    """Collect mosdepth summary rows across samples and references.

    The current alignment-QC workflow writes files such as:
      <outdir>/<batch>/<sample>/02.alignqc/<ref>/<sample>_<ref>.mosdepth.summary.txt

    This helper parses all such files and annotates each row with batch, sample,
    reference and relative file path so that the report shows content, not only
    discovered filenames.
    """
    patterns = [
        "**/02.alignqc/**/*.mosdepth.summary.txt",
        "**/*mosdepth.summary.txt",
    ]

    hits = discover_files(workflow_outdir, patterns) if os.path.exists(workflow_outdir) else []
    rows = []

    for path in sorted(set(hits)):
        parts = path.split(os.sep)

        batch_id = "NA"
        sample_id = "NA"
        reference = "NA"

        try:
            idx = parts.index("02.alignqc")
            reference = parts[idx + 1]
            sample_id = parts[idx - 1]
            batch_id = parts[idx - 2]
        except Exception:
            pass

        try:
            df = pd.read_csv(path, sep="\t")

            if df.empty:
                rows.append({
                    "batch_id": batch_id,
                    "sample_id": sample_id,
                    "reference": reference,
                    "file": relpath(path, workflow_outdir),
                    "status": "empty",
                })
                continue

            # Prefer the total/genome row when present; otherwise keep the first row.
            chosen = df.iloc[0]
            for col in df.columns:
                if col.lower() in {"chrom", "chromosome", "contig"}:
                    total_rows = df[df[col].astype(str).str.lower().isin(["total", "genome"])].copy()
                    if not total_rows.empty:
                        chosen = total_rows.iloc[0]
                    break

            row = chosen.to_dict()
            row = {str(k): v for k, v in row.items()}
            row.update({
                "batch_id": batch_id,
                "sample_id": sample_id,
                "reference": reference,
                "file": relpath(path, workflow_outdir),
                "status": "ok",
            })
            rows.append(row)

        except Exception as exc:
            rows.append({
                "batch_id": batch_id,
                "sample_id": sample_id,
                "reference": reference,
                "file": relpath(path, workflow_outdir),
                "status": f"read_error: {exc}",
            })

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)

    # Move metadata columns first, regardless of mosdepth version/column names.
    front = [c for c in ["batch_id", "sample_id", "reference", "status", "file"] if c in df.columns]
    rest = [c for c in df.columns if c not in front]
    return df[front + rest]


def discover_representative_read_qc_images(cohort_results_dir):
    read_qc_dir = os.path.join(cohort_results_dir, "read_qc")

    patterns = [
        "**/*LengthvsQuality*.png",
        "**/*Yield_By_Length*.png",
        "**/*HistogramReadlength*.png",
        "**/*WeightedHistogramReadlength*.png",
        "**/*Non_weightedHistogramReadlength*.png",
        "**/*NanoPlot*.png",
        "**/*nanoplot*.png",
    ]

    hits = discover_files(read_qc_dir, patterns) if os.path.exists(read_qc_dir) else []
    return sorted(set(hits))[:16]


def discover_representative_alignment_images(workflow_outdir, cohort_results_dir=None):
    patterns = [
        "**/*alfred*.png",
        "**/*coverage*.png",
        "**/*depth*.png",
        "**/*alignment*q*.png",
        "**/*align*q*.png",
    ]

    hits = []

    if workflow_outdir and os.path.exists(workflow_outdir):
        hits.extend(discover_files(workflow_outdir, patterns))

    if cohort_results_dir and os.path.exists(cohort_results_dir):
        for root in [
            os.path.join(cohort_results_dir, "alignment_qc"),
            os.path.join(cohort_results_dir, "align_qc"),
            os.path.join(cohort_results_dir, "alignment"),
        ]:
            if os.path.exists(root):
                hits.extend(discover_files(root, patterns))

    return sorted(set(hits))[:12]


def discover_needlr_trio_outputs(needlr_trio_dir):
    if not needlr_trio_dir or not os.path.exists(needlr_trio_dir):
        return pd.DataFrame()

    rows = []

    family_dirs = sorted(
        d for d in glob.glob(os.path.join(needlr_trio_dir, "*"))
        if os.path.isdir(d)
    )

    for family_dir in family_dirs:
        family_id = os.path.basename(family_dir)

        done_files = glob.glob(
            os.path.join(family_dir, "*.needLR_comparator.done")
        )

        comparator_dir = os.path.join(family_dir, "needLR_comparator")
        vcf_dir = os.path.join(family_dir, "vcfs")

        comparator_files = discover_files(
            comparator_dir,
            ["**/*"]
        ) if os.path.exists(comparator_dir) else []

        comparator_files = [
            x for x in comparator_files
            if os.path.isfile(x)
        ]

        subset_vcfs = discover_files(
            vcf_dir,
            ["*.vcf.gz", "*.vcf.gz.tbi"]
        ) if os.path.exists(vcf_dir) else []

        rows.append({
            "family_id": family_id,
            "family_dir": family_dir,
            "done": len(done_files) > 0,
            "n_subset_vcf_files": len(subset_vcfs),
            "n_comparator_output_files": len(comparator_files),
            "example_comparator_files": "; ".join(
                relpath(x, needlr_trio_dir) for x in comparator_files[:8]
            ),
        })

    return pd.DataFrame(rows)


# =============================================================================
# SECTION RENDERERS
# =============================================================================

def section(title, anchor, body):
    return f"""
    <section id="{html_escape(anchor)}">
      <h2>{html_escape(title)}</h2>
      {body}
    </section>
    """


def render_read_qc_section(
    read_qc_df,
    global_read_qc_df,
    cohort_results_dir,
    report_dir,
):
    if read_qc_df is None or read_qc_df.empty:
        return "<p class='note'>No read QC discovery table was generated.</p>"

    n_samples = len(read_qc_df)
    n_with_files = int((read_qc_df["n_read_qc_files_detected"] > 0).sum())
    total_files = int(read_qc_df["n_read_qc_files_detected"].sum())

    compact_cols = [
        "batch_id",
        "sample_id",
        "sample_dir_exists",
        "n_read_qc_files_detected",
    ]

    compact = read_qc_df[[c for c in compact_cols if c in read_qc_df.columns]].copy()

    html = f"""
    <p>
      Read-level QC file discovery found <b>{total_files:,}</b> files across
      <b>{n_samples:,}</b> samples. Samples with detected read QC files:
      <b>{n_with_files:,}/{n_samples:,}</b>.
    </p>
    """

    if total_files == 0:
        html += """
        <div class="missing">
          No read QC files were detected in the expected sample folders or in
          <code>cohort_results/read_qc</code>. If NanoPlot rules were just
          added, make sure the final report rule depends on the MultiQC output.
        </div>
        """

    html += compact_table_html(compact, max_rows=100)

    if total_files > 0:
        html += collapsed_table_html(
            "Show example read QC file paths by sample",
            read_qc_df,
            max_rows=100,
            drop_columns=None,
        )

    if global_read_qc_df is not None and not global_read_qc_df.empty:
        html += collapsed_table_html(
            "Show all generated read QC files",
            global_read_qc_df,
            max_rows=120,
            drop_columns=None,
        )

    read_qc_images = discover_representative_read_qc_images(cohort_results_dir)

    if read_qc_images:
        html += "<h3>Representative NanoPlot/read QC images</h3>"

        for path in read_qc_images:
            html += png_figure_html(
                path,
                title=relpath(path, cohort_results_dir),
                caption="Automatically discovered read QC image.",
                cohort_results_dir=cohort_results_dir,
            )

    return html


def render_alignment_qc_section(
    alignment_qc_df,
    workflow_outdir,
    cohort_results_dir,
    plot_paths=None,
    tables=None,
    max_rows=20,
):
    if alignment_qc_df is None or alignment_qc_df.empty:
        return "<p class='note'>No alignment QC discovery table was generated.</p>"

    n_samples = len(alignment_qc_df)
    n_with_files = int((alignment_qc_df["n_alignment_qc_files_detected"] > 0).sum())
    total_files = int(alignment_qc_df["n_alignment_qc_files_detected"].sum())

    compact_cols = [
        "batch_id",
        "sample_id",
        "sample_dir_exists",
        "n_alignment_qc_files_detected",
    ]

    compact = alignment_qc_df[
        [c for c in compact_cols if c in alignment_qc_df.columns]
    ].copy()

    html = f"""
    <p>
      Alignment QC file discovery found <b>{total_files:,}</b> files across
      <b>{n_samples:,}</b> samples. Samples with detected alignment QC files:
      <b>{n_with_files:,}/{n_samples:,}</b>.
    </p>
    """

    if total_files == 0:
        html += missing_box(
            "No alignment QC files detected",
            (
                "The report searched the per-sample output folders and common "
                "cohort_results/alignment_qc locations. If alignment QC exists "
                "under a different folder name, add that folder to "
                "discover_sample_alignment_qc_files()."
            )
        )

    html += compact_table_html(compact, max_rows=100)

    html += collapsed_table_html(
        "Show example alignment QC file paths",
        alignment_qc_df,
        max_rows=100,
        drop_columns=None,
    )

    if plot_paths is not None:
        html += """
        <h3>Alignment QC summary plots</h3>
        <details open>
          <summary>Show coverage and mapping example figures</summary>
        """

        html += png_figure_html(
            plot_paths.get("alignment_qc_alfred"),
            "ALFRED mapping and alignment error metrics",
            caption=(
                "Comparison of GRCh38 and CHM13 mapping metrics parsed from "
                "ALFRED outputs. BAM-derived unique mapping is skipped in the "
                "standard report to avoid expensive full-BAM scans."
            ),
            cohort_results_dir=cohort_results_dir,
        )

        html += png_figure_html(
            plot_paths.get("alignment_qc_mosdepth"),
            "mosdepth coverage summary",
            caption=(
                "Coverage summary across samples and references parsed from "
                "mosdepth summary files."
            ),
            cohort_results_dir=cohort_results_dir,
        )

        html += png_figure_html(
            plot_paths.get("alignment_qc_mosdepth_dist"),
            "mosdepth global coverage distribution",
            caption=(
                "Genome-wide coverage distribution parsed from mosdepth "
                "global distribution files."
            ),
            cohort_results_dir=cohort_results_dir,
        )

        html += png_figure_html(
            plot_paths.get("alignment_qc_mosdepth_contig"),
            "mosdepth coverage by contig",
            caption=(
                "Mean coverage by contig for the preserved representative sample "
                "selected from the trio proband when a trio file is provided."
            ),
            cohort_results_dir=cohort_results_dir,
        )

        html += """
        </details>
        """

    if tables is not None:
        alfred_summary = read_table(tables.get("alignment_qc_alfred_summary"))
        mosdepth_summary = read_table(tables.get("alignment_qc_mosdepth_summary"))
        mosdepth_dist = read_table(tables.get("alignment_qc_mosdepth_dist"))

        html += "<h3>Alignment QC summary tables</h3>"

        html += collapsed_table_html(
            "Show ALFRED mapping summary table",
            alfred_summary,
            max_rows=max_rows,
            drop_columns=None,
        )

        html += collapsed_table_html(
            "Show mosdepth coverage summary table",
            mosdepth_summary,
            max_rows=max_rows,
            drop_columns=None,
        )

        html += collapsed_table_html(
            "Show mosdepth global distribution table",
            mosdepth_dist,
            max_rows=max_rows,
            drop_columns=None,
        )

    return html


def render_needlr_trio_section(needlr_trio_dir, max_rows=20):
    trio_df = discover_needlr_trio_outputs(needlr_trio_dir)

    if trio_df is None or trio_df.empty:
        return f"""
        <p>
          This report was generated in <b>trio mode</b>, but no needLR comparator
          family outputs were detected.
        </p>
        {missing_box(
            "No trio needLR outputs found",
            f"Expected trio directory: {needlr_trio_dir}"
        )}
        """

    n_families = len(trio_df)
    n_done = int(trio_df["done"].sum())
    total_comparator_files = int(trio_df["n_comparator_output_files"].sum())

    compact_cols = [
        "family_id",
        "done",
        "n_subset_vcf_files",
        "n_comparator_output_files",
    ]

    compact = trio_df[[c for c in compact_cols if c in trio_df.columns]].copy()

    html = f"""
    <p>
      This run used <b>needLR comparator trio mode</b>. The upstream SV calling,
      cohort merge, and force-genotyping steps are identical to the standard
      cohort workflow, but final needLR annotation is performed per family by
      comparing the proband VCF against maternal and paternal VCFs.
    </p>

    {metric_cards([
        ("Families detected", format_int(n_families)),
        ("Completed comparator runs", f"{n_done:,}/{n_families:,}"),
        ("Comparator output files", format_int(total_comparator_files)),
    ])}

    <h3>Trio comparator output summary</h3>
    {dataframe_to_html(compact, max_rows=100)}

    {collapsed_table_html(
        "Show trio comparator output details",
        trio_df,
        max_rows=max_rows,
        drop_columns=None
    )}

    <table class="mini-table">
      <tr><th>needLR trio directory</th><td>{html_escape(needlr_trio_dir)}</td></tr>
    </table>
    """

    return html


def render_needlr_trio_exploration_section(plot_paths, tables, cohort_results_dir, max_rows=20):
    """
    Render exploratory plots and tables generated from the needLR trio comparator VCF.

    The section is optional and is only displayed in trio mode. Missing files are
    shown as warning boxes so the report remains buildable while the rule is being
    developed or when no candidate table is produced.
    """

    inheritance_counts = read_table(tables.get("trio_exploration_inheritance_counts"))
    annotation_summary = read_table(tables.get("trio_exploration_annotation_summary"))
    high_priority = read_table(tables.get("trio_exploration_high_priority"))
    rare_denovo = read_table(tables.get("trio_exploration_rare_denovo"))
    context_summary = read_table(tables.get("trio_exploration_context_summary"))

    cards = []

    if inheritance_counts is not None and not inheritance_counts.empty:
        if "n_variants" in inheritance_counts.columns:
            total_variants = inheritance_counts["n_variants"].sum()
            cards.append(("Trio comparator SVs", format_int(total_variants)))

        if {"Inheritance", "n_variants"}.issubset(inheritance_counts.columns):
            denovo = inheritance_counts.loc[
                inheritance_counts["Inheritance"].astype(str) == "de_novo",
                "n_variants",
            ]
            if not denovo.empty:
                cards.append(("De novo candidates", format_int(denovo.iloc[0])))

    if high_priority is not None:
        cards.append(("High-priority candidates", format_int(len(high_priority))))

    if rare_denovo is not None:
        cards.append(("Rare annotated de novo", format_int(len(rare_denovo))))

    metric_html = metric_cards(cards) if cards else "<p class='note'>No trio exploration summary metrics available.</p>"

    html = f"""
    <p>
      This section summarizes the trio-specific needLR comparator VCF. Variants are
      stratified by inheritance class, SV type, proband and parental read support,
      control-population frequency, functional annotation, and genomic-context flags.
    </p>

    <p>
      Candidate prioritization is exploratory. Scores are based on inheritance,
      rarity in controls, overlap with genes or disease resources, read support,
      and genomic-context flags. These tables are intended for biological review,
      not as validated clinical classifications.
    </p>

    {metric_html}

    {png_figure_html(
        plot_paths.get("trio_exploration_summary"),
        "Trio inheritance, SVTYPE, read support and annotation burden",
        caption=(
            "Overview of inheritance classes, SVTYPE composition, proband/parental "
            "alt-read fractions, and annotation burden by inheritance class."
        ),
        cohort_results_dir=cohort_results_dir,
    )}

    {png_figure_html(
        plot_paths.get("trio_exploration_population"),
        "Trio population frequency, SV size, candidate priority and genomic context",
        caption=(
            "Control population frequency, SVLEN distribution, exploratory priority "
            "scores, and genomic-context flags for needLR trio comparator variants."
        ),
        cohort_results_dir=cohort_results_dir,
    )}

    <h3>Inheritance and annotation summary</h3>

    <h4>Inheritance class counts</h4>
    {dataframe_to_html(inheritance_counts, max_rows=max_rows)}

    <h4>Annotation summary by inheritance class</h4>
    {dataframe_to_html(annotation_summary, max_rows=max_rows)}

    <h4>Genomic context summary</h4>
    {dataframe_to_html(context_summary, max_rows=max_rows)}

    <h3>Candidate tables</h3>

    {collapsed_table_html(
        "Show high-priority exploratory candidates",
        high_priority,
        max_rows=max_rows,
        drop_columns=None,
    )}

    {collapsed_table_html(
        "Show rare annotated de novo candidates",
        rare_denovo,
        max_rows=max_rows,
        drop_columns=None,
    )}

    <table class="mini-table">
      <tr><th>Exploration directory</th><td>{html_escape(os.path.join(cohort_results_dir, "needlr_trio_exploration"))}</td></tr>
    </table>
    """

    return html



# =============================================================================
# REPORT BUILDING
# =============================================================================

def build_html_report(
    args,
    samples_df,
    sample_paths_df,
    read_qc_df,
    global_read_qc_df,
    alignment_qc_df,
    cohort_summary,
    confirmation_summary,
    plot_paths,
    html_paths,
    tables,
    summary_rows,
):
    cohort_results_dir = args.cohort_results_dir
    workflow_outdir = os.path.dirname(cohort_results_dir.rstrip("/"))
    report_dir = os.path.dirname(args.out_html)

    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    n_samples = samples_df["sample_id"].nunique()
    n_batches = samples_df["batch_id"].nunique()
    carrier_count_label = f"1-{n_samples}"

    overview_metrics = [
        ("Samples", format_int(n_samples)),
        ("Batches", format_int(n_batches)),
        ("Cohort VCF records", format_int(cohort_summary["total_records"])),
        ("PASS records", format_int(cohort_summary["pass_records"])),
        ("Canonical PASS SVs", format_int(cohort_summary["canonical_pass_records"])),
        ("Native and lifted support", format_int(cohort_summary["confirmed_native_and_lifted"])),
        ("Native only", format_int(cohort_summary["native_only"])),
        ("Lifted only", format_int(cohort_summary["lifted_only"])),
    ]

    if args.mode == "trio":
        mode_description = """
        This report was generated in <b>trio mode</b>. The workflow uses the same
        upstream GRCh38/CHM13 SV integration and force-genotyping backbone, but
        replaces standard cohort needLR annotation with per-family needLR
        comparator annotation.
        """
    else:
        mode_description = """
        This report was generated in <b>cohort mode</b>. The workflow uses
        standard needLR cohort annotation and population-frequency summaries.
        """

    overview_body = f"""
    {metric_cards(overview_metrics)}
    <p>
      This report summarizes the GRCh38 structural-variant workflow using the
      integrated no-Truvari confirmation strategy. Cross-reference confirmation
      is derived directly from INFO fields in the integrated GRCh38 cohort VCF:
      <code>NATIVE_GRCH38_SUPP</code>, <code>LIFTED_CHM13_GRCH38_SUPP</code>,
      <code>ANY_NATIVE_GRCH38</code>, and <code>ANY_LIFTED_CHM13_GRCH38</code>.
    </p>

    <p>
      {mode_description}
    </p>

    <table class="mini-table">
      <tr><th>Generated</th><td>{html_escape(now)}</td></tr>
      <tr><th>Mode</th><td>{html_escape(args.mode)}</td></tr>
      <tr><th>Samples file</th><td>{html_escape(args.samples)}</td></tr>
      <tr><th>Cohort results dir</th><td>{html_escape(args.cohort_results_dir)}</td></tr>
      <tr><th>Cohort VCF</th><td>{html_escape(args.cohort_vcf)}</td></tr>
    </table>
    """

    read_qc_body = f"""
    <p>
      This section summarizes NanoPlot and MultiQC read-level QC artifacts.
      Long file-path lists are hidden to keep the report readable.
    </p>
    {render_read_qc_section(
        read_qc_df=read_qc_df,
        global_read_qc_df=global_read_qc_df,
        cohort_results_dir=cohort_results_dir,
        report_dir=report_dir,
    )}
    """

    alignment_qc_body = f"""
    <p>
      This section summarizes alignment-level QC artifacts, including common
      outputs from Alfred, flagstat, idxstats, mosdepth, and coverage summaries
      when present.
    </p>
    {render_alignment_qc_section(alignment_qc_df, workflow_outdir, cohort_results_dir, plot_paths=plot_paths, tables=tables, max_rows=args.max_table_rows)}
    """

    svtype_df = pd.DataFrame(
        [
            {"SVTYPE": k, "count": v}
            for k, v in sorted(cohort_summary["svtype_counts"].items())
        ]
    )

    cohort_body = f"""
    <p>
      The final integrated GRCh38 cohort VCF combines native GRCh38 support and
      CHM13 calls lifted to GRCh38. The summary below is parsed directly from the
      final VCF without bcftools.
    </p>
    {metric_cards([
        ("Total records", format_int(cohort_summary["total_records"])),
        ("PASS records", format_int(cohort_summary["pass_records"])),
        ("Canonical PASS SVs", format_int(cohort_summary["canonical_pass_records"])),
        ("Native supported records", format_int(cohort_summary["native_supported"])),
        ("Lifted supported records", format_int(cohort_summary["lifted_supported"])),
    ])}

    <h3>Canonical SV type counts</h3>
    {dataframe_to_html(svtype_df, max_rows=20)}
    """

    upset_sets = read_table(tables["upset_sets"])

    upset_body = f"""
    <p>
      The UpSet plot summarizes support across native GRCh38 and lifted CHM13
      tool/reference groups. Bar color represents mean cohort sample support.
    </p>
    {png_figure_html(
        plot_paths["upset"],
        "Tool/reference support UpSet",
        cohort_results_dir=cohort_results_dir
    )}

    <h3>Set sizes</h3>
    {dataframe_to_html(upset_sets, max_rows=20)}
    """

    if args.mode == "trio":
        needlr_body = render_needlr_trio_section(
            needlr_trio_dir=args.needlr_trio_dir,
            max_rows=args.max_table_rows,
        )
        needlr_section_title = "needLR trio comparator annotation"
        needlr_toc_label = "needLR trio comparator annotation"
    else:
        needlr_burden = read_table(tables["needlr_burden"])
        needlr_context = read_table(tables["needlr_context"])
        needlr_carriers = read_table(tables["needlr_carriers"])
        popfreq_summary = read_table(tables["needlr_popfreq_summary"])
        carrier_dynamic_summary = read_table(tables["carrier_dynamic_summary"])

        needlr_body = f"""
        <p>
          needLR annotations provide gene, coding, OMIM, genomic-context, cohort
          frequency, and 1KGP population frequency context for the final GRCh38
          genotyped cohort.
        </p>

        {png_figure_html(
            plot_paths["needlr_burden"],
            "needLR annotation burden and support",
            cohort_results_dir=cohort_results_dir
        )}

        <h3>Annotation burden summary</h3>
        {dataframe_to_html(needlr_burden, max_rows=args.max_table_rows)}

        <h3>Genomic context burden</h3>
        {dataframe_to_html(needlr_context, max_rows=args.max_table_rows)}

        <h3>Carrier count distribution</h3>
        {dataframe_to_html(needlr_carriers, max_rows=20)}

        {png_figure_html(
            plot_paths["needlr_popfreq"],
            "needLR control population frequency summary",
            cohort_results_dir=cohort_results_dir
        )}

        <h3>Population frequency summary</h3>
        {dataframe_to_html(popfreq_summary, max_rows=20)}

        {png_figure_html(
            plot_paths["carrier_dynamic_present"],
            "Control population frequency distributions for carrier counts {carrier_count_label}, present variants",
            cohort_results_dir=cohort_results_dir
        )}

        {png_figure_html(
            plot_paths["carrier_dynamic_absent"],
            "Control population frequency distributions for carrier counts {carrier_count_label}, including absent variants",
            cohort_results_dir=cohort_results_dir
        )}

        {png_figure_html(
            plot_paths["carrier_dynamic_summary"],
            "Carrier count {carrier_count_label} population frequency trends",
            cohort_results_dir=cohort_results_dir
        )}

        <h3>Carrier count {carrier_count_label} summary</h3>
        {dataframe_to_html(carrier_dynamic_summary, max_rows=24)}
        """

        needlr_section_title = "needLR annotation and population frequency"
        needlr_toc_label = "needLR annotation and population frequency"

    if args.mode == "trio":
        trio_exploration_body = render_needlr_trio_exploration_section(
            plot_paths=plot_paths,
            tables=tables,
            cohort_results_dir=cohort_results_dir,
            max_rows=args.max_table_rows,
        )
    else:
        trio_exploration_body = ""

    crossref_metrics = read_metrics_tsv(tables["crossref_metrics"])
    crossref_metrics_df = read_table(tables["crossref_metrics"])
    crossref_table = read_table(tables["crossref_table"])

    confirmation_cards = []

    if crossref_metrics:
        for key in [
            "n_grch38_canonical",
            "n_confirmed_by_chm13",
            "n_grch38_only",
            "pct_confirmed_by_chm13",
            "n_integrated_shared_native_and_lifted",
            "n_integrated_lifted_only",
            "n_integrated_native_only",
        ]:
            if key in crossref_metrics:
                value = crossref_metrics[key]

                if "pct" in key:
                    try:
                        value = f"{float(value):.2f}%"
                    except Exception:
                        pass
                else:
                    try:
                        value = f"{int(float(value)):,}"
                    except Exception:
                        pass

                confirmation_cards.append((key, value))

    if not confirmation_cards and confirmation_summary:
        for key, value in confirmation_summary.items():
            if isinstance(value, (int, float, str)):
                confirmation_cards.append((key, value))

    crossref_body = f"""
    <p>
      Confirmation is computed directly from the integrated GRCh38 INFO fields.
      A variant is considered confirmed when it has lifted CHM13-to-GRCh38
      support. This replaces the previous Truvari-based confirmation step.
    </p>

    {metric_cards(confirmation_cards[:8]) if confirmation_cards else "<p class='note'>No cross-reference metrics available.</p>"}

    {png_figure_html(
        plot_paths["crossref_fig1"],
        "Cross-reference confirmation summary",
        cohort_results_dir=cohort_results_dir
    )}

    {png_figure_html(
        plot_paths["crossref_fig2"],
        "Cross-reference confirmation patterns",
        cohort_results_dir=cohort_results_dir
    )}

    {png_figure_html(
        plot_paths["crossref_fig3"],
        "Chromosome-level confirmation",
        cohort_results_dir=cohort_results_dir
    )}

    <h3>Cross-reference metrics</h3>
    {dataframe_to_html(crossref_metrics_df, max_rows=100)}

    {collapsed_table_html(
        "Show confirmation table preview",
        crossref_table,
        max_rows=args.max_table_rows,
        drop_columns=None
    )}
    """

    breakpoint_summary = read_table(tables["svlen_breakpoint_summary"])

    breakpoint_body = f"""
    <p>
      This analysis compares the coordinates and SV lengths reported by native
      GRCh38 and lifted CHM13 support within the same integrated GRCh38 cohort
      record. It is derived from FORMAT-level support information and the
      GRCh38 tool/reference metadata.
    </p>

    {png_figure_html(
        plot_paths["svlen_breakpoint"],
        "Native GRCh38 vs lifted CHM13 support distance",
        cohort_results_dir=cohort_results_dir
    )}

    <h3>Breakpoint and SVLEN distance summary</h3>
    {dataframe_to_html(breakpoint_summary, max_rows=20)}
    """

    output_rows = []

    for label, path in plot_paths.items():
        output_rows.append({
            "category": "figure",
            "name": label,
            "path": relpath(path, cohort_results_dir),
            "exists": os.path.exists(path),
        })

    for label, path in html_paths.items():
        output_rows.append({
            "category": "html",
            "name": label,
            "path": relpath(path, cohort_results_dir),
            "exists": os.path.exists(path),
        })

    for label, path in tables.items():
        output_rows.append({
            "category": "table",
            "name": label,
            "path": relpath(path, cohort_results_dir),
            "exists": os.path.exists(path),
        })

    if args.mode == "trio" and args.needlr_trio_dir:
        output_rows.append({
            "category": "directory",
            "name": "needlr_trio_dir",
            "path": relpath(args.needlr_trio_dir, cohort_results_dir),
            "exists": os.path.exists(args.needlr_trio_dir),
        })

    output_df = pd.DataFrame(output_rows)

    outputs_body = f"""
    <p>
      This table lists the figure, HTML, and table outputs expected by the QC report.
      Large result tables are not fully embedded in the HTML; they are saved as
      TSV files and listed here.
    </p>
    {dataframe_to_html(output_df, max_rows=200)}
    """

    trio_exploration_toc = ""
    if args.mode == "trio":
        trio_exploration_toc = '<li><a href="#trio-exploration">needLR trio exploration</a></li>'

    toc = f"""
    <nav class="toc">
      <h2>Contents</h2>
      <ol>
        <li><a href="#overview">Run overview</a></li>
        <li><a href="#read-qc">Read QC</a></li>
        <li><a href="#alignment-qc">Alignment QC</a></li>
        <li><a href="#cohort">Cohort SV construction</a></li>
        <li><a href="#upset">Tool/reference support UpSet</a></li>
        <li><a href="#needlr">{html_escape(needlr_toc_label)}</a></li>
        {trio_exploration_toc}
        <li><a href="#crossref">GRCh38/CHM13 confirmation</a></li>
        <li><a href="#breakpoint">Breakpoint and SVLEN distance</a></li>
        <li><a href="#outputs">Output files</a></li>
      </ol>
    </nav>
    """

    body_sections = [
        section("Run overview", "overview", overview_body),
        section("Read QC", "read-qc", read_qc_body),
        section("Alignment QC", "alignment-qc", alignment_qc_body),
        section("Cohort SV construction", "cohort", cohort_body),
        section("Tool/reference support UpSet", "upset", upset_body),
        section(needlr_section_title, "needlr", needlr_body),
    ]

    if args.mode == "trio":
        body_sections.append(
            section("needLR trio exploration", "trio-exploration", trio_exploration_body)
        )

    body_sections.extend([
        section("GRCh38/CHM13 confirmation", "crossref", crossref_body),
        section("Breakpoint and SVLEN distance", "breakpoint", breakpoint_body),
        section("Output files", "outputs", outputs_body),
    ])

    body = "\n".join(body_sections)

    css = """
    <style>
      body {
        font-family: Arial, Helvetica, sans-serif;
        margin: 0;
        background: #f6f7f9;
        color: #222;
        line-height: 1.45;
      }

      header {
        background: #1f2937;
        color: white;
        padding: 28px 44px;
      }

      header h1 {
        margin: 0 0 8px 0;
        font-size: 30px;
      }

      header p {
        margin: 0;
        color: #d1d5db;
      }

      main {
        max-width: 1280px;
        margin: 0 auto;
        padding: 28px 32px 60px 32px;
      }

      .toc {
        background: white;
        padding: 20px 26px;
        border-radius: 12px;
        margin-bottom: 28px;
        box-shadow: 0 1px 5px rgba(0,0,0,0.08);
      }

      .toc h2 {
        margin-top: 0;
      }

      .toc a {
        color: #2563eb;
        text-decoration: none;
      }

      section {
        background: white;
        padding: 26px;
        border-radius: 12px;
        margin-bottom: 28px;
        box-shadow: 0 1px 5px rgba(0,0,0,0.08);
      }

      h2 {
        margin-top: 0;
        border-bottom: 2px solid #e5e7eb;
        padding-bottom: 8px;
        color: #111827;
      }

      h3 {
        margin-top: 28px;
        color: #111827;
      }

      h4 {
        margin: 0 0 10px 0;
      }

      code {
        background: #f3f4f6;
        padding: 2px 5px;
        border-radius: 4px;
      }

      .button-link {
        display: inline-block;
        background: #2563eb;
        color: white !important;
        text-decoration: none;
        padding: 10px 14px;
        border-radius: 8px;
        font-weight: bold;
        margin: 12px 0;
      }

      .metric-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(190px, 1fr));
        gap: 14px;
        margin: 18px 0 24px 0;
      }

      .metric-card {
        background: #f9fafb;
        border: 1px solid #e5e7eb;
        border-radius: 10px;
        padding: 16px;
      }

      .metric-value {
        font-size: 24px;
        font-weight: bold;
        color: #111827;
      }

      .metric-label {
        margin-top: 5px;
        color: #4b5563;
        font-size: 13px;
      }

      .figure-block {
        margin: 24px 0;
        padding: 14px;
        border: 1px solid #e5e7eb;
        border-radius: 10px;
        background: #ffffff;
      }

      .figure-block img {
        width: 100%;
        height: auto;
        display: block;
        border: 1px solid #eeeeee;
        border-radius: 6px;
      }

      .caption {
        color: #4b5563;
        font-size: 13px;
      }

      .path {
        color: #6b7280;
        font-size: 12px;
        word-break: break-all;
      }

      .missing {
        margin: 16px 0;
        padding: 14px;
        background: #fff7ed;
        border: 1px solid #fed7aa;
        border-radius: 8px;
        color: #9a3412;
      }

      .note {
        color: #6b7280;
        font-style: italic;
      }

      .table-scroll {
        width: 100%;
        max-width: 100%;
        overflow-x: auto;
        overflow-y: visible;
        margin: 16px 0;
        border: 1px solid #e5e7eb;
        border-radius: 8px;
        background: white;
      }

      .data-table {
        border-collapse: collapse;
        width: max-content;
        min-width: 100%;
        max-width: none;
        font-size: 13px;
        margin: 0;
      }

      .data-table th {
        background: #f3f4f6;
        border: 1px solid #e5e7eb;
        padding: 8px;
        text-align: left;
        position: sticky;
        top: 0;
        white-space: nowrap;
      }

      .data-table td {
        border: 1px solid #e5e7eb;
        padding: 7px;
        vertical-align: top;
        white-space: nowrap;
      }

      .data-table tr:nth-child(even) {
        background: #fafafa;
      }

      .mini-table {
        border-collapse: collapse;
        margin-top: 16px;
        width: 100%;
      }

      .mini-table th {
        width: 210px;
        text-align: left;
        background: #f3f4f6;
      }

      .mini-table th,
      .mini-table td {
        border: 1px solid #e5e7eb;
        padding: 8px;
        word-break: break-all;
      }

      details {
        margin: 16px 0;
        padding: 12px 14px;
        background: #f9fafb;
        border: 1px solid #e5e7eb;
        border-radius: 8px;
      }

      summary {
        cursor: pointer;
        font-weight: bold;
        color: #374151;
      }

      footer {
        color: #6b7280;
        text-align: center;
        padding: 26px;
      }
    </style>
    """

    if args.mode == "trio":
        html_title = "GRCh38 full pipeline QC report - trio mode"
        html_subtitle = (
            "Integrated native GRCh38 and CHM13-to-GRCh38 structural-variant "
            "workflow with needLR trio comparator annotation"
        )
    else:
        html_title = "GRCh38 full pipeline QC report"
        html_subtitle = (
            "Integrated native GRCh38 and CHM13-to-GRCh38 structural-variant workflow"
        )

    html = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>{html_escape(html_title)}</title>
  {css}
</head>
<body>
  <header>
    <h1>{html_escape(html_title)}</h1>
    <p>{html_escape(html_subtitle)}</p>
  </header>
  <main>
    {toc}
    {body}
  </main>
  <footer>
    Generated by build_full_qc_report.py on {html_escape(now)}
  </footer>
</body>
</html>
"""

    return html


# =============================================================================
# MAIN
# =============================================================================

def main():
    args = parse_args()

    os.makedirs(os.path.dirname(args.out_html), exist_ok=True)
    os.makedirs(os.path.dirname(args.out_summary), exist_ok=True)

    cohort_results_dir = args.cohort_results_dir
    workflow_outdir = os.path.dirname(cohort_results_dir.rstrip("/"))

    if args.mode == "trio" and not args.needlr_trio_dir:
        raise ValueError("--needlr-trio-dir is required when --mode trio")

    samples_df = load_samples(args.samples)
    n_samples = samples_df["sample_id"].nunique()
    carrier_count_dir = f"needlr_population_frequency_carriers_1_{n_samples}"
    carrier_count_suffix = f"1_{n_samples}"
    carrier_count_label = f"1-{n_samples}"
    sample_paths_df = classify_sample_paths(samples_df, workflow_outdir)

    read_qc_df = build_read_qc_discovery(
        samples_df=samples_df,
        workflow_outdir=workflow_outdir,
        cohort_results_dir=cohort_results_dir,
    )

    global_read_qc_df = build_global_read_qc_discovery(cohort_results_dir)

    alignment_qc_df = build_alignment_qc_discovery(
        samples_df=samples_df,
        workflow_outdir=workflow_outdir,
        cohort_results_dir=cohort_results_dir,
    )

    cohort_summary = summarize_vcf(args.cohort_vcf)
    confirmation_summary = read_confirmation_summary(args.confirmation_summary)

    html_paths = {
        "read_qc_multiqc": os.path.join(
            cohort_results_dir,
            "read_qc",
            "multiqc",
            "multiqc_report.html"
        ),
    }

    plot_paths = {
        "upset": os.path.join(
            cohort_results_dir,
            "tool_reference_upset_new_cohort",
            "GRCh38_integrated_toolref_upset_mean_sample_frequency.png"
        ),
        "needlr_burden": os.path.join(
            cohort_results_dir,
            "needlr_annotation_plots",
            "needlr_annotation_burden_and_support.png"
        ),
        "needlr_popfreq": os.path.join(
            cohort_results_dir,
            "needlr_annotation_plots",
            "needlr_control_population_frequency_summary.png"
        ),
        "carrier_dynamic_present": os.path.join(
            cohort_results_dir,
            carrier_count_dir,
            f"needlr_popfreq_violin_carrier_counts_{carrier_count_suffix}_present_only.png"
        ),
        "carrier_dynamic_absent": os.path.join(
            cohort_results_dir,
            carrier_count_dir,
            f"needlr_popfreq_violin_carrier_counts_{carrier_count_suffix}_including_absent_variants.png"
        ),
        "carrier_dynamic_summary": os.path.join(
            cohort_results_dir,
            carrier_count_dir,
            f"needlr_popfreq_carrier_counts_{carrier_count_suffix}_summary_lines.png"
        ),
        "crossref_fig1": os.path.join(
            cohort_results_dir,
            "crossref_confirmation_infofield",
            "crossref_fig1_summary.png"
        ),
        "crossref_fig2": os.path.join(
            cohort_results_dir,
            "crossref_confirmation_infofield",
            "crossref_fig2_confirmation_patterns.png"
        ),
        "crossref_fig3": os.path.join(
            cohort_results_dir,
            "crossref_confirmation_infofield",
            "crossref_fig3_chromosome_confirmation.png"
        ),
        "svlen_breakpoint": os.path.join(
            cohort_results_dir,
            "crossref_confirmation_from_integrated_GRCh38",
            "native_vs_lifted_support_distance",
            "native_grch38_vs_lifted_chm13_support_distance.png"
        ),
                "alignment_qc_alfred": os.path.join(
            cohort_results_dir,
            "alignment_qc_plots",
            "alignment_qc_alfred_mapping_unique.png"
        ),
        "alignment_qc_mosdepth": os.path.join(
            cohort_results_dir,
            "alignment_qc_plots",
            "alignment_qc_mosdepth_coverage_summary.png"
        ),
        "alignment_qc_mosdepth_dist": os.path.join(
            cohort_results_dir,
            "alignment_qc_plots",
            "alignment_qc_mosdepth_global_distribution.png"
        ),
        "alignment_qc_mosdepth_contig": os.path.join(
            cohort_results_dir,
            "alignment_qc_plots",
            "alignment_qc_mosdepth_contig_coverage.png"
        ),
    }

    tables = {
        "upset_variant_table": os.path.join(
            cohort_results_dir,
            "tool_reference_upset_new_cohort",
            "toolref_variant_support_table.tsv"
        ),
        "upset_intersections": os.path.join(
            cohort_results_dir,
            "tool_reference_upset_new_cohort",
            "toolref_upset_intersections.tsv"
        ),
        "upset_sets": os.path.join(
            cohort_results_dir,
            "tool_reference_upset_new_cohort",
            "toolref_upset_set_sizes.tsv"
        ),
        "needlr_burden": os.path.join(
            cohort_results_dir,
            "needlr_annotation_plots",
            "needlr_annotation_burden_summary.tsv"
        ),
        "needlr_context": os.path.join(
            cohort_results_dir,
            "needlr_annotation_plots",
            "needlr_genomic_context_burden.tsv"
        ),
        "needlr_carriers": os.path.join(
            cohort_results_dir,
            "needlr_annotation_plots",
            "needlr_carrier_count_distribution.tsv"
        ),
        "needlr_popfreq_summary": os.path.join(
            cohort_results_dir,
            "needlr_annotation_plots",
            "needlr_control_population_frequency_summary.tsv"
        ),
        "carrier_dynamic_summary": os.path.join(
            cohort_results_dir,
            carrier_count_dir,
            f"needlr_popfreq_summary_carrier_counts_{carrier_count_suffix}.tsv"
        ),
        "crossref_table": os.path.join(
            cohort_results_dir,
            "crossref_confirmation_infofield",
            "crossref_infofield_confirmation_table.tsv"
        ),
        "crossref_metrics": os.path.join(
            cohort_results_dir,
            "crossref_confirmation_infofield",
            "crossref_infofield_summary_metrics.tsv"
        ),
        "svlen_breakpoint_summary": os.path.join(
            cohort_results_dir,
            "crossref_confirmation_from_integrated_GRCh38",
            "native_vs_lifted_support_distance",
            "native_grch38_vs_lifted_chm13_support_distance_summary.tsv"
        ),
                "alignment_qc_alfred_summary": os.path.join(
            cohort_results_dir,
            "alignment_qc_plots",
            "alignment_qc_alfred_mapping_summary.tsv"
        ),
        "alignment_qc_mosdepth_summary": os.path.join(
            cohort_results_dir,
            "alignment_qc_plots",
            "alignment_qc_mosdepth_summary.tsv"
        ),
        "alignment_qc_mosdepth_dist": os.path.join(
            cohort_results_dir,
            "alignment_qc_plots",
            "alignment_qc_mosdepth_global_distribution.tsv"
        ),
        "alignment_qc_mosdepth_contig": os.path.join(
            cohort_results_dir,
            "alignment_qc_plots",
            "alignment_qc_mosdepth_contig_coverage.tsv"
        ),
    }

    if args.mode == "trio":
        plot_paths.update({
            "trio_exploration_summary": os.path.join(
                cohort_results_dir,
                "needlr_trio_exploration",
                "needlr_trio_exploration_summary.png"
            ),
            "trio_exploration_population": os.path.join(
                cohort_results_dir,
                "needlr_trio_exploration",
                "needlr_trio_exploration_population_priority.png"
            ),
        })

        tables.update({
            "trio_exploration_variant_table": os.path.join(
                cohort_results_dir,
                "needlr_trio_exploration",
                "needlr_trio_exploration_variant_table.tsv"
            ),
            "trio_exploration_inheritance_counts": os.path.join(
                cohort_results_dir,
                "needlr_trio_exploration",
                "needlr_trio_exploration_inheritance_counts.tsv"
            ),
            "trio_exploration_inheritance_svtype": os.path.join(
                cohort_results_dir,
                "needlr_trio_exploration",
                "needlr_trio_exploration_inheritance_svtype_counts.tsv"
            ),
            "trio_exploration_annotation_summary": os.path.join(
                cohort_results_dir,
                "needlr_trio_exploration",
                "needlr_trio_exploration_annotation_summary.tsv"
            ),
            "trio_exploration_context_summary": os.path.join(
                cohort_results_dir,
                "needlr_trio_exploration",
                "needlr_trio_exploration_context_summary.tsv"
            ),
            "trio_exploration_high_priority": os.path.join(
                cohort_results_dir,
                "needlr_trio_exploration",
                "needlr_trio_exploration_high_priority_candidates.tsv"
            ),
            "trio_exploration_rare_denovo": os.path.join(
                cohort_results_dir,
                "needlr_trio_exploration",
                "needlr_trio_exploration_rare_de_novo_candidates.tsv"
            ),
        })

        # Trio reports use needLR comparator outputs instead of standard cohort needLR plots.
        for key in [
            "needlr_burden",
            "needlr_popfreq",
            "carrier_dynamic_present",
            "carrier_dynamic_absent",
            "carrier_dynamic_summary",
        ]:
            plot_paths.pop(key, None)

        for key in [
            "needlr_burden",
            "needlr_context",
            "needlr_carriers",
            "needlr_popfreq_summary",
            "carrier_dynamic_summary",
        ]:
            tables.pop(key, None)

    summary_rows = []

    summary_rows.append({
        "section": "run",
        "metric": "mode",
        "value": args.mode,
    })

    summary_rows.extend([
        {"section": "overview", "metric": "n_samples", "value": samples_df["sample_id"].nunique()},
        {"section": "overview", "metric": "n_batches", "value": samples_df["batch_id"].nunique()},
        {"section": "cohort", "metric": "total_vcf_records", "value": cohort_summary["total_records"]},
        {"section": "cohort", "metric": "pass_vcf_records", "value": cohort_summary["pass_records"]},
        {"section": "cohort", "metric": "canonical_pass_svs", "value": cohort_summary["canonical_pass_records"]},
        {"section": "cohort", "metric": "native_supported_records", "value": cohort_summary["native_supported"]},
        {"section": "cohort", "metric": "lifted_supported_records", "value": cohort_summary["lifted_supported"]},
        {"section": "cohort", "metric": "native_and_lifted_records", "value": cohort_summary["confirmed_native_and_lifted"]},
        {"section": "cohort", "metric": "native_only_records", "value": cohort_summary["native_only"]},
        {"section": "cohort", "metric": "lifted_only_records", "value": cohort_summary["lifted_only"]},
    ])

    summary_rows.extend([
        {
            "section": "read_qc",
            "metric": "total_read_qc_files_detected_by_sample",
            "value": int(read_qc_df["n_read_qc_files_detected"].sum()) if not read_qc_df.empty else 0,
        },
        {
            "section": "read_qc",
            "metric": "samples_with_read_qc_files",
            "value": int((read_qc_df["n_read_qc_files_detected"] > 0).sum()) if not read_qc_df.empty else 0,
        },
        {
            "section": "read_qc",
            "metric": "total_global_read_qc_files",
            "value": len(global_read_qc_df) if global_read_qc_df is not None else 0,
        },
        {
            "section": "alignment_qc",
            "metric": "total_alignment_qc_files_detected",
            "value": int(alignment_qc_df["n_alignment_qc_files_detected"].sum()) if not alignment_qc_df.empty else 0,
        },
        {
            "section": "alignment_qc",
            "metric": "samples_with_alignment_qc_files",
            "value": int((alignment_qc_df["n_alignment_qc_files_detected"] > 0).sum()) if not alignment_qc_df.empty else 0,
        },
    ])

    if args.mode == "trio":
        trio_df = discover_needlr_trio_outputs(args.needlr_trio_dir)

        summary_rows.extend([
            {
                "section": "needlr_trio",
                "metric": "needlr_trio_dir",
                "value": args.needlr_trio_dir if args.needlr_trio_dir else "NA",
            },
            {
                "section": "needlr_trio",
                "metric": "n_families_detected",
                "value": len(trio_df) if trio_df is not None else 0,
            },
            {
                "section": "needlr_trio",
                "metric": "n_completed_comparator_runs",
                "value": int(trio_df["done"].sum()) if trio_df is not None and not trio_df.empty else 0,
            },
            {
                "section": "needlr_trio",
                "metric": "n_comparator_output_files",
                "value": int(trio_df["n_comparator_output_files"].sum()) if trio_df is not None and not trio_df.empty else 0,
            },
        ])

        trio_inheritance_counts = read_table(tables["trio_exploration_inheritance_counts"])
        trio_high_priority = read_table(tables["trio_exploration_high_priority"])
        trio_rare_denovo = read_table(tables["trio_exploration_rare_denovo"])

        if trio_inheritance_counts is not None and not trio_inheritance_counts.empty:
            if "n_variants" in trio_inheritance_counts.columns:
                summary_rows.append({
                    "section": "needlr_trio_exploration",
                    "metric": "n_trio_exploration_variants",
                    "value": int(trio_inheritance_counts["n_variants"].sum()),
                })

            if {"Inheritance", "n_variants"}.issubset(trio_inheritance_counts.columns):
                for _, row in trio_inheritance_counts.iterrows():
                    summary_rows.append({
                        "section": "needlr_trio_exploration_inheritance",
                        "metric": str(row["Inheritance"]),
                        "value": row["n_variants"],
                    })

        summary_rows.append({
            "section": "needlr_trio_exploration",
            "metric": "n_high_priority_candidates",
            "value": len(trio_high_priority) if trio_high_priority is not None else 0,
        })

        summary_rows.append({
            "section": "needlr_trio_exploration",
            "metric": "n_rare_annotated_de_novo_candidates",
            "value": len(trio_rare_denovo) if trio_rare_denovo is not None else 0,
        })

    for svtype, count in cohort_summary["svtype_counts"].items():
        summary_rows.append({
            "section": "cohort_svtype",
            "metric": svtype,
            "value": count,
        })

    for key, value in confirmation_summary.items():
        if isinstance(value, (int, float, str)):
            summary_rows.append({
                "section": "confirmation_summary_json",
                "metric": key,
                "value": value,
            })

    crossref_metrics = read_metrics_tsv(tables["crossref_metrics"])
    for key, value in crossref_metrics.items():
        summary_rows.append({
            "section": "crossref_infofield",
            "metric": key,
            "value": value,
        })

    for label, path in plot_paths.items():
        summary_rows.append({
            "section": "output_figures",
            "metric": label,
            "value": "present" if os.path.exists(path) else "missing",
        })

    for label, path in html_paths.items():
        summary_rows.append({
            "section": "output_html",
            "metric": label,
            "value": "present" if os.path.exists(path) else "missing",
        })

    for label, path in tables.items():
        summary_rows.append({
            "section": "output_tables",
            "metric": label,
            "value": "present" if os.path.exists(path) else "missing",
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(args.out_summary, sep="\t", index=False)

    html = build_html_report(
        args=args,
        samples_df=samples_df,
        sample_paths_df=sample_paths_df,
        read_qc_df=read_qc_df,
        global_read_qc_df=global_read_qc_df,
        alignment_qc_df=alignment_qc_df,
        cohort_summary=cohort_summary,
        confirmation_summary=confirmation_summary,
        plot_paths=plot_paths,
        html_paths=html_paths,
        tables=tables,
        summary_rows=summary_rows,
    )

    with open(args.out_html, "w") as out:
        out.write(html)

    print("Saved HTML report:")
    print(args.out_html)
    print("Saved summary TSV:")
    print(args.out_summary)


if __name__ == "__main__":
    main()