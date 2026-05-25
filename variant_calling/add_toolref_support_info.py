#!/usr/bin/env python3

import argparse
import gzip
import sys
from collections import OrderedDict


# ==============================================================================
# I/O helpers
# ==============================================================================

def open_text(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_info(info_str):
    info = OrderedDict()

    if info_str in [".", "", None]:
        return info

    for item in info_str.split(";"):
        if item == "":
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True

    return info


def format_info(info):
    parts = []

    for key, value in info.items():
        if value is True:
            parts.append(key)
        elif value is False or value is None:
            continue
        else:
            parts.append(f"{key}={value}")

    return ";".join(parts) if parts else "."


# ==============================================================================
# Metadata
# ==============================================================================

def load_metadata(path):
    rows = []

    with open(path, "r") as fh:
        header = fh.readline().rstrip("\n").split("\t")

        required = {
            "vector_index",
            "set_index",
            "set_label",
            "tool",
            "source",
            "batch",
            "sample",
            "path",
        }

        missing = required - set(header)
        if missing:
            raise RuntimeError(
                f"Metadata file is missing required columns: {', '.join(sorted(missing))}"
            )

        idx = {name: header.index(name) for name in header}

        for line in fh:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            row = {
                "vector_index": int(fields[idx["vector_index"]]),
                "set_index": int(fields[idx["set_index"]]),
                "set_label": fields[idx["set_label"]],
                "tool": fields[idx["tool"]],
                "source": fields[idx["source"]],
                "batch": fields[idx["batch"]],
                "sample": fields[idx["sample"]],
                "path": fields[idx["path"]],
            }

            rows.append(row)

    if not rows:
        raise RuntimeError(f"No rows found in metadata file: {path}")

    metadata_by_vector_index = {}
    set_indices = set()
    samples = []

    for row in rows:
        vi = row["vector_index"]

        if vi in metadata_by_vector_index:
            raise RuntimeError(f"Duplicate vector_index in metadata: {vi}")

        metadata_by_vector_index[vi] = row
        set_indices.add(row["set_index"])

        if row["sample"] not in samples:
            samples.append(row["sample"])

    n_vectors = max(metadata_by_vector_index) + 1
    n_sets = max(set_indices) + 1

    expected_vectors = set(range(n_vectors))
    observed_vectors = set(metadata_by_vector_index)

    if expected_vectors != observed_vectors:
        missing = sorted(expected_vectors - observed_vectors)
        raise RuntimeError(
            f"Metadata vector_index values are not continuous. Missing: {missing[:20]}"
        )

    sample_to_idx = {sample: i for i, sample in enumerate(samples)}

    return rows, metadata_by_vector_index, n_vectors, n_sets, samples, sample_to_idx


# ==============================================================================
# Header annotation
# ==============================================================================

def inject_info_headers(
    header_lines,
    native_info_prefix,
    lifted_info_prefix,
    any_native_flag,
    any_lifted_flag,
):
    new_info_headers = [
        '##INFO=<ID=TOOLREF_SUPP_VEC,Number=1,Type=String,Description="Six-bit vector indicating support by tool/reference group. Set order is defined in the corresponding *_toolref_48way_vcf_metadata.tsv file.">',
        '##INFO=<ID=TOOLREF_SUPP,Number=1,Type=Integer,Description="Number of tool/reference groups supporting this merged variant.">',
        '##INFO=<ID=SAMPLE_SUPP_VEC,Number=1,Type=String,Description="Vector indicating which cohort samples support this merged variant. Sample order is defined by first appearance in the corresponding *_toolref_48way_vcf_metadata.tsv file.">',
        '##INFO=<ID=SAMPLE_SUPP,Number=1,Type=Integer,Description="Number of cohort samples supporting this merged variant.">',
        f'##INFO=<ID={native_info_prefix}_SUPP,Number=1,Type=Integer,Description="Number of native reference tool groups supporting this merged variant.">',
        f'##INFO=<ID={lifted_info_prefix}_SUPP,Number=1,Type=Integer,Description="Number of lifted reference tool groups supporting this merged variant.">',
        f'##INFO=<ID={any_native_flag},Number=0,Type=Flag,Description="Variant has support from at least one native reference tool group.">',
        f'##INFO=<ID={any_lifted_flag},Number=0,Type=Flag,Description="Variant has support from at least one lifted reference tool group.">',
    ]

    existing_ids = set()

    for line in header_lines:
        if line.startswith("##INFO=<ID="):
            try:
                info_id = line.split("ID=", 1)[1].split(",", 1)[0]
                existing_ids.add(info_id)
            except Exception:
                pass

    filtered_new = []
    for line in new_info_headers:
        info_id = line.split("ID=", 1)[1].split(",", 1)[0]
        if info_id not in existing_ids:
            filtered_new.append(line)

    out = []

    inserted = False
    for line in header_lines:
        if line.startswith("#CHROM") and not inserted:
            out.extend(filtered_new)
            inserted = True
        out.append(line)

    if not inserted:
        out.extend(filtered_new)

    return out


# ==============================================================================
# Main annotation logic
# ==============================================================================

def annotate_vcf(
    vcf_path,
    metadata_path,
    out_path,
    native_source,
    lifted_source,
    native_info_prefix,
    lifted_info_prefix,
    any_native_flag,
    any_lifted_flag,
):
    (
        metadata_rows,
        metadata_by_vector_index,
        n_vectors,
        n_sets,
        samples,
        sample_to_idx,
    ) = load_metadata(metadata_path)

    native_set_indices = {
        row["set_index"]
        for row in metadata_rows
        if row["source"] == native_source
    }

    lifted_set_indices = {
        row["set_index"]
        for row in metadata_rows
        if row["source"] == lifted_source
    }

    if not native_set_indices:
        raise RuntimeError(
            f"No metadata rows found with source == {native_source!r}"
        )

    if not lifted_set_indices:
        raise RuntimeError(
            f"No metadata rows found with source == {lifted_source!r}"
        )

    header_lines = []
    body_count = 0

    with open_text(vcf_path, "rt") as inp, open(out_path, "w") as out:
        for line in inp:
            line = line.rstrip("\n")

            if line.startswith("#"):
                header_lines.append(line)
                continue

            if header_lines:
                annotated_header = inject_info_headers(
                    header_lines,
                    native_info_prefix=native_info_prefix,
                    lifted_info_prefix=lifted_info_prefix,
                    any_native_flag=any_native_flag,
                    any_lifted_flag=any_lifted_flag,
                )
                for h in annotated_header:
                    out.write(h + "\n")
                header_lines = []

            fields = line.split("\t")
            if len(fields) < 8:
                continue

            info = parse_info(fields[7])

            supp_vec = info.get("SUPP_VEC")

            if supp_vec is None:
                raise RuntimeError(
                    "Input VCF record does not contain INFO/SUPP_VEC. "
                    "This script expects a SURVIVOR merged VCF."
                )

            supp_vec = str(supp_vec)

            if len(supp_vec) != n_vectors:
                raise RuntimeError(
                    f"SUPP_VEC length {len(supp_vec)} does not match metadata "
                    f"vector count {n_vectors}. Example record: "
                    f"{fields[0]}:{fields[1]}"
                )

            toolref_vec = ["0"] * n_sets
            sample_vec = ["0"] * len(samples)

            for vector_index, bit in enumerate(supp_vec):
                if bit != "1":
                    continue

                meta = metadata_by_vector_index[vector_index]

                set_index = meta["set_index"]
                sample = meta["sample"]

                toolref_vec[set_index] = "1"
                sample_vec[sample_to_idx[sample]] = "1"

            toolref_supp = sum(1 for x in toolref_vec if x == "1")
            sample_supp = sum(1 for x in sample_vec if x == "1")

            native_supp = sum(
                1
                for idx, bit in enumerate(toolref_vec)
                if bit == "1" and idx in native_set_indices
            )

            lifted_supp = sum(
                1
                for idx, bit in enumerate(toolref_vec)
                if bit == "1" and idx in lifted_set_indices
            )

            # Remove previous versions if re-annotating
            for key in [
                "TOOLREF_SUPP_VEC",
                "TOOLREF_SUPP",
                "SAMPLE_SUPP_VEC",
                "SAMPLE_SUPP",
                f"{native_info_prefix}_SUPP",
                f"{lifted_info_prefix}_SUPP",
                any_native_flag,
                any_lifted_flag,
            ]:
                if key in info:
                    del info[key]

            info["TOOLREF_SUPP_VEC"] = "".join(toolref_vec)
            info["TOOLREF_SUPP"] = str(toolref_supp)
            info["SAMPLE_SUPP_VEC"] = "".join(sample_vec)
            info["SAMPLE_SUPP"] = str(sample_supp)
            info[f"{native_info_prefix}_SUPP"] = str(native_supp)
            info[f"{lifted_info_prefix}_SUPP"] = str(lifted_supp)

            if native_supp > 0:
                info[any_native_flag] = True

            if lifted_supp > 0:
                info[any_lifted_flag] = True

            fields[7] = format_info(info)
            out.write("\t".join(fields) + "\n")

            body_count += 1

        # Handle header-only files
        if header_lines:
            annotated_header = inject_info_headers(
                header_lines,
                native_info_prefix=native_info_prefix,
                lifted_info_prefix=lifted_info_prefix,
                any_native_flag=any_native_flag,
                any_lifted_flag=any_lifted_flag,
            )
            for h in annotated_header:
                out.write(h + "\n")

    if body_count == 0:
        raise RuntimeError(f"No variant records were written from {vcf_path}")

    sys.stderr.write(f"Wrote {body_count:,} annotated records to {out_path}\n")
    sys.stderr.write(f"Native source: {native_source}\n")
    sys.stderr.write(f"Lifted source: {lifted_source}\n")
    sys.stderr.write(f"Samples: {','.join(samples)}\n")


# ==============================================================================
# CLI
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Annotate a SURVIVOR 48-way merged VCF with tool/reference support "
            "and sample support using a metadata table."
        )
    )

    parser.add_argument("--vcf", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out", required=True)

    parser.add_argument("--native-source", default="native_grch38")
    parser.add_argument("--lifted-source", default="lifted_chm13_to_grch38")

    parser.add_argument("--native-info-prefix", default="NATIVE_GRCH38")
    parser.add_argument("--lifted-info-prefix", default="LIFTED_CHM13_GRCH38")

    parser.add_argument("--any-native-flag", default="ANY_NATIVE_GRCH38")
    parser.add_argument("--any-lifted-flag", default="ANY_LIFTED_CHM13_GRCH38")

    args = parser.parse_args()

    annotate_vcf(
        vcf_path=args.vcf,
        metadata_path=args.metadata,
        out_path=args.out,
        native_source=args.native_source,
        lifted_source=args.lifted_source,
        native_info_prefix=args.native_info_prefix,
        lifted_info_prefix=args.lifted_info_prefix,
        any_native_flag=args.any_native_flag,
        any_lifted_flag=args.any_lifted_flag,
    )


if __name__ == "__main__":
    main()