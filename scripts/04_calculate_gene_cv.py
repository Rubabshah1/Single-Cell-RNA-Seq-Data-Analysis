#!/usr/bin/env python3
"""Calculate per-gene expression CV from the exported sparse count matrix."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.io import mmread


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def calculate_group_stats(
    raw_counts: sparse.csr_matrix,
    normalized_counts: sparse.csr_matrix,
    gene_names: list[str],
    group_name: str,
) -> list[dict[str, object]]:
    """Return per-gene raw and normalized expression statistics for one group."""
    n_cells = raw_counts.shape[1]
    raw_sum = np.asarray(raw_counts.sum(axis=1)).ravel()
    raw_sum_squares = np.asarray(raw_counts.multiply(raw_counts).sum(axis=1)).ravel()
    norm_sum = np.asarray(normalized_counts.sum(axis=1)).ravel()
    norm_sum_squares = np.asarray(
        normalized_counts.multiply(normalized_counts).sum(axis=1)
    ).ravel()

    def sample_sd(sum_values: np.ndarray, sum_squares: np.ndarray) -> np.ndarray:
        if n_cells < 2:
            return np.full(sum_values.shape, np.nan)
        variance = (sum_squares - np.square(sum_values) / n_cells) / (n_cells - 1)
        return np.sqrt(np.maximum(variance, 0))

    raw_mean = raw_sum / n_cells
    norm_mean = norm_sum / n_cells
    raw_sd = sample_sd(raw_sum, raw_sum_squares)
    norm_sd = sample_sd(norm_sum, norm_sum_squares)
    raw_cv = np.divide(
        raw_sd, raw_mean, out=np.full_like(raw_sd, np.nan), where=raw_mean > 0
    )
    norm_cv = np.divide(
        norm_sd, norm_mean, out=np.full_like(norm_sd, np.nan), where=norm_mean > 0
    )
    pct_detected = np.asarray(normalized_counts.getnnz(axis=1)).ravel() / n_cells * 100

    return [
        {
            "gene_index": index,
            "gene": gene,
            "group": group_name,
            "n_cells": n_cells,
            "mean_normalized": norm_mean[index - 1],
            "sd_normalized": norm_sd[index - 1],
            "cv": norm_cv[index - 1],
            "pct_detected": pct_detected[index - 1],
            "mean_raw": raw_mean[index - 1],
            "sd_raw": raw_sd[index - 1],
            "cv_raw": raw_cv[index - 1],
        }
        for index, gene in enumerate(gene_names, start=1)
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate gene-wise CV on linear counts normalized to a fixed "
            "library size. By default, cells are analyzed together."
        )
    )
    parser.add_argument(
        "--counts",
        type=Path,
        default=PROJECT_ROOT / "outputs" / "counts.mtx",
        help="Gene-by-cell Matrix Market count matrix.",
    )
    parser.add_argument(
        "--genes",
        type=Path,
        default=PROJECT_ROOT / "outputs" / "genes.tsv",
        help="One gene name per matrix row.",
    )
    parser.add_argument(
        "--barcodes",
        type=Path,
        default=PROJECT_ROOT / "outputs" / "barcodes.tsv",
        help="One cell barcode per matrix column.",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=PROJECT_ROOT / "outputs" / "cell_metadata.csv",
        help="Cell metadata CSV, indexed by barcode.",
    )
    parser.add_argument(
        "--group-by",
        help="Optional metadata column; calculate separate CVs within each level.",
    )
    parser.add_argument(
        "--scale-factor",
        type=float,
        default=10_000,
        help="Library-size normalization target (default: 10000).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=PROJECT_ROOT / "outputs" / "gene_cv.csv",
        help="Output CSV path.",
    )
    args = parser.parse_args()
    if not np.isfinite(args.scale_factor) or args.scale_factor <= 0:
        parser.error("--scale-factor must be a positive finite number")
    return args


def main() -> None:
    args = parse_args()
    counts = mmread(args.counts).tocsr().astype(np.float64)
    counts.eliminate_zeros()
    genes = args.genes.read_text(encoding="utf-8").splitlines()
    barcodes = args.barcodes.read_text(encoding="utf-8").splitlines()

    if counts.shape != (len(genes), len(barcodes)):
        raise ValueError(
            "Matrix dimensions do not match genes/barcodes: "
            f"{counts.shape[0]} genes x {counts.shape[1]} cells versus "
            f"{len(genes)} genes x {len(barcodes)} barcodes."
        )
    if len(set(barcodes)) != len(barcodes):
        raise ValueError("Barcodes must be unique to align metadata reliably.")
    if not np.isfinite(counts.data).all() or np.any(counts.data < 0):
        raise ValueError("Counts must be finite and non-negative.")

    library_sizes = np.asarray(counts.sum(axis=0)).ravel()
    zero_library = np.flatnonzero(library_sizes == 0)
    if zero_library.size:
        examples = ", ".join(barcodes[i] for i in zero_library[:5])
        raise ValueError(
            f"{zero_library.size} cells have zero total counts "
            f"(for example: {examples}); filter these cells before calculating CV."
        )

    normalized = counts.multiply(args.scale_factor / library_sizes).tocsr()

    if args.group_by:
        with args.metadata.open(newline="", encoding="utf-8") as metadata_file:
            reader = csv.DictReader(metadata_file)
            fieldnames = reader.fieldnames or []
            if not fieldnames:
                raise ValueError("Metadata CSV has no header row.")
            if args.group_by not in fieldnames:
                raise ValueError(
                    f"Metadata has no column named {args.group_by!r}. "
                    f"Available columns: {', '.join(fieldnames)}"
                )
            barcode_column = fieldnames[0]
            metadata_by_barcode = {}
            for row in reader:
                barcode = row[barcode_column]
                if barcode in metadata_by_barcode:
                    raise ValueError(f"Metadata barcode {barcode!r} is duplicated.")
                metadata_by_barcode[barcode] = row

        missing = [barcode for barcode in barcodes if barcode not in metadata_by_barcode]
        if missing:
            raise ValueError(
                f"Metadata is missing {len(missing)} matrix barcodes "
                f"(for example: {', '.join(missing[:5])})."
            )
        group_values = [metadata_by_barcode[cell][args.group_by] for cell in barcodes]
        if any(value is None or value.strip() == "" for value in group_values):
            raise ValueError(f"Metadata column {args.group_by!r} has missing values.")
    else:
        group_values = ["all_cells"] * len(barcodes)

    cells_by_group: dict[str, list[int]] = {}
    for cell_index, group_name in enumerate(group_values):
        cells_by_group.setdefault(group_name, []).append(cell_index)

    fieldnames = [
        "gene_index",
        "gene",
        "group",
        "n_cells",
        "mean_normalized",
        "sd_normalized",
        "cv",
        "pct_detected",
        "mean_raw",
        "sd_raw",
        "cv_raw",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=fieldnames)
        writer.writeheader()
        for group_name, cell_indices in cells_by_group.items():
            writer.writerows(
                calculate_group_stats(
                    counts[:, cell_indices],
                    normalized[:, cell_indices],
                    genes,
                    group_name,
                )
            )
    print(
        f"Wrote CV for {len(genes)} genes across {len(cells_by_group)} group(s) "
        f"to {args.output}"
    )
    print(
        "CV uses sample SD / mean of linear counts normalized to "
        f"{args.scale_factor:g} counts per cell library."
    )


if __name__ == "__main__":
    main()
