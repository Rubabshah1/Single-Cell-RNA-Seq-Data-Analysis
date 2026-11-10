"""
03_inspect_all_3ca_robust.py

Robust inspector for all 3CA lung studies.

Reads archives directly without extracting:
    <data-root>/<study>/data.tar.gz
    <data-root>/<study>/metadata.tar.gz

Features:
- Case-insensitive archive/member matching
- Automatically detects CSV metadata files
- Reports all metadata columns
- Detects likely patient/donor columns WITHOUT assuming equivalence
- Detects likely tumor/normal/source columns and reports their values
- Detects likely cell-type/subtype columns
- Summarizes annotated cell types other than explicit malignant labels
- Automatically detects Matrix Market (.mtx / .mtx.gz) files
- Previews the first N (default 5) matrix entries (row, col, value)
  WITHOUT loading the full expression matrix
- Correctly reports dimension consistency
- Produces CSV + TXT summaries

Run:
python 02_inspect_3ca_data.py \
  --data-root "/Users/rubabshah/Desktop/lums research/single cell/data/3ca_lung"
"""

import argparse
import gzip
import io
import re
import tarfile
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------
# Matching helpers
# ---------------------------------------------------------------------

def norm(s):
    """Normalize a name for case-insensitive/flexible matching."""
    return re.sub(r"[^a-z0-9]", "", str(s).casefold())


def open_tar(path):
    return tarfile.open(path, mode="r:gz")


def archive_files(path):
    with open_tar(path) as tar:
        return [m.name for m in tar.getmembers() if m.isfile()]


def find_member(tar, filename):
    """Exact basename match, case-insensitive."""
    target = norm(filename)

    matches = [
        m for m in tar.getmembers()
        if m.isfile() and norm(Path(m.name).name) == target
    ]

    return matches[0] if matches else None


def all_members(tar):
    return [m for m in tar.getmembers() if m.isfile()]


# ---------------------------------------------------------------------
# CSV detection
# ---------------------------------------------------------------------

def csv_members(path):
    with open_tar(path) as tar:
        return [
            m for m in all_members(tar)
            if Path(m.name).suffix.casefold() == ".csv"
        ]


def read_csv_member(tar, member):
    f = tar.extractfile(member)
    if f is None:
        raise RuntimeError(f"Could not read {member.name}")
    return pd.read_csv(f, low_memory=False)


def identify_csv_role(member_name, columns):
    """
    Identify whether a CSV looks like cell metadata or sample metadata.
    This is heuristic and is reported explicitly.
    """
    cols = {norm(c) for c in columns}

    cell_score = 0
    sample_score = 0

    cell_indicators = {
        "cellname", "cellid", "barcode", "cellbarcode",
        "celltype", "cellsubtype", "cellannotation",
        "umap1", "umap2"
    }

    sample_indicators = {
        "sample", "sampleid", "sampletype", "technology",
        "platform", "patient", "patientid", "donor", "donorid"
    }

    cell_score = len(cols & cell_indicators)
    sample_score = len(cols & sample_indicators)

    if cell_score > sample_score:
        return "cell-like"
    if sample_score > cell_score:
        return "sample-like"

    # Filename fallback
    name = norm(Path(member_name).stem)
    if "cell" in name:
        return "cell-like"
    if "sample" in name:
        return "sample-like"

    return "unknown"


def load_metadata_from_archives(archives):
    """
    Inspect every CSV and select likely cell/sample metadata.
    Returns:
        cell_df, sample_df, metadata_inventory
    """
    inventory = []
    cell_candidates = []
    sample_candidates = []

    for archive in archives:
        with open_tar(archive) as tar:
            for member in all_members(tar):
                if Path(member.name).suffix.casefold() != ".csv":
                    continue

                try:
                    df = read_csv_member(tar, member)
                    role = identify_csv_role(member.name, df.columns)

                    item = {
                        "archive": archive.name,
                        "file": member.name,
                        "rows": len(df),
                        "columns": len(df.columns),
                        "role": role,
                        "column_names": list(df.columns),
                    }
                    inventory.append(item)

                    if role == "cell-like":
                        cell_candidates.append((df, item))
                    elif role == "sample-like":
                        sample_candidates.append((df, item))

                except Exception as e:
                    inventory.append({
                        "archive": archive.name,
                        "file": member.name,
                        "rows": None,
                        "columns": None,
                        "role": "ERROR",
                        "column_names": [],
                        "error": str(e),
                    })

    # Prefer the largest cell-like CSV and largest sample-like CSV.
    cell_df = None
    cell_info = None
    if cell_candidates:
        cell_df, cell_info = max(cell_candidates, key=lambda x: len(x[0]))

    sample_df = None
    sample_info = None
    if sample_candidates:
        sample_df, sample_info = max(
            sample_candidates, key=lambda x: len(x[0])
        )

    return cell_df, cell_info, sample_df, sample_info, inventory


# ---------------------------------------------------------------------
# Expression matrix detection
# ---------------------------------------------------------------------

MATRIX_EXTENSIONS = (
    ".mtx",
    ".mtx.gz",
)

MATRIX_NAME_HINTS = (
    "matrix",
    "count",
    "counts",
    "expression",
    "expdata",
    "umis",
    "umi",
    "data",
)

# How many matrix data entries to preview by default. Kept as a module
# level constant so it's easy to change in one place.
DEFAULT_MATRIX_PREVIEW_ROWS = 5


def matrix_candidates(path):
    candidates = []

    with open_tar(path) as tar:
        for member in all_members(tar):
            lower = Path(member.name).name.casefold()

            if lower.endswith(MATRIX_EXTENSIONS):
                candidates.append(member)

    return candidates


def matrix_score(member_name):
    name = Path(member_name).name.casefold()

    score = 0
    for hint in MATRIX_NAME_HINTS:
        if hint in name:
            score += 1

    # Prefer uncompressed .mtx if both exist.
    if name.endswith(".mtx"):
        score += 1

    return score


def read_matrix_header(path, n_preview=DEFAULT_MATRIX_PREVIEW_ROWS):
    """
    Reads the header/dimension line of the selected Matrix Market file
    and previews the first `n_preview` data entries (row, col, value)
    that immediately follow it -- WITHOUT loading the full matrix.

    Returns:
        (info_dict_or_None, candidates)

    info_dict contains:
        file, format, n_rows, n_cols, nnz, valid, preview_rows
    where preview_rows is a list of up to `n_preview`
    (row, col, value) string tuples, in file order.
    """
    candidates = matrix_candidates(path)

    if not candidates:
        return None, []

    candidates = sorted(
        candidates,
        key=lambda m: matrix_score(m.name),
        reverse=True
    )

    selected = candidates[0]

    with open_tar(path) as tar:
        f = tar.extractfile(selected)

        if f is None:
            return None, candidates

        if selected.name.casefold().endswith(".gz"):
            # tar.extractfile gives compressed bytes for a .gz member.
            f = gzip.GzipFile(fileobj=io.BytesIO(f.read()))

        raw = f.readline()
        first = raw.decode("utf-8", errors="replace").strip()

        if not first.startswith("%%MatrixMarket"):
            return {
                "file": selected.name,
                "format": first,
                "n_rows": None,
                "n_cols": None,
                "nnz": None,
                "valid": False,
                "preview_rows": [],
            }, candidates

        n_rows = n_cols = nnz = None
        preview_rows = []

        for raw in f:
            line = raw.decode("utf-8", errors="replace").strip()

            if not line or line.startswith("%"):
                continue

            parts = line.split()

            # The first non-comment line after the format header is the
            # dimension line: n_rows n_cols nnz.
            if n_rows is None:
                if len(parts) == 3:
                    try:
                        n_rows, n_cols, nnz = map(int, parts)
                    except ValueError:
                        # Malformed dimension line -- stop trying to parse.
                        break
                continue

            # Every line after that is a data entry: row col value.
            # We only ever read up to n_preview of these, then stop --
            # the rest of the (potentially huge) matrix is never read.
            if len(parts) >= 3 and len(preview_rows) < n_preview:
                preview_rows.append((parts[0], parts[1], parts[2]))

            if len(preview_rows) >= n_preview:
                break

        if n_rows is None:
            return None, candidates

        return {
            "file": selected.name,
            "format": first,
            "n_rows": n_rows,
            "n_cols": n_cols,
            "nnz": nnz,
            "valid": True,
            "preview_rows": preview_rows,
        }, candidates


# ---------------------------------------------------------------------
# Metadata concept detection
# ---------------------------------------------------------------------

PATIENT_PATTERNS = [
    "patient",
    "patientid",
    "patient_id",
    "patientnumber",
    "patientnumberid",
    "donor",
    "donorid",
    "donor_id",
    "subject",
    "subjectid",
    "subject_id",
    "individual",
    "individualid",
]

SOURCE_PATTERNS = [
    "source",
    "sampletype",
    "sample_type",
    "tissue",
    "tissuetype",
    "condition",
    "diseasestatus",
    "disease_status",
    "tumornormal",
    "tumor_normal",
    "tumour",
    "tumor",
]

CELLTYPE_PATTERNS = [
    "celltype",
    "cell_type",
    "cellannotation",
    "cell_annotation",
    "annotation",
    "cellclass",
    "cell_class",
]

SUBTYPE_PATTERNS = [
    "cellsubtype",
    "cell_subtype",
    "subtype",
    "cellsubclass",
    "cell_subclass",
]


def matching_columns(columns, patterns):
    result = []

    for col in columns:
        n = norm(col)

        for pattern in patterns:
            p = norm(pattern)

            if n == p or p in n or n in p:
                result.append(col)
                break

    return list(dict.fromkeys(result))


def format_non_malignant_cell_type_counts(cells, celltype_column):
    """Summarize cell types except labels explicitly identifying malignant cells."""
    malignant_labels = {
        "malignant",
        "malignantcell",
        "malignantcells",
        "tumor",
        "tumour",
        "tumorcell",
        "tumourcell",
        "tumorcells",
        "tumourcells",
        "cancer",
        "cancercell",
        "cancercells",
    }
    labels = cells[celltype_column].astype("string").str.strip()
    normalized_labels = labels.map(norm)
    counts = (
        cells.loc[
            normalized_labels.notna()
            & ~normalized_labels.isin(malignant_labels),
            celltype_column,
        ]
        .dropna()
        .value_counts()
    )
    return "; ".join(
        f"{cell_type} ({count:,})"
        for cell_type, count in counts.items()
    )


def summarize_candidate_columns(df, columns, max_values=15):
    output = {}

    if df is None:
        return output

    for col in columns:
        values = df[col].dropna()

        # Convert to string so mixed types don't cause sorting problems.
        values = values.astype(str).str.strip()
        values = values[values != ""]

        counts = values.value_counts().head(max_values)

        output[col] = {
            "unique": int(values.nunique()),
            "top_values": [(str(k), int(v)) for k, v in counts.items()],
        }

    return output


# ---------------------------------------------------------------------
# Study inspection
# ---------------------------------------------------------------------

def inspect_study(study_dir, n_preview=DEFAULT_MATRIX_PREVIEW_ROWS):

    files = {
        p.name.casefold(): p
        for p in study_dir.iterdir()
        if p.is_file()
    }

    data_archive = files.get("data.tar.gz")
    metadata_archive = files.get("metadata.tar.gz")

    result = {
        "study": study_dir.name,
        "status": "OK",
        "error": "",
        "genes": None,
        "cells": None,
        "samples": None,
        "patients": None,
        "tumor_cells": None,
        "normal_cells": None,
        "non_malignant_cell_type_counts": "",
        "unannotated_cells": None,
        "annotated_cells": None,
        "cell_types": None,
        "cell_subtypes": None,
        "technologies": None,
        "matrix_file": None,
        "matrix_preview": [],
        "matrix_genes": None,
        "matrix_cells": None,
        "nnz": None,
        "density_pct": None,
        "sparsity_pct": None,
        "genes_match": None,
        "cells_match": None,
        "patient_columns": "",
        "source_columns": "",
        "celltype_columns": "",
        "subtype_columns": "",
        "top_cell_types": "",
        "metadata_files": "",
    }

    details = {
        "patient_summary": {},
        "source_summary": {},
        "celltype_summary": {},
        "subtype_summary": {},
        "cell_columns": [],
        "sample_columns": [],
        "metadata_inventory": [],
        "matrix_candidates": [],
    }

    try:
        if data_archive is None:
            raise FileNotFoundError(
                "data.tar.gz not found (case-insensitive search)"
            )

        if metadata_archive is None:
            raise FileNotFoundError(
                "metadata.tar.gz not found (case-insensitive search)"
            )

        archives = [data_archive, metadata_archive]

        # -------------------------
        # Metadata
        # -------------------------
        (
            cells,
            cell_info,
            samples,
            sample_info,
            inventory,
        ) = load_metadata_from_archives(archives)

        details["metadata_inventory"] = inventory

        if cells is None:
            raise FileNotFoundError(
                "Could not identify a cell-like CSV metadata file"
            )

        result["cells"] = len(cells)

        if samples is not None:
            result["samples"] = len(samples)

        details["cell_columns"] = list(cells.columns)

        if samples is not None:
            details["sample_columns"] = list(samples.columns)

        result["metadata_files"] = "; ".join(
            f"{x['archive']}:{x['file']} [{x['role']}]"
            for x in inventory
        )

        # -------------------------
        # Patient detection
        # -------------------------
        patient_cols = matching_columns(
            cells.columns, PATIENT_PATTERNS
        )

        if not patient_cols and samples is not None:
            patient_cols = matching_columns(
                samples.columns, PATIENT_PATTERNS
            )

        result["patient_columns"] = ", ".join(map(str, patient_cols))

        if patient_cols:
            # Use the first detected column only for the headline count.
            # Other detected columns remain visible in the report.
            col = patient_cols[0]
            source_df = cells if col in cells.columns else samples

            values = (
                source_df[col]
                .dropna()
                .astype(str)
                .str.strip()
            )
            values = values[values != ""]

            result["patients"] = int(values.nunique())

            details["patient_summary"] = summarize_candidate_columns(
                source_df, patient_cols
            )

        # -------------------------
        # Tumor / normal detection
        # -------------------------
        source_cols = matching_columns(
            cells.columns, SOURCE_PATTERNS
        )

        if not source_cols and samples is not None:
            source_cols = matching_columns(
                samples.columns, SOURCE_PATTERNS
            )

        result["source_columns"] = ", ".join(map(str, source_cols))

        if source_cols:
            col = source_cols[0]
            source_df = cells if col in cells.columns else samples

            details["source_summary"] = summarize_candidate_columns(
                source_df, source_cols
            )

            # Count tumor/normal only when values can be confidently
            # recognized. Do not label arbitrary values as tumor/normal.
            values = source_df[col].astype("string").str.strip()
            normalized = values.str.casefold()

            tumor_mask = normalized.isin({
                "tumor", "tumour", "tumor tissue",
                "tumour tissue", "cancer", "malignant"
            })

            normal_mask = normalized.isin({
                "normal", "normal tissue",
                "adjacent normal", "nonmalignant",
                "non-malignant"
            })

            # Only count if the source column is in cell metadata.
            if col in cells.columns:
                result["tumor_cells"] = int(tumor_mask.sum())
                result["normal_cells"] = int(normal_mask.sum())

        # -------------------------
        # Cell type
        # -------------------------
        celltype_cols = matching_columns(
            cells.columns, CELLTYPE_PATTERNS
        )

        result["celltype_columns"] = ", ".join(
            map(str, celltype_cols)
        )

        if celltype_cols:
            col = celltype_cols[0]
            details["celltype_summary"] = summarize_candidate_columns(
                cells, celltype_cols
            )

            missing = cells[col].isna()
            result["unannotated_cells"] = int(missing.sum())
            result["annotated_cells"] = int((~missing).sum())
            result["cell_types"] = int(
                cells[col].nunique(dropna=True)
            )

            top = cells[col].value_counts(dropna=True).head(5)

            result["top_cell_types"] = "; ".join(
                f"{idx} ({val:,})"
                for idx, val in top.items()
            )

            result["non_malignant_cell_type_counts"] = (
                format_non_malignant_cell_type_counts(cells, col)
            )

        # -------------------------
        # Cell subtype
        # -------------------------
        subtype_cols = matching_columns(
            cells.columns, SUBTYPE_PATTERNS
        )

        result["subtype_columns"] = ", ".join(
            map(str, subtype_cols)
        )

        if subtype_cols:
            col = subtype_cols[0]

            details["subtype_summary"] = summarize_candidate_columns(
                cells, subtype_cols
            )

            result["cell_subtypes"] = int(
                cells[col].nunique(dropna=True)
            )

        # -------------------------
        # Technology
        # -------------------------
        if samples is not None:
            tech_cols = matching_columns(
                samples.columns,
                ["technology", "platform", "sequencingtechnology"]
            )

            if tech_cols:
                vals = (
                    samples[tech_cols[0]]
                    .dropna()
                    .astype(str)
                    .str.strip()
                    .unique()
                )
                result["technologies"] = ", ".join(sorted(vals))

        # -------------------------
        # Genes
        # -------------------------
        genes = None

        for archive in archives:
            with open_tar(archive) as tar:
                member = find_member(tar, "genes.txt")

                if member is not None:
                    f = tar.extractfile(member)

                    if f is not None:
                        genes = [
                            line.decode(
                                "utf-8", errors="replace"
                            ).strip()
                            for line in f
                            if line.strip()
                        ]
                        break

        if genes is not None:
            result["genes"] = len(genes)

        # -------------------------
        # Matrix (header + first N entries preview)
        # -------------------------
        matrix, candidates = read_matrix_header(
            data_archive, n_preview=n_preview
        )

        if matrix is None:
            # Some datasets may keep the matrix in metadata archive.
            matrix, candidates = read_matrix_header(
                metadata_archive, n_preview=n_preview
            )

        details["matrix_candidates"] = [
            m.name for m in candidates
        ]

        if matrix is None:
            raise FileNotFoundError(
                "No readable Matrix Market (.mtx/.mtx.gz) expression matrix found"
            )

        result["matrix_file"] = matrix["file"]
        result["matrix_preview"] = matrix.get("preview_rows", [])
        result["matrix_genes"] = matrix["n_rows"]
        result["matrix_cells"] = matrix["n_cols"]
        result["nnz"] = matrix["nnz"]

        if matrix["n_rows"] and matrix["n_cols"]:
            total = matrix["n_rows"] * matrix["n_cols"]
            density = matrix["nnz"] / total

            result["density_pct"] = density * 100
            result["sparsity_pct"] = (1 - density) * 100

        result["genes_match"] = (
            result["genes"] is not None
            and result["genes"] == result["matrix_genes"]
        )

        result["cells_match"] = (
            result["cells"] == result["matrix_cells"]
        )

    except Exception as e:
        result["status"] = "ERROR"
        result["error"] = str(e)

    return result, details


# ---------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------

def fmt_int(x):
    if x is None or pd.isna(x):
        return "NA"
    return f"{int(x):,}"


def print_matrix_preview(preview_rows):
    """
    Prints the first N (row, col, value) entries read directly from the
    Matrix Market file, before anything else about the matrix is shown.
    This is a raw preview only -- no matrix is loaded into memory.
    """
    if not preview_rows:
        print("First matrix entries: NONE FOUND")
        return

    print(f"First {len(preview_rows)} matrix entries (row, col, value):")
    for i, (row, col, value) in enumerate(preview_rows, 1):
        print(f"  {i}. row={row}, col={col}, value={value}")


def print_study(result, details, number, total):

    print("-" * 100)
    print(f"[{number}/{total}] {result['study']}")
    print("-" * 100)

    if result["status"] == "ERROR":
        print(f"ERROR: {result['error']}")
        return

    # Matrix entry preview is printed first, ahead of every other
    # matrix/metadata statistic, per request.
    print_matrix_preview(result.get("matrix_preview"))
    print()

    print(f"Genes:             {fmt_int(result['genes'])}")
    print(f"Cells:             {fmt_int(result['cells'])}")
    print(f"Samples:           {fmt_int(result['samples'])}")
    print(f"Patients:          {fmt_int(result['patients'])}")
    print(f"Tumor cells:       {fmt_int(result['tumor_cells'])}")
    print(f"Normal cells:      {fmt_int(result['normal_cells'])}")
    print(
        "Non-malignant cell types (excluding explicit malignant labels): "
        f"{result['non_malignant_cell_type_counts'] or 'NA'}"
    )
    print(f"Unannotated cells: {fmt_int(result['unannotated_cells'])}")
    print(f"Cell types:        {fmt_int(result['cell_types'])}")
    print(f"Cell subtypes:     {fmt_int(result['cell_subtypes'])}")
    print(f"Technology:        {result['technologies'] or 'NA'}")
    print(f"Matrix file:       {result['matrix_file']}")
    print(
        f"Matrix:            "
        f"{fmt_int(result['matrix_genes'])} × "
        f"{fmt_int(result['matrix_cells'])}"
    )
    print(f"Nonzero entries:   {fmt_int(result['nnz'])}")
    print(
        f"Sparsity:          "
        f"{result['sparsity_pct']:.2f}%"
        if result["sparsity_pct"] is not None
        else "Sparsity:          NA"
    )
    print(f"Genes match:       {result['genes_match']}")
    print(f"Cells match:       {result['cells_match']}")
    print(f"Patient columns:   {result['patient_columns'] or 'NONE DETECTED'}")
    print(f"Source columns:    {result['source_columns'] or 'NONE DETECTED'}")
    print(f"Cell-type columns: {result['celltype_columns'] or 'NONE DETECTED'}")
    print(f"Subtype columns:   {result['subtype_columns'] or 'NONE DETECTED'}")
    print(f"Top cell types:    {result['top_cell_types'] or 'NA'}")

    if result["patient_columns"]:
        print("\nPatient-like column values:")
        for col, info in details["patient_summary"].items():
            vals = "; ".join(
                f"{v} ({n:,})"
                for v, n in info["top_values"]
            )
            print(f"  {col}: {info['unique']:,} unique | {vals}")

    if result["source_columns"]:
        print("\nSource/tumor-normal-like column values:")
        for col, info in details["source_summary"].items():
            vals = "; ".join(
                f"{v} ({n:,})"
                for v, n in info["top_values"]
            )
            print(f"  {col}: {info['unique']:,} unique | {vals}")

    if result["matrix_file"]:
        print("\nMatrix candidates:")
        for name in details["matrix_candidates"]:
            print(f"  {name}")


def main():

    parser = argparse.ArgumentParser(
        description="Robustly inspect all 3CA lung studies."
    )

    parser.add_argument(
        "--data-root",
        required=True,
        help="Directory containing the 3CA study folders."
    )

    parser.add_argument(
        "--matrix-preview-rows",
        type=int,
        default=DEFAULT_MATRIX_PREVIEW_ROWS,
        help=(
            "Number of raw (row, col, value) matrix entries to preview "
            "per study, read directly from the .mtx file header "
            f"(default: {DEFAULT_MATRIX_PREVIEW_ROWS})."
        ),
    )

    args = parser.parse_args()

    root = Path(args.data_root).expanduser().resolve()

    if not root.exists():
        raise FileNotFoundError(
            f"Data root does not exist: {root}"
        )

    study_dirs = []

    for p in sorted(root.iterdir()):
        if not p.is_dir():
            continue

        filenames = {
            child.name.casefold()
            for child in p.iterdir()
            if child.is_file()
        }

        if (
            "data.tar.gz" in filenames
            and "metadata.tar.gz" in filenames
        ):
            study_dirs.append(p)

    print("=" * 100)
    print("3CA LUNG — ROBUST ALL-STUDY INSPECTION")
    print("=" * 100)
    print(f"Root: {root}")
    print(f"Studies found: {len(study_dirs)}")
    print()
    print("Case-insensitive file matching: ON")
    print("Archive extraction: OFF")
    print("Full expression matrix loading: OFF")
    print(f"Matrix entry preview per study: {args.matrix_preview_rows}")
    print()

    if not study_dirs:
        print("No study folders containing data.tar.gz + metadata.tar.gz found.")
        return

    results = []
    all_details = {}

    for i, study_dir in enumerate(study_dirs, 1):

        result, details = inspect_study(
            study_dir, n_preview=args.matrix_preview_rows
        )

        results.append(result)
        all_details[result["study"]] = details

        print_study(
            result,
            details,
            i,
            len(study_dirs)
        )

    df = pd.DataFrame(results)

    # -----------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------

    print()
    print("=" * 100)
    print("SUMMARY OF ALL STUDIES")
    print("=" * 100)

    summary_cols = [
        "study",
        "genes",
        "cells",
        "samples",
        "patients",
        "tumor_cells",
        "normal_cells",
        "non_malignant_cell_type_counts",
        "unannotated_cells",
        "cell_types",
        "matrix_genes",
        "matrix_cells",
        "nnz",
        "sparsity_pct",
        "genes_match",
        "cells_match",
        "status",
    ]

    summary = df[summary_cols].copy()

    print(summary.to_string(index=False))

    # -----------------------------------------------------------------
    # Correct dimension mismatch calculation
    # -----------------------------------------------------------------

    successful = df[df["status"] == "OK"].copy()

    if len(successful):

        gene_match = successful["genes_match"].fillna(False).astype(bool)
        cell_match = successful["cells_match"].fillna(False).astype(bool)

        dimension_mismatches = successful[
            ~(gene_match & cell_match)
        ]

        print()
        print("=" * 100)
        print("OVERALL SUMMARY")
        print("=" * 100)

        print(
            f"Successfully inspected: "
            f"{len(successful)}/{len(df)} studies"
        )

        print(
            f"Total cells across studies: "
            f"{successful['cells'].sum():,.0f}"
        )

        print(
            f"Total samples across studies: "
            f"{successful['samples'].sum():,.0f}"
        )

        print(
            f"Total tumor cells detected: "
            f"{successful['tumor_cells'].fillna(0).sum():,.0f}"
        )

        print(
            f"Total normal cells detected: "
            f"{successful['normal_cells'].fillna(0).sum():,.0f}"
        )

        print(
            f"Total unannotated cells: "
            f"{successful['unannotated_cells'].fillna(0).sum():,.0f}"
        )

        print(
            f"Median genes/study: "
            f"{successful['genes'].median():,.0f}"
        )

        print(
            f"Median cells/study: "
            f"{successful['cells'].median():,.0f}"
        )

        print(
            f"Median sparsity: "
            f"{successful['sparsity_pct'].median():.2f}%"
        )

        print(
            f"Studies with true matrix/metadata dimension mismatch: "
            f"{len(dimension_mismatches)}"
        )

        if len(dimension_mismatches):
            print("\nDimension mismatches:")
            print(
                dimension_mismatches[
                    [
                        "study",
                        "genes",
                        "matrix_genes",
                        "cells",
                        "matrix_cells",
                    ]
                ].to_string(index=False)
            )

    errors = df[df["status"] == "ERROR"]

    if len(errors):
        print()
        print("Studies with errors:")
        print(
            errors[["study", "error"]]
            .to_string(index=False)
        )

    # -----------------------------------------------------------------
    # Save machine-readable CSV
    # -----------------------------------------------------------------

    summary_path = root / "3ca_lung_inspection_summary.csv"
    report_path = root / "3ca_lung_inspection_report.txt"

    # matrix_preview holds tuples, which are awkward in a flat CSV --
    # store a compact string version there and keep the structured
    # version only in the in-memory results / text report.
    df_csv = df.copy()
    df_csv["matrix_preview"] = df_csv["matrix_preview"].apply(
        lambda rows: "; ".join(
            f"({r},{c},{v})" for r, c, v in rows
        ) if rows else ""
    )
    df_csv.to_csv(summary_path, index=False)

    # -----------------------------------------------------------------
    # Detailed text report
    # -----------------------------------------------------------------

    with open(report_path, "w", encoding="utf-8") as f:

        f.write("3CA LUNG — ROBUST ALL-STUDY INSPECTION REPORT\n")
        f.write("=" * 100 + "\n")
        f.write(f"Root: {root}\n")
        f.write(f"Studies found: {len(study_dirs)}\n")
        f.write("Archive extraction: OFF\n")
        f.write("Full matrix loading: OFF\n")
        f.write(f"Matrix entry preview per study: {args.matrix_preview_rows}\n\n")

        for result in results:

            f.write("\n" + "-" * 100 + "\n")
            f.write(f"STUDY: {result['study']}\n")
            f.write("-" * 100 + "\n")

            f.write("\nFIRST MATRIX ENTRIES (row, col, value):\n")
            preview_rows = result.get("matrix_preview") or []
            if preview_rows:
                for i, (row, col, value) in enumerate(preview_rows, 1):
                    f.write(f"  {i}. row={row}, col={col}, value={value}\n")
            else:
                f.write("  NONE FOUND\n")

            for key, value in result.items():
                if key == "matrix_preview":
                    continue
                f.write(f"{key}: {value}\n")

            if result["non_malignant_cell_type_counts"]:
                f.write(
                    "\nANNOTATED CELL TYPES EXCLUDING EXPLICIT MALIGNANT LABELS "
                    "(cell type: cell count):\n"
                )
                for entry in result["non_malignant_cell_type_counts"].split("; "):
                    f.write(f"  {entry}\n")

            details = all_details[result["study"]]

            f.write("\nCELL METADATA COLUMNS:\n")
            f.write(", ".join(map(str, details["cell_columns"])) + "\n")

            f.write("\nSAMPLE METADATA COLUMNS:\n")
            f.write(", ".join(map(str, details["sample_columns"])) + "\n")

            f.write("\nMETADATA FILE INVENTORY:\n")
            for item in details["metadata_inventory"]:
                f.write(
                    f"  {item['archive']} :: {item['file']} "
                    f"[{item['role']}] "
                    f"rows={item['rows']} "
                    f"columns={item['columns']}\n"
                )

            f.write("\nMATRIX CANDIDATES:\n")
            for name in details["matrix_candidates"]:
                f.write(f"  {name}\n")

            f.write("\n")

        f.write("\nOVERALL SUMMARY\n")
        f.write("=" * 100 + "\n")
        f.write(df.drop(columns=["matrix_preview"]).to_string(index=False))
        f.write("\n")

    print()
    print("=" * 100)
    print("FILES CREATED")
    print("=" * 100)
    print(f"CSV:  {summary_path}")
    print(f"TXT:  {report_path}")
    print()
    print("Inspection complete. No archive was extracted.")


if __name__ == "__main__":
    main()