
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# Project structure
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Actual input data for one study
STUDY_NAME = "Bischoff_et_al._2021"

DATA_DIR = PROJECT_ROOT / "data" / "3ca_lung" / STUDY_NAME / 'Data_Bischoff2021_Lung'
OUTPUT_DIR = PROJECT_ROOT / "outputs"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

MATRIX_PATH = DATA_DIR / "Exp_data_UMIcounts.mtx"
CELLS_PATH = DATA_DIR / "Cells.csv"
SAMPLES_PATH = DATA_DIR / "Samples.csv"
GENES_PATH = DATA_DIR / "Genes.txt"

# Load metadata and gene names
cells = pd.read_csv(CELLS_PATH)
samples = pd.read_csv(SAMPLES_PATH)
genes = pd.read_csv(
    GENES_PATH, header=None, names=["gene"]
)["gene"].astype(str)

# Load genes x cells sparse matrix, then transpose
adata = sc.read_mtx(str(MATRIX_PATH)).transpose()
adata.X = adata.X.tocsr()

print("Matrix shape (cells x genes):", adata.shape)
print("Cell metadata shape:", cells.shape)
print("Sample metadata shape:", samples.shape)
print("Number of genes:", len(genes))

# Verify dimensions
assert adata.n_vars == len(genes), (
    f"Matrix has {adata.n_vars} genes; genes.txt has {len(genes)}"
)
assert adata.n_obs == len(cells), (
    f"Matrix has {adata.n_obs} cells; cells.csv has {len(cells)}"
)

# Assign gene names
adata.var_names = genes.to_numpy()
adata.var_names_make_unique()

# Check the actual cell identifier column
print("Cell metadata columns:", cells.columns.tolist())

# Your example uses cell_name as the barcode identifier
assert "cell_name" in cells.columns

# Ensure the cell metadata is aligned to the matrix column order.
# This assumes cells.csv is already in the same order as the matrix.
# Verify this against the study's barcode file or documentation.
adata.obs_names = cells["cell_name"].astype(str).to_numpy()
adata.obs = cells.set_index("cell_name").copy()

# Preserve raw UMI counts
adata.layers["counts"] = adata.X.copy()

print(adata)
print("Samples:", samples["sample"].nunique())
print("Patients:", samples["patient"].nunique())
print("Cell types:")
print(adata.obs["cell_type"].value_counts().head(10))


# QC metrics calculation
import matplotlib.pyplot as plt

# Identify mitochondrial and ribosomal genes
gene_names = adata.var_names.astype(str).str.upper()

adata.var["mt"] = gene_names.str.startswith("MT-")
adata.var["ribo"] = gene_names.str.startswith(("RPS", "RPL"))

print("Mitochondrial genes:", int(adata.var["mt"].sum()))
print("Ribosomal genes:", int(adata.var["ribo"].sum()))

# Calculate cell-level and gene-level QC metrics
sc.pp.calculate_qc_metrics(
    adata,
    layer="counts",
    qc_vars=["mt", "ribo"],
    percent_top=None,
    log1p=False,
    inplace=True
)

# Inspect cell-level metrics
qc_columns = [
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mt",
    "pct_counts_ribo"
]

print("\nCell-level QC summary:")
print(adata.obs[qc_columns].describe())

# Inspect gene-level metrics
print("\nGene-level QC summary:")
print(
    adata.var[
        ["total_counts", "n_cells_by_counts", "mt", "ribo"]
    ].describe()
)