from pathlib import Path
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# Project paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent
STUDY_NAME = "Bischoff_et_al._2021"

DATA_DIR = PROJECT_ROOT / "data" / "3ca_lung" / STUDY_NAME
OUTPUT_DIR = PROJECT_ROOT / "outputs" / STUDY_NAME
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Assumes adata has already been loaded and metadata attached.
# Preserve the original raw counts in adata.layers["counts"].

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