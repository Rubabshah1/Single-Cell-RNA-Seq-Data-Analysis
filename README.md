# Single-cell variability analysis project

This repository contains a single-cell RNA-seq workflow for cancer datasets, with two separate fetch pipelines:

1. a 3CA scRNA fetcher for patient/tissue datasets, and
2. a cell-line dataset fetcher for cancer cell line atlases.

The project then inspects the archived metadata and matrix structure and runs a variability-analysis workflow on the resulting count matrix.

---

## What this project does

### 1. scRNA data acquisition
The scRNA fetcher in [scripts/01a_fetch_scrna.py](scripts/01a_fetch_scrna.py) downloads 3CA tissue/patient datasets.

It:
- maps an organ name to a 3CA category slug,
- scrapes the category page for study rows,
- finds the data and metadata download links,
- downloads the paired archives for each study,
- saves them under a study-specific folder in the output directory.

Example usage:
- `python scripts/01a_fetch_scrna.py --organ lung`
- `python scripts/01a_fetch_scrna.py --organ breast --study "Chen"`
- `python scripts/01a_fetch_scrna.py --organ lung --list-3ca`

### 2. Cell-line data acquisition
The cell-line fetcher in [scripts/01b_fetch_cellline.py](scripts/01b_fetch_cellline.py) downloads cancer cell-line datasets that are not consistently available through 3CA.

It supports:
- Broad Single Cell Portal downloads via SCP study accessions and auth codes,
- GEO-based downloads for registry entries such as breast cancer cell-line atlases.

Example usage:
- `python scripts/01b_fetch_cellline.py --list`
- `python scripts/01b_fetch_cellline.py --dataset breast_gambardella2022`
- `python scripts/01b_fetch_cellline.py --dataset pancancer_kinker2020 --scp-auth-code XXXXXX`

### 3. 3CA archive inspection
The inspection script in [scripts/02_inspect_3ca_data.py](scripts/02_inspect_3ca_data.py) reads `.tar.gz` archives directly without extracting them.

It checks for:
- metadata CSV files,
- likely patient/donor columns,
- source / tumor-vs-normal annotations,
- cell-type and subtype columns,
- counts cell-type annotations other than explicit malignant labels,
- Matrix Market expression files,
- dimension consistency between metadata and matrix files.

The results are written to:
- [data/3ca_lung/3ca_lung_inspection_summary.csv](data/3ca_lung/3ca_lung_inspection_summary.csv)
- [data/3ca_lung/3ca_lung_inspection_report.txt](data/3ca_lung/3ca_lung_inspection_report.txt)

### 3. Variability benchmarking
The R script [02_run_scVar.R](02_run_scVar.R) loads the exported count matrix and metadata, defines a grouping variable, and runs several variability metrics using scVar.

It includes a representative comparison across metrics such as:
- `SD`
- `CV`
- `scran`
- `Seurat_vst`

It then:
- combines results into a single gene-by-metric table,
- writes [outputs/variability_scores.csv](outputs/variability_scores.csv),
- computes inter-metric Spearman correlations,
- evaluates a ribosomal-gene stability check.

### 4. Per-gene coefficient of variation
The standalone pipeline in [scripts/04_calculate_gene_cv.py](scripts/04_calculate_gene_cv.py) computes per-gene expression statistics from the exported sparse matrix:

```bash
python scripts/04_calculate_gene_cv.py
```

By default, it calculates CV across all cells and writes [outputs/gene_cv.csv](outputs/gene_cv.csv). For a within-cell-type calculation, pass the metadata column:

```bash
python scripts/04_calculate_gene_cv.py --group-by cell_type --output outputs/gene_cv_by_cell_type.csv
```

CV is calculated as sample standard deviation divided by mean on linear counts normalized per cell to 10,000 total counts. The output also includes raw-count statistics and the percent of cells detecting each gene. CV is undefined for a gene with mean expression of zero (reported as `NaN`). Raw counts have additional library-size dependence, so use the normalized `cv` column for comparisons. Grouping by cell type (or another suitable metadata field) avoids combining distinct cell populations into one variability estimate.

The script requires Python with NumPy and SciPy installed.

---

## Repository structure

- [scripts/01a_fetch_scrna.py](scripts/01a_fetch_scrna.py)  
  downloads 3CA patient/tissue scRNA datasets.

- [scripts/01b_fetch_cellline.py](scripts/01b_fetch_cellline.py)  
  downloads cancer cell-line datasets from SCP and GEO.

- [scripts/02_inspect_3ca_data.py](scripts/02_inspect_3ca_data.py)  
  robust archive inspection and reporting tool.

- [02_run_scVar.R](02_run_scVar.R)  
  variability benchmarking workflow in R.

- [data/3ca_lung](data/3ca_lung)  
  lung study archives and inspection output.

- [outputs](outputs)  
  processed matrix and metadata outputs used by the R analysis.

- [README.md](README.md)  
  project overview and usage notes.

---

## Current generated outputs

The repository currently includes or is expected to generate these outputs:

- [outputs/counts.mtx](outputs/counts.mtx) — sparse count matrix
- [outputs/genes.tsv](outputs/genes.tsv) — gene names
- [outputs/barcodes.tsv](outputs/barcodes.tsv) — cell barcodes
- [outputs/cell_metadata.csv](outputs/cell_metadata.csv) — cell-level metadata
- [outputs/variability_scores.csv](outputs/variability_scores.csv) — variability score table across metrics
- [outputs/gene_cv.csv](outputs/gene_cv.csv) — per-gene raw and library-size-normalized CV statistics
- [data/3ca_lung/3ca_lung_inspection_summary.csv](data/3ca_lung/3ca_lung_inspection_summary.csv) — aggregate summary across lung studies
- [data/3ca_lung/3ca_lung_inspection_report.txt](data/3ca_lung/3ca_lung_inspection_report.txt) — detailed per-study report

---

## Current status

The codebase has already completed the following work:

- 3CA scRNA fetch pipeline implemented in [scripts/01a_fetch_scrna.py](scripts/01a_fetch_scrna.py).
- cell-line fetch pipeline implemented in [scripts/01b_fetch_cellline.py](scripts/01b_fetch_cellline.py).
- 3CA lung archive inspection implemented and run for the available lung studies.
- summary and report files generated for the 3CA lung dataset collection.
- scVar-based variability workflow prepared and written in R.

---

## Planned next steps

- use the fetch script to pull additional 3CA categories beyond lung,
- inspect additional studies in the same format,
- refine the grouping variable in [02_run_scVar.R](02_run_scVar.R) for tumor-vs-normal or cell-state comparisons,
- add plotting and downstream comparisons for the variability score matrix,
- expand the benchmarking with additional control-gene checks.

---

## Notes

This project is a research workflow rather than a packaged software tool. The scripts are intentionally simple and modular so they can be reused for new 3CA categories and different single-cell variability analyses.
