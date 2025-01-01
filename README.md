# Breast Cancer ImmunoProfiling Analysis

This project provides a comprehensive analysis pipeline for immunoprofiling of breast cancer using single-cell RNA sequencing (scRNA-Seq) data. The analysis involves preprocessing, immuno-cell profiling, and biomarker identification. We use Seurat for data processing, clustering, and differential gene expression analysis.

## Goals of the Analysis:
- To identify immune cell populations in breast cancer samples.
- To profile immune-related markers and their expression across cancer subtypes.
- To discover novel biomarkers for breast cancer immunotherapy.

## Files:
- `immunoProfiling_analysis.R`: The R script that contains the analysis pipeline.
- `immune_marker_differential_expression.csv`: Differential expression results between immune cell populations.
- `breast_cancer_immunoprofiling_seurat_obj.rds`: Saved Seurat object for future use.

## How to Run the Analysis:
1. Download raw scRNA-Seq data from NCBI SRA (use the `fastq-dump` command).
2. Run the R script (`immunoProfiling_analysis.R`) to perform preprocessing, clustering, and differential expression analysis.

## Requirements:
- R (v4.1.1 or later)
- Seurat package
- Other required libraries: `ggplot2`, `dplyr`, `cowplot`

## License:
MIT License
