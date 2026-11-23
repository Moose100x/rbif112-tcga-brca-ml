# TCGA-BRCA RNA-seq: Exploratory Analysis and Classification (RBIF 112)

An R Markdown analysis of TCGA breast cancer RNA-seq data: 226 samples (113 primary tumor and 113 solid tissue normal) using STAR gene counts for protein-coding genes.

Course project, Brandeis University M.S. in Bioinformatics, RBIF 112.

## Contents

| File | Purpose |
|---|---|
| `data_acquisition.R` | Queries and downloads TCGA-BRCA STAR counts from the GDC with TCGAbiolinks |
| `tcga_brca_ml.Rmd` | Main analysis, which knits to an HTML report |

## Analysis

- **Preparation:** subset to protein-coding genes, pick tumor and normal samples, compare normalization methods (including quantile normalization with EDASeq), and log2-transform.
- **QC and exploration:** library sizes, expression distributions, missing-value checks, and the most and least variable genes.
- **Dimensionality reduction:** PCA on high-variance genes, variance explained, the number of PCs needed for 50% and 75% of variance, and 2D/3D plots colored by phenotype.
- **Association testing:** gene-level Spearman correlations with age, using one tumor sample per patient.
- **Classification:** logistic regression, support vector machines (`e1071`), random forests (`randomForest`), and a neural network (`neuralnet`), compared with cross-validation, along with feature selection.

## Reproduce

1. Run `data_acquisition.R` to download TCGA-BRCA STAR counts. Save the prepared `SummarizedExperiment` as `TCGA_BRCA_SummarizedExperiment.rds` in this folder.
2. Knit `tcga_brca_ml.Rmd` in RStudio. Figures are written to `HW5_figures/`.

**R packages:** TCGAbiolinks, EDASeq, SummarizedExperiment, matrixStats, dplyr, scatterplot3d, randomForest, e1071, neuralnet.

The data is public TCGA data from the NCI Genomic Data Commons. It is not included because of its size.
