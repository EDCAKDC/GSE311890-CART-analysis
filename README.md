# GSE311890 CAR-T Single-Cell Analysis

Reproducible analysis code for the **GSE311890** CAR-T single-cell study, with a focus on long-term CAR-T persistence, transcriptional state, TCR clonotype dynamics, regulon activity, and external reference comparison.

> **Associated publication:** *Nature Medicine* (2026)  
> **DOI:** [10.1038/s41591-026-04578-1](https://doi.org/10.1038/s41591-026-04578-1)

## Overview

This repository contains analysis scripts used to study CAR-T cells at peak expansion and at long-term follow-up (year 9.3). The workflow includes:

- scRNA-seq preprocessing and Seurat object preparation
- Harmony integration, dimensional reduction, clustering, and visualization
- CAR-positive versus CAR-negative T-cell comparisons
- TCR clonotype and repertoire-related analyses
- cell-cycle-aware differential expression
- pathway and signature analysis
- pySCENIC regulon analysis
- external CITE-seq similarity analysis

The repository is intended to provide **transparent, reproducible computational analysis code**. Large raw and processed data objects are not duplicated here.

## Data Availability

Raw and processed data are available through GEO under accession **GSE311890**.

Large input files are excluded from this repository because of file-size and data-source restrictions. If processed GEO objects have already been downloaded, the initial object-construction sections of the scripts can be skipped and downstream analysis can begin from the processed `.rds` objects.

## Repository Structure

```text
GSE311890-CART-analysis/
├── README.md
├── scripts/
│   ├── 1_9y.R
│   ├── 2_peak.R
│   ├── 3_9y and peak car.R
│   ├── 4_9yclean.R
│   └── 5_PT1_PT2_CITEseq_G1_similarity.R
└── pyscenic_out_G1_C4/
    ├── C4_regulon_DEG.tsv
    └── regulons.csv
```

## Analysis Modules

### 1. Year 9.3 processing — `1_9y.R`

Processes the year 9.3 T-cell / CAR-T single-cell object, including object preparation, normalization, dimensional reduction, annotation checking, and preparation for downstream analyses.

### 2. Peak expansion processing — `2_peak.R`

Processes the peak-expansion T-cell / CAR-T object and prepares the data for comparison with long-term CAR-positive cells.

### 3. Peak vs. year 9.3 CAR-T analysis — `3_9y and peak car.R`

Merges peak and year 9.3 CAR-positive cells and performs:

- normalization and dimensional reduction
- Harmony integration
- clustering and visualization
- cell-cycle scoring
- differential expression analysis
- clonotype-related analysis
- signature scoring
- CAR-T state visualization

### 4. Year 9.3 annotation refinement — `4_9yclean.R`

Refines the long-term CAR-T annotation and performs:

- removal/reannotation of the cytotoxic NK-like γδ T-cell cluster
- G1-specific CAR-positive vs. CAR-negative comparison
- differential expression analysis
- over-representation analysis
- signature scoring
- heatmap visualization

### 5. External CITE-seq similarity analysis — `5_PT1_PT2_CITEseq_G1_similarity.R`

Compares the GSE311890 CAR-positive G1 reference with an external PT1Y9 CITE-seq reference using logistic-regression-based prediction and probability/logit heatmaps.

This script requires:

```r
obj2 <- readRDS("PT1Y9_filtered_with_AUCell.RDS")
```

The external object is not included because of size and source restrictions. It should be obtained from the original data source.

## pySCENIC Regulon Analysis

The `pyscenic_out_G1_C4/` directory contains selected outputs used for the C4 regulon analysis:

```text
pyscenic_out_G1_C4/
├── C4_regulon_DEG.tsv
└── regulons.csv
```

- `C4_regulon_DEG.tsv`: differential regulon-activity results
- `regulons.csv`: pySCENIC regulon output used to recover target genes

Large matrices, database files, and intermediate pySCENIC files are intentionally excluded.

## Suggested Run Order

```text
1_9y.R
2_peak.R
3_9y and peak car.R
4_9yclean.R
5_PT1_PT2_CITEseq_G1_similarity.R
```

## Main R Dependencies

```text
Seurat
dplyr
ggplot2
openxlsx
clusterProfiler
ReactomePA
ComplexHeatmap
pheatmap
glmnet
DropletUtils
AnnotationDbi
org.Hs.eg.db
readr
stringr
tidyr
purrr
igraph
tidygraph
ggraph
```

## Files Not Included

Examples of excluded large files:

```text
*.rds
*.RDS
*.h5
*.loom
*.gz
large matrix files
raw Cell Ranger output folders
pySCENIC database files
```

File paths in the scripts may need to be adapted to the local working directory.

## Citation

If this repository is useful for reproducing or extending the associated analysis, please cite the related publication using DOI **10.1038/s41591-026-04578-1** and the GEO accession **GSE311890**.
