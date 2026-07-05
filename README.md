# GSE311890 CART Analysis

This repository contains analysis scripts for the GSE311890 CAR-T single-cell project, including year 9.3 T-cell/CAR-T processing, peak vs year 9.3 CAR+ comparison, G1-specific analysis, pySCENIC regulon analysis, external CITE-seq similarity analysis, and downstream visualization.

## Data availability

Large input files are not included in this repository because of file size limitations.

Raw and processed data files are available from GEO under accession **GSE311890**.

Some scripts contain an initial section for constructing Seurat objects from raw or intermediate files. If the processed GEO files have already been downloaded, this object-construction section can be skipped, and the downstream analysis can be started directly from the processed `.rds` objects.

## Repository structure

~~~text
GSE311890-CART-analysis/
├── README.md
├── scripts/
│   ├── 1_9y.R
│   ├── 2_peak.R
│   ├── 3_9y and peak car.R
│   ├── 4_9yclean.R
│   └── 5_PT1_PT2_CITEseq_G1_similarity.R
│
└── pyscenic_out_G1_C4/
    ├── C4_regulon_DEG.tsv
    └── regulons.csv
~~~

## Script description

### `1_9y.R`

Processes the year 9.3 T-cell / CAR-T single-cell object.

This script includes Seurat object preparation, normalization, dimensional reduction, annotation checking, and preparation of downstream analysis objects.

The initial object-construction section can be skipped if the processed GEO object has already been downloaded.

### `2_peak.R`

Processes the peak expansion single-cell object.

This script prepares the peak CAR+ / T-cell data for downstream comparison with the year 9.3 CAR+ population.

The initial object-construction section can be skipped if the processed GEO object has already been downloaded.

### `3_9y and peak car.R`

Merges peak and year 9.3 CAR+ cells.

This script performs normalization, dimensional reduction, Harmony integration, clustering, cell-cycle scoring, differential expression analysis, clonotype-related analysis, signature scoring, and visualization of CAR+ cell states.

The initial object-construction section can be skipped if the processed GEO object has already been downloaded.

### `4_9yclean.R`

Cleans and updates the year 9.3 CAR-T annotation.

This script includes reannotation/removal of the cytotoxic NK-like gamma-delta T-cell cluster, G1-specific CART vs CAR-negative comparison, DEG analysis, ORA enrichment analysis, signature scoring, and heatmap visualization.

The initial object-construction section can be skipped if the processed GEO object has already been downloaded.

### `5_PT1_PT2_CITEseq_G1_similarity.R.R`

Performs external CITE-seq similarity analysis.

This script compares the GSE311890 CAR+ G1 reference object with an external PT1Y9 CITE-seq reference object using logistic regression-based prediction and probability/logit heatmaps.

This script requires the following external file:

~~~r
obj2 = readRDS("PT1Y9_filtered_with_AUCell.RDS")
~~~

`PT1Y9_filtered_with_AUCell.RDS` is not included in this repository because of its large file size and data-source restrictions. This object was obtained from an external Zenodo CITE-seq dataset. Please obtain it from the original source or contact the repository owner if needed.

## pySCENIC regulon analysis

The folder `pyscenic_out_G1_C4/` contains key files used for the C4 regulon network analysis.

~~~text
pyscenic_out_G1_C4/
├── C4_regulon_DEG.tsv
└── regulons.csv
~~~

`C4_regulon_DEG.tsv` contains differential regulon activity results for the C4 comparison.

`regulons.csv` is the pySCENIC regulon output used to extract target genes for the top regulons.

The C4 comparison refers to the G1 CAR+ vs CAR-negative regulon activity comparison.

Large pySCENIC input matrices, database files, and intermediate files are not included in this repository.

## Files not included

The following file types are not included because of file size limitations:

~~~text
*.rds
*.RDS
*.h5
*.loom
*.gz
large matrix files
raw Cell Ranger output folders
pySCENIC database files
~~~

Users should download the required input files from GEO, Zenodo, or the original data source before running the scripts.

## Suggested run order

A typical analysis order is:

~~~text
1_9y.R
2_peak.R
3_9y and peak car.R
4_9yclean.R
5_PT1_PT2_CITEseq_G1_similarity.R
~~~

The pySCENIC C4 regulon network analysis requires:

~~~text
C4_regulon_DEG.tsv
regulons.csv
~~~

These files should be placed in:

~~~text
pyscenic_out_G1_C4/
~~~

## Main R dependencies

The main R packages used in this repository include:

~~~r
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
~~~

## Notes

File paths in the scripts may need to be adjusted according to the local working directory.

Large input objects are intentionally excluded from this repository. This repository is intended to provide reproducible analysis code and selected downstream regulon outputs, rather than to store all raw and processed data files.
