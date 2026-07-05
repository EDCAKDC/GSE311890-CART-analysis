# GSE311890 CART Analysis

This repository contains the analysis scripts used for the GSE311890 CAR-T single-cell analysis. The analysis includes year 9.3 CAR-T/T-cell processing, peak vs year 9.3 CAR+ comparison, G1-specific analysis, pySCENIC regulon analysis, external CITE-seq similarity analysis, and downstream visualization.

## Data availability

Large input files are not included in this repository because of file size limitations.

Raw and processed data files are available from GEO under accession **GSE311890**.

Some scripts contain an initial section for constructing Seurat objects from raw or intermediate files. If the processed GEO files have already been downloaded, this object-construction section can be skipped, and the downstream analysis can be started directly from the processed `.rds` objects.

## Repository structure

```text
GSE311890-CART-analysis/
├── README.md
├── scripts/
│   ├── 1_9y.R
│   ├── 2_peak.R
│   ├── 3_9y_peak_car.R
│   ├── 4_9yclean.R
│   └── 5_PT1_PT2_CITEseq_G1_similarity.R
│
└── pyscenic_out_G1_C4/
    ├── C4_regulon_DEG.tsv
    └── regulons.csv
