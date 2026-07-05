library(Seurat)
library(scRepertoire)
library(harmony)
library(stringr)
library(openxlsx)
library(dplyr)
library(scales)
library(ggplot2)
## Read rds files
datap_T_new = readRDS("datap_T_new.rds")
data9_car_new = readRDS("data9_car_new.rds")

## Collect CAR+ in peak
datap_car_new = subset(datap_T_new, subset = celltype_simple == "CART")

## Tag sample source
datap_car_new$sample = "peak"
data9_car_new$sample = "y9"

## Merge two CAR+ objects
obj = merge(datap_car_new, y = data9_car_new, add.cell.ids = c("peak","y9"), project = "CARpos")

## Normalization / HVGs / PCA
obj$sample = sub("_.*$", "", colnames(obj))
DefaultAssay(obj) = "RNA"
obj = obj |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA(npcs = 50)

## Harmony and clustering
obj = RunHarmony(obj, group.by.vars = "sample")
obj = RunUMAP(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
obj = FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj = FindClusters(obj, resolution = 0.33)## 0.47 5 clusters

## Standard plots
DimPlot(obj, label = TRUE)
FeaturePlot(obj, features = c("CD3D","CD14","CD8A","CD4","FCGR3A","NCAM1","CTL019","GZMB","MKI67"))

## Counts by cluster × sample
table(obj$seurat_clusters, obj$sample)

## Merge RNA layers and Normalization
DefaultAssay(obj) = "RNA"
obj = JoinLayers(obj, assay = "RNA")
obj = NormalizeData(obj, assay = "RNA")

## Cell cycle scoring
DefaultLayer(obj[["RNA"]]) = "data"
obj = CellCycleScoring(
  obj,
  s.features   = Seurat::cc.genes.updated.2019$s.genes,
  g2m.features = Seurat::cc.genes.updated.2019$g2m.genes,
  set.ident    = FALSE,
  assay        = "RNA"
)

## Plot stacked bar chart
df = obj@meta.data %>%
  group_by(seurat_clusters, Phase) %>%
  summarise(Count = n()) %>%
  group_by(seurat_clusters) %>%
  mutate(Proportion = Count / sum(Count))
ggplot(df, aes(x = seurat_clusters, y = Proportion, fill = Phase)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count),
            position = position_stack(vjust = 0.5),
            size = 3,
            color = "black") +
  theme_classic() +
  labs(x = "Cluster", y = "Proportion of cells") +
  scale_fill_manual(values = c("G1" = "#66c2a5", "S" = "#8da0cb", "G2M" = "#fc8d62")) +
  ggtitle("Cell cycle phase distribution across clusters")

## DEG analysis
Idents(obj) = "seurat_clusters"
all_deg = FindAllMarkers(
  obj,
  assay = "RNA",         
  slot  = "data",        
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
deg_to_save = all_deg[, c("gene", "cluster", "avg_log2FC", "p_val_adj")]
write.xlsx(deg_to_save, file = "9p_car_Log2FC_Padj.xlsx", rowNames = FALSE)

## DEG analysis in sample
Idents(obj) = "sample"
all_deg = FindAllMarkers(
  obj,
  assay = "RNA",         
  slot  = "data",     
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
deg_to_save = all_deg[, c("gene", "cluster", "avg_log2FC", "p_val_adj")]
write.xlsx(deg_to_save, file = "9p_car_Log2FC_Padj_sample.xlsx", rowNames = FALSE)

Idents(obj_g1) = "sample"
all_deg = FindAllMarkers(
  obj_g1,
  assay = "RNA",         
  slot  = "data",     
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
deg_to_save = all_deg[, c("gene", "cluster", "avg_log2FC", "p_val_adj")]
write.xlsx(deg_to_save, file = "9p_g1_car_Log2FC_Padj_sample_new.xlsx", rowNames = FALSE)
## UMAP in Phase and sample
DimPlot(obj, reduction = "umap", group.by = "Phase")
DimPlot(obj, reduction = "umap", group.by = "sample")

## Clusters × sample
tab = table(obj$seurat_clusters, obj$sample)

## Write counts, column-normalized props, row-normalized props, and per-sample totals
openxlsx::write.xlsx(
  list(
    counts = as.data.frame.matrix(tab),                                         # absolute counts
    prop_by_timepoint = as.data.frame.matrix(round(prop.table(tab, 2), 4)),     # column-norm (sample)
    prop_by_cluster   = as.data.frame.matrix(round(prop.table(tab, 1), 4)),     # row-norm (cluster)
    totals_timepoint  = data.frame(sample = colnames(tab),
                                   total  = colSums(tab),
                                   row.names = NULL)
  ),
  file = "y9_peak_cluster_summary.xlsx",
  overwrite = TRUE,
  rowNames = FALSE
)

## Collect dominant clonotype in y9 
top_clone = names(sort(table(obj$CTaa[obj$sample == "y9"]), decreasing = TRUE))[1]

## Cells of that clonotype in y9 and peak
cells_y9 = colnames(obj)[obj$sample == "y9" & obj$CTaa == top_clone & !is.na(obj$CTaa)]
cells_peak = colnames(obj)[obj$sample == "peak" & obj$CTaa == top_clone & !is.na(obj$CTaa)]

## UMAP: y9=blue, peak=red, others=grey
p = DimPlot(
  obj, reduction = "umap",
  cells.highlight = list(y9_top = cells_y9, peak_top = cells_peak),
  cols.highlight  = c("blue", "red"),
  cols = "grey",                 
  sizes.highlight = 1.2
)
p + scale_color_manual(
  values = c(y9_top = "blue", peak_top = "red", Unselected = "grey"),
  breaks = c("y9_top", "peak_top", "Unselected"),
  labels = c("y9_top", "peak_top", "others")
) +
  guides(colour = guide_legend(title = NULL)) +
  ggtitle(NULL) +
  theme(plot.title = element_blank())

## Density plots
genes_fig1 = c(
  "CD3E", "CD4", "CD8A", "CD8B",
  "CCR7", "SELL", "IL7R", "CD27",
  "CTLA4", "LAG3", "PDCD1", "TIGIT",
  "PCNA", "MKI67", "CD69", "CD38",
  "GZMA", "GZMB", "GZMH", "GZMK",
  "LAMP1", "IFNG", "PRF1","NKG7"
)
genes_fig2 = c(
  "IRF1", "IRF7", "STAT1", "ITGAE", "BHLHE40",
  "ITGA1", "CXCR6", "ZNF683", "PRDM1",
  "S1PR1", "KLF2", "RUNX3", "TBX21", "EOMES",
  "IFIT1", "IFIT3", "ISG15", "OAS1", "OAS3",
  "MX1",    "IL15",
  "CCL3", "CCL4", "CCL5", "CXCL13", 
  "CCR5", "CXCR3", "CXCR4", "CXCR5", "IL15RA"
)
genes_fig3 = c(
  "FXYD2","HMOX1","GPR183","TIGIT","HLA-DRA",
  "TMEM155","FCMR","MTSS1","TTN","DENND2D",
  "NUCB2","EOMES","CAV1","LYAR","ANTXR2",
  "PASK","CEP128","CMTM7","ABHD17B","LIMD2","GNA12"
)
genes_fig1 = genes_fig1[genes_fig1 %in% rownames(obj)]
genes_fig2 = genes_fig2[genes_fig2 %in% rownames(obj)]
genes_fig3 = genes_fig3[genes_fig3 %in% rownames(obj)]
plot_density(
  obj,
  features = genes_fig1,
  combine = TRUE 
)
plot_density(
  obj,
  features = genes_fig2,
  combine = TRUE 
)
plot_density(
  obj,
  features = genes_fig3,
  combine = TRUE  
)

## Save results
saveRDS(obj, "car_9p_new.rds")























