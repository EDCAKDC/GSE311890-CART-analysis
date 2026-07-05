library(Seurat)
library(scRepertoire)
library(harmony)
library(stringr)
library(openxlsx)
library(dplyr)
library(scales)
library(ggplot2)
## Subset peak sample from the seurat
datap = subset(seurat, subset = orig.ident == "p")

## Remove doublets
library(SingleCellExperiment)
library(scDblFinder)
set.seed(123)
scep = as.SingleCellExperiment(datap)
scep = scDblFinder(scep)
datap$scDblFinder.class = scep$scDblFinder.class
data9$scDblFinder.score = scep$scDblFinder.score
datap_clean = subset(datap, subset = scDblFinder.class == "singlet")

## Re-normalize, dimensional reduction, clustering
datap_clean = NormalizeData(datap_clean)
datap_clean = FindVariableFeatures(datap_clean)
datap_clean = ScaleData(datap_clean)
datap_clean = RunPCA(datap_clean)
datap_clean = RunUMAP(datap_clean, dims = 1:30)
datap_clean = FindNeighbors(datap_clean, dims = 1:30)
datap_clean = FindClusters(datap_clean, resolution = 0.3)

## Standard plots
DimPlot(datap_clean, label = TRUE)
FeaturePlot(datap_clean, features = c("CD3D", "CD14", "CD8A", "CD4", "FCGR3A", "NCAM1", "CAR", "GZMB", "MKI67"))
VlnPlot(datap_clean, features = c("CD3D", "CD14", "CD8A", "CD4", "FCGR3A", "NCAM1", "CAR", "GZMB", "MKI67"), pt.size = 0, ncol = 3)

## Annotation
datap_clean$celltype_simple = "CAR Neg Lymphocytes"
datap_clean$celltype_simple[datap_clean$seurat_clusters %in% c(3, 13)] = "monocytes"
datap_clean$celltype_simple[datap_clean$seurat_clusters %in% c(1, 9, 14)] = "CART"

## Check counts of celltype and plot
table(datap_clean$celltype_simple)
DimPlot(
  datap_clean,
  group.by = "celltype_simple", 
  label = FALSE,                 
  repel = TRUE
) +
  ggtitle(NULL) +              
  theme(plot.title = element_blank()) 

## Remove monocytes
datap_T_new = subset(datap_clean, subset = !(seurat_clusters %in% c(3,13)))

## Re-run standard workflow
datap_T_new = NormalizeData(datap_T_new)
datap_T_new = FindVariableFeatures(datap_T_new)
datap_T_new = ScaleData(datap_T_new)
datap_T_new = RunPCA(datap_T_new, npcs = 30)
datap_T_new = RunUMAP(datap_T_new, dims = 1:30)
datap_T_new = FindNeighbors(datap_T_new, dims = 1:30)
datap_T_new = FindClusters(datap_T_new, resolution = 0.3)

## Plot clustering results
DimPlot(datap_T_new, label = TRUE)

## Remove Unknown Cells (10) and plot
datap_T_new = subset(datap_T_new, subset = !(seurat_clusters %in% c(10)))
datap_T_new = FindClusters(datap_T_new, resolution = 0.3)
DimPlot(datap_T_new, label = TRUE)
FeaturePlot(datap_T_new, features = c("CD3D", "CD14", "CD8A", "CD4", "FCGR3A", "NCAM1", "CAR", "GZMB", "MKI67"))

## Re-Annotation
datap_T_new$celltype_simple = "CAR Neg Lymphocytes"
datap_T_new$celltype_simple[datap_T_new$seurat_clusters %in% c(11)] = "NK Cells"
datap_T_new$celltype_simple[datap_T_new$seurat_clusters %in% c(0, 13)] = "CART"

## Re-Check counts of celltype and plot
table(datap_T_new$celltype_simple)
DimPlot(
  datap_T_new,
  group.by = "celltype_simple",  
  label = FALSE,                 
  repel = TRUE
) +
  ggtitle(NULL) +                
  theme(plot.title = element_blank()) 

## Save results
saveRDS(datap_T_new, file = "datap_T_new.rds")

## Remive NK cells
datap_T_noNK = subset(datap_T_new, subset = celltype_simple != "NK Cells")
table(datap_T_noNK$celltype_simple)

## Keep productive cells only
anno_p = read.csv("pall_contig_annotations.csv", stringsAsFactors = FALSE)
anno_p_prod = subset(anno_p, is_cell == "true" & productive == "true")
bc_p_prod = paste0("p_", unique(anno_p_prod$barcode))
datap_prod = subset(datap_T_noNK, cells = bc_p_prod)

## Clonotypes and proportions
count_prop = function(df) {
  df %>%
    filter(!is.na(CTaa), CTaa != "") %>%
    count(CTaa, name = "count") %>%
    arrange(desc(count)) %>%
    mutate(prop = count / sum(count))
}

## Save results
tbl_neg = count_prop(datap_prod@meta.data %>% filter(celltype_simple == "CAR Neg Lymphocytes"))
tbl_pos = count_prop(datap_prod@meta.data %>% filter(celltype_simple == "CART"))
tbl_all = count_prop(datap_prod@meta.data)
write.xlsx(
  list(
    CAR_Neg_Lymphocytes = tbl_neg,
    CART = tbl_pos,
    ALL = tbl_all
  ),
  file = "CTaa_counts_peak_productive.xlsx",
  rowNames = FALSE
)
















