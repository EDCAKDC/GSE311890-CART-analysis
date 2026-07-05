library(Seurat)
library(scRepertoire)
library(harmony)
library(stringr)
library(openxlsx)
library(dplyr)
library(scales)
library(ggplot2)
samples = c("9", "p")                                    # sample prefixes
matrixzip = paste0(samples, "raw_feature_bc_matrix.zip") # count matrix zip
tcrcsv = paste0(samples, "all_contig_annotations.csv")   # TCR annotation

## Unzip raw_feature_bc_matrix
for (i in seq_along(samples)) {
  unzip(matrixzip[i], exdir = paste0(samples[i], "_matrix"), overwrite = FALSE)
}

## Read 10x counts and build Seurat objects
seurat.list = lapply(samples, function(s) {
  path = file.path(paste0(s, "_matrix"), "raw_feature_bc_matrix")
  counts = Read10X(data.dir = path, gene.column = 2)
  obj = CreateSeuratObject(counts, project = s)
  obj$sample = s
  obj
})
names(seurat.list) = samples

## Merge samples
seurat = merge(
  seurat.list[[1]],
  y            = seurat.list[-1],
  add.cell.ids = samples,
  project      = "Combined"
)

## QC filtering
seurat[["percent.mt"]] = PercentageFeatureSet(seurat, pattern = "^MT-")  

seurat = subset(
  seurat,
  subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10
)

## Read and integrate TCR
tcr.list = lapply(tcrcsv, function(f) {
  df = read.csv(f)
  subset(df, is_cell == "true")
})
names(tcr.list) = samples

tcr.comb = combineTCR(tcr.list, samples = samples, ID = samples)
barcode_suffix = str_extract(colnames(seurat)[1], "-[0-9]+$")

for (i in seq_along(tcr.comb)) {
  df = tcr.comb[[i]]
  
  df$barcode = gsub("__", "", df$barcode)  # remove double underscores
  df$barcode = gsub(paste0("^(", samples[i], "_)+"),
                    paste0(samples[i], "_"),
                    df$barcode)
  df$barcode = gsub("(-1)+$", "-1", df$barcode)  # fix duplicated "-1"
  tcr.comb[[i]] = df
}

## Add TCR to seurat
seurat = combineExpression(tcr.comb, seurat, cloneCall = "aa")

## Normalization / HVGs / PCA
seurat = NormalizeData(seurat) |>
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA(npcs = 50)

## Harmony and clustering
seurat = RunHarmony(seurat, group.by.vars = "sample")
seurat = RunUMAP(seurat, reduction = "harmony", dims = 1:30)
seurat = FindNeighbors(seurat, reduction = "harmony", dims = 1:30)
seurat = FindClusters(seurat, resolution = 0.06)

## Save object and plot 
Idents(seurat) = "seurat_clusters"
DimPlot(seurat, label = TRUE)
FeaturePlot(seurat, features = c("CD3D", "CD14", "CD8A", "CD4", "FCGR3A", "NCAM1", "CAR", "GZMB", "MKI67"))
saveRDS(seurat, "seurat_9p_TCR_harmony.rds")

## Subset y9 sample from the seurat
data9 = subset(seurat, subset = orig.ident == "9")

## Remove doublets
library(SingleCellExperiment)
library(scDblFinder)
set.seed(123)
sce9 = as.SingleCellExperiment(data9)
sce9 = scDblFinder(sce9)
data9$scDblFinder.class = sce9$scDblFinder.class
data9$scDblFinder.score = sce9$scDblFinder.score
data9_clean = subset(data9, subset = scDblFinder.class == "singlet")

## Re-normalize, dimensional reduction, clustering
data9_clean = NormalizeData(data9_clean)
data9_clean = FindVariableFeatures(data9_clean)
data9_clean = ScaleData(data9_clean)
data9_clean = RunPCA(data9_clean)
data9_clean = RunUMAP(data9_clean, dims = 1:30)
data9_clean = FindNeighbors(data9_clean, dims = 1:30)
data9_clean = FindClusters(data9_clean, resolution = 0.3)

## Save object and plot
DimPlot(data9_clean, label = TRUE)
FeaturePlot(data9_clean, features = c("CD3D", "CD14", "CD8A", "CD4", "FCGR3A", "NCAM1", "CAR", "GZMB", "MKI67"))
saveRDS(data9_clean, file = "data9_clean_new.rds")

## Density plots for TCR genes
library(Nebulosa)
library(patchwork)
genes_tcr = c(
  # T cell markers
  "CD3E", "CD4", "CD8A", "CD8B",
  
  # αβ TCR
  "TRAC", "TRBC1", "TRBC2",
  
  # γδ TCR
  "TRDC", "TRGC1", "TRGC2",'CD14','CAR'
)
genes_tcr = genes_tcr[genes_tcr %in% rownames(data9_clean)]
plot_density(
  data9_clean,
  features = genes_tcr,
  combine = TRUE
)

## DEG analysis
all_deg = FindAllMarkers(
  data9_clean,
  assay = "RNA",
  slot  = "data",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
deg_to_save = all_deg[, c("gene", "cluster", "avg_log2FC", "p_val_adj")]
write.xlsx(deg_to_save, file = "y9_Log2FC_Padj.xlsx", rowNames = FALSE)

## Annotation
data9_clean$celltype_simple = "CAR Neg Lymphocytes"
data9_clean$celltype_simple[data9_clean$seurat_clusters %in% c(1, 6, 8, 14)] = "monocytes"
data9_clean$celltype_simple[data9_clean$seurat_clusters %in% c(10, 7)] = "CART"

## Check counts of celltype and plot
table(data9_clean$celltype_simple)
DimPlot(
  data9_clean,
  group.by = "celltype_simple",
  label = FALSE,
  repel = TRUE
) +
  ggtitle(NULL) +
  theme(plot.title = element_blank())

## Remove monocytes
data9_T_new = subset(data9_clean, subset = !(seurat_clusters %in% c(1, 6, 8, 14)))

## Re-run standard workflow
data9_T_new = NormalizeData(data9_T_new)
data9_T_new = FindVariableFeatures(data9_T_new)
data9_T_new = ScaleData(data9_T_new)
data9_T_new = RunPCA(data9_T_new, npcs = 30)
data9_T_new = RunUMAP(data9_T_new, dims = 1:30)
data9_T_new = FindNeighbors(data9_T_new, dims = 1:30)
data9_T_new = FindClusters(data9_T_new, resolution = 0.06)
data9_T_new = FindClusters(data9_T_new, resolution = 0.08)
data9_T_new = FindClusters(data9_T_new, resolution = 0.1)
data9_T_new = FindClusters(data9_T_new, resolution = 1.2)
data9_T_new = FindClusters(data9_T_new, resolution = 0.3)
## Plot clustering results
DimPlot(data9_T_new, label = TRUE)
FeaturePlot(data9_T_new, features = c("CD3D", "CD14", "CD8A", "CD4", "FCGR3A", "NCAM1", "CAR", "GZMB", "MKI67"))
FeaturePlot(data9_T_new, features = c("TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2", "TRGV9", "TRDV1", "CAR"))
VlnPlot(
  data9_T_new,
  features = c("CD3D", "CD14", "CD8A", "CD4", "FCGR3A", "NCAM1", "CAR", "GZMB", "MKI67"),
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 3
)
p_mt <- VlnPlot(
  data9_T_new,
  features = c("MT-ND1", "MT-CO1", "MT-CO2"),
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 3
)

p_mt
VlnPlot(
  data9_T_new,
  features = c("CD3E", "CD4", "CD8A", "CD8B",
               "CCR7", "SELL", "IL7R", "CD27",
               "CTLA4", "LAG3", "PDCD1", "TIGIT",
               "PCNA", "MKI67", "CD69", "CD38",
               "GZMA", "GZMB", "GZMH", "GZMK",
               "LAMP1", "IFNG", "PRF1", "NKG7"),
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 3
)
VlnPlot(
  data9_T_new,
  features = c("IRF1", "IRF7", "STAT1", "ITGAE", "BHLHE40",
               "ITGA1", "CXCR6", "ZNF683", "PRDM1",
               "S1PR1", "KLF2", "RUNX3", "TBX21", "EOMES",
               "IFIT1", "IFIT3", "ISG15", "OAS1", "OAS3",
               "MX1", "IL15",
               "CCL3", "CCL4", "CCL5", "CXCL13", 
               "CCR5", "CXCR3", "CXCR4", "CXCR5", "IL15RA"),
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 3
)
genes_fig1 = c(
  "CD3E", "CD4", "CD8A", "CD8B",
  "CCR7", "SELL", "IL7R", "CD27",
  "CTLA4", "LAG3", "PDCD1", "TIGIT",
  "PCNA", "MKI67", "CD69", "CD38",
  "GZMA", "GZMB", "GZMH", "GZMK",
  "LAMP1", "IFNG", "PRF1", "NKG7"
)
genes_fig2 = c(
  "IRF1", "IRF7", "STAT1", "ITGAE", "BHLHE40",
  "ITGA1", "CXCR6", "ZNF683", "PRDM1",
  "S1PR1", "KLF2", "RUNX3", "TBX21", "EOMES",
  "IFIT1", "IFIT3", "ISG15", "OAS1", "OAS3",
  "MX1", "IL15",
  "CCL3", "CCL4", "CCL5", "CXCL13", 
  "CCR5", "CXCR3", "CXCR4", "CXCR5", "IL15RA"
)










## Remove Unknown Cells (13 and 14)
data9_T_new = subset(data9_T_new, subset = !(seurat_clusters %in% c(13, 14)))

## Density plots for T cell markers
genes_tcr1 = c("CD3E", "CD4", "CD8A", "CD8B", "CAR")
genes_fig1 = genes_tcr1[genes_tcr1 %in% rownames(data9_T_new)]
plot_density(
  data9_T_new,
  features = genes_fig1,
  combine = TRUE   
)

## Density plots for TCR genes
genes_tcr2 = c("TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2", "TRGV9", "TRDV1", "CAR")
genes_fig2 = genes_tcr2[genes_tcr2 %in% rownames(data9_T_new)]
plot_density(
  data9_T_new,
  features = genes_fig2,
  combine = TRUE  
)

## Feature plots 
FeaturePlot(data9_T_new, features = c("CD3D", "CD14", "CD8A", "CD4", "FCGR3A", "NCAM1", "CAR", "GZMB", "MKI67"))
FeaturePlot(data9_T_new, features = c("TRAC", "TRBC1", "TRBC2", "TRDC",  "TRGC1", "TRGC2", "TRGV9", "TRDV1", "CAR"))

## DEG analysis
all_deg = FindAllMarkers(
  data9_T_new,
  assay = "RNA",         
  slot  = "data",       
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
deg_to_save = all_deg[, c("gene", "cluster", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
deg_to_save$pctw = deg_to_save$pct.1 - deg_to_save$pct.2
write.xlsx(deg_to_save, file = "T_9y_deg_with_pct.xlsx", rowNames = FALSE)

## Annotation
data9_T_new$celltype_simple[data9_T_new$seurat_clusters %in% c(1,2,3,6, 8, 0)] = "CAR Neg Lymphocytes"
data9_T_new$celltype_simple[data9_T_new$seurat_clusters %in% c(4, 7)] = "CART"
data9_T_new$celltype_simple[data9_T_new$seurat_clusters %in% c(5)] = "NK Cells"
## Check counts of celltypes and plot
table(data9_T_new$celltype_simple)
DimPlot(
  data9_T_new,
  group.by = "celltype_simple", 
  label = FALSE,                 
  repel = TRUE
) +
  ggtitle(NULL) +                
  theme(plot.title = element_blank())  

## Re-Clustering
data9_T_new = FindClusters(data9_T_new, resolution = 0.3)
DimPlot(data9_T_new, label = TRUE)

## Re-Annotation
data9_T_new$celltype_simple = "CAR Neg Lymphocytes"
data9_T_new$celltype_simple[data9_T_new$seurat_clusters %in% c(6)] = "NK Cells"
data9_T_new$celltype_simple[data9_T_new$seurat_clusters %in% c(5, 9)] = "CART"
DimPlot(
  data9_T_new,
  group.by = "celltype_simple",  
  label = FALSE,                 
  repel = TRUE
) +
  ggtitle(NULL) +                
  theme(plot.title = element_blank())  

## Re-Density plots for T cell markers
genes_fig1 = genes_tcr[genes_tcr %in% rownames(data9_T_new)]
plot_density(
  data9_T_new,
  features = genes_fig1,
  combine = TRUE  
)

## RE-DEG analysis
all_deg = FindAllMarkers(
  data9_T_new,
  assay = "RNA",         
  slot  = "data",    
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
deg_to_save = all_deg[, c("gene", "cluster", "avg_log2FC", "p_val_adj")]
write.xlsx(deg_to_save, file = "y9_T_Log2FC_Padj.xlsx", rowNames = FALSE)

## Remove NK Cells and save results
data9_T_new = subset(data9_T_new, 
                      subset = !(seurat_clusters %in% c( 6)))
saveRDS(data9_T_new, file = "data9_T_new.rds")

## Subset CAR clusters (5 and 9)
data9_car_new = subset(data9_T_new, subset = seurat_clusters %in% c(5, 9))

## Standard workflow
data9_car_new = NormalizeData(data9_car_new)
data9_car_new = FindVariableFeatures(data9_car_new)
data9_car_new = ScaleData(data9_car_new)
data9_car_new = RunPCA(data9_car_new)
data9_car_new = RunUMAP(data9_car_new, dims = 1:30)
data9_car_new = FindNeighbors(data9_car_new, dims = 1:30)
data9_car_new = FindClusters(data9_car_new, resolution = 0.3)  

## Plot clustering results
DimPlot(data9_car_new, label = TRUE)
FeaturePlot(data9_car_new, features = c("CD3D", "CD14", "CD8A", "CD4", "FCGR3A", "NCAM1", "CAR", "GZMB", "MKI67"))

## DEG analysis
all_deg = FindAllMarkers(
  data9_car_new,
  assay = "RNA",        
  slot  = "data",        
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
deg_to_save = all_deg[, c("gene", "cluster", "avg_log2FC", "p_val_adj")]
write.xlsx(deg_to_save, file = "y9_car_Log2FC_Padj.xlsx", rowNames = FALSE)

## Subset CAR clusters (0–4)
data9_car_new = subset(data9_car_new, subset = seurat_clusters %in% c(0, 1, 2, 3, 4))

## Standard workflow
data9_car_new = NormalizeData(data9_car_new)
data9_car_new = FindVariableFeatures(data9_car_new)
data9_car_new = ScaleData(data9_car_new)
data9_car_new = RunPCA(data9_car_new)
data9_car_new = RunUMAP(data9_car_new, dims = 1:30)
data9_car_new = FindNeighbors(data9_car_new, dims = 1:30)
data9_car_new = FindClusters(data9_car_new, resolution = 0.1)

## Plot clustering results
DimPlot(data9_car_new, label = TRUE)
FeaturePlot(data9_car_new, features = c("CD3D", "CD14", "CD8A", "CD4", "FCGR3A", "NCAM1", "CAR", "GZMB", "MKI67"))
table(data9_car_new$celltype_simple)
## Cell cycle scoring
data9_car_new = CellCycleScoring(
  data9_car_new,
  s.features   = Seurat::cc.genes.updated.2019$s.genes,
  g2m.features = Seurat::cc.genes.updated.2019$g2m.genes,
  set.ident    = FALSE,
  assay        = "RNA"
)

## Plot by cell cycle phase
DimPlot(data9_car_new, reduction = "umap", group.by = "Phase")

## Cluster proportions 
prop.table(table(data9_car_new$seurat_clusters))

## Bar plot 
df = data9_car_new@meta.data %>%
  count(seurat_clusters, name = "n") %>%
  mutate(
    cluster = factor(seurat_clusters),
    pct = n / sum(n),
    pct_lab = percent(pct, accuracy = 0.1)
  )
ggplot(df, aes(x = cluster, y = pct)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = n), color = "white", fontface = "bold", size = 4, vjust = 1.2) +
  geom_text(aes(y = pct + 0.03, label = pct_lab), fontface = "bold", size = 3.6) +
  scale_y_continuous(labels = percent, limits = c(0, max(df$pct) * 1.15)) +
  labs(x = "Cluster", y = "Proportion") +
  theme_classic()

## Density plot
genes_tcr = c("TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2", "TRGV9", "TRDV1",'CAR')
genes_fig1 = c(
  "CD3E", "CD4", "CD8A", "CD8B",
  "CCR7", "SELL", "IL7R", "CD27",
  "CTLA4", "LAG3", "PDCD1", "TIGIT",
  "PCNA", "MKI67", "CD69", "CD38",
  "GZMA", "GZMB", "GZMH", "GZMK",
  "LAMP1", "IFNG", "PRF1", "NKG7"
)
genes_fig2 = c(
  "IRF1", "IRF7", "STAT1", "ITGAE", "BHLHE40",
  "ITGA1", "CXCR6", "ZNF683", "PRDM1",
  "S1PR1", "KLF2", "RUNX3", "TBX21", "EOMES",
  "IFIT1", "IFIT3", "ISG15", "OAS1", "OAS3",
  "MX1", "IL15",
  "CCL3", "CCL4", "CCL5", "CXCL13", 
  "CCR5", "CXCR3", "CXCR4", "CXCR5", "IL15RA"
)
genes_tcr = genes_tcr[genes_tcr %in% rownames(data9_car_new)]
genes_fig1 = genes_fig1[genes_fig1 %in% rownames(data9_car_new)]
genes_fig2 = genes_fig2[genes_fig2 %in% rownames(data9_car_new)]
plot_density(
  data9_car_new,
  features = genes_tcr,
  combine = TRUE 
)
plot_density(
  data9_car_new,
  features = genes_fig1,
  combine = TRUE  
)
plot_density(
  data9_car_new,
  features = genes_fig2,
  combine = TRUE  
)

## Plot UMAP for top 10 clonotypes
top10 = data9_car_new@meta.data %>%
  filter(!is.na(CTaa), CTaa != "") %>%
  count(CTaa) %>%
  arrange(desc(n)) %>%
  slice_head(n = 10) %>%
  pull(CTaa)
data9_car_new$CTaa_top10 = ifelse(
  data9_car_new$CTaa %in% top10,
  paste0("clonotype", match(data9_car_new$CTaa, top10)),
  "Others"
)
data9_car_new$CTaa_top10 = factor(
  data9_car_new$CTaa_top10,
  levels = c(paste0("clonotype", 1:10), "Others") # Set clonotypes order
)
cols = c(setNames(hue_pal()(10), paste0("clonotype", 1:10)), Others = "grey")
DimPlot(data9_car_new, reduction = "umap", group.by = "CTaa_top10", cols = cols) +
  theme(plot.title = element_blank())

## Plot UMAP for dominant clonotype
dominant = data9_car_new@meta.data %>%
  filter(!is.na(CTaa), CTaa != "") %>%
  count(CTaa, sort = TRUE) %>%
  slice_head(n = 1) %>%
  pull(CTaa)
data9_car_new$CTaa_dom = ifelse(data9_car_new$CTaa == dominant,
                                "dominant clonotype", "Others")
data9_car_new$CTaa_dom = factor(data9_car_new$CTaa_dom,
                                levels = c("dominant clonotype", "Others"))

## Replace NA with Others
data9_car_new$CTaa_dom[is.na(data9_car_new$CTaa_dom)] = "Others"

## Convert to factor again
data9_car_new$CTaa_dom = factor(
  data9_car_new$CTaa_dom,
  levels = c("dominant clonotype", "Others")
)

## Set colors and plot
cols = c("dominant clonotype" = "#8491B4FF", "Others" = "grey")
DimPlot(data9_car_new, reduction = "umap",
        group.by = "CTaa_dom", cols = cols) +
  theme(plot.title = element_blank())

## DEG analysis
all_deg = FindAllMarkers(
  data9_car_new,
  assay = "RNA",         
  slot  = "data",       
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
deg_to_save = all_deg[, c("gene", "cluster", "avg_log2FC", "p_val_adj")]
write.xlsx(deg_to_save, file = "y9_car_Log2FC_Padj.xlsx", rowNames = FALSE)


## DEG analysis: Dominant clonotype vs Others
Idents(data9_car_new) = data9_car_new$CTaa_dom
deg_dom = FindMarkers(
  data9_car_new,
  ident.1 = "dominant clonotype",
  ident.2 = "Others",
  assay   = "RNA",
  slot    = "data",       
  only.pos = FALSE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
deg_to_save = deg_dom[, c("avg_log2FC", "p_val_adj")]
deg_to_save = tibble::rownames_to_column(as.data.frame(deg_to_save), "gene")
write.xlsx(deg_to_save, file = "y9_dominant_vs_others_DEG.xlsx", rowNames = FALSE)

library(openxlsx)
library(dplyr)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
## Read DEG table 
df = read.xlsx("y9_dominant_vs_others_DEG.xlsx")

## Build geneList (named vector: names = SYMBOL, values = log2FC)
geneList = df |>
  distinct(gene, .keep_all = TRUE) |>
  arrange(desc(avg_log2FC))
geneList = setNames(geneList$avg_log2FC, geneList$gene)
geneList = sort(geneList, decreasing = TRUE)

## Hallmark gene sets 
hallmark_sym = msigdbr(species = "Homo sapiens", category = "H") |>
  dplyr::select(gs_name, gene_symbol) |>
  unique()

## GSEA analysis
egsea = GSEA(
  geneList,
  TERM2GENE = hallmark_sym,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

## Save GSEA results
write.xlsx(as.data.frame(egsea), "9y_GSEA_dominant_vs_others_Hallmark.xlsx",rowNames = FALSE)

## Plots
dotplot(egsea, showCategory = 20)
clusterProfiler::barplot(
  egsea,
  showCategory = 20,
  title = "GSEA - Dominant vs Others",
  font.size = 12
)

## Save results
saveRDS(data9_car_new, file = "data9_car_new.rds")


