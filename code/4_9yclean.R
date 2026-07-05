library(Seurat)
library(scRepertoire)
library(harmony)
library(stringr)
library(openxlsx)
library(dplyr)
library(scales)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)
library(tibble)

## Read updated rds files
data9_car_new = readRDS("data9_car_new.rds")
data9_T_new = readRDS("data9_T_new.rds")

## Re-annotate cluster 3 in Y9 T object
cells_cluster3 = colnames(data9_car_new)[data9_car_new$seurat_clusters == 3]
cells_cluster3 = intersect(cells_cluster3, colnames(data9_T_new))

data9_T_new$celltype_simple[cells_cluster3] = "Cytotoxic NK like Gamma Delta T cells"

table(data9_T_new$celltype_simple)

DimPlot(data9_T_new, group.by = "celltype_simple")

FeaturePlot(
  data9_car_new,
  features = c("CD3D", "CD14", "CD8A", "CD4", "FCGR3A", "NCAM1", "CAR", "GZMB", "MKI67")
)

## Remove cluster 3 from Y9 CART object
data9_car_new = subset(
  data9_car_new,
  subset = !(seurat_clusters %in% c(3))
)

DimPlot(data9_car_new, group.by = "seurat_clusters")

saveRDS(data9_T_new, file = "data9_T_new.rds")
saveRDS(data9_car_new, file = "data9_car_new.rds")


## Read peak/y9 CAR+ object
obj = readRDS("car_9p_new.rds")

obj = subset(
  obj,
  subset = seurat_clusters %in% c(0, 1, 2, 3)
)

## Re-run normalization and clustering
DefaultAssay(obj) = "RNA"

obj = obj |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA(npcs = 50)

obj = RunHarmony(obj, group.by.vars = "sample")
obj = RunUMAP(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
obj = FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj = FindClusters(obj, resolution = 0.33)

DimPlot(obj, label = TRUE)
DimPlot(obj, group.by = "sample")

## Cell cycle scoring
DefaultLayer(obj[["RNA"]]) = "data"

obj = CellCycleScoring(
  obj,
  s.features = Seurat::cc.genes.updated.2019$s.genes,
  g2m.features = Seurat::cc.genes.updated.2019$g2m.genes,
  set.ident = FALSE,
  assay = "RNA"
)

## Keep G1 cells only
obj_G1 = subset(
  obj,
  subset = Phase == "G1" & sample %in% c("y9", "peak")
)

table(obj_G1$sample)


## Read objects for Y9 G1 comparison
data9_T_new = readRDS("data9_T_new.rds")
data9_car_new = readRDS("data9_car_new.rds")

table(data9_T_new$Phase, useNA = "ifany")
table(data9_car_new$Phase, useNA = "ifany")
table(data9_T_new$celltype_simple, useNA = "ifany")

## Collect CART G1 cells
cart_g1 = subset(
  data9_car_new,
  subset = Phase == "G1"
)

cart_g1$celltype_simple = "CART"

## Collect CAR-negative G1 cells
carneg_g1 = subset(
  data9_T_new,
  subset = Phase == "G1" & celltype_simple == "CAR Neg Lymphocytes"
)

carneg_g1$celltype_simple = "CAR Neg Lymphocytes"

ncol(cart_g1)
ncol(carneg_g1)

table(cart_g1$Phase, useNA = "ifany")
table(carneg_g1$Phase, useNA = "ifany")

## Merge two G1 groups
obj_g1 = merge(
  cart_g1,
  y = carneg_g1,
  add.cell.ids = c("CART", "CARneg"),
  project = "Y9_CART_vs_CARneg_G1"
)

table(obj_g1$celltype_simple, useNA = "ifany")
table(obj_g1$Phase, useNA = "ifany")

## Join RNA layers and normalize
DefaultAssay(obj_g1) = "RNA"

obj_g1 = NormalizeData(obj_g1)
obj_g1 = JoinLayers(obj_g1)
obj_g1 = NormalizeData(obj_g1)

Idents(obj_g1) = "celltype_simple"

## CART G1 vs CAR-negative G1 DEG
deg = FindMarkers(
  obj_g1,
  ident.1 = "CART",
  ident.2 = "CAR Neg Lymphocytes",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.05,
  return.thresh = 1
)

deg = deg |>
  rownames_to_column("gene")

head(deg)
summary(deg$avg_log2FC)
summary(deg$p_val_adj)

## Split up/down genes
lfc_cut = 0
padj_cut = 0.05

up_tbl = deg |>
  filter(avg_log2FC > lfc_cut, p_val_adj < padj_cut)

down_tbl = deg |>
  filter(avg_log2FC < -lfc_cut, p_val_adj < padj_cut)

## Save DEG table
out_xlsx = "CART_vs_CARneg_G1_DEG_separate_cellcycle.xlsx"

wb = createWorkbook()

addWorksheet(wb, "Upregulated")
writeData(wb, "Upregulated", up_tbl)

addWorksheet(wb, "Downregulated")
writeData(wb, "Downregulated", down_tbl)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

nrow(up_tbl)
nrow(down_tbl)


## ORA function
enrich_from_xlsx = function(xlsx_file, prefix, split_up_down = TRUE, pcut = 0.05, show_n = 20) {
  
  up = tryCatch(
    read.xlsx(xlsx_file, sheet = "Upregulated"),
    error = function(e) NULL
  )
  
  down = tryCatch(
    read.xlsx(xlsx_file, sheet = "Downregulated"),
    error = function(e) NULL
  )
  
  get_genes = function(df) {
    if (is.null(df) || nrow(df) == 0) {
      return(character(0))
    }
    
    if (!"gene" %in% names(df)) {
      names(df)[1] = "gene"
    }
    
    unique(na.omit(df$gene))
  }
  
  genes_up = get_genes(up)
  genes_down = get_genes(down)
  
  do_one_enrich = function(gene_symbols, tag) {
    
    if (length(gene_symbols) < 5) {
      message(prefix, "_", tag, ": too few genes")
      return(invisible(NULL))
    }
    
    map = bitr(
      gene_symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    entrez = unique(na.omit(map$ENTREZID))
    
    if (length(entrez) < 5) {
      message(prefix, "_", tag, ": too few mapped genes")
      return(invisible(NULL))
    }
    
    ## GO enrichment
    for (ont in c("BP", "CC", "MF")) {
      
      ego = enrichGO(
        gene = entrez,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = pcut,
        readable = TRUE
      )
      
      df = as.data.frame(ego)
      
      if (nrow(df) > 0) {
        
        p = barplot(
          ego,
          showCategory = show_n,
          title = paste0(prefix, " ", tag, " GO ", ont),
          font.size = 10
        )
        
        ggsave(
          paste0(prefix, "_", tag, "_GO_", ont, ".pdf"),
          p,
          width = 25,
          height = 16,
          dpi = 300
        )
        
        write.xlsx(
          df,
          paste0(prefix, "_", tag, "_GO_", ont, ".xlsx"),
          rowNames = FALSE
        )
      }
    }
    
    ## Reactome enrichment
    rea = enrichPathway(
      gene = entrez,
      organism = "human",
      pvalueCutoff = pcut,
      readable = TRUE
    )
    
    dfR = as.data.frame(rea)
    
    if (nrow(dfR) > 0) {
      
      p2 = barplot(
        rea,
        showCategory = show_n,
        title = paste0(prefix, " ", tag, " Reactome"),
        font.size = 12,
        label_format = 40
      )
      
      ggsave(
        paste0(prefix, "_", tag, "_Reactome.pdf"),
        p2,
        width = 25,
        height = 16,
        dpi = 300
      )
      
      write.xlsx(
        dfR,
        paste0(prefix, "_", tag, "_Reactome.xlsx"),
        rowNames = FALSE
      )
    }
  }
  
  if (split_up_down) {
    
    if (length(genes_up) > 0) {
      do_one_enrich(genes_up, "Up")
    }
    
    if (length(genes_down) > 0) {
      do_one_enrich(genes_down, "Down")
    }
    
  } else {
    
    genes_all = unique(c(genes_up, genes_down))
    do_one_enrich(genes_all, "AllSig")
  }
}

## Run ORA
enrich_from_xlsx(
  out_xlsx,
  prefix = "CAR_G1_separate_cellcycle",
  split_up_down = TRUE
)




library(Seurat)
library(ggplot2)

## Define late CAR-T signature genes
late_cart_genes = c(
  "FXYD2", "HMOX1", "GPR183", "TIGIT", "HLA-DRA", "TMEM155", "FCMR", "MTSS1", "TTN",
  "DENND2D", "NUCB2", "EOMES", "CAV1", "LYAR", "ANTXR2", "PASK", "CEP128", "CMTM7",
  "ABHD17B", "LIMD2", "GNA12"
)

## Match genes in object
obj_genes_upper = toupper(rownames(obj))
want_upper = toupper(late_cart_genes)

hit_idx = match(want_upper, obj_genes_upper)
sig_genes = rownames(obj)[na.omit(hit_idx)]
missing = late_cart_genes[is.na(hit_idx)]

sig_genes
missing

## Add module score
obj = AddModuleScore(
  object = obj,
  features = list(sig_genes),
  name = "Late_CART_Score",
  assay = DefaultAssay(obj)
)

score_col = "Late_CART_Score1"

summary(obj@meta.data[[score_col]])

## Z-score normalization
obj$Late_CART_Score_z = as.numeric(scale(obj[[score_col]][, 1]))

summary(obj$Late_CART_Score_z)

## Prepare plotting data
df = obj@meta.data[, c("Late_CART_Score_z", "sample")]

table(df$sample)

## Wilcoxon test
res = wilcox.test(
  Late_CART_Score_z ~ sample,
  data = df
)

res

pval = res$p.value

p_txt = ifelse(
  pval < 2.2e-16,
  "p < 2.2e-16",
  paste0("p = ", formatC(pval, format = "e", digits = 2))
)

p_txt

## Violin plot
p = VlnPlot(
  obj,
  features = "Late_CART_Score_z",
  group.by = "sample"
) +
  theme_minimal() +
  ylab("Late CAR-T signature score") +
  xlab(NULL) +
  ggtitle("Late CAR-T signature score by sample")

p = p +
  annotate(
    "text",
    x = 1.5,
    y = Inf,
    vjust = 1.4,
    label = p_txt,
    size = 5
  )

p

library(Seurat)
library(ggplot2)
library(dplyr)

## Read Y9 T object
data9_T_new = readRDS("data9_T_new.rds")

## Keep CART and CAR-negative cells
obj = subset(
  data9_T_new,
  subset = celltype_simple %in% c("CART", "CAR Neg Lymphocytes")
)

DefaultAssay(obj) = "RNA"

obj = obj |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA(npcs = 50)

obj = RunUMAP(obj, dims = 1:30)
obj = FindNeighbors(obj, dims = 1:30)
obj = FindClusters(obj, resolution = 0.3)

table(obj$celltype_simple)

DimPlot(obj, group.by = "celltype_simple")
DimPlot(obj, label = TRUE)

## Define late CAR-T signature genes
late_cart_genes = c(
  "FXYD2", "HMOX1", "GPR183", "TIGIT", "HLA-DRA", "TMEM155", "FCMR", "MTSS1", "TTN",
  "DENND2D", "NUCB2", "EOMES", "CAV1", "LYAR", "ANTXR2", "PASK", "CEP128", "CMTM7",
  "ABHD17B", "LIMD2", "GNA12"
)

## Match genes in object
obj_genes_upper = toupper(rownames(obj))
want_upper = toupper(late_cart_genes)

hit_idx = match(want_upper, obj_genes_upper)

sig_genes = rownames(obj)[na.omit(hit_idx)]
missing = late_cart_genes[is.na(hit_idx)]

sig_genes
missing

## Add module score
obj = AddModuleScore(
  object = obj,
  features = list(sig_genes),
  name = "Late_CART_Score",
  assay = DefaultAssay(obj)
)

score_col = "Late_CART_Score1"

## Z-score normalization
obj$Late_CART_Score_z = as.numeric(scale(obj[[score_col]][, 1]))

summary(obj[[score_col]][, 1])
summary(obj$Late_CART_Score_z)
range(obj$Late_CART_Score_z, na.rm = TRUE)

## Set group order
obj$celltype_simple = factor(
  obj$celltype_simple,
  levels = c("CAR Neg Lymphocytes", "CART")
)

## Prepare plotting data
df = obj@meta.data[, c("Late_CART_Score_z", "celltype_simple")]

table(df$celltype_simple)

## Wilcoxon test
res = wilcox.test(
  Late_CART_Score_z ~ celltype_simple,
  data = df
)

res

pval = res$p.value

p_txt = ifelse(
  pval < 2.2e-16,
  "p < 2.2e-16",
  paste0("p = ", formatC(pval, format = "e", digits = 2))
)

p_txt

## Violin plot
p = VlnPlot(
  obj,
  features = "Late_CART_Score_z",
  group.by = "celltype_simple",
  pt.size = 0
) +
  theme_minimal() +
  ylab("Late CAR-T signature score") +
  xlab(NULL) +
  ggtitle("Late CAR-T signature score in Y9 CAR+ vs CAR-")

p = p +
  annotate(
    "text",
    x = 1.5,
    y = max(df$Late_CART_Score_z, na.rm = TRUE) * 1.05,
    label = p_txt,
    size = 5
  )

p


library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

## Read Y9 T object
data9_T_new = readRDS("data9_T_new.rds")

DefaultAssay(data9_T_new) = "RNA"

data9_T_new = JoinLayers(data9_T_new, assay = "RNA")

if (!"data" %in% Layers(data9_T_new[["RNA"]])) {
  data9_T_new = NormalizeData(data9_T_new, assay = "RNA")
}

DefaultLayer(data9_T_new[["RNA"]]) = "data"

## Define heatmap groups
data9_T_new$heatmap_group = case_when(
  data9_T_new$celltype_simple == "NK Cells" ~ "NK",
  data9_T_new$celltype_simple == "CAR Neg Lymphocytes" ~ "Car negative T-cells",
  data9_T_new$celltype_simple == "CART" ~ "CAR+ T",
  data9_T_new$celltype_simple == "Cytotoxic NK like Gamma Delta T cells" ~ "Cytotoxic NK-like γδ CAR-T cells",
  TRUE ~ NA_character_
)

group_order = c(
  "NK",
  "Car negative T-cells",
  "CAR+ T",
  "Cytotoxic NK-like γδ CAR-T cells"
)

data9_T_new = subset(
  data9_T_new,
  subset = heatmap_group %in% group_order
)

data9_T_new$heatmap_group = factor(
  data9_T_new$heatmap_group,
  levels = group_order
)

table(data9_T_new$heatmap_group)

## Define marker genes
genes_use = c(
  "FCGR3A",
  "NCAM1",
  "KLRD1",
  "CD5",
  "CD3E",
  "TRBC1",
  "TRBC2",
  "TRAC",
  "TRGV9",
  "TRDV1",
  "CAR"
)

genes_present = genes_use[genes_use %in% rownames(data9_T_new)]
missing_genes = setdiff(genes_use, genes_present)

genes_present
missing_genes

stopifnot(length(genes_present) > 0)

## Average expression by group
avg = AverageExpression(
  data9_T_new,
  features = genes_present,
  group.by = "heatmap_group",
  assays = "RNA",
  slot = "data"
)$RNA

avg = avg[, group_order, drop = FALSE]
avg = avg[genes_present, , drop = FALSE]

## Row z-score
scale_rows = function(m) {
  m = as.matrix(m)
  m = t(scale(t(m)))
  m[is.na(m)] = 0
  m
}

mat = scale_rows(avg)

mat[mat > 2] = 2
mat[mat < -2] = -2

## Save matrix
write.csv(
  avg,
  "Y9_lymphocyte_marker_average_expression.csv",
  row.names = TRUE
)

write.csv(
  mat,
  "Y9_lymphocyte_marker_average_expression_zscore.csv",
  row.names = TRUE
)

## Plot heatmap
col_fun = colorRamp2(
  c(-2, 0, 2),
  c("#4B4BA8", "white", "#F06449")
)

ht = Heatmap(
  mat,
  name = "Expr (z)",
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  show_row_names = TRUE,
  show_column_names = TRUE,
  
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  
  column_title = "Heatmap of gene expression defining clusters of lymphocytes (Y9.3)",
  column_title_gp = gpar(fontsize = 13, fontface = "bold"),
  
  rect_gp = gpar(col = "white", lwd = 0.5),
  
  heatmap_legend_param = list(
    title = "Expr (z)",
    at = c(-2, -1, 0, 1, 2),
    legend_height = unit(3, "cm")
  )
)

draw(ht)

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr)

## Read CAR+ object
obj = readRDS("car_9p_new.rds")

DefaultAssay(obj) = "RNA"

obj = JoinLayers(obj, assay = "RNA")

if (!"data" %in% Layers(obj[["RNA"]])) {
  obj = NormalizeData(obj, assay = "RNA")
}

DefaultLayer(obj[["RNA"]]) = "data"

## Keep target clusters
obj$cluster_label = paste0("Cluster ", obj$seurat_clusters)

target_clusters = paste0("Cluster ", 0:3)

obj = subset(
  obj,
  subset = cluster_label %in% target_clusters
)

obj$cluster_label = factor(
  obj$cluster_label,
  levels = target_clusters
)

Idents(obj) = "cluster_label"

table(obj$cluster_label)
table(obj$cluster_label, obj$sample)

## Rename cluster labels
cluster_name_map = c(
  "Cluster 0" = "Cytotoxic CD8+ CAR T",
  "Cluster 1" = "Memory-like CAR T",
  "Cluster 2" = "IFN-Activated-Exhausted CAR T",
  "Cluster 3" = "Proliferating CAR T"
)

## Define marker genes
genes_use = c(
  "GNLY", "CD8A", "KLRD1", "GZMH", "GZMB", "CCL5", "NKG7",
  "ZFP36L2", "TCF7", "BTG1",
  "IFI6", "IL10", "IFIH1", "TIGIT", "CD38",
  "HAVCR2", "LAG3", "CD69", "IFI35",
  "STMN1", "CLSPN", "CCNA2", "HIST1H1B",
  "MKI67", "TYMS", "PCNA", "UBE2C"
)

genes_present = genes_use[genes_use %in% rownames(obj)]
missing_genes = setdiff(genes_use, genes_present)

genes_present
missing_genes

stopifnot(length(genes_present) > 0)

## Average expression by cluster
avg = AverageExpression(
  obj,
  features = genes_present,
  group.by = "cluster_label",
  assays = "RNA",
  slot = "data"
)$RNA

avg = avg[, target_clusters, drop = FALSE]
avg = avg[genes_present, , drop = FALSE]

colnames(avg) = cluster_name_map[colnames(avg)]

## Row z-score
scale_rows = function(m) {
  m = as.matrix(m)
  m = t(scale(t(m)))
  m[is.na(m)] = 0
  m
}

mat = scale_rows(avg)

mat[mat > 2] = 2
mat[mat < -2] = -2

## Save matrix
write.csv(
  avg,
  "Y9_peak_CART_cluster_average_expression.csv",
  row.names = TRUE
)

write.csv(
  mat,
  "Y9_peak_CART_cluster_average_expression_zscore.csv",
  row.names = TRUE
)

## Plot heatmap
col_fun = colorRamp2(
  c(-2, 0, 2),
  c("#4B4BA8", "white", "#F06449")
)

ht = Heatmap(
  mat,
  name = "Expr (z)",
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  show_row_names = TRUE,
  show_column_names = TRUE,
  
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  
  rect_gp = gpar(col = "white", lwd = 0.6),
  
  heatmap_legend_param = list(
    title = "Expr (z)",
    at = c(-2, -1, 0, 1, 2),
    legend_height = unit(3, "cm")
  )
)

draw(ht)

library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(tidyr)

## Read CAR+ object
obj = readRDS("car_9p_new.rds")

obj = subset(
  obj,
  subset = seurat_clusters %in% c(0, 1, 2, 3)
)

## Keep peak and y9
obj_sub = subset(
  obj,
  subset = sample %in% c("peak", "y9")
)

obj_sub$sample = factor(
  obj_sub$sample,
  levels = c("peak", "y9")
)

table(obj_sub$sample, useNA = "ifany")

## Define target genes
genes_use = c("DNMT3A", "TET2")

genes_present = intersect(genes_use, rownames(obj_sub))
missing_genes = setdiff(genes_use, genes_present)

genes_present
missing_genes

if (length(genes_present) == 0) {
  stop("Target genes are not found in the object.")
}

## Extract RNA counts
DefaultAssay(obj_sub) = "RNA"

if (inherits(obj_sub[["RNA"]], "Assay5")) {
  obj_sub = JoinLayers(obj_sub)
}

expr_mat = GetAssayData(
  obj_sub,
  assay = "RNA",
  layer = "counts"
)[genes_present, , drop = FALSE]

## Convert to long format
df = as.data.frame(t(as.matrix(expr_mat))) |>
  rownames_to_column("cell") |>
  cbind(obj_sub@meta.data[, "sample", drop = FALSE]) |>
  pivot_longer(
    cols = all_of(genes_present),
    names_to = "gene",
    values_to = "count"
  )

df$sample = factor(
  df$sample,
  levels = c("peak", "y9")
)

## Summary statistics
summary_df = df |>
  group_by(gene, sample) |>
  summarise(
    n_cells = n(),
    mean_count = mean(count, na.rm = TRUE),
    median_count = median(count, na.rm = TRUE),
    pct_expr = mean(count > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  )

summary_df

write.csv(
  summary_df,
  "DNMT3A_TET2_counts_summary_peak_vs_y9.csv",
  row.names = FALSE
)

## Wilcoxon test
stat_df = df |>
  group_by(gene) |>
  summarise(
    p_value = wilcox.test(count ~ sample)$p.value,
    .groups = "drop"
  ) |>
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  )

stat_df

write.csv(
  stat_df,
  "DNMT3A_TET2_counts_wilcox_peak_vs_y9.csv",
  row.names = FALSE
)

## Violin plot
p = ggplot(
  df,
  aes(x = sample, y = count, fill = sample)
) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = "white"
  ) +
  facet_wrap(
    ~gene,
    scales = "free_y"
  ) +
  stat_compare_means(
    comparisons = list(c("peak", "y9")),
    method = "wilcox.test",
    label = "p.format"
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = "RNA counts",
    title = "DNMT3A and TET2 counts in peak vs y9 CAR+ cells"
  )

p

ggsave(
  "DNMT3A_TET2_counts_violin_peak_vs_y9.pdf",
  p,
  width = 8,
  height = 5
)

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(tibble)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggplot2)
library(grid)

## Set working directory
setwd("pyscenic_out_G1_C4")

## Read regulon DEG result
deg = read_tsv(
  "C4_regulon_DEG.tsv",
  show_col_types = FALSE
)

colnames(deg)
head(deg)

## Select top regulons enriched in CAR+ T
top3_deg = deg |>
  filter(fdr < 0.05, delta_mean > 0) |>
  arrange(fdr) |>
  slice_head(n = 3)

top3_deg

top3_regs = top3_deg |>
  pull(regulon)

top3_tfs = gsub(
  "\\(\\+\\)",
  "",
  top3_regs
)

top3_regs
top3_tfs

## Read pySCENIC regulon table
regs = read_csv(
  "regulons.csv",
  skip = 2,
  show_col_types = FALSE
)

colnames(regs) = c(
  "TF",
  "MotifID",
  "AUC",
  "NES",
  "MotifSimilarityQvalue",
  "OrthologousIdentity",
  "Annotation",
  "Context",
  "TargetGenes",
  "RankAtMax"
)

head(regs, 3)

## Keep top TFs
regs_top3 = regs |>
  filter(TF %in% top3_tfs)

table(regs_top3$TF)

## Choose best motif for each TF
best_regs = regs_top3 |>
  group_by(TF) |>
  arrange(desc(NES), .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup()

best_regs |>
  select(TF, MotifID, AUC, NES, TargetGenes) |>
  print(n = 20, width = Inf)

missing_tfs = setdiff(
  top3_tfs,
  unique(best_regs$TF)
)

missing_tfs

## Parse target genes
parse_targets = function(tf, s) {
  
  parts = stringr::str_match_all(
    s,
    "\\('([^']+)',\\s*([0-9eE.+-]+)\\)"
  )[[1]]
  
  if (nrow(parts) == 0) {
    return(
      tibble(
        TF = character(),
        target = character(),
        score = numeric()
      )
    )
  }
  
  tibble(
    TF = tf,
    target = parts[, 2],
    score = as.numeric(parts[, 3])
  )
}

edges_scored = map2_dfr(
  best_regs$TF,
  best_regs$TargetGenes,
  parse_targets
)

dim(edges_scored)
head(edges_scored, 20)
table(edges_scored$TF)

## Build network table
edges_all = edges_scored |>
  distinct(TF, target, .keep_all = TRUE)

shared_targets_all = edges_all |>
  count(target, sort = TRUE) |>
  filter(n >= 2)

nodes_all = tibble(
  name = unique(c(edges_all$TF, edges_all$target))
) |>
  mutate(
    node_type = case_when(
      name %in% edges_all$TF ~ "TF",
      name %in% shared_targets_all$target ~ "SharedTarget",
      TRUE ~ "Target"
    )
  )

dim(edges_all)
dim(nodes_all)
table(nodes_all$node_type)

## Export network files
write.csv(
  edges_all,
  "C4_top3_alltargets_network_edges.csv",
  row.names = FALSE
)

write.csv(
  nodes_all,
  "C4_top3_alltargets_network_nodes.csv",
  row.names = FALSE
)

write.csv(
  shared_targets_all,
  "C4_top3_alltargets_shared_targets.csv",
  row.names = FALSE
)

edges_export = edges_all |>
  rename(from = TF, to = target) |>
  left_join(
    nodes_all |>
      rename(from = name, from_type = node_type),
    by = "from"
  ) |>
  left_join(
    nodes_all |>
      rename(to = name, to_type = node_type),
    by = "to"
  )

write.csv(
  edges_export,
  "C4_top3_alltargets_network_edges_annotated.csv",
  row.names = FALSE
)

## Build graph
g_all = tbl_graph(
  nodes = nodes_all,
  edges = edges_all |>
    select(from = TF, to = target, score),
  directed = TRUE
)

## Plot network
p_all = ggraph(
  g_all,
  layout = "fr",
  niter = 5000
) +
  geom_edge_link(
    aes(width = score),
    alpha = 0.22,
    colour = "grey45",
    arrow = arrow(length = unit(1.8, "mm")),
    end_cap = circle(2.8, "mm"),
    show.legend = FALSE
  ) +
  geom_node_point(
    aes(color = node_type, size = node_type)
  ) +
  geom_node_text(
    aes(label = name, filter = node_type == "Target"),
    repel = TRUE,
    size = 3.2,
    bg.color = "white",
    bg.r = 0.10,
    point.padding = unit(0.22, "lines"),
    box.padding = unit(0.30, "lines"),
    max.overlaps = Inf
  ) +
  geom_node_text(
    aes(label = name, filter = node_type == "SharedTarget"),
    repel = TRUE,
    size = 3.5,
    fontface = "bold",
    bg.color = "white",
    bg.r = 0.12,
    point.padding = unit(0.28, "lines"),
    box.padding = unit(0.35, "lines"),
    max.overlaps = Inf
  ) +
  geom_node_text(
    aes(label = name, filter = node_type == "TF"),
    repel = TRUE,
    size = 5.2,
    fontface = "bold",
    bg.color = "white",
    bg.r = 0.16,
    point.padding = unit(0.40, "lines"),
    box.padding = unit(0.50, "lines"),
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(
      "TF" = "limegreen",
      "Target" = "orchid",
      "SharedTarget" = "deepskyblue"
    )
  ) +
  scale_size_manual(
    values = c(
      "TF" = 8.5,
      "Target" = 3.2,
      "SharedTarget" = 3.8
    )
  ) +
  scale_edge_width(
    range = c(0.2, 1.2)
  ) +
  theme_void() +
  labs(
    title = "Top 3 regulons enriched in CAR+ T"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 20,
      face = "bold"
    )
  )

p_all

ggsave(
  "C4_top3_alltargets_network.pdf",
  p_all,
  width = 12,
  height = 10
)