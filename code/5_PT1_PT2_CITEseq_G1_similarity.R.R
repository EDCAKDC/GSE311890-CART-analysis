library(Seurat)
library(Matrix)
library(glmnet)
library(dplyr)
library(pheatmap)

set.seed(123)

## Read input objects
obj = readRDS("car_9p_new.rds")

obj = subset(
  obj,
  subset = seurat_clusters %in% c(0, 1, 2, 3)
)

obj_G1 = subset(
  obj,
  subset = Phase == "G1"
)

obj2 = readRDS("PT1Y9_filtered_with_AUCell.RDS")

obj2_G1 = subset(
  obj2,
  subset = Phase == "G1"
)

ref_obj = obj_G1
qry_obj = obj2_G1

## Prepare Seurat object
prepare_seurat_for_lr = function(obj, assay = "RNA") {
  
  DefaultAssay(obj) = assay
  
  if (inherits(obj[[assay]], "Assay5")) {
    obj = JoinLayers(obj)
  }
  
  dat = tryCatch(
    LayerData(obj, assay = assay, layer = "data"),
    error = function(e) NULL
  )
  
  if (is.null(dat) || nrow(dat) == 0 || ncol(dat) == 0) {
    obj = NormalizeData(obj)
  }
  
  obj
}

## Extract expression matrix
get_lr_matrix = function(obj, assay = "RNA", layer = "data") {
  
  DefaultAssay(obj) = assay
  
  LayerData(
    obj,
    assay = assay,
    layer = layer
  )
}

## Define labels
keep_labels = c(
  "CAR T: CD4+",
  "Non-CAR T: CD4+",
  "Non-CAR T: CD8+"
)

sample_order = c("peak", "y9")

ref_obj = prepare_seurat_for_lr(ref_obj, assay = "RNA")
qry_obj = prepare_seurat_for_lr(qry_obj, assay = "RNA")

table(ref_obj$sample, useNA = "ifany")
table(qry_obj$tcell.type.labels, useNA = "ifany")

## Keep target groups
ref_obj = subset(
  ref_obj,
  subset = !is.na(sample) & sample %in% sample_order
)

qry_obj = subset(
  qry_obj,
  subset = !is.na(tcell.type.labels) & tcell.type.labels %in% keep_labels
)

ref_obj$sample = factor(
  ref_obj$sample,
  levels = sample_order
)

qry_obj$tcell.type.labels = factor(
  qry_obj$tcell.type.labels,
  levels = keep_labels
)

table(ref_obj$sample, useNA = "ifany")
table(qry_obj$tcell.type.labels, useNA = "ifany")

ncol(ref_obj)
ncol(qry_obj)

## Build matrices
ref_mat = get_lr_matrix(ref_obj, assay = "RNA", layer = "data")
qry_mat = get_lr_matrix(qry_obj, assay = "RNA", layer = "data")

common_genes = intersect(
  rownames(ref_mat),
  rownames(qry_mat)
)

## Remove noisy genes
common_genes = common_genes[!grepl("^MT-", common_genes)]
common_genes = common_genes[!grepl("^RPL|^RPS", common_genes)]
common_genes = common_genes[!grepl("^TRAV|^TRAJ|^TRBV|^TRBJ|^TRBC|^TRAC", common_genes)]
common_genes = common_genes[!grepl("^IGH|^IGK|^IGL", common_genes)]

## Use reference HVGs
ref_obj = FindVariableFeatures(
  ref_obj,
  selection.method = "vst",
  nfeatures = 3000
)

features_use = intersect(
  VariableFeatures(ref_obj),
  common_genes
)

length(features_use)

ref_mat = ref_mat[features_use, , drop = FALSE]
qry_mat = qry_mat[features_use, , drop = FALSE]

x_ref = t(as.matrix(ref_mat))
x_qry = t(as.matrix(qry_mat))

y_ref = ref_obj$sample

## Keep valid training cells
keep_cells = colnames(ref_obj)[
  !is.na(y_ref) & y_ref %in% sample_order
]

x_ref_sub = x_ref[keep_cells, , drop = FALSE]

y_ref_sub = y_ref[
  match(keep_cells, colnames(ref_obj))
]

y_ref_sub = factor(
  y_ref_sub,
  levels = sample_order
)

table(y_ref_sub)

## Train logistic regression
cvfit = cv.glmnet(
  x = x_ref_sub,
  y = y_ref_sub,
  family = "multinomial",
  alpha = 0,
  type.measure = "class",
  nfolds = 5
)

pdf("obj_G1_reference_logreg_cv_curve.pdf", width = 7, height = 5)
plot(cvfit)
dev.off()

## Predict query cells
pred_prob = predict(
  cvfit,
  newx = x_qry,
  s = "lambda.min",
  type = "response"
)

if (is.list(pred_prob)) {
  pred_prob_mat = do.call(cbind, pred_prob)
} else if (length(dim(pred_prob)) == 3) {
  pred_prob_mat = pred_prob[, , 1]
} else {
  pred_prob_mat = pred_prob
}

pred_prob_mat = as.matrix(pred_prob_mat)
rownames(pred_prob_mat) = colnames(qry_obj)

pred_label = colnames(pred_prob_mat)[
  max.col(pred_prob_mat, ties.method = "first")
]

pred_score = apply(
  pred_prob_mat,
  1,
  max
)

qry_obj$predicted_sample = pred_label
qry_obj$predicted_score = pred_score

for (lab in colnames(pred_prob_mat)) {
  qry_obj[[paste0("prob_", lab)]] = pred_prob_mat[, lab]
}

table(qry_obj$predicted_sample, useNA = "ifany")

## Mean probability by query label
all_prob_cols = grep(
  "^prob_",
  colnames(qry_obj@meta.data),
  value = TRUE
)

prob_map = setNames(
  all_prob_cols,
  sub("^prob_", "", all_prob_cols)
)

meta_df = as.data.frame(qry_obj@meta.data)

celltype_prob_summary = meta_df |>
  mutate(
    tcell.type.labels = as.character(tcell.type.labels)
  ) |>
  filter(
    !is.na(tcell.type.labels),
    tcell.type.labels %in% keep_labels
  ) |>
  select(
    tcell.type.labels,
    all_of(prob_map[sample_order])
  ) |>
  group_by(tcell.type.labels) |>
  summarise(
    across(everything(), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

celltype_prob_df = as.data.frame(celltype_prob_summary)

rownames(celltype_prob_df) = celltype_prob_df$tcell.type.labels
celltype_prob_df$tcell.type.labels = NULL

colnames(celltype_prob_df) = sub(
  "^prob_",
  "",
  colnames(celltype_prob_df)
)

celltype_prob_df = celltype_prob_df[
  keep_labels,
  sample_order,
  drop = FALSE
]

celltype_prob_df

write.csv(
  celltype_prob_df,
  "obj2_G1_selected_tcell_labels_mean_probability.csv",
  row.names = TRUE
)

## Convert probability to logit
eps = 1e-6

celltype_logit_df = log(
  (celltype_prob_df + eps) / (1 - celltype_prob_df + eps)
)

celltype_logit_df

write.csv(
  celltype_logit_df,
  "obj2_G1_selected_tcell_labels_mean_logit.csv",
  row.names = TRUE
)

## Plot clipped logit heatmap
mat_plot = as.matrix(celltype_logit_df)

mat_plot[mat_plot > 5] = 5
mat_plot[mat_plot < -5] = -5

bk = seq(-5, 5, length.out = 101)

pdf(
  "obj2_G1_selected_tcell_labels_logit_heatmap_clipped.pdf",
  width = 5,
  height = 4
)

pheatmap(
  mat = mat_plot,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  border_color = "grey80",
  main = "Predicted similarity (G1 logit)",
  angle_col = 0,
  breaks = bk
)

dev.off()

## Plot raw logit heatmap
pdf(
  "obj2_G1_selected_tcell_labels_logit_heatmap_raw.pdf",
  width = 5,
  height = 4
)

pheatmap(
  mat = as.matrix(celltype_logit_df),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  border_color = "grey80",
  main = "Predicted similarity (G1 logit)",
  angle_col = 0
)

dev.off()

library(Seurat)
library(Matrix)
library(glmnet)
library(dplyr)
library(pheatmap)

set.seed(123)

## Read barcode annotation
barcode_metadata = read.table(
  "CARTcell_barcodes.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

## Read query object
g1_qry = readRDS("G1_four_groups_with_logreg_predictions.rds")

## Add barcode ID
g1_meta = g1_qry@meta.data
g1_meta$CellID_current = rownames(g1_meta)

g1_meta$CellID = g1_meta$CellID_current
g1_meta$CellID = gsub("CARPALL_G1_", "", g1_meta$CellID, fixed = TRUE)
g1_meta$CellID = gsub("YEAR9_G1_", "", g1_meta$CellID, fixed = TRUE)

nrow(g1_meta)
nrow(barcode_metadata)
sum(g1_meta$CellID %in% barcode_metadata$CellID)

## Merge author annotation
g1_meta2 = g1_meta |>
  left_join(
    barcode_metadata[, c("CellID", "CellType1", "CellType2", "CARTcell")],
    by = "CellID"
  )

rownames(g1_meta2) = g1_meta2$CellID_current
g1_qry@meta.data = g1_meta2

table(is.na(g1_qry$CellType1))
table(g1_qry$CellType2, useNA = "ifany")[1:20]

head(
  g1_qry@meta.data[, c("CellID_current", "CellID", "CellType1", "CellType2", "CARTcell")],
  10
)

## Keep query cells with CellType2
qry_obj = subset(
  g1_qry,
  subset = !is.na(CellType2)
)

ncol(ref_obj)
ncol(qry_obj)

## Prepare Seurat object
prepare_seurat_for_lr = function(obj, assay = "RNA") {
  
  DefaultAssay(obj) = assay
  
  if (inherits(obj[[assay]], "Assay5")) {
    obj = JoinLayers(obj)
  }
  
  dat = tryCatch(
    LayerData(obj, assay = assay, layer = "data"),
    error = function(e) NULL
  )
  
  if (is.null(dat) || nrow(dat) == 0 || ncol(dat) == 0) {
    obj = NormalizeData(obj)
  }
  
  obj
}

ref_obj = prepare_seurat_for_lr(ref_obj, assay = "RNA")
qry_obj = prepare_seurat_for_lr(qry_obj, assay = "RNA")

## Reference labels
label_col = "sample"

table(ref_obj[[label_col]][, 1], useNA = "ifany")

## Build matrices
ref_mat = GetAssayData(
  ref_obj,
  assay = "RNA",
  layer = "data"
)

qry_mat = GetAssayData(
  qry_obj,
  assay = "RNA",
  layer = "data"
)

common_genes = intersect(
  rownames(ref_mat),
  rownames(qry_mat)
)

## Remove noisy genes
common_genes = common_genes[!grepl("^MT-", common_genes)]
common_genes = common_genes[!grepl("^RPL|^RPS", common_genes)]
common_genes = common_genes[!grepl("^TRAV|^TRAJ|^TRBV|^TRBJ|^TRBC|^TRAC", common_genes)]
common_genes = common_genes[!grepl("^IGH|^IGK|^IGL", common_genes)]

## Use reference HVGs
ref_obj = FindVariableFeatures(
  ref_obj,
  selection.method = "vst",
  nfeatures = 3000
)

features_use = intersect(
  VariableFeatures(ref_obj),
  common_genes
)

length(features_use)

ref_mat = ref_mat[features_use, , drop = FALSE]
qry_mat = qry_mat[features_use, , drop = FALSE]

x_ref = t(as.matrix(ref_mat))
x_qry = t(as.matrix(qry_mat))

y_ref = ref_obj[[label_col]][, 1]

## Downsample reference cells
cells_by_class = split(
  colnames(ref_obj),
  y_ref
)

target_n = min(
  min(sapply(cells_by_class, length)),
  5000
)

train_cells = unlist(
  lapply(cells_by_class, function(v) {
    sample(v, size = min(length(v), target_n))
  })
)

x_ref_sub = x_ref[train_cells, , drop = FALSE]

y_ref_sub = y_ref[
  match(train_cells, colnames(ref_obj))
]

table(y_ref_sub)

## Train logistic regression
cvfit = cv.glmnet(
  x = x_ref_sub,
  y = y_ref_sub,
  family = "multinomial",
  alpha = 0,
  type.measure = "class",
  nfolds = 5
)

pdf("G1_sample_reference_logreg_cv_curve.pdf", width = 7, height = 5)
plot(cvfit)
dev.off()

## Predict query cells
pred_prob = predict(
  cvfit,
  newx = x_qry,
  s = "lambda.min",
  type = "response"
)

if (is.list(pred_prob)) {
  pred_prob_mat = do.call(cbind, pred_prob)
} else if (length(dim(pred_prob)) == 3) {
  pred_prob_mat = pred_prob[, , 1]
} else {
  pred_prob_mat = pred_prob
}

pred_prob_mat = as.matrix(pred_prob_mat)
rownames(pred_prob_mat) = colnames(qry_obj)

pred_label = colnames(pred_prob_mat)[
  max.col(pred_prob_mat, ties.method = "first")
]

pred_score = apply(
  pred_prob_mat,
  1,
  max
)

qry_obj$predicted_sample = pred_label
qry_obj$predicted_score = pred_score

for (lab in colnames(pred_prob_mat)) {
  qry_obj[[paste0("prob_", lab)]] = pred_prob_mat[, lab]
}

table(qry_obj$predicted_sample, useNA = "ifany")

## Mean probability by CellType2
prob_cols_now = paste0(
  "prob_",
  sort(unique(as.character(ref_obj$sample)))
)

prob_cols_now

celltype_prob_summary = qry_obj@meta.data |>
  mutate(
    CellID = rownames(qry_obj@meta.data)
  ) |>
  select(
    CellID,
    CellType2,
    all_of(prob_cols_now)
  ) |>
  group_by(CellType2) |>
  summarise(
    across(all_of(prob_cols_now), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

celltype_prob_df = as.data.frame(celltype_prob_summary)

rownames(celltype_prob_df) = celltype_prob_df$CellType2
celltype_prob_df$CellType2 = NULL

colnames(celltype_prob_df) = sub(
  "^prob_",
  "",
  colnames(celltype_prob_df)
)

celltype_prob_df = celltype_prob_df[, c("peak", "y9"), drop = FALSE]

celltype_prob_df

write.csv(
  celltype_prob_df,
  "G1_CellType2_mean_probability_by_sample_reference.csv",
  row.names = TRUE
)

## Convert probability to logit
eps = 1e-6

celltype_logit_df = log(
  (celltype_prob_df + eps) / (1 - celltype_prob_df + eps)
)

write.csv(
  celltype_logit_df,
  "G1_CellType2_mean_logit_by_sample_reference.csv",
  row.names = TRUE
)

## Set row order
desired_order = c(
  "Product:TCF7+ Tcm",
  "Late:DN TEM",
  "Early:DN gdT",
  "Early:CD8 CTL",
  "Early:CD8 TEM",
  "Mid:CD8 TEM",
  "Early:DN TEM",
  "Late:Proliferating CAR T",
  "Late:CD8 TEM",
  "Late:CD8 CTL",
  "Late:CD4 TEM",
  "Late:CD4 CTL",
  "Late:CD8 gdT",
  "Late:CD4 Naive/TCM",
  "Late:CD8 Naive/TCM",
  "Late:DN gdT",
  "Late:TCF7+ Tcm",
  "Mid:DN TEM",
  "Mid:CD8 CTL",
  "Mid:CD4 TEM",
  "Mid:CD8 gdT",
  "Mid:DN gdT",
  "Mid:CD4 CTL",
  "Mid:CD4 Naive/TCM",
  "Mid:CD8 Naive/TCM",
  "Mid:Proliferating CAR T",
  "Mid:TCF7+ Tcm",
  "Early:CD4 TEM",
  "Early:CD4 CTL",
  "Early:CD8 gdT",
  "Early:CD8 Naive/TCM",
  "Early:Proliferating CAR T",
  "Early:TCF7+ Tcm"
)

desired_order = desired_order[
  desired_order %in% rownames(celltype_logit_df)
]

remaining_rows = setdiff(
  rownames(celltype_logit_df),
  desired_order
)

final_row_order = c(
  desired_order,
  remaining_rows
)

celltype_logit_df = celltype_logit_df[
  final_row_order,
  ,
  drop = FALSE
]

## Plot heatmap
pdf(
  "G1_CellType2_by_sample_reference_logit_heatmap.pdf",
  width = 6,
  height = 10
)

pheatmap(
  mat = as.matrix(celltype_logit_df),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  border_color = "grey80",
  main = "Predicted similarity (logit, G1)",
  angle_col = 0
)

dev.off()