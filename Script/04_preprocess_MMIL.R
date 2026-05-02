############################################################
# Stage-based pilot dataset for MMIL / MCEM
# Dataset: GSE131907 lung adenocarcinoma
############################################################

proj_dir <- "/Users/rongrong/Desktop/JHU/Academics/Course/大三/大三下/Monte Carlo Methods/Project/"
setwd(proj_dir)

library(data.table)
library(dplyr)
library(Matrix)

############################################################
# 1. Load pilot expression subset and metadata
############################################################

expr_file <- file.path(proj_dir, "GSE131907_pilot_expression_subset_dt.rds")
meta_file <- file.path(proj_dir, "GSE131907_pilot_sampled_cell_metadata_with_clinical.csv")

expr_dt <- readRDS(expr_file)

meta <- read.csv(
  meta_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Expression subset dimensions:\n")
print(dim(expr_dt))

cat("\nMetadata dimensions:\n")
print(dim(meta))

cat("\nMetadata columns:\n")
print(colnames(meta))

############################################################
# 2. Convert expression data.table to matrix
############################################################

# First column is gene name.
gene_names <- expr_dt[[1]]

# Remaining columns are cell IDs.
cell_ids <- colnames(expr_dt)[-1]

# Convert expression values to matrix: genes x cells
expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
rownames(expr_mat) <- gene_names
colnames(expr_mat) <- cell_ids

cat("\nExpression matrix dimensions, genes x cells:\n")
print(dim(expr_mat))

############################################################
# 3. Match metadata to expression columns
############################################################

meta <- meta[match(cell_ids, meta$Index), ]

cat("\nNumber of unmatched cells after metadata match:\n")
print(sum(is.na(meta$Index)))

cat("\nCheck first few matched IDs:\n")
print(head(data.frame(expr_cell = cell_ids, meta_cell = meta$Index), 10))

############################################################
# 4. Create broad stage and Early/Advanced labels
############################################################

meta <- meta %>%
  mutate(
    Stage_broad = case_when(
      Stage %in% c("IA", "IA2", "IA3", "IB") ~ "I",
      Stage %in% c("IIA", "IIB") ~ "II",
      Stage %in% c("IIIA") ~ "III",
      Stage %in% c("IV") ~ "IV",
      TRUE ~ NA_character_
    ),
    Stage_group = case_when(
      Stage_broad %in% c("I", "II") ~ "Early",
      Stage_broad %in% c("III", "IV") ~ "Advanced",
      TRUE ~ NA_character_
    )
  )

meta$Stage_broad <- factor(meta$Stage_broad, levels = c("I", "II", "III", "IV"))
meta$Stage_group <- factor(meta$Stage_group, levels = c("Early", "Advanced"))

cat("\nCell counts by detailed Stage:\n")
print(table(meta$Stage, useNA = "ifany"))

cat("\nCell counts by broad Stage:\n")
print(table(meta$Stage_broad, useNA = "ifany"))

cat("\nCell counts by Stage_group:\n")
print(table(meta$Stage_group, useNA = "ifany"))

cat("\nPatient counts by Stage_group:\n")
print(
  meta %>%
    distinct(Patient, Stage_group) %>%
    count(Stage_group, name = "n_patients")
)

############################################################
# 5. Remove cells with missing stage group
############################################################

keep_cells <- !is.na(meta$Stage_group)

expr_mat <- expr_mat[, keep_cells]
meta <- meta[keep_cells, ]

cat("\nAfter filtering missing stage labels:\n")
cat("Expression matrix dimensions:\n")
print(dim(expr_mat))
cat("Metadata dimensions:\n")
print(dim(meta))

############################################################
# 6. Basic gene filtering
############################################################

# Keep genes expressed in at least 1% of pilot cells.
min_cells <- ceiling(0.01 * ncol(expr_mat))

gene_detected_cells <- rowSums(expr_mat > 0)
keep_genes <- gene_detected_cells >= min_cells

expr_mat_filt <- expr_mat[keep_genes, ]

cat("\nNumber of genes before filtering:\n")
print(nrow(expr_mat))

cat("\nNumber of genes after filtering:\n")
print(nrow(expr_mat_filt))

############################################################
# 7. Normalize and log-transform
############################################################

# Library-size normalization to counts per 10,000, then log1p.
lib_size <- colSums(expr_mat_filt)

expr_norm <- t(t(expr_mat_filt) / lib_size) * 10000
expr_log <- log1p(expr_norm)

cat("\nSummary of library sizes after gene filtering:\n")
print(summary(lib_size))

############################################################
# 8. Select highly variable genes
############################################################

gene_means <- rowMeans(expr_log)
gene_vars <- apply(expr_log, 1, var)

hvg_df <- data.frame(
  gene = rownames(expr_log),
  mean = gene_means,
  variance = gene_vars,
  dispersion = gene_vars / (gene_means + 1e-8)
)

# Keep top 2000 highly variable genes, or fewer if not enough genes.
n_hvg <- min(2000, nrow(hvg_df))

hvg_genes <- hvg_df %>%
  arrange(desc(dispersion)) %>%
  slice(1:n_hvg) %>%
  pull(gene)

expr_hvg <- expr_log[hvg_genes, ]

cat("\nNumber of HVGs selected:\n")
print(length(hvg_genes))

############################################################
# 9. PCA: cells x PCs
############################################################

# prcomp expects observations as rows, features as columns.
# So transpose: cells x genes.
set.seed(2026)

pca_fit <- prcomp(
  t(expr_hvg),
  center = TRUE,
  scale. = TRUE,
  rank. = 20
)

pc_mat <- pca_fit$x

cat("\nPC matrix dimensions, cells x PCs:\n")
print(dim(pc_mat))

cat("\nVariance explained by first 10 PCs:\n")
var_explained <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
print(round(var_explained[1:10], 4))

############################################################
# 10. Build modeling dataframe
############################################################

model_df <- data.frame(
  Cell = rownames(pc_mat),
  Patient = meta$Patient,
  Sample = meta$Sample,
  Tissue = meta$Tissue,
  Sample_Origin = meta$Sample_Origin,
  Stage = meta$Stage,
  Stage_broad = meta$Stage_broad,
  Stage_group = meta$Stage_group,
  Cell_type = meta$Cell_type,
  Cell_type_refined = meta$Cell_type.refined,
  Cell_subtype = meta$Cell_subtype,
  pc_mat,
  stringsAsFactors = FALSE
)

cat("\nModel dataframe dimensions:\n")
print(dim(model_df))

cat("\nFirst few rows:\n")
print(head(model_df[, 1:12]))

cat("\nPatient counts by Stage_group in model_df:\n")
print(
  model_df %>%
    distinct(Patient, Stage_group) %>%
    count(Stage_group, name = "n_patients")
)

cat("\nCell counts by Stage_group and Cell_type:\n")
print(table(model_df$Stage_group, model_df$Cell_type))

############################################################
# 11. Save stage-based pilot modeling data
############################################################

write.csv(
  model_df,
  file = file.path(proj_dir, "GSE131907_stage_pilot_PCs_metadata.csv"),
  row.names = FALSE
)

saveRDS(
  pca_fit,
  file = file.path(proj_dir, "GSE131907_stage_pilot_pca_fit.rds")
)

saveRDS(
  expr_hvg,
  file = file.path(proj_dir, "GSE131907_stage_pilot_log_norm_HVG_expression.rds")
)

cat("\nDone. Saved:\n")
cat("GSE131907_stage_pilot_PCs_metadata.csv\n")
cat("GSE131907_stage_pilot_pca_fit.rds\n")
cat("GSE131907_stage_pilot_log_norm_HVG_expression.rds\n")