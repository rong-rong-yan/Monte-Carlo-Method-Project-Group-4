############################################################
# 02_feature_extractors.R
# Feature extraction module for PCA, ICA, and NMF.
#
# Important:
#   Feature extraction is fit only on train cells.
#   Test cells are transformed using the train-fitted feature model.
#
# Supported feature methods:
#   1. PCA
#   2. ICA
#   3. NMF
############################################################

load_pkg("dplyr")

############################################################
# Basic preprocessing
############################################################

filter_genes_train_only <- function(
  train_expr,
  test_expr,
  min_cell_fraction = 0.01
) {
  min_cells <- ceiling(min_cell_fraction * ncol(train_expr))
  
  gene_detected_cells <- rowSums(train_expr > 0)
  keep_genes <- gene_detected_cells >= min_cells
  
  cat("\nGene filtering:\n")
  cat("Minimum detected train cells required:", min_cells, "\n")
  cat("Genes before filtering:", nrow(train_expr), "\n")
  cat("Genes after filtering:", sum(keep_genes), "\n")
  
  train_filt <- train_expr[keep_genes, , drop = FALSE]
  test_filt <- test_expr[keep_genes, , drop = FALSE]
  
  return(
    list(
      train_expr = train_filt,
      test_expr = test_filt,
      keep_genes = rownames(train_expr)[keep_genes]
    )
  )
}

normalize_log_train_test <- function(
  train_expr,
  test_expr,
  scale_factor = 10000
) {
  train_norm <- col_normalize_counts(
    train_expr,
    scale_factor = scale_factor
  )
  
  test_norm <- col_normalize_counts(
    test_expr,
    scale_factor = scale_factor
  )
  
  train_log <- safe_log1p(train_norm)
  test_log <- safe_log1p(test_norm)
  
  return(
    list(
      train_log = train_log,
      test_log = test_log
    )
  )
}

select_hvg_train_only <- function(
  train_log,
  test_log,
  n_hvg = 2000
) {
  gene_means <- rowMeans(train_log)
  gene_vars <- apply(train_log, 1, var)
  
  hvg_df <- data.frame(
    gene = rownames(train_log),
    mean = gene_means,
    variance = gene_vars,
    dispersion = gene_vars / (gene_means + 1e-8),
    stringsAsFactors = FALSE
  )
  
  hvg_df <- hvg_df %>%
    dplyr::filter(is.finite(dispersion)) %>%
    dplyr::arrange(dplyr::desc(dispersion))
  
  n_select <- min(n_hvg, nrow(hvg_df))
  hvg_genes <- hvg_df$gene[seq_len(n_select)]
  
  cat("\nHVG selection:\n")
  cat("Requested HVGs:", n_hvg, "\n")
  cat("Selected HVGs:", length(hvg_genes), "\n")
  
  train_hvg <- train_log[hvg_genes, , drop = FALSE]
  test_hvg <- test_log[hvg_genes, , drop = FALSE]
  
  return(
    list(
      train_hvg = train_hvg,
      test_hvg = test_hvg,
      hvg_genes = hvg_genes,
      hvg_stats = hvg_df
    )
  )
}

preprocess_for_feature_extraction <- function(
  train_expr,
  test_expr,
  min_cell_fraction = 0.01,
  n_hvg = 2000,
  scale_factor = 10000
) {
  filtered <- filter_genes_train_only(
    train_expr = train_expr,
    test_expr = test_expr,
    min_cell_fraction = min_cell_fraction
  )
  
  normalized <- normalize_log_train_test(
    train_expr = filtered$train_expr,
    test_expr = filtered$test_expr,
    scale_factor = scale_factor
  )
  
  hvg <- select_hvg_train_only(
    train_log = normalized$train_log,
    test_log = normalized$test_log,
    n_hvg = n_hvg
  )
  
  return(
    list(
      train_hvg = hvg$train_hvg,
      test_hvg = hvg$test_hvg,
      keep_genes = filtered$keep_genes,
      hvg_genes = hvg$hvg_genes,
      hvg_stats = hvg$hvg_stats
    )
  )
}

############################################################
# Scaling helper
############################################################

scale_train_test_by_train <- function(train_x, test_x, eps = 1e-8) {
  train_center <- colMeans(train_x)
  train_scale <- apply(train_x, 2, sd)
  train_scale[train_scale < eps] <- 1
  
  train_scaled <- sweep(train_x, 2, train_center, "-")
  train_scaled <- sweep(train_scaled, 2, train_scale, "/")
  
  test_scaled <- sweep(test_x, 2, train_center, "-")
  test_scaled <- sweep(test_scaled, 2, train_scale, "/")
  
  return(
    list(
      train_scaled = train_scaled,
      test_scaled = test_scaled,
      center = train_center,
      scale = train_scale
    )
  )
}

############################################################
# PCA feature extraction
############################################################

fit_transform_PCA <- function(
  train_hvg,
  test_hvg,
  n_components = 20
) {
  cat("\nFitting PCA feature extractor:\n")
  cat("Requested number of components:", n_components, "\n")
  
  # prcomp expects observations x features.
  # Here observations are cells and features are HVG genes.
  train_x <- t(train_hvg)
  test_x <- t(test_hvg)
  
  max_components <- min(nrow(train_x), ncol(train_x))
  n_components_eff <- min(n_components, max_components)
  
  if (n_components_eff < n_components) {
    warning(
      paste0(
        "Requested PCA components = ", n_components,
        ", but maximum possible components = ", max_components,
        ". Using n_components = ", n_components_eff, "."
      )
    )
  }
  
  pca_fit <- prcomp(
    train_x,
    center = TRUE,
    scale. = TRUE,
    rank. = n_components_eff
  )
  
  train_features <- pca_fit$x[, seq_len(n_components_eff), drop = FALSE]
  
  # Project test cells into the PCA space fitted on train cells.
  test_scaled <- scale(
    test_x,
    center = pca_fit$center,
    scale = pca_fit$scale
  )
  
  test_features <- test_scaled %*% pca_fit$rotation[, seq_len(n_components_eff), drop = FALSE]
  
  feature_names <- paste0("PC", seq_len(n_components_eff))
  colnames(train_features) <- feature_names
  colnames(test_features) <- feature_names
  
  rownames(train_features) <- colnames(train_hvg)
  rownames(test_features) <- colnames(test_hvg)
  
  return(
    list(
      x_train = as.matrix(train_features),
      x_test = as.matrix(test_features),
      feature_names = feature_names,
      feature_model = pca_fit,
      preprocessing = list(
        center = pca_fit$center,
        scale = pca_fit$scale,
        rotation = pca_fit$rotation
      )
    )
  )
}

############################################################
# ICA feature extraction
############################################################

fit_transform_ICA <- function(
  train_hvg,
  test_hvg,
  n_components = 20,
  seed = 2026
) {
  require_pkg("fastICA")
  
  set.seed(seed)
  
  cat("\nFitting ICA feature extractor:\n")
  cat("Requested number of components:", n_components, "\n")
  
  # fastICA expects observations x features.
  # Here observations are cells and features are HVG genes.
  train_x <- t(train_hvg)
  test_x <- t(test_hvg)
  
  max_components <- min(nrow(train_x), ncol(train_x))
  n_components_eff <- min(n_components, max_components)
  
  if (n_components_eff < n_components) {
    warning(
      paste0(
        "Requested ICA components = ", n_components,
        ", but maximum possible components = ", max_components,
        ". Using n_components = ", n_components_eff, "."
      )
    )
  }
  
  scaled <- scale_train_test_by_train(
    train_x = train_x,
    test_x = test_x
  )
  
  ica_fit <- fastICA::fastICA(
    X = scaled$train_scaled,
    n.comp = n_components_eff,
    alg.typ = "parallel",
    fun = "logcosh",
    alpha = 1,
    method = "C",
    row.norm = FALSE,
    maxit = 200,
    tol = 1e-04,
    verbose = FALSE
  )
  
  train_features <- ica_fit$S
  
  # Project test data using the train-fitted ICA transformation.
  test_features <- scaled$test_scaled %*% ica_fit$K %*% ica_fit$W
  
  feature_names <- paste0("IC", seq_len(n_components_eff))
  colnames(train_features) <- feature_names
  colnames(test_features) <- feature_names
  
  rownames(train_features) <- colnames(train_hvg)
  rownames(test_features) <- colnames(test_hvg)
  
  return(
    list(
      x_train = as.matrix(train_features),
      x_test = as.matrix(test_features),
      feature_names = feature_names,
      feature_model = ica_fit,
      preprocessing = list(
        center = scaled$center,
        scale = scaled$scale
      )
    )
  )
}

############################################################
# NMF feature extraction
############################################################

project_nmf_test <- function(basis_mat, test_hvg) {
  require_pkg("nnls")
  
  # basis_mat: genes x rank
  # test_hvg: genes x test cells
  n_test <- ncol(test_hvg)
  rank <- ncol(basis_mat)
  
  coef_mat <- matrix(
    NA_real_,
    nrow = rank,
    ncol = n_test
  )
  
  colnames(coef_mat) <- colnames(test_hvg)
  rownames(coef_mat) <- paste0("NMF", seq_len(rank))
  
  for (j in seq_len(n_test)) {
    fit_j <- nnls::nnls(
      A = basis_mat,
      b = test_hvg[, j]
    )
    coef_mat[, j] <- fit_j$x
  }
  
  return(coef_mat)
}

fit_transform_NMF <- function(
  train_hvg,
  test_hvg,
  n_components = 20,
  seed = 2026,
  nmf_method = "brunet",
  nrun = 1
) {
  require_pkg("NMF")
  require_pkg("nnls")
  
  set.seed(seed)
  
  cat("\nFitting NMF feature extractor:\n")
  cat("Requested number of components:", n_components, "\n")
  
  max_components <- min(nrow(train_hvg), ncol(train_hvg))
  n_components_eff <- min(n_components, max_components)
  
  if (n_components_eff < n_components) {
    warning(
      paste0(
        "Requested NMF components = ", n_components,
        ", but maximum possible components = ", max_components,
        ". Using n_components = ", n_components_eff, "."
      )
    )
  }
  
  # NMF requires a non-negative matrix.
  train_hvg <- pmax(train_hvg, 0)
  test_hvg <- pmax(test_hvg, 0)
  
  # NMF expects features x samples.
  # Here features are genes and samples are cells.
  nmf_fit <- NMF::nmf(
    x = train_hvg,
    rank = n_components_eff,
    method = nmf_method,
    nrun = nrun,
    seed = seed,
    .options = "v"
  )
  
  basis_mat <- NMF::basis(nmf_fit)
  train_coef <- NMF::coef(nmf_fit)
  
  test_coef <- project_nmf_test(
    basis_mat = basis_mat,
    test_hvg = test_hvg
  )
  
  train_features <- t(train_coef)
  test_features <- t(test_coef)
  
  feature_names <- paste0("NMF", seq_len(n_components_eff))
  colnames(train_features) <- feature_names
  colnames(test_features) <- feature_names
  
  rownames(train_features) <- colnames(train_hvg)
  rownames(test_features) <- colnames(test_hvg)
  
  return(
    list(
      x_train = as.matrix(train_features),
      x_test = as.matrix(test_features),
      feature_names = feature_names,
      feature_model = nmf_fit,
      preprocessing = list(
        basis = basis_mat
      )
    )
  )
}

############################################################
# Main feature extraction wrapper
############################################################

extract_features <- function(
  train_expr,
  test_expr,
  feature_method = c("PCA", "ICA", "NMF"),
  n_components = 20,
  min_cell_fraction = 0.01,
  n_hvg = 2000,
  scale_factor = 10000,
  seed = 2026,
  nmf_method = "brunet",
  nmf_nrun = 1
) {
  feature_method <- match.arg(feature_method)
  
  preprocessed <- preprocess_for_feature_extraction(
    train_expr = train_expr,
    test_expr = test_expr,
    min_cell_fraction = min_cell_fraction,
    n_hvg = n_hvg,
    scale_factor = scale_factor
  )
  
  if (feature_method == "PCA") {
    features <- fit_transform_PCA(
      train_hvg = preprocessed$train_hvg,
      test_hvg = preprocessed$test_hvg,
      n_components = n_components
    )
  }
  
  if (feature_method == "ICA") {
    features <- fit_transform_ICA(
      train_hvg = preprocessed$train_hvg,
      test_hvg = preprocessed$test_hvg,
      n_components = n_components,
      seed = seed
    )
  }
  
  if (feature_method == "NMF") {
    features <- fit_transform_NMF(
      train_hvg = preprocessed$train_hvg,
      test_hvg = preprocessed$test_hvg,
      n_components = n_components,
      seed = seed,
      nmf_method = nmf_method,
      nrun = nmf_nrun
    )
  }
  
  features$hvg_genes <- preprocessed$hvg_genes
  features$keep_genes <- preprocessed$keep_genes
  features$hvg_stats <- preprocessed$hvg_stats
  features$feature_method <- feature_method
  features$n_components <- length(features$feature_names)
  features$n_components_requested <- n_components
  
  return(features)
}

############################################################
# Build modeling data frames
############################################################

build_feature_modeling_dfs <- function(
  train_meta,
  test_meta,
  feature_result
) {
  x_train <- feature_result$x_train
  x_test <- feature_result$x_test
  
  if (!all(rownames(x_train) == train_meta$Cell)) {
    stop("Row names of x_train do not match train_meta$Cell.")
  }
  
  if (!all(rownames(x_test) == test_meta$Cell)) {
    stop("Row names of x_test do not match test_meta$Cell.")
  }
  
  train_feature_df <- data.frame(
    train_meta,
    x_train,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  test_feature_df <- data.frame(
    test_meta,
    x_test,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  model_df <- rbind(train_feature_df, test_feature_df)
  
  return(
    list(
      train_df = train_feature_df,
      test_df = test_feature_df,
      model_df = model_df,
      feature_cols = feature_result$feature_names
    )
  )
}