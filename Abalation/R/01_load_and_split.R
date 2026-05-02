############################################################
# 01_load_and_split.R
# Load expression matrix and metadata, align cells, create labels,
# and split train/test at patient level.
#
# Supports:
#   1. Binary task:
#        Stage_group = Early / Advanced
#        z_obs = 0 / 1
#
#   2. Categorical task:
#        Stage_cat = I / II_III / IV      if three_stage
#        Stage_cat = I / II / III / IV    if four_stage
#        y_cat = 0, 1, 2, ... for model backends such as XGBoost
############################################################

############################################################
# Expected config fields:
#   config$proj_dir
#   config$in_dir
#   config$out_dir
#   config$expr_file
#   config$meta_file
#   config$seed
#   config$test_fraction
#
# Optional task fields:
#   config$task_type = "binary" or "categorical"
#   config$categorical_scheme = "three_stage" or "four_stage"
############################################################

load_pkg("data.table")
load_pkg("dplyr")

############################################################
# Load expression matrix
############################################################

load_expression_matrix <- function(expr_file) {
  log_msg("Loading expression data from: ", expr_file)
  
  expr_dt <- readRDS(expr_file)
  
  print_dim(expr_dt, "Expression data.table")
  
  gene_names <- expr_dt[[1]]
  cell_ids <- colnames(expr_dt)[-1]
  
  expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
  rownames(expr_mat) <- gene_names
  colnames(expr_mat) <- cell_ids
  
  print_dim(expr_mat, "Expression matrix, genes x cells")
  
  return(expr_mat)
}

############################################################
# Load metadata
############################################################

load_metadata <- function(meta_file) {
  log_msg("Loading metadata from: ", meta_file)
  
  meta <- read.csv(
    meta_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  print_dim(meta, "Metadata")
  
  cat("\nMetadata columns:\n")
  print(colnames(meta))
  
  return(meta)
}

############################################################
# Match metadata rows to expression matrix columns
############################################################

match_metadata_to_expression <- function(
  meta,
  expr_mat,
  cell_id_col = "Index"
) {
  check_required_cols(meta, c(cell_id_col), object_name = "metadata")
  
  cell_ids <- colnames(expr_mat)
  
  matched_index <- match(cell_ids, meta[[cell_id_col]])
  n_unmatched <- sum(is.na(matched_index))
  
  cat("\nNumber of unmatched cells after metadata matching:\n")
  print(n_unmatched)
  
  if (n_unmatched > 0) {
    stop(
      paste0(
        "Some expression matrix cell IDs were not found in metadata$",
        cell_id_col,
        ". Please check cell IDs."
      )
    )
  }
  
  meta_matched <- meta[matched_index, , drop = FALSE]
  
  cat("\nFirst few matched cell IDs:\n")
  print(
    head(
      data.frame(
        expr_cell = cell_ids,
        meta_cell = meta_matched[[cell_id_col]]
      ),
      10
    )
  )
  
  meta_matched$Cell <- meta_matched[[cell_id_col]]
  
  return(meta_matched)
}

############################################################
# Filter cells with missing task label
############################################################

filter_missing_task_label_cells <- function(
  expr_mat,
  meta,
  task_type = "binary"
) {
  task_type <- normalize_task_type(task_type)
  label_col <- get_task_label_col(task_type)
  
  check_required_cols(meta, c(label_col), object_name = "metadata")
  
  keep_cells <- !is.na(meta[[label_col]])
  
  expr_mat_filtered <- expr_mat[, keep_cells, drop = FALSE]
  meta_filtered <- meta[keep_cells, , drop = FALSE]
  
  cat("\nAfter filtering cells with missing task label:\n")
  cat("Task type:", task_type, "\n")
  cat("Label column:", label_col, "\n")
  
  print_dim(expr_mat_filtered, "Filtered expression matrix")
  print_dim(meta_filtered, "Filtered metadata")
  
  return(
    list(
      expr_mat = expr_mat_filtered,
      meta = meta_filtered
    )
  )
}

############################################################
# Add observed labels for binary and categorical tasks
############################################################

add_observed_labels <- function(
  meta,
  categorical_scheme = "three_stage"
) {
  categorical_scheme <- normalize_categorical_scheme(categorical_scheme)
  class_levels <- get_categorical_class_levels(categorical_scheme)
  
  ############################################################
  # Binary observed label
  ############################################################
  
  meta$z_obs <- create_binary_label(meta$Stage_group)
  
  ############################################################
  # Categorical observed label
  ############################################################
  
  meta$y_cat <- create_categorical_label(
    stage_cat = meta$Stage_cat,
    class_levels = class_levels
  )
  
  meta$Stage_cat <- factor(
    meta$Stage_cat,
    levels = class_levels
  )
  
  attr(meta, "class_levels") <- class_levels
  
  cat("\nBinary label check: Stage_group by z_obs\n")
  print(table(meta$Stage_group, meta$z_obs, useNA = "ifany"))
  
  cat("\nCategorical label check: Stage_cat by y_cat\n")
  print(table(meta$Stage_cat, meta$y_cat, useNA = "ifany"))
  
  cat("\nPatient counts by Stage_group:\n")
  print(
    meta %>%
      dplyr::distinct(Patient, Stage_group, z_obs) %>%
      dplyr::count(Stage_group, z_obs, name = "n_patients")
  )
  
  cat("\nPatient counts by Stage_cat:\n")
  print(
    meta %>%
      dplyr::distinct(Patient, Stage_cat, y_cat) %>%
      dplyr::count(Stage_cat, y_cat, name = "n_patients")
  )
  
  return(meta)
}

############################################################
# Add train/test split by patient
############################################################

add_train_test_split <- function(
  meta,
  task_type = "binary",
  patient_col = "Patient",
  test_fraction = 0.25,
  seed = 2026
) {
  task_type <- normalize_task_type(task_type)
  label_col <- get_task_label_col(task_type)
  
  meta$split <- make_patient_split(
    meta = meta,
    label_col = label_col,
    patient_col = patient_col,
    test_fraction = test_fraction,
    seed = seed
  )
  
  cat("\nPatient split counts:\n")
  print(
    meta %>%
      dplyr::distinct(
        .data[[patient_col]],
        .data[[label_col]],
        split
      ) %>%
      dplyr::count(split, .data[[label_col]], name = "n_patients")
  )
  
  cat("\nCell split counts:\n")
  print(table(meta$split, meta[[label_col]], useNA = "ifany"))
  
  return(meta)
}

############################################################
# Main data loading and split function
############################################################

load_and_split_data <- function(config) {
  required_config <- c(
    "proj_dir",
    "in_dir",
    "out_dir",
    "expr_file",
    "meta_file",
    "seed",
    "test_fraction"
  )
  
  missing_config <- setdiff(required_config, names(config))
  
  if (length(missing_config) > 0) {
    stop(
      paste(
        "Missing config fields:",
        paste(missing_config, collapse = ", ")
      )
    )
  }
  
  ############################################################
  # Resolve task settings
  ############################################################
  
  task_type <- normalize_task_type(config$task_type %||% "binary")
  
  categorical_scheme <- normalize_categorical_scheme(
    config$categorical_scheme %||% "three_stage"
  )
  
  class_levels <- get_categorical_class_levels(categorical_scheme)
  
  cat("\nData loading task settings:\n")
  cat("task_type:", task_type, "\n")
  cat("categorical_scheme:", categorical_scheme, "\n")
  cat("categorical class levels:", paste(class_levels, collapse = ", "), "\n")
  
  ensure_dir(config$out_dir)
  
  ############################################################
  # 1. Load expression and metadata
  ############################################################
  
  expr_mat <- load_expression_matrix(config$expr_file)
  meta <- load_metadata(config$meta_file)
  
  check_required_cols(
    meta,
    c("Index", "Patient", "Stage"),
    object_name = "metadata"
  )
  
  ############################################################
  # 2. Align metadata to expression columns
  ############################################################
  
  meta <- match_metadata_to_expression(
    meta = meta,
    expr_mat = expr_mat,
    cell_id_col = "Index"
  )
  
  ############################################################
  # 3. Create binary and categorical labels
  ############################################################
  
  meta <- create_stage_labels(
    meta = meta,
    categorical_scheme = categorical_scheme
  )
  
  cat("\nCell counts by detailed Stage:\n")
  print(table(meta$Stage, useNA = "ifany"))
  
  cat("\nCell counts by broad Stage:\n")
  print(table(meta$Stage_broad, useNA = "ifany"))
  
  cat("\nCell counts by binary Stage_group:\n")
  print(table(meta$Stage_group, useNA = "ifany"))
  
  cat("\nCell counts by categorical Stage_cat:\n")
  print(table(meta$Stage_cat, useNA = "ifany"))
  
  ############################################################
  # 4. Filter cells with missing label for selected task
  ############################################################
  
  filtered <- filter_missing_task_label_cells(
    expr_mat = expr_mat,
    meta = meta,
    task_type = task_type
  )
  
  expr_mat <- filtered$expr_mat
  meta <- filtered$meta
  
  ############################################################
  # 5. Add observed labels for both binary and categorical tasks
  ############################################################
  
  meta <- add_observed_labels(
    meta = meta,
    categorical_scheme = categorical_scheme
  )
  
  ############################################################
  # 6. Patient-level train/test split by selected task label
  ############################################################
  
  meta <- add_train_test_split(
    meta = meta,
    task_type = task_type,
    patient_col = "Patient",
    test_fraction = config$test_fraction,
    seed = config$seed
  )
  
  ############################################################
  # 7. Build train/test expression and metadata objects
  ############################################################
  
  train_idx <- meta$split == "train"
  test_idx <- meta$split == "test"
  
  train_expr <- expr_mat[, train_idx, drop = FALSE]
  test_expr <- expr_mat[, test_idx, drop = FALSE]
  
  train_meta <- meta[train_idx, , drop = FALSE]
  test_meta <- meta[test_idx, , drop = FALSE]
  
  ############################################################
  # 8. Build split info
  ############################################################
  
  split_info_cols <- c(
    "Cell",
    "Patient",
    "Stage",
    "Stage_broad",
    "Stage_group",
    "z_obs",
    "Stage_cat",
    "y_cat",
    "split"
  )
  
  split_info_cols <- intersect(split_info_cols, colnames(meta))
  
  split_info <- meta[, split_info_cols, drop = FALSE]
  
  ############################################################
  # 9. Return data object
  ############################################################
  
  return(
    list(
      task_type = task_type,
      categorical_scheme = categorical_scheme,
      class_levels = class_levels,
      
      expr_mat = expr_mat,
      meta = meta,
      
      train_expr = train_expr,
      test_expr = test_expr,
      
      train_meta = train_meta,
      test_meta = test_meta,
      
      split_info = split_info
    )
  )
}