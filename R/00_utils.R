############################################################
# 00_utils.R
# General utility functions for MMIL / MCEM ablation pipeline
#
# Supports:
#   1. Binary task:
#        Early vs Advanced
#   2. Categorical task:
#        three_stage: I / II_III / IV
#        four_stage:  I / II / III / IV
############################################################

############################################################
# Basic helper
############################################################

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

############################################################
# Package helpers
############################################################

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0(
        "Package '", pkg, "' is required but not installed. ",
        "Please install it before running this script."
      )
    )
  }
}

load_pkg <- function(pkg) {
  require_pkg(pkg)
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

############################################################
# File and directory helpers
############################################################

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  invisible(path)
}

safe_write_csv <- function(x, file, row.names = FALSE) {
  ensure_dir(dirname(file))
  write.csv(x, file = file, row.names = row.names)
  invisible(file)
}

log_msg <- function(...) {
  cat("\n")
  cat(paste0(...))
  cat("\n")
}

print_dim <- function(x, name = "Object") {
  cat("\n", name, " dimensions:\n", sep = "")
  print(dim(x))
}

############################################################
# Data checks
############################################################

check_required_cols <- function(df, required_cols, object_name = "data frame") {
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing required columns in ", object_name, ": ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  invisible(TRUE)
}

stop_if_any_na <- function(x, message) {
  if (any(is.na(x))) {
    stop(message)
  }
  invisible(TRUE)
}

############################################################
# Task helpers
############################################################

normalize_task_type <- function(task_type = "binary") {
  task_type <- tolower(task_type)
  
  if (!(task_type %in% c("binary", "categorical"))) {
    stop("task_type must be either 'binary' or 'categorical'.")
  }
  
  return(task_type)
}

normalize_categorical_scheme <- function(categorical_scheme = "three_stage") {
  categorical_scheme <- tolower(categorical_scheme)
  
  if (!(categorical_scheme %in% c("three_stage", "four_stage"))) {
    stop("categorical_scheme must be either 'three_stage' or 'four_stage'.")
  }
  
  return(categorical_scheme)
}

get_categorical_class_levels <- function(categorical_scheme = "three_stage") {
  categorical_scheme <- normalize_categorical_scheme(categorical_scheme)
  
  if (categorical_scheme == "three_stage") {
    return(c("I", "II_III", "IV"))
  }
  
  if (categorical_scheme == "four_stage") {
    return(c("I", "II", "III", "IV"))
  }
}

get_task_label_col <- function(task_type = "binary") {
  task_type <- normalize_task_type(task_type)
  
  if (task_type == "binary") {
    return("Stage_group")
  }
  
  if (task_type == "categorical") {
    return("Stage_cat")
  }
}

############################################################
# Stage label helpers
############################################################

create_stage_labels <- function(
  meta,
  categorical_scheme = "three_stage"
) {
  check_required_cols(meta, c("Stage"), object_name = "metadata")
  
  categorical_scheme <- normalize_categorical_scheme(categorical_scheme)
  
  ############################################################
  # Broad stage: I / II / III / IV
  ############################################################
  
  meta$Stage_broad <- dplyr::case_when(
    meta$Stage %in% c("IA", "IA2", "IA3", "IB") ~ "I",
    meta$Stage %in% c("IIA", "IIB") ~ "II",
    meta$Stage %in% c("IIIA") ~ "III",
    meta$Stage %in% c("IV") ~ "IV",
    TRUE ~ NA_character_
  )
  
  meta$Stage_broad <- factor(
    meta$Stage_broad,
    levels = c("I", "II", "III", "IV")
  )
  
  ############################################################
  # Binary label: Early / Advanced
  ############################################################
  
  meta$Stage_group <- dplyr::case_when(
    meta$Stage_broad %in% c("I", "II") ~ "Early",
    meta$Stage_broad %in% c("III", "IV") ~ "Advanced",
    TRUE ~ NA_character_
  )
  
  meta$Stage_group <- factor(
    meta$Stage_group,
    levels = c("Early", "Advanced")
  )
  
  ############################################################
  # Categorical label
  ############################################################
  
  if (categorical_scheme == "three_stage") {
    meta$Stage_cat <- dplyr::case_when(
      meta$Stage_broad == "I" ~ "I",
      meta$Stage_broad %in% c("II", "III") ~ "II_III",
      meta$Stage_broad == "IV" ~ "IV",
      TRUE ~ NA_character_
    )
    
    meta$Stage_cat <- factor(
      meta$Stage_cat,
      levels = c("I", "II_III", "IV")
    )
  }
  
  if (categorical_scheme == "four_stage") {
    meta$Stage_cat <- as.character(meta$Stage_broad)
    
    meta$Stage_cat <- factor(
      meta$Stage_cat,
      levels = c("I", "II", "III", "IV")
    )
  }
  
  return(meta)
}

create_binary_label <- function(stage_group) {
  ifelse(stage_group == "Advanced", 1, 0)
}

create_categorical_label <- function(stage_cat, class_levels) {
  stage_cat <- factor(stage_cat, levels = class_levels)
  as.integer(stage_cat) - 1
}

make_onehot <- function(y, class_levels) {
  y <- factor(y, levels = class_levels)
  
  out <- matrix(
    0,
    nrow = length(y),
    ncol = length(class_levels)
  )
  
  colnames(out) <- class_levels
  
  for (k in seq_along(class_levels)) {
    out[, k] <- as.numeric(y == class_levels[k])
  }
  
  return(out)
}

############################################################
# Cell type helper
############################################################

detect_celltype_column <- function(meta) {
  possible_celltype_cols <- c(
    "Cell_type",
    "CellType",
    "cell_type",
    "Major_cell_type",
    "MajorCellType",
    "Annotation",
    "annotation",
    "Cell_subtype",
    "Cell_subtype_1",
    "Cell_type.refined"
  )
  
  celltype_col <- possible_celltype_cols[
    possible_celltype_cols %in% colnames(meta)
  ][1]
  
  if (is.na(celltype_col)) {
    stop(
      paste0(
        "Could not automatically detect a cell type column. ",
        "Please check colnames(meta) and set celltype_col manually."
      )
    )
  }
  
  return(celltype_col)
}

############################################################
# Matrix helpers
############################################################

row_normalize <- function(mat, eps = 1e-12) {
  mat <- pmax(mat, eps)
  mat / rowSums(mat)
}

col_normalize_counts <- function(mat, scale_factor = 10000) {
  lib_size <- colSums(mat)
  
  if (any(lib_size <= 0)) {
    stop("Some cells have zero total counts after gene filtering.")
  }
  
  t(t(mat) / lib_size) * scale_factor
}

safe_log1p <- function(mat) {
  log1p(mat)
}

############################################################
# Train/test split helpers
############################################################

make_patient_split <- function(
  meta,
  label_col = "Stage_group",
  patient_col = "Patient",
  test_fraction = 0.25,
  seed = 2026
) {
  check_required_cols(
    meta,
    c(patient_col, label_col),
    object_name = "metadata"
  )

  set.seed(seed)

  patient_labels <- meta %>%
    dplyr::filter(!is.na(.data[[label_col]])) %>%
    dplyr::distinct(
      .data[[patient_col]],
      .data[[label_col]]
    )

  colnames(patient_labels) <- c("Patient", "Label")

  test_patients <- c()

  for (lab in unique(patient_labels$Label)) {
    pts <- patient_labels$Patient[patient_labels$Label == lab]
    n_pts <- length(pts)

    if (n_pts < 2) {
      warning(
        sprintf(
          "Only %d patient(s) for class '%s'; keeping all in train so glmnet can see this class.",
          n_pts,
          lab
        )
      )
      next
    }

    n_test <- ceiling(test_fraction * n_pts)
    n_test <- max(1, n_test)
    n_test <- min(n_test, n_pts - 1)

    test_patients <- c(
      test_patients,
      sample(pts, size = n_test)
    )
  }

  split <- ifelse(
    meta[[patient_col]] %in% test_patients,
    "test",
    "train"
  )

  train_labels <- unique(meta[[label_col]][split == "train"])
  missing_train_labels <- setdiff(unique(patient_labels$Label), train_labels)

  if (length(missing_train_labels) > 0) {
    stop(
      paste0(
        "Train split is missing class(es): ",
        paste(missing_train_labels, collapse = ", "),
        ". glmnet multinomial cannot be run reliably."
      )
    )
  }

  return(split)
}

############################################################
# Binary metrics
############################################################

log_loss_binary <- function(y, p, eps = 1e-8) {
  p <- pmin(pmax(p, eps), 1 - eps)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

compute_auc_safe <- function(y, score) {
  require_pkg("pROC")
  
  if (length(unique(y)) < 2) {
    return(NA_real_)
  }
  
  if (length(unique(score)) < 2) {
    return(NA_real_)
  }
  
  auc_obj <- pROC::auc(
    response = y,
    predictor = score,
    quiet = TRUE
  )
  
  as.numeric(auc_obj)
}

############################################################
# Multiclass metrics
############################################################

multiclass_log_loss <- function(
  y_true,
  prob_mat,
  class_levels,
  eps = 1e-8
) {
  y_true <- factor(y_true, levels = class_levels)
  
  prob_mat <- as.matrix(prob_mat)
  prob_mat <- prob_mat[, class_levels, drop = FALSE]
  prob_mat <- row_normalize(prob_mat)
  prob_mat <- pmin(pmax(prob_mat, eps), 1 - eps)
  
  idx <- cbind(
    seq_along(y_true),
    as.integer(y_true)
  )
  
  -mean(log(prob_mat[idx]))
}

multiclass_accuracy <- function(
  y_true,
  prob_mat,
  class_levels
) {
  y_true <- factor(y_true, levels = class_levels)
  
  prob_mat <- as.matrix(prob_mat)
  prob_mat <- prob_mat[, class_levels, drop = FALSE]
  
  y_pred <- class_levels[
    max.col(prob_mat, ties.method = "first")
  ]
  
  mean(y_pred == as.character(y_true))
}

macro_f1_score <- function(
  y_true,
  y_pred,
  class_levels
) {
  y_true <- factor(y_true, levels = class_levels)
  y_pred <- factor(y_pred, levels = class_levels)
  
  f1_values <- sapply(
    class_levels,
    function(cls) {
      tp <- sum(y_true == cls & y_pred == cls)
      fp <- sum(y_true != cls & y_pred == cls)
      fn <- sum(y_true == cls & y_pred != cls)
      
      precision <- ifelse(tp + fp == 0, NA_real_, tp / (tp + fp))
      recall <- ifelse(tp + fn == 0, NA_real_, tp / (tp + fn))
      
      if (is.na(precision) || is.na(recall) || precision + recall == 0) {
        return(NA_real_)
      }
      
      2 * precision * recall / (precision + recall)
    }
  )
  
  mean(f1_values, na.rm = TRUE)
}

balanced_accuracy_multiclass <- function(
  y_true,
  y_pred,
  class_levels
) {
  y_true <- factor(y_true, levels = class_levels)
  y_pred <- factor(y_pred, levels = class_levels)
  
  recall_values <- sapply(
    class_levels,
    function(cls) {
      denom <- sum(y_true == cls)
      if (denom == 0) {
        return(NA_real_)
      }
      sum(y_true == cls & y_pred == cls) / denom
    }
  )
  
  mean(recall_values, na.rm = TRUE)
}

############################################################
# Soft-label expansion helpers
############################################################

expand_binary_soft_labels <- function(x, q, eps = 1e-8) {
  q <- pmin(pmax(q, 0), 1)
  n <- nrow(x)
  
  x_expanded <- x[rep(seq_len(n), times = 2), , drop = FALSE]
  y_expanded <- rep(c(0, 1), each = n)
  w_expanded <- c(1 - q, q)
  
  keep <- w_expanded > eps
  
  x_expanded <- x_expanded[keep, , drop = FALSE]
  y_expanded <- y_expanded[keep]
  w_expanded <- w_expanded[keep]
  
  w_expanded <- w_expanded / mean(w_expanded)
  
  out <- list(
    x = x_expanded,
    y = y_expanded,
    weight = w_expanded
  )
  
  return(out)
}

expand_multiclass_soft_labels <- function(
  x,
  q_mat,
  class_levels,
  eps = 1e-8
) {
  q_mat <- as.matrix(q_mat)
  q_mat <- row_normalize(q_mat)
  
  if (is.null(colnames(q_mat))) {
    colnames(q_mat) <- class_levels
  }
  
  q_mat <- q_mat[, class_levels, drop = FALSE]
  
  n <- nrow(x)
  K <- ncol(q_mat)
  
  if (K != length(class_levels)) {
    stop("Number of columns in q_mat must match length(class_levels).")
  }
  
  x_expanded <- x[rep(seq_len(n), times = K), , drop = FALSE]
  y_expanded <- rep(class_levels, each = n)
  w_expanded <- as.vector(q_mat)
  
  keep <- w_expanded > eps
  
  x_expanded <- x_expanded[keep, , drop = FALSE]
  y_expanded <- y_expanded[keep]
  w_expanded <- w_expanded[keep]
  
  w_expanded <- w_expanded / mean(w_expanded)
  
  out <- list(
    x = x_expanded,
    y = y_expanded,
    weight = w_expanded
  )
  
  return(out)
}

############################################################
# Experiment naming helper
############################################################

make_experiment_name <- function(config) {
  task_type <- config$task_type %||% "binary"
  categorical_scheme <- config$categorical_scheme %||% NULL
  
  paste(
    task_type,
    if (!is.null(categorical_scheme) && task_type == "categorical") {
      categorical_scheme
    } else {
      NULL
    },
    config$feature_method,
    paste0("k", config$n_components),
    config$model_backend,
    if (!is.null(config$glmnet_alpha)) {
      paste0("alpha", config$glmnet_alpha)
    } else {
      NULL
    },
    paste0("seed", config$seed),
    sep = "_"
  )
}