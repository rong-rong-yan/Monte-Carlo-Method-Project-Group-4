############################################################
# 05_patient_level_evaluation.R
# Patient-level aggregation evaluation for binary and categorical
# MMIL / MCEM experiments.
#
# Supported task types:
#   1. binary
#   2. categorical
############################################################

load_pkg("dplyr")
load_pkg("tidyr")
load_pkg("ggplot2")
load_pkg("pROC")
load_pkg("readr")

############################################################
# Shared helper
############################################################

sanitize_class_name_eval <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

get_aggregation_cols_binary <- function() {
  c(
    "mean_prob",
    "median_prob",
    "q90_prob",
    "q99_prob",
    "prop_gt_0.5",
    "prop_gt_0.9"
  )
}

get_aggregation_cols_categorical <- function() {
  c(
    "mean_prob",
    "median_prob",
    "q90_prob",
    "q99_prob"
  )
}

############################################################
# Binary evaluation: load and check prediction file
############################################################

load_binary_prediction_file <- function(pred_file) {
  cat("\nLoading binary cell-level prediction file:\n")
  cat(pred_file, "\n")
  
  cell_pred <- read.csv(
    pred_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat("\nLoaded cell-level prediction file dimensions:\n")
  print(dim(cell_pred))
  
  cat("\nColumn names:\n")
  print(colnames(cell_pred))
  
  required_cols <- c(
    "Patient",
    "split",
    "Stage_group",
    "naive_prob",
    "em_mmil_prob",
    "mcem_mmil_prob"
  )
  
  check_required_cols(
    cell_pred,
    required_cols,
    object_name = "binary cell-level prediction file"
  )
  
  if (!("Cell_type" %in% colnames(cell_pred))) {
    cell_pred$Cell_type <- "Unknown"
  }
  
  if (!("Cell" %in% colnames(cell_pred))) {
    cell_pred$Cell <- paste0("Cell_", seq_len(nrow(cell_pred)))
  }
  
  cat("\nCell counts by split and Stage_group:\n")
  print(table(cell_pred$split, cell_pred$Stage_group))
  
  cat("\nPatient counts by split and Stage_group:\n")
  print(
    cell_pred %>%
      dplyr::distinct(Patient, split, Stage_group) %>%
      dplyr::count(split, Stage_group, name = "n_patients")
  )
  
  return(cell_pred)
}

############################################################
# Binary evaluation: prepare long probability table
############################################################

prepare_binary_probability_long <- function(cell_pred) {
  cell_pred <- cell_pred %>%
    dplyr::mutate(
      y_patient = ifelse(Stage_group == "Advanced", 1, 0)
    )
  
  cat("\nCheck binary patient-level label:\n")
  print(table(cell_pred$Stage_group, cell_pred$y_patient))
  
  prob_long <- cell_pred %>%
    dplyr::select(
      Cell,
      Patient,
      split,
      Stage_group,
      y_patient,
      Cell_type,
      naive_prob,
      em_mmil_prob,
      mcem_mmil_prob
    ) %>%
    tidyr::pivot_longer(
      cols = c(
        naive_prob,
        em_mmil_prob,
        mcem_mmil_prob
      ),
      names_to = "model",
      values_to = "cell_prob"
    ) %>%
    dplyr::mutate(
      model = dplyr::case_when(
        model == "naive_prob" ~ "Naive inherited-label classifier",
        model == "em_mmil_prob" ~ "Deterministic EM-MMIL",
        model == "mcem_mmil_prob" ~ "MCEM-MMIL",
        TRUE ~ model
      )
    )
  
  cat("\nBinary probability long format dimensions:\n")
  print(dim(prob_long))
  
  return(prob_long)
}

############################################################
# Binary evaluation: patient-level aggregation
############################################################

aggregate_binary_patient_scores <- function(prob_long) {
  patient_scores <- prob_long %>%
    dplyr::group_by(
      split,
      Patient,
      Stage_group,
      y_patient,
      model
    ) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      
      mean_prob = mean(cell_prob, na.rm = TRUE),
      median_prob = median(cell_prob, na.rm = TRUE),
      q90_prob = as.numeric(
        quantile(cell_prob, probs = 0.90, na.rm = TRUE)
      ),
      q99_prob = as.numeric(
        quantile(cell_prob, probs = 0.99, na.rm = TRUE)
      ),
      
      prop_gt_0.5 = mean(cell_prob > 0.5, na.rm = TRUE),
      prop_gt_0.9 = mean(cell_prob > 0.9, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  cat("\nBinary patient-level scores dimensions:\n")
  print(dim(patient_scores))
  
  cat("\nFirst few binary patient-level scores:\n")
  print(head(patient_scores, 10))
  
  return(patient_scores)
}

############################################################
# Binary evaluation: AUROC
############################################################

compute_binary_patient_auc <- function(patient_scores) {
  aggregation_cols <- get_aggregation_cols_binary()
  
  auc_results <- patient_scores %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(aggregation_cols),
      names_to = "aggregation",
      values_to = "patient_score"
    ) %>%
    dplyr::group_by(
      split,
      model,
      aggregation
    ) %>%
    dplyr::summarise(
      n_patients = dplyr::n(),
      n_early = sum(y_patient == 0),
      n_advanced = sum(y_patient == 1),
      auroc = compute_auc_safe(y_patient, patient_score),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      split,
      model,
      dplyr::desc(auroc)
    )
  
  cat("\nBinary patient-level AUROC results:\n")
  print(auc_results, n = 100)
  
  return(auc_results)
}

############################################################
# Binary evaluation: log loss
############################################################

compute_binary_patient_logloss <- function(patient_scores) {
  patient_logloss <- patient_scores %>%
    dplyr::group_by(
      split,
      model
    ) %>%
    dplyr::summarise(
      n_patients = dplyr::n(),
      logloss_mean_prob = log_loss_binary(y_patient, mean_prob),
      logloss_median_prob = log_loss_binary(y_patient, median_prob),
      logloss_q90_prob = log_loss_binary(y_patient, q90_prob),
      logloss_q99_prob = log_loss_binary(y_patient, q99_prob),
      .groups = "drop"
    )
  
  cat("\nBinary patient-level log loss results:\n")
  print(patient_logloss)
  
  return(patient_logloss)
}

############################################################
# Binary evaluation: best test AUROC
############################################################

get_binary_best_test_auc <- function(auc_results) {
  best_test_auc <- auc_results %>%
    dplyr::filter(split == "test") %>%
    dplyr::group_by(model) %>%
    dplyr::arrange(
      dplyr::desc(auroc),
      .by_group = TRUE
    ) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  cat("\nBest binary test AUROC aggregation per model:\n")
  print(best_test_auc)
  
  return(best_test_auc)
}

############################################################
# Binary evaluation: cell-type enrichment
############################################################

compute_binary_celltype_enrichment <- function(prob_long) {
  celltype_enrichment <- prob_long %>%
    dplyr::group_by(
      split,
      model,
      Stage_group,
      Cell_type
    ) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      mean_prob = mean(cell_prob, na.rm = TRUE),
      prop_gt_0.5 = mean(cell_prob > 0.5, na.rm = TRUE),
      prop_gt_0.9 = mean(cell_prob > 0.9, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      split,
      model,
      Stage_group,
      dplyr::desc(mean_prob)
    )
  
  cat("\nBinary cell-type summary for biological interpretation:\n")
  print(celltype_enrichment, n = 80)
  
  return(celltype_enrichment)
}

############################################################
# Binary plots
############################################################

plot_binary_patient_scores <- function(patient_scores, out_dir) {
  aggregation_cols <- get_aggregation_cols_binary()
  
  patient_scores_long <- patient_scores %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(aggregation_cols),
      names_to = "aggregation",
      values_to = "patient_score"
    )
  
  p_scores <- ggplot2::ggplot(
    patient_scores_long %>% dplyr::filter(split == "test"),
    ggplot2::aes(
      x = Stage_group,
      y = patient_score,
      fill = Stage_group
    )
  ) +
    ggplot2::geom_boxplot(outlier.size = 1) +
    ggplot2::geom_jitter(
      width = 0.15,
      size = 1.5,
      alpha = 0.8
    ) +
    ggplot2::facet_grid(
      model ~ aggregation,
      scales = "free_y"
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Binary patient-level aggregated scores on test patients",
      x = "Patient-level stage group",
      y = "Aggregated predicted probability"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
    )
  
  plot_file <- file.path(
    out_dir,
    "plot_patient_level_binary_MMIL_scores.png"
  )
  
  ggplot2::ggsave(
    filename = plot_file,
    plot = p_scores,
    width = 14,
    height = 8,
    dpi = 300
  )
  
  return(plot_file)
}

plot_binary_patient_auc <- function(auc_results, out_dir) {
  p_auc <- ggplot2::ggplot(
    auc_results %>% dplyr::filter(split == "test"),
    ggplot2::aes(
      x = aggregation,
      y = auroc,
      fill = model
    )
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Binary patient-level test AUROC by aggregation method",
      x = "Aggregation method",
      y = "Patient-level AUROC",
      fill = "Model"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
    )
  
  plot_file <- file.path(
    out_dir,
    "plot_patient_level_binary_MMIL_AUROC.png"
  )
  
  ggplot2::ggsave(
    filename = plot_file,
    plot = p_auc,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  return(plot_file)
}

plot_binary_celltype_enrichment <- function(celltype_enrichment, out_dir) {
  p_celltype <- celltype_enrichment %>%
    dplyr::filter(split == "test") %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = reorder(Cell_type, prop_gt_0.5),
        y = prop_gt_0.5,
        fill = Stage_group
      )
    ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ model) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Binary proportion of cells with predicted probability > 0.5",
      x = "Cell type",
      y = "Proportion of cells",
      fill = "Stage group"
    )
  
  plot_file <- file.path(
    out_dir,
    "plot_celltype_binary_MMIL_prop_gt_0.5.png"
  )
  
  ggplot2::ggsave(
    filename = plot_file,
    plot = p_celltype,
    width = 12,
    height = 7,
    dpi = 300
  )
  
  return(plot_file)
}

############################################################
# Main binary evaluation function
############################################################

evaluate_binary_predictions <- function(
  pred_file,
  out_dir,
  save_plots = TRUE
) {
  ensure_dir(out_dir)
  
  cell_pred <- load_binary_prediction_file(pred_file)
  prob_long <- prepare_binary_probability_long(cell_pred)
  patient_scores <- aggregate_binary_patient_scores(prob_long)
  auc_results <- compute_binary_patient_auc(patient_scores)
  patient_logloss <- compute_binary_patient_logloss(patient_scores)
  best_test_auc <- get_binary_best_test_auc(auc_results)
  celltype_enrichment <- compute_binary_celltype_enrichment(prob_long)
  
  patient_scores_file <- file.path(
    out_dir,
    "patient_level_binary_MMIL_aggregation_scores.csv"
  )
  
  auc_results_file <- file.path(
    out_dir,
    "patient_level_binary_MMIL_AUROC_results.csv"
  )
  
  patient_logloss_file <- file.path(
    out_dir,
    "patient_level_binary_MMIL_logloss_results.csv"
  )
  
  best_test_auc_file <- file.path(
    out_dir,
    "patient_level_binary_MMIL_best_test_AUROC.csv"
  )
  
  celltype_enrichment_file <- file.path(
    out_dir,
    "celltype_enrichment_binary_MMIL_probabilities.csv"
  )
  
  write.csv(patient_scores, patient_scores_file, row.names = FALSE)
  write.csv(auc_results, auc_results_file, row.names = FALSE)
  write.csv(patient_logloss, patient_logloss_file, row.names = FALSE)
  write.csv(best_test_auc, best_test_auc_file, row.names = FALSE)
  write.csv(celltype_enrichment, celltype_enrichment_file, row.names = FALSE)
  
  plot_files <- list()
  
  if (save_plots) {
    plot_files$patient_scores <- plot_binary_patient_scores(
      patient_scores = patient_scores,
      out_dir = out_dir
    )
    
    plot_files$auc <- plot_binary_patient_auc(
      auc_results = auc_results,
      out_dir = out_dir
    )
    
    plot_files$celltype <- plot_binary_celltype_enrichment(
      celltype_enrichment = celltype_enrichment,
      out_dir = out_dir
    )
  }
  
  cat("\nBinary patient-level evaluation done.\n")
  
  return(
    list(
      task_type = "binary",
      cell_pred = cell_pred,
      prob_long = prob_long,
      patient_scores = patient_scores,
      auc_results = auc_results,
      patient_logloss = patient_logloss,
      best_test_auc = best_test_auc,
      celltype_enrichment = celltype_enrichment,
      files = list(
        patient_scores_file = patient_scores_file,
        auc_results_file = auc_results_file,
        patient_logloss_file = patient_logloss_file,
        best_test_auc_file = best_test_auc_file,
        celltype_enrichment_file = celltype_enrichment_file,
        plot_files = plot_files
      )
    )
  )
}

############################################################
# Categorical evaluation: load and check prediction file
############################################################

get_categorical_prob_columns <- function(class_levels, prefix) {
  paste0(prefix, "_", sanitize_class_name_eval(class_levels))
}

load_categorical_prediction_file <- function(
  pred_file,
  class_levels
) {
  cat("\nLoading categorical cell-level prediction file:\n")
  cat(pred_file, "\n")
  
  cell_pred <- read.csv(
    pred_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat("\nLoaded categorical prediction file dimensions:\n")
  print(dim(cell_pred))
  
  cat("\nColumn names:\n")
  print(colnames(cell_pred))
  
  naive_cols <- get_categorical_prob_columns(class_levels, "naive_prob")
  em_cols <- get_categorical_prob_columns(class_levels, "cat_em_prob")
  mcem_cols <- get_categorical_prob_columns(class_levels, "cat_mcem_prob")
  
  required_cols <- c(
    "Patient",
    "split",
    "Stage_cat",
    naive_cols,
    em_cols,
    mcem_cols
  )
  
  check_required_cols(
    cell_pred,
    required_cols,
    object_name = "categorical cell-level prediction file"
  )
  
  if (!("Cell_type" %in% colnames(cell_pred))) {
    cell_pred$Cell_type <- "Unknown"
  }
  
  if (!("Cell" %in% colnames(cell_pred))) {
    cell_pred$Cell <- paste0("Cell_", seq_len(nrow(cell_pred)))
  }
  
  cell_pred$Stage_cat <- factor(
    cell_pred$Stage_cat,
    levels = class_levels
  )
  
  cat("\nCell counts by split and Stage_cat:\n")
  print(table(cell_pred$split, cell_pred$Stage_cat))
  
  cat("\nPatient counts by split and Stage_cat:\n")
  print(
    cell_pred %>%
      dplyr::distinct(Patient, split, Stage_cat) %>%
      dplyr::count(split, Stage_cat, name = "n_patients")
  )
  
  return(cell_pred)
}

############################################################
# Categorical evaluation: prepare long probability table
############################################################

prepare_categorical_probability_long <- function(
  cell_pred,
  class_levels
) {
  model_specs <- list(
    list(
      model = "Naive inherited-label classifier",
      prefix = "naive_prob"
    ),
    list(
      model = "Categorical EM-MMIL",
      prefix = "cat_em_prob"
    ),
    list(
      model = "Categorical MCEM-MMIL",
      prefix = "cat_mcem_prob"
    )
  )
  
  long_list <- lapply(
    model_specs,
    function(spec) {
      prob_cols <- get_categorical_prob_columns(
        class_levels = class_levels,
        prefix = spec$prefix
      )
      
      tmp <- cell_pred %>%
        dplyr::select(
          Cell,
          Patient,
          split,
          Stage_cat,
          Cell_type,
          dplyr::all_of(prob_cols)
        )
      
      colnames(tmp)[
        match(prob_cols, colnames(tmp))
      ] <- class_levels
      
      tmp_long <- tmp %>%
        tidyr::pivot_longer(
          cols = dplyr::all_of(class_levels),
          names_to = "class",
          values_to = "cell_prob"
        ) %>%
        dplyr::mutate(
          model = spec$model,
          y_true = as.character(Stage_cat)
        )
      
      tmp_long
    }
  )
  
  prob_long <- dplyr::bind_rows(long_list)
  
  prob_long$class <- factor(
    prob_long$class,
    levels = class_levels
  )
  
  prob_long$Stage_cat <- factor(
    prob_long$Stage_cat,
    levels = class_levels
  )
  
  cat("\nCategorical probability long format dimensions:\n")
  print(dim(prob_long))
  
  return(prob_long)
}

############################################################
# Categorical evaluation: patient-level aggregation
############################################################

aggregate_categorical_patient_scores <- function(prob_long) {
  patient_scores <- prob_long %>%
    dplyr::group_by(
      split,
      Patient,
      Stage_cat,
      y_true,
      model,
      class
    ) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      mean_prob = mean(cell_prob, na.rm = TRUE),
      median_prob = median(cell_prob, na.rm = TRUE),
      q90_prob = as.numeric(
        quantile(cell_prob, probs = 0.90, na.rm = TRUE)
      ),
      q99_prob = as.numeric(
        quantile(cell_prob, probs = 0.99, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  
  cat("\nCategorical patient-level scores dimensions:\n")
  print(dim(patient_scores))
  
  cat("\nFirst few categorical patient-level scores:\n")
  print(head(patient_scores, 10))
  
  return(patient_scores)
}

############################################################
# Categorical evaluation: convert patient scores to wide format
############################################################

categorical_patient_scores_wide <- function(
  patient_scores,
  class_levels,
  aggregation
) {
  wide <- patient_scores %>%
    dplyr::select(
      split,
      Patient,
      Stage_cat,
      y_true,
      model,
      class,
      patient_score = dplyr::all_of(aggregation)
    ) %>%
    tidyr::pivot_wider(
      names_from = class,
      values_from = patient_score
    )
  
  missing_classes <- setdiff(class_levels, colnames(wide))
  
  for (cls in missing_classes) {
    wide[[cls]] <- 0
  }
  
  wide <- wide %>%
    dplyr::select(
      split,
      Patient,
      Stage_cat,
      y_true,
      model,
      dplyr::all_of(class_levels)
    )
  
  prob_mat <- as.matrix(wide[, class_levels, drop = FALSE])
  prob_mat <- row_normalize(prob_mat)
  wide[, class_levels] <- prob_mat
  
  return(wide)
}

############################################################
# Categorical evaluation: metrics
############################################################

compute_categorical_patient_metrics <- function(
  patient_scores,
  class_levels
) {
  aggregation_cols <- get_aggregation_cols_categorical()
  
  metric_list <- lapply(
    aggregation_cols,
    function(aggregation) {
      wide <- categorical_patient_scores_wide(
        patient_scores = patient_scores,
        class_levels = class_levels,
        aggregation = aggregation
      )
      
      metrics <- wide %>%
        dplyr::group_by(split, model) %>%
        dplyr::group_modify(
          function(.x, .y) {
            prob_mat <- as.matrix(.x[, class_levels, drop = FALSE])
            prob_mat <- row_normalize(prob_mat)
            
            y_true <- as.character(.x$Stage_cat)
            
            y_pred <- class_levels[
              max.col(prob_mat, ties.method = "first")
            ]
            
            data.frame(
              n_patients = nrow(.x),
              multiclass_logloss = multiclass_log_loss(
                y_true = y_true,
                prob_mat = prob_mat,
                class_levels = class_levels
              ),
              accuracy = mean(y_pred == y_true),
              balanced_accuracy = balanced_accuracy_multiclass(
                y_true = y_true,
                y_pred = y_pred,
                class_levels = class_levels
              ),
              macro_f1 = macro_f1_score(
                y_true = y_true,
                y_pred = y_pred,
                class_levels = class_levels
              ),
              stringsAsFactors = FALSE
            )
          }
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          aggregation = aggregation
        )
      
      metrics
    }
  )
  
  metrics_results <- dplyr::bind_rows(metric_list) %>%
    dplyr::arrange(
      split,
      model,
      multiclass_logloss
    )
  
  cat("\nCategorical patient-level metrics:\n")
  print(metrics_results, n = 100)
  
  return(metrics_results)
}

############################################################
# Categorical evaluation: best test metrics
############################################################

get_categorical_best_test_metrics <- function(metrics_results) {
  best_test_metrics <- metrics_results %>%
    dplyr::filter(split == "test") %>%
    dplyr::group_by(model) %>%
    dplyr::arrange(
      dplyr::desc(accuracy),
      multiclass_logloss,
      .by_group = TRUE
    ) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  cat("\nBest categorical test metrics per model:\n")
  print(best_test_metrics)
  
  return(best_test_metrics)
}

############################################################
# Categorical evaluation: confusion matrix
############################################################

compute_categorical_confusion <- function(
  patient_scores,
  class_levels
) {
  aggregation_cols <- get_aggregation_cols_categorical()
  
  confusion_list <- lapply(
    aggregation_cols,
    function(aggregation) {
      wide <- categorical_patient_scores_wide(
        patient_scores = patient_scores,
        class_levels = class_levels,
        aggregation = aggregation
      )
      
      prob_mat <- as.matrix(wide[, class_levels, drop = FALSE])
      prob_mat <- row_normalize(prob_mat)
      
      wide$predicted_class <- class_levels[
        max.col(prob_mat, ties.method = "first")
      ]
      
      wide %>%
        dplyr::count(
          split,
          model,
          aggregation = aggregation,
          true_class = Stage_cat,
          predicted_class,
          name = "n_patients"
        )
    }
  )
  
  confusion_results <- dplyr::bind_rows(confusion_list)
  
  cat("\nCategorical confusion matrix rows:\n")
  print(confusion_results, n = 100)
  
  return(confusion_results)
}

############################################################
# Categorical evaluation: cell-type enrichment
############################################################

compute_categorical_celltype_enrichment <- function(prob_long) {
  celltype_enrichment <- prob_long %>%
    dplyr::group_by(
      split,
      model,
      Stage_cat,
      Cell_type,
      class
    ) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      mean_prob = mean(cell_prob, na.rm = TRUE),
      q90_prob = as.numeric(
        quantile(cell_prob, probs = 0.90, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      split,
      model,
      Stage_cat,
      class,
      dplyr::desc(mean_prob)
    )
  
  cat("\nCategorical cell-type summary for biological interpretation:\n")
  print(celltype_enrichment, n = 100)
  
  return(celltype_enrichment)
}

############################################################
# Categorical plots
############################################################

plot_categorical_metrics <- function(metrics_results, out_dir) {
  p_acc <- ggplot2::ggplot(
    metrics_results %>% dplyr::filter(split == "test"),
    ggplot2::aes(
      x = aggregation,
      y = accuracy,
      fill = model
    )
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Categorical patient-level test accuracy",
      x = "Aggregation method",
      y = "Accuracy",
      fill = "Model"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
    )
  
  acc_file <- file.path(
    out_dir,
    "plot_patient_level_categorical_accuracy.png"
  )
  
  ggplot2::ggsave(
    filename = acc_file,
    plot = p_acc,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  p_logloss <- ggplot2::ggplot(
    metrics_results %>% dplyr::filter(split == "test"),
    ggplot2::aes(
      x = aggregation,
      y = multiclass_logloss,
      fill = model
    )
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Categorical patient-level test multiclass log loss",
      x = "Aggregation method",
      y = "Multiclass log loss",
      fill = "Model"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
    )
  
  logloss_file <- file.path(
    out_dir,
    "plot_patient_level_categorical_logloss.png"
  )
  
  ggplot2::ggsave(
    filename = logloss_file,
    plot = p_logloss,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  return(
    list(
      accuracy = acc_file,
      logloss = logloss_file
    )
  )
}

############################################################
# Main categorical evaluation function
############################################################

evaluate_categorical_predictions <- function(
  pred_file,
  out_dir,
  class_levels,
  save_plots = TRUE
) {
  ensure_dir(out_dir)
  
  cell_pred <- load_categorical_prediction_file(
    pred_file = pred_file,
    class_levels = class_levels
  )
  
  prob_long <- prepare_categorical_probability_long(
    cell_pred = cell_pred,
    class_levels = class_levels
  )
  
  patient_scores <- aggregate_categorical_patient_scores(prob_long)
  
  metrics_results <- compute_categorical_patient_metrics(
    patient_scores = patient_scores,
    class_levels = class_levels
  )
  
  best_test_metrics <- get_categorical_best_test_metrics(metrics_results)
  
  confusion_results <- compute_categorical_confusion(
    patient_scores = patient_scores,
    class_levels = class_levels
  )
  
  celltype_enrichment <- compute_categorical_celltype_enrichment(prob_long)
  
  patient_scores_file <- file.path(
    out_dir,
    "patient_level_categorical_MMIL_aggregation_scores.csv"
  )
  
  metrics_results_file <- file.path(
    out_dir,
    "patient_level_categorical_MMIL_metrics.csv"
  )
  
  best_test_metrics_file <- file.path(
    out_dir,
    "patient_level_categorical_MMIL_best_test_metrics.csv"
  )
  
  confusion_results_file <- file.path(
    out_dir,
    "patient_level_categorical_MMIL_confusion_matrix.csv"
  )
  
  celltype_enrichment_file <- file.path(
    out_dir,
    "celltype_enrichment_categorical_MMIL_probabilities.csv"
  )
  
  write.csv(patient_scores, patient_scores_file, row.names = FALSE)
  write.csv(metrics_results, metrics_results_file, row.names = FALSE)
  write.csv(best_test_metrics, best_test_metrics_file, row.names = FALSE)
  write.csv(confusion_results, confusion_results_file, row.names = FALSE)
  write.csv(celltype_enrichment, celltype_enrichment_file, row.names = FALSE)
  
  plot_files <- list()
  
  if (save_plots) {
    plot_files$metrics <- plot_categorical_metrics(
      metrics_results = metrics_results,
      out_dir = out_dir
    )
  }
  
  cat("\nCategorical patient-level evaluation done.\n")
  
  return(
    list(
      task_type = "categorical",
      class_levels = class_levels,
      cell_pred = cell_pred,
      prob_long = prob_long,
      patient_scores = patient_scores,
      metrics_results = metrics_results,
      best_test_metrics = best_test_metrics,
      confusion_results = confusion_results,
      celltype_enrichment = celltype_enrichment,
      files = list(
        patient_scores_file = patient_scores_file,
        metrics_results_file = metrics_results_file,
        best_test_metrics_file = best_test_metrics_file,
        confusion_results_file = confusion_results_file,
        celltype_enrichment_file = celltype_enrichment_file,
        plot_files = plot_files
      )
    )
  )
}

############################################################
# Unified evaluation dispatcher
############################################################

evaluate_predictions <- function(
  pred_file,
  out_dir,
  task_type = "binary",
  class_levels = NULL,
  save_plots = TRUE
) {
  task_type <- normalize_task_type(task_type)
  
  if (task_type == "binary") {
    return(
      evaluate_binary_predictions(
        pred_file = pred_file,
        out_dir = out_dir,
        save_plots = save_plots
      )
    )
  }
  
  if (task_type == "categorical") {
    if (is.null(class_levels)) {
      stop("class_levels must be provided for categorical evaluation.")
    }
    
    return(
      evaluate_categorical_predictions(
        pred_file = pred_file,
        out_dir = out_dir,
        class_levels = class_levels,
        save_plots = save_plots
      )
    )
  }
  
  stop("Unsupported task_type.")
}