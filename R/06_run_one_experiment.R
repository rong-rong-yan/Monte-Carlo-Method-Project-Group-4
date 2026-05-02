############################################################
# 06_run_one_experiment.R
# Run one MMIL / MCEM ablation experiment.
#
# Supports:
#   1. binary
#   2. categorical
#
# Output structure:
#   Output/
#   ├── binary/
#   │   └── experiments/
#   └── categorical/
#       └── experiments/
#
# Usage from command line:
#   Rscript R/06_run_one_experiment.R config/single_experiment.yml
#
# Resource-sharing extension:
#   run_one_experiment() now optionally accepts precomputed_data,
#   precomputed_features, and precomputed_model_dfs. When provided,
#   the corresponding pipeline stage is SKIPPED and the precomputed
#   object is reused. This lets the grid runner share data loading
#   and feature extraction across experiments that have identical
#   data + feature configurations (e.g. same task / feature_method /
#   n_components / seed but different model backends).
############################################################

############################################################
# Source helper
############################################################

source_pipeline_scripts <- function(r_dir = "R") {
  scripts <- c(
    "00_utils.R",
    "01_load_and_split.R",
    "02_feature_extractors.R",
    "03_model_backends.R",
    "04_mmil_wrappers.R",
    "05_patient_level_evaluation.R"
  )
  
  for (script in scripts) {
    script_path <- file.path(r_dir, script)
    
    if (!file.exists(script_path)) {
      stop(paste0("Cannot find required script: ", script_path))
    }
    
    source(script_path)
  }
  
  invisible(TRUE)
}

############################################################
# Small helper functions
############################################################

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

is_absolute_path <- function(path) {
  grepl("^([A-Za-z]:[\\/]|[\\/])", path)
}

resolve_path <- function(path, base_dir) {
  if (is.null(path)) {
    return(NULL)
  }
  
  if (is_absolute_path(path)) {
    return(path)
  }
  
  file.path(base_dir, path)
}

############################################################
# Read and normalize one experiment config
############################################################

read_experiment_config <- function(config_file) {
  load_pkg("yaml")
  
  if (!file.exists(config_file)) {
    stop(paste0("Config file not found: ", config_file))
  }
  
  yaml::read_yaml(config_file)
}

normalize_experiment_config <- function(raw_config) {
  project <- raw_config$project %||% list()
  experiment <- raw_config$experiment %||% list()
  backend <- raw_config$backend %||% list()
  
  ############################################################
  # Project paths
  ############################################################
  
  proj_dir <- project$proj_dir %||% getwd()
  in_dir <- project$in_dir %||% file.path(proj_dir, "Input")
  out_dir <- project$out_dir %||% file.path(proj_dir, "Output")
  
  expr_file <- project$expr_file %||% "GSE131907_pilot_expression_subset_dt.rds"
  meta_file <- project$meta_file %||% "GSE131907_pilot_sampled_cell_metadata_with_clinical.csv"
  
  expr_file <- resolve_path(expr_file, in_dir)
  meta_file <- resolve_path(meta_file, in_dir)
  
  ############################################################
  # Task settings
  ############################################################
  
  task_type <- normalize_task_type(
    experiment$task_type %||% "binary"
  )
  
  categorical_scheme <- normalize_categorical_scheme(
    experiment$categorical_scheme %||% "three_stage"
  )
  
  class_levels <- get_categorical_class_levels(
    categorical_scheme = categorical_scheme
  )
  
  ############################################################
  # General experiment settings
  ############################################################
  
  seed <- experiment$seed %||% 2026
  test_fraction <- experiment$test_fraction %||% 0.25
  
  feature_method <- experiment$feature_method %||% "PCA"
  n_components <- experiment$n_components %||% 20
  
  min_cell_fraction <- experiment$min_cell_fraction %||% 0.01
  n_hvg <- experiment$n_hvg %||% 2000
  scale_factor <- experiment$scale_factor %||% 10000
  
  nmf_method <- experiment$nmf_method %||% "brunet"
  nmf_nrun <- experiment$nmf_nrun %||% 1
  
  ############################################################
  # MMIL settings
  ############################################################
  
  rho <- experiment$rho %||% 0.70
  label_strength <- experiment$label_strength %||% 0.70
  
  max_iter <- experiment$max_iter %||% 30
  tol <- experiment$tol %||% 1e-4
  mcem_samples <- experiment$mcem_samples %||% 30
  
  save_plots <- experiment$save_plots %||% TRUE
  
  ############################################################
  # Backend settings
  ############################################################
  
  model_backend <- backend$backend %||% "glmnet"
  glmnet_alpha <- backend$glmnet_alpha %||% NA_real_
  
  ############################################################
  # Experiment ID
  ############################################################
  
  tmp_config <- list(
    task_type = task_type,
    categorical_scheme = categorical_scheme,
    feature_method = feature_method,
    n_components = n_components,
    model_backend = model_backend,
    glmnet_alpha = glmnet_alpha,
    seed = seed
  )
  
  experiment_id <- experiment$experiment_id %||% make_experiment_name(tmp_config)
  
  ############################################################
  # Output directories
  ############################################################
  
  task_out_dir <- file.path(
    out_dir,
    task_type
  )
  
  experiment_out_dir <- file.path(
    task_out_dir,
    "experiments",
    experiment_id
  )
  
  ensure_dir(experiment_out_dir)
  
  ############################################################
  # Backend config
  ############################################################
  
  backend_config <- backend
  backend_config$backend <- model_backend
  
  if (model_backend == "glmnet" && is.null(backend_config$glmnet_alpha)) {
    backend_config$glmnet_alpha <- 0
  }
  
  ############################################################
  # Normalized config object
  ############################################################
  
  normalized <- list(
    experiment_id = experiment_id,
    
    task_type = task_type,
    categorical_scheme = categorical_scheme,
    class_levels = class_levels,
    
    proj_dir = proj_dir,
    in_dir = in_dir,
    out_dir = out_dir,
    task_out_dir = task_out_dir,
    experiment_out_dir = experiment_out_dir,
    
    expr_file = expr_file,
    meta_file = meta_file,
    
    seed = seed,
    test_fraction = test_fraction,
    
    feature_method = feature_method,
    n_components = n_components,
    min_cell_fraction = min_cell_fraction,
    n_hvg = n_hvg,
    scale_factor = scale_factor,
    nmf_method = nmf_method,
    nmf_nrun = nmf_nrun,
    
    rho = rho,
    label_strength = label_strength,
    max_iter = max_iter,
    tol = tol,
    mcem_samples = mcem_samples,
    
    save_plots = save_plots,
    
    model_backend = model_backend,
    glmnet_alpha = glmnet_alpha,
    backend_config = backend_config
  )
  
  return(normalized)
}

############################################################
# Save config used
############################################################

save_config_used <- function(raw_config, normalized_config) {
  load_pkg("yaml")
  
  out_dir <- normalized_config$experiment_out_dir
  
  raw_file <- file.path(out_dir, "config_used_raw.yml")
  normalized_file <- file.path(out_dir, "config_used_normalized.yml")
  
  yaml::write_yaml(raw_config, raw_file)
  yaml::write_yaml(normalized_config, normalized_file)
  
  invisible(
    list(
      raw_file = raw_file,
      normalized_file = normalized_file
    )
  )
}

############################################################
# Validation helpers for precomputed objects
############################################################

validate_precomputed_features <- function(cfg, precomputed_features) {
  if (is.null(precomputed_features)) {
    return(invisible(TRUE))
  }
  
  if (!is.null(precomputed_features$feature_method)) {
    if (precomputed_features$feature_method != cfg$feature_method) {
      stop(
        sprintf(
          paste0(
            "precomputed_features$feature_method ('%s') does not ",
            "match cfg$feature_method ('%s')."
          ),
          precomputed_features$feature_method,
          cfg$feature_method
        )
      )
    }
  }
  
  if (!is.null(precomputed_features$n_components)) {
    if (precomputed_features$n_components != cfg$n_components) {
      stop(
        sprintf(
          paste0(
            "precomputed_features$n_components (%d) does not ",
            "match cfg$n_components (%d)."
          ),
          precomputed_features$n_components,
          cfg$n_components
        )
      )
    }
  }
  
  invisible(TRUE)
}

validate_precomputed_data <- function(cfg, precomputed_data) {
  if (is.null(precomputed_data)) {
    return(invisible(TRUE))
  }
  
  if (!is.null(precomputed_data$task_type)) {
    if (precomputed_data$task_type != cfg$task_type) {
      stop(
        sprintf(
          "precomputed_data$task_type ('%s') does not match cfg$task_type ('%s').",
          precomputed_data$task_type,
          cfg$task_type
        )
      )
    }
  }
  
  if (!is.null(precomputed_data$categorical_scheme) &&
      cfg$task_type == "categorical") {
    if (precomputed_data$categorical_scheme != cfg$categorical_scheme) {
      stop(
        sprintf(
          paste0(
            "precomputed_data$categorical_scheme ('%s') does not ",
            "match cfg$categorical_scheme ('%s')."
          ),
          precomputed_data$categorical_scheme,
          cfg$categorical_scheme
        )
      )
    }
  }
  
  invisible(TRUE)
}

############################################################
# Build experiment summary table
############################################################

build_experiment_summary <- function(config, eval_result) {
  ############################################################
  # Binary summary
  ############################################################
  
  if (eval_result$task_type == "binary") {
    summary_df <- eval_result$best_test_auc
    
    summary_df$metric_name <- "auroc"
    summary_df$metric_value <- summary_df$auroc
    
    summary_df$experiment_id <- config$experiment_id
    summary_df$task_type <- config$task_type
    summary_df$categorical_scheme <- NA_character_
    summary_df$feature_method <- config$feature_method
    summary_df$n_components <- config$n_components
    summary_df$model_backend <- config$model_backend
    summary_df$glmnet_alpha <- config$glmnet_alpha
    summary_df$seed <- config$seed
    summary_df$rho <- config$rho
    summary_df$label_strength <- NA_real_
    summary_df$mcem_samples <- config$mcem_samples
    
    summary_df$multiclass_logloss <- NA_real_
    summary_df$accuracy <- NA_real_
    summary_df$balanced_accuracy <- NA_real_
    summary_df$macro_f1 <- NA_real_
    
    summary_df <- summary_df[, c(
      "experiment_id",
      "task_type",
      "categorical_scheme",
      "feature_method",
      "n_components",
      "model_backend",
      "glmnet_alpha",
      "seed",
      "rho",
      "label_strength",
      "mcem_samples",
      "model",
      "aggregation",
      "n_patients",
      "n_early",
      "n_advanced",
      "auroc",
      "multiclass_logloss",
      "accuracy",
      "balanced_accuracy",
      "macro_f1",
      "metric_name",
      "metric_value"
    )]
    
    return(summary_df)
  }
  
  ############################################################
  # Categorical summary
  ############################################################
  
  if (eval_result$task_type == "categorical") {
    summary_df <- eval_result$best_test_metrics
    
    summary_df$metric_name <- "accuracy"
    summary_df$metric_value <- summary_df$accuracy
    
    summary_df$experiment_id <- config$experiment_id
    summary_df$task_type <- config$task_type
    summary_df$categorical_scheme <- config$categorical_scheme
    summary_df$feature_method <- config$feature_method
    summary_df$n_components <- config$n_components
    summary_df$model_backend <- config$model_backend
    summary_df$glmnet_alpha <- config$glmnet_alpha
    summary_df$seed <- config$seed
    summary_df$rho <- NA_real_
    summary_df$label_strength <- config$label_strength
    summary_df$mcem_samples <- config$mcem_samples
    
    summary_df$n_early <- NA_integer_
    summary_df$n_advanced <- NA_integer_
    summary_df$auroc <- NA_real_
    
    summary_df <- summary_df[, c(
      "experiment_id",
      "task_type",
      "categorical_scheme",
      "feature_method",
      "n_components",
      "model_backend",
      "glmnet_alpha",
      "seed",
      "rho",
      "label_strength",
      "mcem_samples",
      "model",
      "aggregation",
      "n_patients",
      "n_early",
      "n_advanced",
      "auroc",
      "multiclass_logloss",
      "accuracy",
      "balanced_accuracy",
      "macro_f1",
      "metric_name",
      "metric_value"
    )]
    
    return(summary_df)
  }
  
  stop("Unsupported eval_result$task_type.")
}

############################################################
# Helper: prepare shared (data + features + model_dfs) once
############################################################

prepare_shared_data_features_and_dfs <- function(config) {
  ############################################################
  # Used by 07_run_grid.R to build the shared objects ONCE
  # for a group of experiments that have identical data /
  # feature configurations.
  #
  # Accepts either a raw YAML config or a normalized config.
  ############################################################
  
  if (is.character(config)) {
    raw_config <- read_experiment_config(config)
    cfg <- normalize_experiment_config(raw_config)
  } else if (!is.null(config$experiment_id) && !is.null(config$task_out_dir)) {
    # Already normalized.
    cfg <- config
  } else {
    raw_config <- config
    cfg <- normalize_experiment_config(raw_config)
  }
  
  cat("\n[group prepare] Loading data and extracting features...\n")
  cat("  Task type:        ", cfg$task_type, "\n", sep = "")
  
  if (cfg$task_type == "categorical") {
    cat("  Categorical:      ", cfg$categorical_scheme, "\n", sep = "")
  }
  
  cat("  Feature method:   ", cfg$feature_method, "\n", sep = "")
  cat("  N components:     ", cfg$n_components, "\n", sep = "")
  cat("  Seed:             ", cfg$seed, "\n", sep = "")
  
  data_obj <- load_and_split_data(cfg)
  
  feature_obj <- extract_features(
    train_expr        = data_obj$train_expr,
    test_expr         = data_obj$test_expr,
    feature_method    = cfg$feature_method,
    n_components      = cfg$n_components,
    min_cell_fraction = cfg$min_cell_fraction,
    n_hvg             = cfg$n_hvg,
    scale_factor      = cfg$scale_factor,
    seed              = cfg$seed,
    nmf_method        = cfg$nmf_method,
    nmf_nrun          = cfg$nmf_nrun
  )
  
  model_dfs <- build_feature_modeling_dfs(
    train_meta     = data_obj$train_meta,
    test_meta      = data_obj$test_meta,
    feature_result = feature_obj
  )
  
  list(
    data      = data_obj,
    features  = feature_obj,
    model_dfs = model_dfs
  )
}

############################################################
# Main function: run one experiment
############################################################

run_one_experiment <- function(
  config,
  r_dir = "R",
  source_scripts = FALSE,
  precomputed_data = NULL,
  precomputed_features = NULL,
  precomputed_model_dfs = NULL
) {
  if (source_scripts) {
    source_pipeline_scripts(r_dir)
  }
  
  load_pkg("yaml")
  
  if (is.character(config)) {
    raw_config <- read_experiment_config(config)
  } else {
    raw_config <- config
  }
  
  cfg <- normalize_experiment_config(raw_config)
  
  cat("\n============================================================\n")
  cat("Running experiment:\n")
  cat(cfg$experiment_id, "\n")
  cat("============================================================\n")
  
  cat("\nExperiment settings:\n")
  cat("Task type:", cfg$task_type, "\n")
  
  if (cfg$task_type == "categorical") {
    cat("Categorical scheme:", cfg$categorical_scheme, "\n")
    cat("Class levels:", paste(cfg$class_levels, collapse = ", "), "\n")
  }
  
  cat("Feature method:", cfg$feature_method, "\n")
  cat("Number of components:", cfg$n_components, "\n")
  cat("Model backend:", cfg$model_backend, "\n")
  
  if (cfg$model_backend == "glmnet") {
    cat("glmnet alpha:", cfg$backend_config$glmnet_alpha, "\n")
  }
  
  cat("Seed:", cfg$seed, "\n")
  cat("Task output directory:", cfg$task_out_dir, "\n")
  cat("Experiment output directory:", cfg$experiment_out_dir, "\n")
  
  save_config_used(
    raw_config = raw_config,
    normalized_config = cfg
  )
  
  ############################################################
  # Validate precomputed args (if any) against the config.
  ############################################################
  
  validate_precomputed_data(cfg, precomputed_data)
  validate_precomputed_features(cfg, precomputed_features)
  
  ############################################################
  # 1. Load data and split by patient (or reuse precomputed)
  ############################################################
  
  if (is.null(precomputed_data)) {
    data_obj <- load_and_split_data(cfg)
  } else {
    cat("\n[reuse] Using precomputed data object (skipping load_and_split_data).\n")
    data_obj <- precomputed_data
  }
  
  ############################################################
  # 2. Feature extraction (or reuse precomputed)
  ############################################################
  
  if (is.null(precomputed_features)) {
    feature_obj <- extract_features(
      train_expr        = data_obj$train_expr,
      test_expr         = data_obj$test_expr,
      feature_method    = cfg$feature_method,
      n_components      = cfg$n_components,
      min_cell_fraction = cfg$min_cell_fraction,
      n_hvg             = cfg$n_hvg,
      scale_factor      = cfg$scale_factor,
      seed              = cfg$seed,
      nmf_method        = cfg$nmf_method,
      nmf_nrun          = cfg$nmf_nrun
    )
  } else {
    cat("\n[reuse] Using precomputed feature object (skipping extract_features).\n")
    feature_obj <- precomputed_features
  }
  
  ############################################################
  # 3. Build modeling data frames (or reuse precomputed)
  ############################################################
  
  if (is.null(precomputed_model_dfs)) {
    model_dfs <- build_feature_modeling_dfs(
      train_meta     = data_obj$train_meta,
      test_meta      = data_obj$test_meta,
      feature_result = feature_obj
    )
  } else {
    cat("\n[reuse] Using precomputed model_dfs (skipping build_feature_modeling_dfs).\n")
    model_dfs <- precomputed_model_dfs
  }
  
  ############################################################
  # 4. Run MMIL according to task type
  ############################################################
  
  if (cfg$task_type == "binary") {
    mmil_result <- run_binary_mmil(
      train_df       = model_dfs$train_df,
      test_df        = model_dfs$test_df,
      feature_cols   = model_dfs$feature_cols,
      backend_config = cfg$backend_config,
      rho            = cfg$rho,
      max_iter       = cfg$max_iter,
      tol            = cfg$tol,
      mcem_samples   = cfg$mcem_samples,
      seed           = cfg$seed
    )
  }
  
  if (cfg$task_type == "categorical") {
    mmil_result <- run_categorical_mmil(
      train_df        = model_dfs$train_df,
      test_df         = model_dfs$test_df,
      feature_cols    = model_dfs$feature_cols,
      backend_config  = cfg$backend_config,
      class_levels    = cfg$class_levels,
      label_strength  = cfg$label_strength,
      max_iter        = cfg$max_iter,
      tol             = cfg$tol,
      mcem_samples    = cfg$mcem_samples,
      seed            = cfg$seed
    )
  }
  
  ############################################################
  # 5. Save cell-level predictions and logs
  ############################################################
  
  saved <- save_mmil_outputs(
    mmil_result   = mmil_result,
    out_dir       = cfg$experiment_out_dir,
    pred_filename = "cell_predictions.csv"
  )
  
  ############################################################
  # 6. Patient-level evaluation
  ############################################################
  
  eval_result <- evaluate_predictions(
    pred_file    = saved$pred_file,
    out_dir      = cfg$experiment_out_dir,
    task_type    = cfg$task_type,
    class_levels = cfg$class_levels,
    save_plots   = cfg$save_plots
  )
  
  ############################################################
  # 7. Save experiment-level summary
  ############################################################
  
  summary_df <- build_experiment_summary(
    config      = cfg,
    eval_result = eval_result
  )
  
  summary_file <- file.path(
    cfg$experiment_out_dir,
    "experiment_summary.csv"
  )
  
  write.csv(
    summary_df,
    file = summary_file,
    row.names = FALSE
  )
  
  cat("\nSaved experiment summary:\n")
  cat(summary_file, "\n")
  
  cat("\nExperiment completed:\n")
  cat(cfg$experiment_id, "\n")
  
  return(
    list(
      config       = cfg,
      data         = data_obj,
      features     = feature_obj,
      model_dfs    = model_dfs,
      mmil_result  = mmil_result,
      eval_result  = eval_result,
      summary      = summary_df,
      files = list(
        prediction_file    = saved$pred_file,
        summary_file       = summary_file,
        task_out_dir       = cfg$task_out_dir,
        experiment_out_dir = cfg$experiment_out_dir
      )
    )
  )
}

############################################################
# Command-line entry point
############################################################

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    stop(
      paste0(
        "Please provide a config file.\n",
        "Example:\n",
        "  Rscript R/06_run_one_experiment.R config/single_experiment.yml"
      )
    )
  }
  
  config_file <- args[1]
  
  run_one_experiment(
    config = config_file,
    r_dir = "R",
    source_scripts = TRUE
  )
}