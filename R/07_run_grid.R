############################################################
# 07_run_grid.R
# Run a grid of binary and categorical MMIL / MCEM ablation experiments.
#
# Key fix in this version:
#   Resume / skip no longer depends only on the unstable EXP001 / EXP002
#   prefix. If an old completed directory has the same experiment content
#   after removing the EXP###_ prefix, it will be recognized and skipped.
#
# Example treated as the same completed experiment:
#   EXP025_categorical_three_stage_PCAk10_glmnet_ridge_seed2026
#   EXP037_categorical_three_stage_PCAk10_glmnet_ridge_seed2026
############################################################

############################################################
# Source required scripts
############################################################

source("R/06_run_one_experiment.R")
source_pipeline_scripts("R")

load_pkg("yaml")
load_pkg("dplyr")

############################################################
# Small helper
############################################################

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

############################################################
# Read grid config
############################################################

read_grid_config <- function(grid_config_file) {
  if (!file.exists(grid_config_file)) {
    stop(paste0("Grid config file not found: ", grid_config_file))
  }
  
  yaml::read_yaml(grid_config_file)
}

############################################################
# Normalize task config
############################################################

normalize_task_grid_item <- function(task_item) {
  if (is.character(task_item)) {
    task_type <- normalize_task_type(task_item)
    
    return(
      list(
        task_type = task_type,
        categorical_scheme = "three_stage"
      )
    )
  }
  
  if (is.list(task_item)) {
    task_type <- normalize_task_type(
      task_item$task_type %||% "binary"
    )
    
    categorical_scheme <- normalize_categorical_scheme(
      task_item$categorical_scheme %||% "three_stage"
    )
    
    return(
      list(
        task_type = task_type,
        categorical_scheme = categorical_scheme
      )
    )
  }
  
  stop("Each task entry must be either a string or a list.")
}

############################################################
# Build one experiment config from grid elements
############################################################

build_single_experiment_config <- function(
  grid_config,
  experiment_id,
  task_config,
  feature_method,
  n_components,
  model_config,
  seed
) {
  project <- grid_config$project %||% list()
  defaults <- grid_config$defaults %||% list()
  
  task_type <- normalize_task_type(
    task_config$task_type %||% "binary"
  )
  
  categorical_scheme <- normalize_categorical_scheme(
    task_config$categorical_scheme %||%
      defaults$categorical_scheme %||%
      "three_stage"
  )
  
  experiment <- list(
    experiment_id = experiment_id,
    
    task_type = task_type,
    categorical_scheme = categorical_scheme,
    
    feature_method = feature_method,
    n_components = n_components,
    seed = seed,
    
    test_fraction = defaults$test_fraction %||% 0.25,
    
    min_cell_fraction = defaults$min_cell_fraction %||% 0.01,
    n_hvg = defaults$n_hvg %||% 2000,
    scale_factor = defaults$scale_factor %||% 10000,
    
    nmf_method = defaults$nmf_method %||% "brunet",
    nmf_nrun = defaults$nmf_nrun %||% 1,
    
    rho = defaults$rho %||% 0.70,
    label_strength = defaults$label_strength %||% 0.70,
    
    max_iter = defaults$max_iter %||% 30,
    tol = defaults$tol %||% 1e-4,
    mcem_samples = defaults$mcem_samples %||% 30,
    
    save_plots = defaults$save_plots %||% TRUE
  )
  
  backend <- model_config
  
  single_config <- list(
    project = project,
    experiment = experiment,
    backend = backend
  )
  
  return(single_config)
}

############################################################
# Expand grid config
############################################################

expand_experiment_grid <- function(grid_config) {
  raw_tasks <- grid_config$tasks %||% list("binary")
  
  task_configs <- lapply(
    raw_tasks,
    normalize_task_grid_item
  )
  
  features <- grid_config$features %||% list()
  
  methods <- features$methods %||% c("PCA", "ICA", "NMF")
  n_components_values <- features$n_components %||% c(20)
  
  models <- grid_config$models %||% list(
    list(
      model_name = "glmnet_ridge",
      backend = "glmnet",
      glmnet_alpha = 0
    )
  )
  
  seeds <- grid_config$seeds %||% c(2026)
  
  experiment_configs <- list()
  experiment_counter <- 1
  
  for (task_config in task_configs) {
    for (feature_method in methods) {
      for (n_components in n_components_values) {
        for (model_config in models) {
          for (seed in seeds) {
            
            task_type <- task_config$task_type
            categorical_scheme <- task_config$categorical_scheme
            model_name <- model_config$model_name %||% model_config$backend
            
            task_name <- if (task_type == "categorical") {
              paste0(task_type, "_", categorical_scheme)
            } else {
              task_type
            }
            
            experiment_id <- paste0(
              "EXP",
              sprintf("%03d", experiment_counter),
              "_",
              task_name,
              "_",
              feature_method,
              "k",
              n_components,
              "_",
              model_name,
              "_seed",
              seed
            )
            
            single_config <- build_single_experiment_config(
              grid_config = grid_config,
              experiment_id = experiment_id,
              task_config = task_config,
              feature_method = feature_method,
              n_components = n_components,
              model_config = model_config,
              seed = seed
            )
            
            experiment_configs[[experiment_counter]] <- single_config
            experiment_counter <- experiment_counter + 1
          }
        }
      }
    }
  }
  
  return(experiment_configs)
}

############################################################
# Save grid expansion table
############################################################

make_grid_table <- function(experiment_configs) {
  rows <- lapply(
    seq_along(experiment_configs),
    function(i) {
      cfg <- experiment_configs[[i]]
      
      data.frame(
        experiment_id = cfg$experiment$experiment_id,
        
        task_type = cfg$experiment$task_type,
        categorical_scheme = cfg$experiment$categorical_scheme,
        
        feature_method = cfg$experiment$feature_method,
        n_components = cfg$experiment$n_components,
        
        model_name = cfg$backend$model_name %||% NA,
        backend = cfg$backend$backend,
        glmnet_alpha = cfg$backend$glmnet_alpha %||% NA,
        
        seed = cfg$experiment$seed,
        
        rho = cfg$experiment$rho,
        label_strength = cfg$experiment$label_strength,
        mcem_samples = cfg$experiment$mcem_samples,
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  dplyr::bind_rows(rows)
}

############################################################
# Resume helpers
############################################################

get_grid_out_dir <- function(grid_config) {
  project <- grid_config$project %||% list()
  proj_dir <- project$proj_dir %||% getwd()
  out_dir <- project$out_dir %||% file.path(proj_dir, "Output")
  return(out_dir)
}

get_experiment_out_dir_for_cfg <- function(grid_config, cfg) {
  out_dir <- get_grid_out_dir(grid_config)
  task_type <- normalize_task_type(cfg$experiment$task_type)
  experiment_id <- cfg$experiment$experiment_id
  
  file.path(
    out_dir,
    task_type,
    "experiments",
    experiment_id
  )
}

is_experiment_complete <- function(experiment_out_dir) {
  summary_file <- file.path(experiment_out_dir, "experiment_summary.csv")
  
  if (!file.exists(summary_file)) {
    return(FALSE)
  }
  
  info <- file.info(summary_file)
  
  if (is.na(info$size) || info$size == 0) {
    return(FALSE)
  }
  
  ok <- tryCatch(
    {
      df <- read.csv(
        summary_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      nrow(df) > 0
    },
    error = function(e) FALSE,
    warning = function(w) FALSE
  )
  
  return(ok)
}

load_existing_summary_safe <- function(experiment_out_dir, cfg = NULL) {
  summary_file <- file.path(experiment_out_dir, "experiment_summary.csv")
  
  out <- tryCatch(
    read.csv(
      summary_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  
  # If summary is loaded from an old EXP-numbered directory,
  # rewrite experiment_id to the current grid ID so combined summaries
  # match the current expanded grid table.
  if (!is.null(out) && !is.null(cfg) && "experiment_id" %in% colnames(out)) {
    out$experiment_id <- cfg$experiment$experiment_id
  }
  
  return(out)
}

strip_experiment_counter_prefix <- function(experiment_id) {
  sub("^EXP[0-9]+_", "", experiment_id)
}

find_existing_experiment_dir <- function(grid_config, cfg, force_rerun = FALSE) {
  expected_dir <- get_experiment_out_dir_for_cfg(grid_config, cfg)
  
  if (force_rerun) {
    return(
      list(
        complete = FALSE,
        experiment_dir = expected_dir,
        expected_dir = expected_dir,
        match_type = "force_rerun"
      )
    )
  }
  
  # 1. First try the exact current directory.
  if (is_experiment_complete(expected_dir)) {
    return(
      list(
        complete = TRUE,
        experiment_dir = expected_dir,
        expected_dir = expected_dir,
        match_type = "exact"
      )
    )
  }
  
  # 2. Then try old directories with different EXP### prefix.
  out_dir <- get_grid_out_dir(grid_config)
  task_type <- normalize_task_type(cfg$experiment$task_type)
  experiments_root <- file.path(out_dir, task_type, "experiments")
  
  if (!dir.exists(experiments_root)) {
    return(
      list(
        complete = FALSE,
        experiment_dir = expected_dir,
        expected_dir = expected_dir,
        match_type = "none"
      )
    )
  }
  
  target_key <- strip_experiment_counter_prefix(cfg$experiment$experiment_id)
  candidate_dirs <- list.dirs(
    experiments_root,
    recursive = FALSE,
    full.names = TRUE
  )
  candidate_keys <- strip_experiment_counter_prefix(basename(candidate_dirs))
  matched_dirs <- candidate_dirs[candidate_keys == target_key]
  
  if (length(matched_dirs) == 0) {
    return(
      list(
        complete = FALSE,
        experiment_dir = expected_dir,
        expected_dir = expected_dir,
        match_type = "none"
      )
    )
  }
  
  completed_dirs <- matched_dirs[
    vapply(matched_dirs, is_experiment_complete, logical(1))
  ]
  
  if (length(completed_dirs) == 0) {
    return(
      list(
        complete = FALSE,
        experiment_dir = expected_dir,
        expected_dir = expected_dir,
        match_type = "matched_but_incomplete"
      )
    )
  }
  
  # If multiple completed dirs match, use the newest one.
  mtimes <- file.info(
    file.path(completed_dirs, "experiment_summary.csv")
  )$mtime
  best_dir <- completed_dirs[which.max(mtimes)]
  
  return(
    list(
      complete = TRUE,
      experiment_dir = best_dir,
      expected_dir = expected_dir,
      match_type = "suffix"
    )
  )
}

############################################################
# Feature signature & grouping
############################################################

make_feature_group_signature <- function(cfg) {
  e <- cfg$experiment
  
  task_part <- if (identical(e$task_type, "categorical")) {
    paste0(e$task_type, ":", e$categorical_scheme)
  } else {
    e$task_type
  }
  
  paste(
    task_part,
    e$seed,
    e$test_fraction,
    e$feature_method,
    e$n_components,
    e$min_cell_fraction,
    e$n_hvg,
    e$scale_factor,
    e$nmf_method,
    e$nmf_nrun,
    sep = "|"
  )
}

group_experiments_by_features <- function(experiment_configs) {
  signatures <- vapply(
    experiment_configs,
    make_feature_group_signature,
    character(1)
  )
  
  unique_sigs <- unique(signatures)
  
  groups <- lapply(
    seq_along(unique_sigs),
    function(g_idx) {
      sig <- unique_sigs[g_idx]
      idx <- which(signatures == sig)
      
      list(
        group_index = g_idx,
        signature   = sig,
        indices     = idx,
        experiments = experiment_configs[idx]
      )
    }
  )
  
  groups
}

############################################################
# Save split summaries / failures / skipped logs
############################################################

save_task_split_summaries <- function(all_summary_df, out_dir) {
  global_summary_dir <- file.path(out_dir, "summary")
  binary_summary_dir <- file.path(out_dir, "binary", "summary")
  categorical_summary_dir <- file.path(out_dir, "categorical", "summary")
  
  ensure_dir(global_summary_dir)
  ensure_dir(binary_summary_dir)
  ensure_dir(categorical_summary_dir)
  
  all_summary_file <- file.path(global_summary_dir, "all_experiment_metrics.csv")
  write.csv(all_summary_df, file = all_summary_file, row.names = FALSE)
  
  binary_summary_df <- all_summary_df %>% dplyr::filter(task_type == "binary")
  binary_summary_file <- file.path(binary_summary_dir, "all_binary_experiment_metrics.csv")
  write.csv(binary_summary_df, file = binary_summary_file, row.names = FALSE)
  
  categorical_summary_df <- all_summary_df %>% dplyr::filter(task_type == "categorical")
  categorical_summary_file <- file.path(categorical_summary_dir, "all_categorical_experiment_metrics.csv")
  write.csv(categorical_summary_df, file = categorical_summary_file, row.names = FALSE)
  
  cat("\nSaved combined experiment metrics:\n")
  cat(all_summary_file, "\n")
  cat(binary_summary_file, "\n")
  cat(categorical_summary_file, "\n")
  
  return(
    list(
      all_summary_file = all_summary_file,
      binary_summary_file = binary_summary_file,
      categorical_summary_file = categorical_summary_file
    )
  )
}

save_task_split_failures <- function(failed_experiments, out_dir) {
  if (length(failed_experiments) == 0) {
    return(NULL)
  }
  
  global_summary_dir <- file.path(out_dir, "summary")
  binary_summary_dir <- file.path(out_dir, "binary", "summary")
  categorical_summary_dir <- file.path(out_dir, "categorical", "summary")
  
  ensure_dir(global_summary_dir)
  ensure_dir(binary_summary_dir)
  ensure_dir(categorical_summary_dir)
  
  failed_df <- dplyr::bind_rows(failed_experiments)
  
  failed_file <- file.path(global_summary_dir, "failed_experiments.csv")
  write.csv(failed_df, file = failed_file, row.names = FALSE)
  
  failed_binary_file <- file.path(binary_summary_dir, "failed_binary_experiments.csv")
  write.csv(failed_df %>% dplyr::filter(task_type == "binary"), file = failed_binary_file, row.names = FALSE)
  
  failed_categorical_file <- file.path(categorical_summary_dir, "failed_categorical_experiments.csv")
  write.csv(failed_df %>% dplyr::filter(task_type == "categorical"), file = failed_categorical_file, row.names = FALSE)
  
  cat("\nSome experiments failed. Saved failure logs:\n")
  cat(failed_file, "\n")
  cat(failed_binary_file, "\n")
  cat(failed_categorical_file, "\n")
  
  return(
    list(
      failed_file = failed_file,
      failed_binary_file = failed_binary_file,
      failed_categorical_file = failed_categorical_file
    )
  )
}

save_task_split_skipped <- function(skipped_experiments, out_dir) {
  if (length(skipped_experiments) == 0) {
    return(NULL)
  }
  
  global_summary_dir <- file.path(out_dir, "summary")
  binary_summary_dir <- file.path(out_dir, "binary", "summary")
  categorical_summary_dir <- file.path(out_dir, "categorical", "summary")
  
  ensure_dir(global_summary_dir)
  ensure_dir(binary_summary_dir)
  ensure_dir(categorical_summary_dir)
  
  skipped_df <- dplyr::bind_rows(skipped_experiments)
  
  skipped_file <- file.path(global_summary_dir, "skipped_experiments.csv")
  write.csv(skipped_df, file = skipped_file, row.names = FALSE)
  
  skipped_binary_file <- file.path(binary_summary_dir, "skipped_binary_experiments.csv")
  write.csv(skipped_df %>% dplyr::filter(task_type == "binary"), file = skipped_binary_file, row.names = FALSE)
  
  skipped_categorical_file <- file.path(categorical_summary_dir, "skipped_categorical_experiments.csv")
  write.csv(skipped_df %>% dplyr::filter(task_type == "categorical"), file = skipped_categorical_file, row.names = FALSE)
  
  cat("\nSaved skipped experiment logs:\n")
  cat(skipped_file, "\n")
  cat(skipped_binary_file, "\n")
  cat(skipped_categorical_file, "\n")
  
  return(
    list(
      skipped_file = skipped_file,
      skipped_binary_file = skipped_binary_file,
      skipped_categorical_file = skipped_categorical_file
    )
  )
}

############################################################
# Pre-flight scan
############################################################

scan_experiments_for_resume <- function(
  experiment_configs,
  grid_config,
  force_rerun = FALSE
) {
  n <- length(experiment_configs)
  
  status_df <- data.frame(
    index = integer(n),
    experiment_id = character(n),
    task_type = character(n),
    categorical_scheme = character(n),
    feature_method = character(n),
    n_components = numeric(n),
    model_name = character(n),
    backend = character(n),
    seed = numeric(n),
    expected_experiment_dir = character(n),
    experiment_dir = character(n),
    resume_match_type = character(n),
    already_complete = logical(n),
    will_run = logical(n),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(experiment_configs)) {
    cfg <- experiment_configs[[i]]
    found <- find_existing_experiment_dir(
      grid_config = grid_config,
      cfg = cfg,
      force_rerun = force_rerun
    )
    
    status_df$index[i] <- i
    status_df$experiment_id[i] <- cfg$experiment$experiment_id
    status_df$task_type[i] <- cfg$experiment$task_type
    status_df$categorical_scheme[i] <- cfg$experiment$categorical_scheme
    status_df$feature_method[i] <- cfg$experiment$feature_method
    status_df$n_components[i] <- cfg$experiment$n_components
    status_df$model_name[i] <- cfg$backend$model_name %||% NA
    status_df$backend[i] <- cfg$backend$backend
    status_df$seed[i] <- cfg$experiment$seed
    status_df$expected_experiment_dir[i] <- found$expected_dir
    status_df$experiment_dir[i] <- found$experiment_dir
    status_df$resume_match_type[i] <- found$match_type
    status_df$already_complete[i] <- found$complete
    status_df$will_run[i] <- !found$complete
  }
  
  status_df
}

############################################################
# Helpers used by run loop
############################################################

record_skipped_experiment <- function(cfg, expected_dir, loaded_summary_ok, match_type = NA) {
  data.frame(
    experiment_id      = cfg$experiment$experiment_id,
    task_type          = cfg$experiment$task_type,
    categorical_scheme = cfg$experiment$categorical_scheme,
    feature_method     = cfg$experiment$feature_method,
    n_components       = cfg$experiment$n_components,
    model_name         = cfg$backend$model_name %||% NA,
    backend            = cfg$backend$backend,
    seed               = cfg$experiment$seed,
    experiment_dir     = expected_dir,
    resume_match_type  = match_type,
    loaded_summary     = loaded_summary_ok,
    stringsAsFactors   = FALSE
  )
}

record_failed_experiment <- function(cfg, error_message) {
  data.frame(
    experiment_id      = cfg$experiment$experiment_id,
    task_type          = cfg$experiment$task_type,
    categorical_scheme = cfg$experiment$categorical_scheme,
    feature_method     = cfg$experiment$feature_method,
    n_components       = cfg$experiment$n_components,
    model_name         = cfg$backend$model_name %||% NA,
    backend            = cfg$backend$backend,
    seed               = cfg$experiment$seed,
    error_message      = error_message,
    stringsAsFactors   = FALSE
  )
}

############################################################
# Main grid runner
############################################################

run_experiment_grid <- function(grid_config_file, force_rerun = NULL) {
  grid_config <- read_grid_config(grid_config_file)
  
  if (is.null(force_rerun)) {
    force_rerun <- isTRUE(grid_config$force_rerun)
  } else {
    force_rerun <- isTRUE(force_rerun)
  }
  
  out_dir <- get_grid_out_dir(grid_config)
  
  global_summary_dir <- file.path(out_dir, "summary")
  binary_summary_dir <- file.path(out_dir, "binary", "summary")
  categorical_summary_dir <- file.path(out_dir, "categorical", "summary")
  
  ensure_dir(global_summary_dir)
  ensure_dir(binary_summary_dir)
  ensure_dir(categorical_summary_dir)
  
  experiment_configs <- expand_experiment_grid(grid_config)
  grid_table <- make_grid_table(experiment_configs)
  
  ############################################################
  # Save expanded grid tables
  ############################################################
  
  grid_table_file <- file.path(global_summary_dir, "expanded_experiment_grid.csv")
  write.csv(grid_table, file = grid_table_file, row.names = FALSE)
  
  binary_grid_table_file <- file.path(binary_summary_dir, "expanded_binary_experiment_grid.csv")
  categorical_grid_table_file <- file.path(categorical_summary_dir, "expanded_categorical_experiment_grid.csv")
  
  write.csv(grid_table %>% dplyr::filter(task_type == "binary"), file = binary_grid_table_file, row.names = FALSE)
  write.csv(grid_table %>% dplyr::filter(task_type == "categorical"), file = categorical_grid_table_file, row.names = FALSE)
  
  cat("\nExpanded experiment grid:\n")
  print(grid_table)
  
  cat("\nSaved expanded grid tables:\n")
  cat(grid_table_file, "\n")
  cat(binary_grid_table_file, "\n")
  cat(categorical_grid_table_file, "\n")
  
  ############################################################
  # Pre-flight resume scan
  ############################################################
  
  resume_status <- scan_experiments_for_resume(
    experiment_configs = experiment_configs,
    grid_config = grid_config,
    force_rerun = force_rerun
  )
  
  resume_status_file <- file.path(global_summary_dir, "resume_status_preflight.csv")
  write.csv(resume_status, file = resume_status_file, row.names = FALSE)
  
  n_total <- nrow(resume_status)
  n_skip  <- sum(resume_status$already_complete)
  n_run   <- sum(resume_status$will_run)
  
  cat("\n############################################################\n")
  cat("Resume / skip pre-flight\n")
  cat("############################################################\n")
  cat("force_rerun:        ", force_rerun, "\n", sep = "")
  cat("Total experiments:  ", n_total, "\n", sep = "")
  cat("Already complete:   ", n_skip, " (will be skipped)\n", sep = "")
  cat("Pending:            ", n_run, " (will be run)\n", sep = "")
  cat("Pre-flight log:     ", resume_status_file, "\n", sep = "")
  cat("############################################################\n")
  
  ############################################################
  # Group experiments by feature signature
  ############################################################
  
  groups <- group_experiments_by_features(experiment_configs)
  
  cat("\n[grouping] ", length(experiment_configs), " experiments grouped into ", length(groups), " feature group(s).\n", sep = "")
  cat("[grouping] Feature extraction will run ", length(groups), "x instead of ", length(experiment_configs), "x.\n", sep = "")
  
  all_summaries <- list()
  failed_experiments <- list()
  skipped_experiments <- list()
  
  ############################################################
  # Loop over groups
  ############################################################
  
  for (g_idx in seq_along(groups)) {
    group <- groups[[g_idx]]
    n_in_group <- length(group$experiments)
    group_status <- resume_status[group$indices, , drop = FALSE]
    
    n_pending_in_group <- sum(group_status$will_run)
    n_skip_in_group <- sum(group_status$already_complete)
    
    cat("\n############################################################\n")
    cat("Group ", g_idx, " / ", length(groups), " | signature: ", group$signature, "\n", sep = "")
    cat("  Experiments in group: ", n_in_group, " (skip: ", n_skip_in_group, ", run: ", n_pending_in_group, ")\n", sep = "")
    cat("############################################################\n")
    
    ##########################################################
    # Case A: all complete
    ##########################################################
    
    if (n_pending_in_group == 0) {
      cat("[group skip] All experiments in group already complete.\n")
      
      for (j in seq_len(n_in_group)) {
        cfg_j <- group$experiments[[j]]
        existing_dir_j <- group_status$experiment_dir[j]
        match_type_j <- group_status$resume_match_type[j]
        
        existing <- load_existing_summary_safe(existing_dir_j, cfg_j)
        
        if (!is.null(existing)) {
          all_summaries[[length(all_summaries) + 1]] <- existing
        }
        
        skipped_experiments[[length(skipped_experiments) + 1]] <-
          record_skipped_experiment(
            cfg = cfg_j,
            expected_dir = existing_dir_j,
            loaded_summary_ok = !is.null(existing),
            match_type = match_type_j
          )
      }
      
      next
    }
    
    ##########################################################
    # Case B: at least one pending, compute shared data/features
    ##########################################################
    
    cat("\n[group prepare] Computing shared data + features...\n")
    
    shared <- tryCatch(
      {
        prepare_shared_data_features_and_dfs(
          config = group$experiments[[1]]
        )
      },
      error = function(e) {
        cat("\n[group failed] Shared data/feature preparation failed:\n  ", e$message, "\n", sep = "")
        return(NULL)
      }
    )
    
    if (is.null(shared)) {
      for (j in seq_len(n_in_group)) {
        cfg_j <- group$experiments[[j]]
        existing_dir_j <- group_status$experiment_dir[j]
        match_type_j <- group_status$resume_match_type[j]
        
        if (group_status$already_complete[j]) {
          existing <- load_existing_summary_safe(existing_dir_j, cfg_j)
          if (!is.null(existing)) {
            all_summaries[[length(all_summaries) + 1]] <- existing
          }
          skipped_experiments[[length(skipped_experiments) + 1]] <-
            record_skipped_experiment(
              cfg = cfg_j,
              expected_dir = existing_dir_j,
              loaded_summary_ok = !is.null(existing),
              match_type = match_type_j
            )
        } else {
          failed_experiments[[length(failed_experiments) + 1]] <-
            record_failed_experiment(
              cfg = cfg_j,
              error_message = "Group feature preparation failed; experiment not attempted."
            )
        }
      }
      next
    }
    
    ##########################################################
    # Run / skip experiments in group
    ##########################################################
    
    for (j in seq_len(n_in_group)) {
      cfg_j <- group$experiments[[j]]
      existing_dir_j <- group_status$experiment_dir[j]
      match_type_j <- group_status$resume_match_type[j]
      experiment_id <- cfg_j$experiment$experiment_id
      
      if (group_status$already_complete[j]) {
        cat("\n  [skip] ", experiment_id, " — already complete (", match_type_j, ").\n", sep = "")
        
        existing <- load_existing_summary_safe(existing_dir_j, cfg_j)
        
        if (!is.null(existing)) {
          all_summaries[[length(all_summaries) + 1]] <- existing
        }
        
        skipped_experiments[[length(skipped_experiments) + 1]] <-
          record_skipped_experiment(
            cfg = cfg_j,
            expected_dir = existing_dir_j,
            loaded_summary_ok = !is.null(existing),
            match_type = match_type_j
          )
        
        next
      }
      
      cat("\n  [run] ", experiment_id, " (group ", g_idx, " / ", length(groups), ", item ", j, " / ", n_in_group, ")\n", sep = "")
      cat("        backend: ", cfg_j$backend$backend,
          if (cfg_j$backend$backend == "glmnet") paste0(", alpha = ", cfg_j$backend$glmnet_alpha %||% NA) else "",
          "\n", sep = "")
      
      result_j <- tryCatch(
        {
          run_one_experiment(
            config                = cfg_j,
            r_dir                 = "R",
            source_scripts        = FALSE,
            precomputed_data      = shared$data,
            precomputed_features  = shared$features,
            precomputed_model_dfs = shared$model_dfs
          )
        },
        error = function(e) {
          cat("\nExperiment failed:\n")
          cat(experiment_id, "\n")
          cat("Error message:\n")
          print(e$message)
          
          failed_experiments[[length(failed_experiments) + 1]] <<-
            record_failed_experiment(
              cfg = cfg_j,
              error_message = e$message
            )
          
          return(NULL)
        }
      )
      
      if (!is.null(result_j)) {
        all_summaries[[length(all_summaries) + 1]] <- result_j$summary
      }
    }
    
    rm(shared)
    invisible(gc(verbose = FALSE))
  }
  
  ############################################################
  # Save combined summaries
  ############################################################
  
  summary_files <- NULL
  
  if (length(all_summaries) > 0) {
    all_summary_df <- dplyr::bind_rows(all_summaries)
    
    summary_files <- save_task_split_summaries(
      all_summary_df = all_summary_df,
      out_dir = out_dir
    )
    
    cat("\nTop results by metric_value:\n")
    print(
      all_summary_df %>%
        dplyr::arrange(dplyr::desc(metric_value)) %>%
        head(20)
    )
    
    if ("auroc" %in% colnames(all_summary_df)) {
      cat("\nTop binary results by AUROC:\n")
      print(
        all_summary_df %>%
          dplyr::filter(task_type == "binary") %>%
          dplyr::arrange(dplyr::desc(auroc)) %>%
          head(20)
      )
    }
    
    if ("accuracy" %in% colnames(all_summary_df)) {
      cat("\nTop categorical results by accuracy:\n")
      print(
        all_summary_df %>%
          dplyr::filter(task_type == "categorical") %>%
          dplyr::arrange(dplyr::desc(accuracy)) %>%
          head(20)
      )
    }
    
  } else {
    cat("\nNo successful experiments to summarize.\n")
  }
  
  failure_files <- save_task_split_failures(
    failed_experiments = failed_experiments,
    out_dir = out_dir
  )
  
  skipped_files <- save_task_split_skipped(
    skipped_experiments = skipped_experiments,
    out_dir = out_dir
  )
  
  ############################################################
  # Final summary banner
  ############################################################
  
  n_run_attempted <- n_run
  n_run_failed    <- length(failed_experiments)
  n_run_succeeded <- n_run_attempted - n_run_failed
  n_skipped       <- length(skipped_experiments)
  n_groups        <- length(groups)
  
  cat("\n############################################################\n")
  cat("Grid run completed.\n")
  cat("############################################################\n")
  cat("Total experiments in grid:     ", n_total,         "\n", sep = "")
  cat("Feature groups:                ", n_groups,        "\n", sep = "")
  cat("Skipped already complete:      ", n_skipped,       "\n", sep = "")
  cat("Attempted this run:            ", n_run_attempted, "\n", sep = "")
  cat("  succeeded:                   ", n_run_succeeded, "\n", sep = "")
  cat("  failed:                      ", n_run_failed,    "\n", sep = "")
  cat("############################################################\n")
  
  return(
    list(
      grid_table         = grid_table,
      resume_status      = resume_status,
      groups             = groups,
      summaries          = all_summaries,
      failed             = failed_experiments,
      skipped            = skipped_experiments,
      summary_files      = summary_files,
      failure_files      = failure_files,
      skipped_files      = skipped_files,
      counts = list(
        total     = n_total,
        groups    = n_groups,
        skipped   = n_skipped,
        attempted = n_run_attempted,
        succeeded = n_run_succeeded,
        failed    = n_run_failed
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
        "Please provide a grid config file.\n",
        "Example:\n",
        "  Rscript R/07_run_grid.R config/experiment_grid.yml\n",
        "  Rscript R/07_run_grid.R config/experiment_grid.yml --force-rerun"
      )
    )
  }
  
  grid_config_file <- args[1]
  force_rerun_cli <- "--force-rerun" %in% args
  
  invisible(
    run_experiment_grid(
      grid_config_file = grid_config_file,
      force_rerun = if (force_rerun_cli) TRUE else NULL
    )
  )
}
