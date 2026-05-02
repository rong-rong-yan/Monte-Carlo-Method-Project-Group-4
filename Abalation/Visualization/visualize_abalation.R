############################################################
# visualize_ablation.R
# Standalone visualization of MMIL / MCEM-MMIL ablation results.
#
# Independent of the main pipeline. Reads the combined
# experiment metrics CSV and produces a focused set of
# figures plus supporting CSVs.
#
# Framing:
#   The ablation grid (3 feature methods x 3 dims x 4 backends
#   = 36 configs per task) is a ROBUSTNESS test for the claim
#   that MCEM-MMIL adds value over Naive and Deterministic EM.
#   The visualizations are organized around that claim, not
#   around finding a single best configuration.
#
# Usage:
#   Rscript visualize_ablation.R \
#     "<.../Output/summary/all_experiment_metrics.csv>" \
#     [output_dir]
#
#   First arg : path to all_experiment_metrics.csv
#   Second arg: output directory (created if missing).
#               If omitted, defaults to "<csv_dir>/../viz".
#
# Outputs (per task -- "binary" and "categorical"):
#
#   01_main_comparison_<task>.png
#       Distribution of primary metric across all 36 configs,
#       split by model. The headline picture: do MCEM/EM beat
#       Naive across the grid?
#
#   02_feature_method_<task>.png
#       Best metric per (feature_method, n_components),
#       holding model = MCEM-MMIL fixed. Sanity check that
#       MCEM's gain is not driven by one feature method.
#
#   03_model_backend_<task>.png
#       Best metric per backend, faceted by feature_method,
#       MCEM only. Sanity check across backends.
#
#   04_dimensionality_<task>.png
#       Metric vs n_components, mean +/- range across backends,
#       MCEM only. Sanity check that MCEM doesn't depend on a
#       specific dimensionality.
#
#   06_win_counts_<task>.png + .csv
#       Across the 36 configs, how often does each model rank
#       1st / 2nd / 3rd? Direct evidence of robustness.
#
#   07_paired_differences_<task>.png + .csv
#       Per-config paired differences: MCEM - Naive and
#       MCEM - EM. Histograms with zero line and annotated
#       summary stats. "Is MCEM's improvement reliably positive?"
#
#   99_appendix_top_configs_<task>.csv
#       Top 10 (model, aggregation, config) rows by primary
#       metric. Kept as appendix material; not the main story.
############################################################

############################################################
# Dependencies
############################################################

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "readr", "scales")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0(
        "Package '", pkg, "' is required. ",
        "Install with: install.packages(\"", pkg, "\")"
      )
    )
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(scales)
})

############################################################
# Constants
############################################################

MODEL_ORDER <- c(
  "Naive inherited-label classifier",
  "Deterministic EM-MMIL",
  "MCEM-MMIL"
)

# Map long names from the pipeline to short labels for plotting.
# Both binary names (em_mmil_prob etc.) and categorical names
# (Categorical EM-MMIL etc.) appear in the data, so handle both.
MODEL_SHORT <- c(
  "Naive inherited-label classifier" = "Naive",
  "Deterministic EM-MMIL"            = "EM-MMIL",
  "MCEM-MMIL"                        = "MCEM-MMIL",
  "Categorical EM-MMIL"              = "EM-MMIL",
  "Categorical MCEM-MMIL"            = "MCEM-MMIL"
)

MODEL_COLORS <- c(
  "Naive"     = "#888888",
  "EM-MMIL"   = "#4C72B0",
  "MCEM-MMIL" = "#C44E52"
)

FEATURE_COLORS <- c(
  "PCA" = "#1F77B4",
  "ICA" = "#2CA02C",
  "NMF" = "#FF7F0E"
)

# Rank colors for the win-count plot.
RANK_COLORS <- c(
  "1st" = "#2CA02C",
  "2nd" = "#FFBB33",
  "3rd" = "#D62728"
)

############################################################
# Argument parsing
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(
    paste0(
      "Usage:\n",
      "  Rscript visualize_ablation.R <path/to/all_experiment_metrics.csv> [output_dir]"
    )
  )
}

input_csv <- args[1]

if (!file.exists(input_csv)) {
  stop("Input CSV not found: ", input_csv)
}

output_dir <- if (length(args) >= 2) {
  args[2]
} else {
  file.path(dirname(input_csv), "..", "viz")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Input CSV:   ", normalizePath(input_csv),  "\n", sep = "")
cat("Output dir:  ", normalizePath(output_dir, mustWork = FALSE), "\n", sep = "")

############################################################
# Load and validate
############################################################

df <- readr::read_csv(input_csv, show_col_types = FALSE)

required_cols <- c(
  "experiment_id", "task_type", "feature_method", "n_components",
  "model_backend", "glmnet_alpha", "model", "aggregation",
  "auroc", "accuracy", "balanced_accuracy", "macro_f1",
  "multiclass_logloss", "metric_value"
)

missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(
    "Input CSV missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

cat("\nLoaded ", nrow(df), " rows across ",
    length(unique(df$experiment_id)), " unique experiments.\n", sep = "")

cat("Task breakdown:\n")
print(
  df %>%
    distinct(experiment_id, task_type) %>%
    count(task_type, name = "n_experiments")
)

############################################################
# Helpers
############################################################

shorten_model <- function(model_name) {
  out <- MODEL_SHORT[as.character(model_name)]
  out[is.na(out)] <- as.character(model_name)[is.na(out)]
  out
}

# Build a "model_label" pretty string from backend + alpha.
# For glmnet, alpha tells us ridge / elasticnet / lasso.
make_model_label <- function(model_backend, glmnet_alpha) {
  ifelse(
    model_backend == "glmnet",
    dplyr::case_when(
      is.na(glmnet_alpha)              ~ "glmnet",
      abs(glmnet_alpha - 0)   < 1e-9   ~ "glmnet (ridge)",
      abs(glmnet_alpha - 0.5) < 1e-9   ~ "glmnet (elasticnet)",
      abs(glmnet_alpha - 1)   < 1e-9   ~ "glmnet (lasso)",
      TRUE ~ paste0("glmnet (alpha=", round(glmnet_alpha, 2), ")")
    ),
    as.character(model_backend)
  )
}

############################################################
# Per-task primary metric
############################################################

# Each task has a different "headline" metric:
#   binary       => auroc
#   categorical  => balanced_accuracy
prepare_task_metric <- function(df_task, task_type) {
  if (task_type == "binary") {
    df_task <- df_task %>%
      mutate(
        primary_metric = auroc,
        primary_label  = "Patient-level test AUROC"
      )
  } else {
    df_task <- df_task %>%
      mutate(
        primary_metric = balanced_accuracy,
        primary_label  = "Patient-level test balanced accuracy"
      )
  }
  
  df_task %>%
    mutate(
      model_short = factor(
        shorten_model(model),
        levels = c("Naive", "EM-MMIL", "MCEM-MMIL")
      ),
      model_label = make_model_label(model_backend, glmnet_alpha),
      n_components = as.integer(n_components)
    )
}

# One row per (experiment_id, model_short), keeping the best
# aggregation per model. This is the natural unit for comparison.
collapse_best_aggregation <- function(df_task) {
  df_task %>%
    group_by(experiment_id, model_short) %>%
    arrange(desc(primary_metric), .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
}

# Reshape to one row per CONFIG (experiment_id) with three columns
# Naive / EM-MMIL / MCEM-MMIL holding their best-aggregation
# primary metric. This is the unit for win counts and paired diffs.
build_config_wide <- function(df_task) {
  df_best <- collapse_best_aggregation(df_task)
  
  config_meta <- df_best %>%
    distinct(
      experiment_id, feature_method, n_components,
      model_backend, glmnet_alpha
    ) %>%
    mutate(model_label = make_model_label(model_backend, glmnet_alpha))
  
  scores_wide <- df_best %>%
    select(experiment_id, model_short, primary_metric) %>%
    pivot_wider(
      names_from  = model_short,
      values_from = primary_metric
    )
  
  config_meta %>%
    left_join(scores_wide, by = "experiment_id")
}

############################################################
# Plot 1 — Main comparison: Naive vs EM vs MCEM
############################################################

plot_main_comparison <- function(df_task, task_type, output_dir) {
  df_best <- collapse_best_aggregation(df_task)
  metric_label <- df_task$primary_label[1]
  
  summary_stats <- df_best %>%
    group_by(model_short) %>%
    summarise(
      median = median(primary_metric, na.rm = TRUE),
      n      = sum(!is.na(primary_metric)),
      .groups = "drop"
    )
  
  p <- ggplot(df_best, aes(x = model_short, y = primary_metric, fill = model_short)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.5) +
    geom_jitter(
      width = 0.18, height = 0, size = 1.6, alpha = 0.7,
      aes(color = model_short)
    ) +
    geom_text(
      data = summary_stats,
      aes(x = model_short, y = median,
          label = sprintf("median = %.3f", median)),
      inherit.aes = FALSE,
      vjust = -0.8, size = 3.3, fontface = "bold"
    ) +
    scale_fill_manual(values = MODEL_COLORS, guide = "none") +
    scale_color_manual(values = MODEL_COLORS, guide = "none") +
    labs(
      title    = paste0("Main comparison (", task_type, "): Naive vs MMIL methods"),
      subtitle = paste0(
        "Each point = one (feature_method, n_components, model_backend) configuration. ",
        "Best aggregation per model."
      ),
      x = NULL,
      y = metric_label
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))
  
  out_file <- file.path(output_dir, paste0("01_main_comparison_", task_type, ".png"))
  ggsave(out_file, p, width = 7, height = 5.5, dpi = 200)
  cat("  Saved: ", out_file, "\n", sep = "")
  invisible(p)
}

############################################################
# Plot 2 — Feature method x dimensionality (MCEM only)
############################################################

plot_feature_method <- function(df_task, task_type, output_dir) {
  df_best <- collapse_best_aggregation(df_task) %>%
    filter(model_short == "MCEM-MMIL")
  
  if (nrow(df_best) == 0) {
    cat("  Skipping plot 2 for ", task_type, " (no MCEM-MMIL rows).\n", sep = "")
    return(invisible(NULL))
  }
  
  metric_label <- df_task$primary_label[1]
  
  df_agg <- df_best %>%
    group_by(feature_method, n_components) %>%
    summarise(best_metric = max(primary_metric, na.rm = TRUE), .groups = "drop") %>%
    mutate(n_components = factor(n_components, levels = sort(unique(n_components))))
  
  p <- ggplot(df_agg, aes(x = n_components, y = best_metric, fill = feature_method)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(
      aes(label = sprintf("%.3f", best_metric)),
      position = position_dodge(width = 0.8),
      vjust = -0.4, size = 3
    ) +
    scale_fill_manual(values = FEATURE_COLORS) +
    labs(
      title    = paste0("Feature method comparison (", task_type, ", MCEM-MMIL)"),
      subtitle = "Best score across model backends, per (feature_method, n_components)",
      x = "Number of components",
      y = metric_label,
      fill = "Feature method"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          legend.position = "top")
  
  out_file <- file.path(output_dir, paste0("02_feature_method_", task_type, ".png"))
  ggsave(out_file, p, width = 8, height = 5.5, dpi = 200)
  cat("  Saved: ", out_file, "\n", sep = "")
  invisible(p)
}

############################################################
# Plot 3 — Model backend comparison
############################################################

plot_model_backend <- function(df_task, task_type, output_dir) {
  df_best <- collapse_best_aggregation(df_task) %>%
    filter(model_short == "MCEM-MMIL")
  
  if (nrow(df_best) == 0) {
    cat("  Skipping plot 3 for ", task_type, " (no MCEM-MMIL rows).\n", sep = "")
    return(invisible(NULL))
  }
  
  metric_label <- df_task$primary_label[1]
  
  df_agg <- df_best %>%
    group_by(feature_method, model_label) %>%
    summarise(best_metric = max(primary_metric, na.rm = TRUE), .groups = "drop")
  
  backend_order <- df_agg %>%
    group_by(model_label) %>%
    summarise(m = mean(best_metric, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(m)) %>%
    pull(model_label)
  
  df_agg$model_label <- factor(df_agg$model_label, levels = backend_order)
  
  p <- ggplot(df_agg, aes(x = model_label, y = best_metric, fill = model_label)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = sprintf("%.3f", best_metric)), vjust = -0.4, size = 3) +
    facet_wrap(~ feature_method) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(
      title    = paste0("Model backend comparison (", task_type, ", MCEM-MMIL)"),
      subtitle = "Best score per backend within each feature method (across n_components)",
      x = NULL,
      y = metric_label
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      strip.background = element_rect(fill = "gray92"),
      strip.text = element_text(face = "bold")
    )
  
  out_file <- file.path(output_dir, paste0("03_model_backend_", task_type, ".png"))
  ggsave(out_file, p, width = 9, height = 5.5, dpi = 200)
  cat("  Saved: ", out_file, "\n", sep = "")
  invisible(p)
}

############################################################
# Plot 4 — Dimensionality trend
############################################################

plot_dimensionality <- function(df_task, task_type, output_dir) {
  df_best <- collapse_best_aggregation(df_task) %>%
    filter(model_short == "MCEM-MMIL")
  
  if (nrow(df_best) == 0) {
    cat("  Skipping plot 4 for ", task_type, " (no MCEM-MMIL rows).\n", sep = "")
    return(invisible(NULL))
  }
  
  metric_label <- df_task$primary_label[1]
  
  df_agg <- df_best %>%
    group_by(feature_method, n_components) %>%
    summarise(
      mean_metric = mean(primary_metric, na.rm = TRUE),
      min_metric  = min(primary_metric, na.rm = TRUE),
      max_metric  = max(primary_metric, na.rm = TRUE),
      .groups = "drop"
    )
  
  p <- ggplot(df_agg,
             aes(x = n_components, y = mean_metric,
                 color = feature_method, fill = feature_method)) +
    geom_ribbon(aes(ymin = min_metric, ymax = max_metric),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.5) +
    scale_color_manual(values = FEATURE_COLORS) +
    scale_fill_manual(values = FEATURE_COLORS, guide = "none") +
    scale_x_continuous(breaks = sort(unique(df_agg$n_components))) +
    labs(
      title    = paste0("Dimensionality trend (", task_type, ", MCEM-MMIL)"),
      subtitle = "Line = mean across model backends; ribbon = min-max range.",
      x = "Number of components",
      y = metric_label,
      color = "Feature method"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          legend.position = "top")
  
  out_file <- file.path(output_dir, paste0("04_dimensionality_", task_type, ".png"))
  ggsave(out_file, p, width = 8, height = 5.5, dpi = 200)
  cat("  Saved: ", out_file, "\n", sep = "")
  invisible(p)
}

############################################################
# Plot 6 — Win counts (NEW: robustness-focused)
############################################################

# For each config (one experiment_id), rank Naive / EM / MCEM
# by their best-aggregation primary metric and tally how often
# each model lands in 1st, 2nd, 3rd place across all configs.
# Stacked bar chart per task makes the robustness pattern obvious.

plot_win_counts <- function(df_task, task_type, output_dir) {
  wide <- build_config_wide(df_task)
  
  needed <- c("Naive", "EM-MMIL", "MCEM-MMIL")
  if (!all(needed %in% colnames(wide))) {
    cat("  Skipping plot 6 for ", task_type,
        " (missing one of Naive/EM/MCEM columns).\n", sep = "")
    return(invisible(NULL))
  }
  
  # Per-config ranking. ties -> "min" so duplicate wins count for everyone tied.
  ranks_long <- wide %>%
    rowwise() %>%
    mutate(
      r_naive = rank(-c(Naive, `EM-MMIL`, `MCEM-MMIL`), ties.method = "min")[1],
      r_em    = rank(-c(Naive, `EM-MMIL`, `MCEM-MMIL`), ties.method = "min")[2],
      r_mcem  = rank(-c(Naive, `EM-MMIL`, `MCEM-MMIL`), ties.method = "min")[3]
    ) %>%
    ungroup() %>%
    select(experiment_id, r_naive, r_em, r_mcem) %>%
    pivot_longer(
      cols = starts_with("r_"),
      names_to = "model_key", values_to = "rank_pos"
    ) %>%
    mutate(
      model_short = case_when(
        model_key == "r_naive" ~ "Naive",
        model_key == "r_em"    ~ "EM-MMIL",
        model_key == "r_mcem"  ~ "MCEM-MMIL"
      ),
      rank_label = case_when(
        rank_pos == 1 ~ "1st",
        rank_pos == 2 ~ "2nd",
        rank_pos == 3 ~ "3rd",
        TRUE          ~ as.character(rank_pos)
      )
    )
  
  win_counts <- ranks_long %>%
    count(model_short, rank_label, name = "n_configs") %>%
    mutate(
      model_short = factor(model_short, levels = c("Naive", "EM-MMIL", "MCEM-MMIL")),
      rank_label  = factor(rank_label,  levels = c("1st", "2nd", "3rd"))
    )
  
  total_configs <- length(unique(ranks_long$experiment_id))
  
  # Save the underlying numbers as CSV for the paper.
  win_counts_wide <- win_counts %>%
    pivot_wider(
      names_from = rank_label, values_from = n_configs, values_fill = 0
    ) %>%
    mutate(
      total_configs = total_configs,
      pct_first  = round(100 * `1st` / total_configs, 1)
    )
  
  csv_file <- file.path(output_dir, paste0("06_win_counts_", task_type, ".csv"))
  readr::write_csv(win_counts_wide, csv_file)
  cat("  Saved: ", csv_file, "\n", sep = "")
  
  p <- ggplot(win_counts,
              aes(x = model_short, y = n_configs, fill = rank_label)) +
    geom_col(width = 0.6, color = "white") +
    geom_text(
      aes(label = n_configs),
      position = position_stack(vjust = 0.5),
      size = 3.6, fontface = "bold", color = "white"
    ) +
    scale_fill_manual(values = RANK_COLORS) +
    labs(
      title    = paste0("Win counts across configurations (", task_type, ")"),
      subtitle = paste0(
        "Out of ", total_configs,
        " configurations, how often does each model rank 1st / 2nd / 3rd? ",
        "Higher 1st-place share = more robust."
      ),
      x = NULL,
      y = "Number of configurations",
      fill = "Rank within config"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          legend.position = "top")
  
  out_file <- file.path(output_dir, paste0("06_win_counts_", task_type, ".png"))
  ggsave(out_file, p, width = 7, height = 5.5, dpi = 200)
  cat("  Saved: ", out_file, "\n", sep = "")
  
  invisible(list(plot = p, csv = win_counts_wide))
}

############################################################
# Plot 7 — Paired differences (NEW: robustness-focused)
############################################################

# For each config, compute MCEM - Naive and MCEM - EM. Plot
# the distribution of these differences as a histogram with a
# vertical zero line. Annotate mean / median / fraction > 0.
#
# This is the cleanest visual answer to "does MCEM reliably help?"
# If both distributions sit clearly above zero, the answer is yes.

plot_paired_differences <- function(df_task, task_type, output_dir) {
  wide <- build_config_wide(df_task)
  
  needed <- c("Naive", "EM-MMIL", "MCEM-MMIL")
  if (!all(needed %in% colnames(wide))) {
    cat("  Skipping plot 7 for ", task_type,
        " (missing one of Naive/EM/MCEM columns).\n", sep = "")
    return(invisible(NULL))
  }
  
  metric_label <- df_task$primary_label[1]
  
  diff_long <- wide %>%
    transmute(
      experiment_id,
      `MCEM - Naive`   = `MCEM-MMIL` - Naive,
      `MCEM - EM-MMIL` = `MCEM-MMIL` - `EM-MMIL`
    ) %>%
    pivot_longer(
      cols = c(`MCEM - Naive`, `MCEM - EM-MMIL`),
      names_to = "comparison",
      values_to = "delta"
    ) %>%
    filter(!is.na(delta))
  
  comparison_levels <- c("MCEM - Naive", "MCEM - EM-MMIL")
  diff_long$comparison <- factor(diff_long$comparison, levels = comparison_levels)
  
  # Per-comparison summary stats for annotation.
  diff_stats <- diff_long %>%
    group_by(comparison) %>%
    summarise(
      n_configs   = n(),
      mean_delta  = mean(delta),
      median_delta = median(delta),
      pct_positive = 100 * mean(delta > 0),
      .groups = "drop"
    ) %>%
    mutate(
      stats_label = sprintf(
        "n = %d\nmean delta = %+.4f\nmedian delta = %+.4f\nfraction > 0 = %.0f%%",
        n_configs, mean_delta, median_delta, pct_positive
      )
    )
  
  # Save numbers as CSV for the paper.
  csv_out <- diff_long %>%
    left_join(
      diff_stats %>% select(comparison, mean_delta, median_delta, pct_positive),
      by = "comparison"
    )
  csv_file <- file.path(output_dir, paste0("07_paired_differences_", task_type, ".csv"))
  readr::write_csv(csv_out, csv_file)
  cat("  Saved: ", csv_file, "\n", sep = "")
  
  # Compute reasonable y position for the annotation box.
  bin_count <- 20
  hist_max_per_panel <- diff_long %>%
    group_by(comparison) %>%
    summarise(
      hist_max = max(
        hist(delta, breaks = bin_count, plot = FALSE)$counts,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  diff_stats <- diff_stats %>%
    left_join(hist_max_per_panel, by = "comparison") %>%
    mutate(
      x_pos = max(diff_long$delta, na.rm = TRUE),
      y_pos = hist_max * 0.95
    )
  
  p <- ggplot(diff_long, aes(x = delta)) +
    geom_histogram(
      bins = bin_count,
      fill = "#4C72B0", color = "white", alpha = 0.85
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(
      data = diff_stats,
      aes(xintercept = mean_delta),
      color = "#C44E52", linewidth = 0.8
    ) +
    geom_text(
      data = diff_stats,
      aes(x = x_pos, y = y_pos, label = stats_label),
      hjust = 1, vjust = 1, size = 3.1, lineheight = 0.95,
      fontface = "plain"
    ) +
    facet_wrap(~ comparison, ncol = 1, scales = "free_y") +
    labs(
      title    = paste0("Paired differences (", task_type,
                        "): MCEM-MMIL vs alternatives"),
      subtitle = paste0(
        "Each observation = one configuration's (MCEM - alternative) score. ",
        "Black dashed = zero; red = mean."
      ),
      x = paste0("Difference in ", tolower(metric_label)),
      y = "Number of configurations"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "gray92"),
      strip.text = element_text(face = "bold")
    )
  
  out_file <- file.path(output_dir, paste0("07_paired_differences_", task_type, ".png"))
  ggsave(out_file, p, width = 8, height = 6, dpi = 200)
  cat("  Saved: ", out_file, "\n", sep = "")
  
  cat("\n  Paired-difference summary for ", task_type, ":\n", sep = "")
  print(as.data.frame(diff_stats %>% select(-x_pos, -y_pos, -hist_max, -stats_label)),
        row.names = FALSE)
  
  invisible(list(plot = p, stats = diff_stats))
}

############################################################
# Appendix CSV — top-10 configurations (kept, demoted)
############################################################

write_appendix_top_configs <- function(df_task, task_type, output_dir, top_n = 10) {
  df_best <- collapse_best_aggregation(df_task)
  
  cols_to_show <- c(
    "experiment_id", "feature_method", "n_components",
    "model_backend", "glmnet_alpha", "model_short", "aggregation",
    "primary_metric"
  )
  
  if (task_type == "binary") {
    extras <- c("auroc")
  } else {
    extras <- c("accuracy", "balanced_accuracy", "macro_f1", "multiclass_logloss")
  }
  
  cols_to_show <- unique(c(cols_to_show, intersect(extras, colnames(df_best))))
  
  out <- df_best %>%
    arrange(desc(primary_metric)) %>%
    head(top_n) %>%
    select(any_of(cols_to_show)) %>%
    mutate(rank = row_number()) %>%
    relocate(rank)
  
  out_file <- file.path(output_dir,
                        paste0("99_appendix_top_configs_", task_type, ".csv"))
  readr::write_csv(out, out_file)
  cat("  Saved: ", out_file, "\n", sep = "")
  
  invisible(out)
}

############################################################
# Driver
############################################################

run_visualization_for_task <- function(df, task_type, output_dir) {
  df_task <- df %>% filter(task_type == !!task_type)
  
  if (nrow(df_task) == 0) {
    cat("\nNo rows for task = ", task_type, ", skipping.\n", sep = "")
    return(invisible(NULL))
  }
  
  cat("\n----------------------------------------\n")
  cat("Task: ", task_type, " (", nrow(df_task), " rows)\n", sep = "")
  cat("----------------------------------------\n")
  
  df_task <- prepare_task_metric(df_task, task_type)
  
  # Sanity-check: do we have all three models in the data?
  models_present <- unique(as.character(df_task$model_short))
  cat("Models present in data: ",
      paste(models_present, collapse = ", "), "\n", sep = "")
  
  plot_main_comparison(df_task, task_type, output_dir)
  plot_feature_method(df_task, task_type, output_dir)
  plot_model_backend(df_task, task_type, output_dir)
  plot_dimensionality(df_task, task_type, output_dir)
  plot_win_counts(df_task, task_type, output_dir)
  plot_paired_differences(df_task, task_type, output_dir)
  write_appendix_top_configs(df_task, task_type, output_dir, top_n = 10)
}

for (tt in c("binary", "categorical")) {
  run_visualization_for_task(df, tt, output_dir)
}

cat("\nAll visualizations written to:\n  ",
    normalizePath(output_dir), "\n", sep = "")