############################################################
# Patient-level aggregation evaluation for binary MMIL / MCEM
# Dataset: GSE131907 stage pilot data
#
# Goal:
# Compare naive baseline, deterministic EM-MMIL, and MCEM-MMIL
# using patient-level aggregation metrics, closer to the MMIL paper.
############################################################

proj_dir <- "/Users/rongrong/Desktop/JHU/Academics/Course/大三/大三下/Monte Carlo Methods/Project/"
setwd(proj_dir)

library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)
library(readr)

############################################################
# 1. Load cell-level prediction results
############################################################

pred_file <- file.path(proj_dir, "stage_binary_MMIL_MCEM_cell_predictions.csv")

cell_pred <- read.csv(
  pred_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Loaded cell-level prediction file:\n")
print(dim(cell_pred))

cat("\nColumn names:\n")
print(colnames(cell_pred))

cat("\nCell counts by split and Stage_group:\n")
print(table(cell_pred$split, cell_pred$Stage_group))

cat("\nPatient counts by split and Stage_group:\n")
print(
  cell_pred %>%
    distinct(Patient, split, Stage_group) %>%
    count(split, Stage_group, name = "n_patients")
)

############################################################
# 2. Check required columns
############################################################

required_cols <- c(
  "Patient",
  "split",
  "Stage_group",
  "naive_prob",
  "em_mmil_prob",
  "mcem_mmil_prob"
)

missing_cols <- setdiff(required_cols, colnames(cell_pred))

if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

############################################################
# 3. Create binary patient-level label
############################################################
# Early = 0
# Advanced = 1

cell_pred <- cell_pred %>%
  mutate(
    y_patient = ifelse(Stage_group == "Advanced", 1, 0)
  )

cat("\nCheck binary label:\n")
print(table(cell_pred$Stage_group, cell_pred$y_patient))

############################################################
# 4. Reshape probabilities into long format
############################################################

prob_long <- cell_pred %>%
  select(
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
  pivot_longer(
    cols = c(naive_prob, em_mmil_prob, mcem_mmil_prob),
    names_to = "model",
    values_to = "cell_prob"
  ) %>%
  mutate(
    model = case_when(
      model == "naive_prob" ~ "Naive inherited-label classifier",
      model == "em_mmil_prob" ~ "Deterministic EM-MMIL",
      model == "mcem_mmil_prob" ~ "MCEM-MMIL",
      TRUE ~ model
    )
  )

cat("\nProbability long format dimensions:\n")
print(dim(prob_long))

############################################################
# 5. Patient-level aggregation
############################################################
# For each patient and each model, aggregate cell probabilities
# into patient-level risk scores.

patient_scores <- prob_long %>%
  group_by(split, Patient, Stage_group, y_patient, model) %>%
  summarise(
    n_cells = n(),
    
    mean_prob = mean(cell_prob, na.rm = TRUE),
    median_prob = median(cell_prob, na.rm = TRUE),
    q90_prob = as.numeric(quantile(cell_prob, probs = 0.90, na.rm = TRUE)),
    q99_prob = as.numeric(quantile(cell_prob, probs = 0.99, na.rm = TRUE)),
    
    prop_gt_0.5 = mean(cell_prob > 0.5, na.rm = TRUE),
    prop_gt_0.9 = mean(cell_prob > 0.9, na.rm = TRUE),
    
    .groups = "drop"
  )

cat("\nPatient-level scores dimensions:\n")
print(dim(patient_scores))

cat("\nFirst few patient-level scores:\n")
print(head(patient_scores, 10))

write.csv(
  patient_scores,
  file = file.path(proj_dir, "patient_level_binary_MMIL_aggregation_scores.csv"),
  row.names = FALSE
)

############################################################
# 6. Patient-level AUROC by model and aggregation method
############################################################

aggregation_cols <- c(
  "mean_prob",
  "median_prob",
  "q90_prob",
  "q99_prob",
  "prop_gt_0.5",
  "prop_gt_0.9"
)

compute_auc_safe <- function(y, score) {
  
  # Need both classes represented.
  if (length(unique(y)) < 2) {
    return(NA_real_)
  }
  
  # Need non-constant score.
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

auc_results <- patient_scores %>%
  pivot_longer(
    cols = all_of(aggregation_cols),
    names_to = "aggregation",
    values_to = "patient_score"
  ) %>%
  group_by(split, model, aggregation) %>%
  summarise(
    n_patients = n(),
    n_early = sum(y_patient == 0),
    n_advanced = sum(y_patient == 1),
    auroc = compute_auc_safe(y_patient, patient_score),
    .groups = "drop"
  ) %>%
  arrange(split, model, desc(auroc))

cat("\nPatient-level AUROC results:\n")
print(auc_results, n = 100)

write.csv(
  auc_results,
  file = file.path(proj_dir, "patient_level_binary_MMIL_AUROC_results.csv"),
  row.names = FALSE
)

############################################################
# 7. Patient-level log loss using mean probability
############################################################
# The MMIL paper emphasizes AUROC-style aggregation, but log loss is
# useful for probability calibration. We use mean probability as the
# patient-level probability.

log_loss <- function(y, p, eps = 1e-8) {
  p <- pmin(pmax(p, eps), 1 - eps)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

patient_logloss <- patient_scores %>%
  group_by(split, model) %>%
  summarise(
    n_patients = n(),
    logloss_mean_prob = log_loss(y_patient, mean_prob),
    logloss_median_prob = log_loss(y_patient, median_prob),
    logloss_q90_prob = log_loss(y_patient, q90_prob),
    logloss_q99_prob = log_loss(y_patient, q99_prob),
    .groups = "drop"
  )

cat("\nPatient-level log loss results:\n")
print(patient_logloss)

write.csv(
  patient_logloss,
  file = file.path(proj_dir, "patient_level_binary_MMIL_logloss_results.csv"),
  row.names = FALSE
)

############################################################
# 8. Best aggregation method per model on test set
############################################################

best_test_auc <- auc_results %>%
  filter(split == "test") %>%
  group_by(model) %>%
  arrange(desc(auroc), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

cat("\nBest test AUROC aggregation per model:\n")
print(best_test_auc)

write.csv(
  best_test_auc,
  file = file.path(proj_dir, "patient_level_binary_MMIL_best_test_AUROC.csv"),
  row.names = FALSE
)

############################################################
# 9. Cell-type enrichment among high-probability cells
############################################################
# This uses Cell_type as biological interpretation, not gold-standard
# stage-associated labels.
#
# For each model, look at cells with probability > 0.5 or > 0.9.

celltype_enrichment <- prob_long %>%
  group_by(split, model, Stage_group, Cell_type) %>%
  summarise(
    n_cells = n(),
    mean_prob = mean(cell_prob, na.rm = TRUE),
    prop_gt_0.5 = mean(cell_prob > 0.5, na.rm = TRUE),
    prop_gt_0.9 = mean(cell_prob > 0.9, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(split, model, Stage_group, desc(mean_prob))

cat("\nCell-type summary for biological interpretation:\n")
print(celltype_enrichment, n = 80)

write.csv(
  celltype_enrichment,
  file = file.path(proj_dir, "celltype_enrichment_binary_MMIL_probabilities.csv"),
  row.names = FALSE
)

############################################################
# 10. Plot patient-level scores
############################################################

patient_scores_long <- patient_scores %>%
  pivot_longer(
    cols = all_of(aggregation_cols),
    names_to = "aggregation",
    values_to = "patient_score"
  )

p_scores <- ggplot(
  patient_scores_long %>% filter(split == "test"),
  aes(
    x = Stage_group,
    y = patient_score,
    fill = Stage_group
  )
) +
  geom_boxplot(outlier.size = 1) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.8) +
  facet_grid(model ~ aggregation, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Patient-level aggregated scores on test patients",
    x = "Patient-level stage group",
    y = "Aggregated predicted probability"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(proj_dir, "plot_patient_level_binary_MMIL_scores.png"),
  plot = p_scores,
  width = 14,
  height = 8,
  dpi = 300
)

############################################################
# 11. Plot AUROC comparison
############################################################

p_auc <- ggplot(
  auc_results %>% filter(split == "test"),
  aes(
    x = aggregation,
    y = auroc,
    fill = model
  )
) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(
    title = "Patient-level test AUROC by aggregation method",
    x = "Aggregation method",
    y = "Patient-level AUROC",
    fill = "Model"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(proj_dir, "plot_patient_level_binary_MMIL_AUROC.png"),
  plot = p_auc,
  width = 10,
  height = 6,
  dpi = 300
)

############################################################
# 12. Plot cell-type high-probability proportions
############################################################

p_celltype <- celltype_enrichment %>%
  filter(split == "test") %>%
  ggplot(
    aes(
      x = reorder(Cell_type, prop_gt_0.5),
      y = prop_gt_0.5,
      fill = Stage_group
    )
  ) +
  geom_col(position = "dodge") +
  coord_flip() +
  facet_wrap(~ model) +
  theme_bw() +
  labs(
    title = "Proportion of cells with predicted probability > 0.5",
    x = "Cell type",
    y = "Proportion of cells",
    fill = "Stage group"
  )

ggsave(
  filename = file.path(proj_dir, "plot_celltype_binary_MMIL_prop_gt_0.5.png"),
  plot = p_celltype,
  width = 12,
  height = 7,
  dpi = 300
)

############################################################
# 13. Final summary
############################################################

cat("\nSaved files:\n")
cat("patient_level_binary_MMIL_aggregation_scores.csv\n")
cat("patient_level_binary_MMIL_AUROC_results.csv\n")
cat("patient_level_binary_MMIL_logloss_results.csv\n")
cat("patient_level_binary_MMIL_best_test_AUROC.csv\n")
cat("celltype_enrichment_binary_MMIL_probabilities.csv\n")
cat("plot_patient_level_binary_MMIL_scores.png\n")
cat("plot_patient_level_binary_MMIL_AUROC.png\n")
cat("plot_celltype_binary_MMIL_prop_gt_0.5.png\n")

cat("\nDone.\n")