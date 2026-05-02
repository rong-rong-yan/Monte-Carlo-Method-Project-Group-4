############################################################
# Categorical MCEM-MMIL
# Dataset: GSE131907 stage pilot data
#
# Classes:
#   Stage I
#   Stage II/III
#   Stage IV
#
# This script adds Monte Carlo sampling to the categorical
# deterministic EM-MMIL baseline.
############################################################

proj_dir <- "/Users/rongrong/Desktop/JHU/Academics/Course/大三/大三下/Monte Carlo Methods/Project/"
setwd(proj_dir)

library(dplyr)
library(ggplot2)
library(glmnet)
library(tidyr)

set.seed(2026)

############################################################
# 1. Load modeling-ready data
############################################################

model_file <- file.path(proj_dir, "GSE131907_stage_pilot_PCs_metadata.csv")

df <- read.csv(
  model_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Loaded modeling data:\n")
print(dim(df))

pc_cols <- grep("^PC", colnames(df), value = TRUE)

cat("\nPC columns:\n")
print(pc_cols)

############################################################
# 2. Create categorical stage label:
#    I, II/III, IV
############################################################

df <- df %>%
  mutate(
    Stage_cat = case_when(
      Stage_broad == "I" ~ "I",
      Stage_broad %in% c("II", "III") ~ "II_III",
      Stage_broad == "IV" ~ "IV",
      TRUE ~ NA_character_
    )
  )

df$Stage_cat <- factor(df$Stage_cat, levels = c("I", "II_III", "IV"))

cat("\nCell counts by Stage_cat:\n")
print(table(df$Stage_cat, useNA = "ifany"))

cat("\nPatient counts by Stage_cat:\n")
print(
  df %>%
    distinct(Patient, Stage_cat) %>%
    count(Stage_cat, name = "n_patients")
)

df <- df %>% filter(!is.na(Stage_cat))

############################################################
# 3. Train/test split by patient
############################################################
# Use the same fixed seed and same logic as the categorical EM script.
# This keeps the comparison fair.

patients <- df %>%
  distinct(Patient, Stage_cat)

test_patients <- patients %>%
  group_by(Stage_cat) %>%
  summarise(
    test_patient = list(sample(Patient, size = ceiling(0.25 * n()))),
    .groups = "drop"
  ) %>%
  pull(test_patient) %>%
  unlist()

df$split <- ifelse(df$Patient %in% test_patients, "test", "train")

cat("\nPatient split counts:\n")
print(
  df %>%
    distinct(Patient, Stage_cat, split) %>%
    count(split, Stage_cat, name = "n_patients")
)

cat("\nCell split counts:\n")
print(table(df$split, df$Stage_cat))

train_df <- df %>% filter(split == "train")
test_df  <- df %>% filter(split == "test")

x_train <- as.matrix(train_df[, pc_cols])
x_test  <- as.matrix(test_df[, pc_cols])

y_train_cat <- train_df$Stage_cat
y_test_cat  <- test_df$Stage_cat

class_levels <- levels(df$Stage_cat)
K <- length(class_levels)

cat("\nClass levels:\n")
print(class_levels)

############################################################
# 4. Helper functions
############################################################

row_normalize <- function(mat, eps = 1e-12) {
  mat <- pmax(mat, eps)
  mat / rowSums(mat)
}

make_onehot <- function(y, class_levels) {
  y <- factor(y, levels = class_levels)
  out <- matrix(0, nrow = length(y), ncol = length(class_levels))
  colnames(out) <- class_levels
  
  for (k in seq_along(class_levels)) {
    out[, k] <- as.numeric(y == class_levels[k])
  }
  
  return(out)
}

multiclass_log_loss <- function(y_true, prob_mat, class_levels, eps = 1e-8) {
  y_true <- factor(y_true, levels = class_levels)
  prob_mat <- pmin(pmax(prob_mat, eps), 1 - eps)
  idx <- cbind(seq_along(y_true), as.integer(y_true))
  -mean(log(prob_mat[idx]))
}

multiclass_accuracy <- function(y_true, prob_mat, class_levels) {
  y_pred <- class_levels[max.col(prob_mat, ties.method = "first")]
  mean(y_pred == as.character(y_true))
}

fit_multinom_hard <- function(x, y, alpha = 0, lambda = NULL) {
  
  y <- factor(y, levels = class_levels)
  
  if (is.null(lambda)) {
    cv_fit <- cv.glmnet(
      x = x,
      y = y,
      family = "multinomial",
      alpha = alpha,
      nfolds = 5,
      type.measure = "deviance",
      standardize = TRUE,
      maxit = 1000000
    )
    
    fit <- glmnet(
      x = x,
      y = y,
      family = "multinomial",
      alpha = alpha,
      lambda = cv_fit$lambda.1se,
      standardize = TRUE,
      maxit = 1000000
    )
    
    return(list(fit = fit, lambda = cv_fit$lambda.1se, cv_fit = cv_fit))
    
  } else {
    fit <- glmnet(
      x = x,
      y = y,
      family = "multinomial",
      alpha = alpha,
      lambda = lambda,
      standardize = TRUE,
      maxit = 1000000
    )
    
    return(list(fit = fit, lambda = lambda, cv_fit = NULL))
  }
}

predict_multinom_prob <- function(fit_obj, x, class_levels) {
  pred <- predict(
    fit_obj$fit,
    newx = x,
    type = "response",
    s = fit_obj$lambda
  )
  
  if (length(dim(pred)) == 3) {
    prob_mat <- pred[, , 1]
  } else {
    prob_mat <- pred
  }
  
  prob_mat <- as.matrix(prob_mat)
  prob_mat <- prob_mat[, class_levels, drop = FALSE]
  prob_mat <- row_normalize(prob_mat)
  
  return(prob_mat)
}

fit_multinom_soft <- function(x, q_mat, alpha = 0, lambda = NULL) {
  
  q_mat <- row_normalize(q_mat)
  
  n <- nrow(x)
  K <- ncol(q_mat)
  
  x_expanded <- x[rep(seq_len(n), times = K), , drop = FALSE]
  y_expanded <- rep(colnames(q_mat), each = n)
  w_expanded <- as.vector(q_mat)
  
  y_expanded <- factor(y_expanded, levels = class_levels)
  
  keep <- w_expanded > 1e-6
  
  x_expanded <- x_expanded[keep, , drop = FALSE]
  y_expanded <- y_expanded[keep]
  w_expanded <- w_expanded[keep]
  
  # Normalize weights for numerical stability.
  w_expanded <- w_expanded / mean(w_expanded)
  
  if (is.null(lambda)) {
    lambda_use <- NULL
  } else {
    lambda_use <- max(lambda, 0.01)
  }
  
  if (is.null(lambda_use)) {
    cv_fit <- cv.glmnet(
      x = x_expanded,
      y = y_expanded,
      weights = w_expanded,
      family = "multinomial",
      alpha = alpha,
      nfolds = 5,
      type.measure = "deviance",
      standardize = TRUE,
      maxit = 1000000
    )
    
    fit <- glmnet(
      x = x_expanded,
      y = y_expanded,
      weights = w_expanded,
      family = "multinomial",
      alpha = alpha,
      lambda = cv_fit$lambda.1se,
      standardize = TRUE,
      maxit = 1000000
    )
    
    return(list(fit = fit, lambda = cv_fit$lambda.1se, cv_fit = cv_fit))
    
  } else {
    fit <- glmnet(
      x = x_expanded,
      y = y_expanded,
      weights = w_expanded,
      family = "multinomial",
      alpha = alpha,
      lambda = lambda_use,
      standardize = TRUE,
      maxit = 1000000
    )
    
    return(list(fit = fit, lambda = lambda_use, cv_fit = NULL))
  }
}

############################################################
# 5. Categorical sampling helper
############################################################
# For each cell, sample one class from its probability vector.
# Repeating this M times gives Monte Carlo samples of latent labels.

sample_categorical_onehot <- function(q_mat, class_levels) {
  
  n <- nrow(q_mat)
  K <- ncol(q_mat)
  
  sampled_mat <- matrix(0, nrow = n, ncol = K)
  colnames(sampled_mat) <- class_levels
  
  for (i in seq_len(n)) {
    sampled_class <- sample(
      x = seq_len(K),
      size = 1,
      prob = q_mat[i, ]
    )
    sampled_mat[i, sampled_class] <- 1
  }
  
  return(sampled_mat)
}

############################################################
# 6. Naive multinomial baseline
############################################################

cat("\n==============================\n")
cat("Fitting naive multinomial baseline\n")
cat("==============================\n")

naive_fit <- fit_multinom_hard(
  x = x_train,
  y = y_train_cat,
  alpha = 0
)

naive_train_prob <- predict_multinom_prob(naive_fit, x_train, class_levels)
naive_test_prob  <- predict_multinom_prob(naive_fit, x_test, class_levels)

cat("\nNaive train log loss:\n")
print(multiclass_log_loss(y_train_cat, naive_train_prob, class_levels))

cat("\nNaive test log loss:\n")
print(multiclass_log_loss(y_test_cat, naive_test_prob, class_levels))

cat("\nNaive train accuracy:\n")
print(multiclass_accuracy(y_train_cat, naive_train_prob, class_levels))

cat("\nNaive test accuracy:\n")
print(multiclass_accuracy(y_test_cat, naive_test_prob, class_levels))

cat("\nNaive test confusion matrix:\n")
naive_test_pred <- class_levels[max.col(naive_test_prob, ties.method = "first")]
print(table(True = y_test_cat, Predicted = naive_test_pred))

############################################################
# 7. Categorical MCEM-MMIL
############################################################
# Latent cell label:
#   Y_i in {I, II_III, IV}
#
# MCEM procedure:
#   1. q_ik is current posterior probability of latent class k for cell i.
#   2. Monte Carlo E-step: sample latent one-hot labels from q_i.
#   3. Average sampled labels over M Monte Carlo samples.
#   4. M-step: fit multinomial classifier using averaged sampled labels.
#   5. Update q using classifier probabilities and weak patient-stage prior.
############################################################

cat("\n==============================\n")
cat("Fitting categorical MCEM-MMIL\n")
cat("==============================\n")

max_iter <- 40
tol <- 1e-4
mcem_samples <- 30

label_strength <- 0.70

cat("\nLabel strength:\n")
print(label_strength)

cat("\nMonte Carlo samples per iteration:\n")
print(mcem_samples)

# Initial q matrix:
# observed class gets label_strength,
# other classes share remaining probability.
q_mat <- matrix(
  (1 - label_strength) / (K - 1),
  nrow = nrow(train_df),
  ncol = K
)
colnames(q_mat) <- class_levels

obs_idx <- match(as.character(y_train_cat), class_levels)

for (i in seq_len(nrow(q_mat))) {
  q_mat[i, obs_idx[i]] <- label_strength
}

q_mat <- row_normalize(q_mat)

cat("\nInitial average q by observed class:\n")
print(
  aggregate(
    q_mat,
    by = list(observed_stage = y_train_cat),
    FUN = mean
  )
)

# Patient-stage prior used in every E-step.
label_prior <- matrix(
  (1 - label_strength) / (K - 1),
  nrow = nrow(train_df),
  ncol = K
)
colnames(label_prior) <- class_levels

for (i in seq_len(nrow(label_prior))) {
  label_prior[i, obs_idx[i]] <- label_strength
}

label_prior <- row_normalize(label_prior)

mcem_log <- data.frame(
  iter = integer(),
  delta = numeric(),
  q_entropy = numeric(),
  mc_entropy = numeric(),
  train_log_loss_inherited = numeric(),
  stringsAsFactors = FALSE
)

mcem_fit <- NULL

for (iter in 1:max_iter) {
  
  old_q <- q_mat
  
  ##########################################################
  # Monte Carlo E-step:
  # sample latent categorical labels from current q.
  ##########################################################
  
  sampled_sum <- matrix(
    0,
    nrow = nrow(q_mat),
    ncol = ncol(q_mat)
  )
  colnames(sampled_sum) <- class_levels
  
  for (m in seq_len(mcem_samples)) {
    sampled_sum <- sampled_sum + sample_categorical_onehot(q_mat, class_levels)
  }
  
  q_mc_avg <- sampled_sum / mcem_samples
  q_mc_avg <- row_normalize(q_mc_avg)
  
  ##########################################################
  # M-step:
  # fit multinomial classifier using Monte Carlo averaged labels.
  ##########################################################
  
  mcem_fit <- fit_multinom_soft(
    x = x_train,
    q_mat = q_mc_avg,
    alpha = 0,
    lambda = max(naive_fit$lambda, 0.01)
  )
  
  p_hat <- predict_multinom_prob(mcem_fit, x_train, class_levels)
  
  ##########################################################
  # Deterministic posterior update for next iteration:
  # combine classifier probabilities with weak patient-stage prior.
  ##########################################################
  
  q_new <- p_hat * label_prior
  q_new <- row_normalize(q_new)
  
  q_mat <- q_new
  
  delta <- mean(abs(q_mat - old_q))
  
  q_entropy <- -mean(rowSums(q_mat * log(pmax(q_mat, 1e-12))))
  mc_entropy <- -mean(rowSums(q_mc_avg * log(pmax(q_mc_avg, 1e-12))))
  
  inherited_log_loss <- multiclass_log_loss(
    y_true = y_train_cat,
    prob_mat = p_hat,
    class_levels = class_levels
  )
  
  mcem_log <- rbind(
    mcem_log,
    data.frame(
      iter = iter,
      delta = delta,
      q_entropy = q_entropy,
      mc_entropy = mc_entropy,
      train_log_loss_inherited = inherited_log_loss
    )
  )
  
  cat(
    "Categorical MCEM iter", iter,
    "| delta =", round(delta, 6),
    "| q entropy =", round(q_entropy, 4),
    "| MC entropy =", round(mc_entropy, 4),
    "| inherited train log loss =", round(inherited_log_loss, 4),
    "\n"
  )
  
  if (delta < tol) {
    cat("Categorical MCEM converged.\n")
    break
  }
}

############################################################
# 8. Evaluate categorical MCEM-MMIL
############################################################

mcem_train_prob <- predict_multinom_prob(mcem_fit, x_train, class_levels)
mcem_test_prob  <- predict_multinom_prob(mcem_fit, x_test, class_levels)

cat("\nCategorical MCEM train log loss against inherited labels:\n")
print(multiclass_log_loss(y_train_cat, mcem_train_prob, class_levels))

cat("\nCategorical MCEM test log loss against inherited labels:\n")
print(multiclass_log_loss(y_test_cat, mcem_test_prob, class_levels))

cat("\nCategorical MCEM train accuracy against inherited labels:\n")
print(multiclass_accuracy(y_train_cat, mcem_train_prob, class_levels))

cat("\nCategorical MCEM test accuracy against inherited labels:\n")
print(multiclass_accuracy(y_test_cat, mcem_test_prob, class_levels))

cat("\nCategorical MCEM test confusion matrix:\n")
mcem_test_pred <- class_levels[max.col(mcem_test_prob, ties.method = "first")]
print(table(True = y_test_cat, Predicted = mcem_test_pred))

############################################################
# 9. Add predictions back to metadata
############################################################

train_results <- train_df

for (k in class_levels) {
  train_results[[paste0("naive_prob_", k)]] <- naive_train_prob[, k]
  train_results[[paste0("cat_mcem_prob_", k)]] <- mcem_train_prob[, k]
  train_results[[paste0("cat_mcem_soft_label_", k)]] <- q_mat[, k]
}

train_results$naive_pred <- class_levels[max.col(naive_train_prob, ties.method = "first")]
train_results$cat_mcem_pred <- class_levels[max.col(mcem_train_prob, ties.method = "first")]
train_results$split <- "train"

test_results <- test_df

for (k in class_levels) {
  test_results[[paste0("naive_prob_", k)]] <- naive_test_prob[, k]
  test_results[[paste0("cat_mcem_prob_", k)]] <- mcem_test_prob[, k]
  test_results[[paste0("cat_mcem_soft_label_", k)]] <- NA_real_
}

test_results$naive_pred <- class_levels[max.col(naive_test_prob, ties.method = "first")]
test_results$cat_mcem_pred <- class_levels[max.col(mcem_test_prob, ties.method = "first")]
test_results$split <- "test"

all_results <- bind_rows(train_results, test_results)

############################################################
# 10. Summarize by observed stage and cell type
############################################################

summary_by_celltype <- all_results %>%
  group_by(split, Stage_cat, Cell_type) %>%
  summarise(
    n_cells = n(),
    mean_naive_prob_I = mean(naive_prob_I, na.rm = TRUE),
    mean_naive_prob_II_III = mean(naive_prob_II_III, na.rm = TRUE),
    mean_naive_prob_IV = mean(naive_prob_IV, na.rm = TRUE),
    mean_cat_mcem_prob_I = mean(cat_mcem_prob_I, na.rm = TRUE),
    mean_cat_mcem_prob_II_III = mean(cat_mcem_prob_II_III, na.rm = TRUE),
    mean_cat_mcem_prob_IV = mean(cat_mcem_prob_IV, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(split, Stage_cat, desc(mean_cat_mcem_prob_IV))

cat("\nSummary by cell type:\n")
print(summary_by_celltype, n = 60)

############################################################
# 11. Save results
############################################################

write.csv(
  all_results,
  file = file.path(proj_dir, "stage_categorical_MCEM_MMIL_cell_predictions.csv"),
  row.names = FALSE
)

write.csv(
  summary_by_celltype,
  file = file.path(proj_dir, "stage_categorical_MCEM_MMIL_summary_by_celltype.csv"),
  row.names = FALSE
)

write.csv(
  mcem_log,
  file = file.path(proj_dir, "stage_categorical_MCEM_MMIL_log.csv"),
  row.names = FALSE
)

############################################################
# 12. Plots
############################################################

# Convergence plot
p_conv <- ggplot(mcem_log, aes(x = iter, y = delta)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(
    title = "Categorical MCEM-MMIL convergence",
    x = "Iteration",
    y = "Mean absolute change in posterior probabilities"
  )

ggsave(
  filename = file.path(proj_dir, "plot_stage_categorical_MCEM_convergence.png"),
  plot = p_conv,
  width = 7,
  height = 5,
  dpi = 300
)

# Entropy plot
entropy_long <- mcem_log %>%
  select(iter, q_entropy, mc_entropy) %>%
  pivot_longer(
    cols = c(q_entropy, mc_entropy),
    names_to = "entropy_type",
    values_to = "entropy"
  )

p_entropy <- ggplot(entropy_long, aes(x = iter, y = entropy, color = entropy_type)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(
    title = "Categorical MCEM entropy over iterations",
    x = "Iteration",
    y = "Entropy",
    color = "Entropy type"
  )

ggsave(
  filename = file.path(proj_dir, "plot_stage_categorical_MCEM_entropy.png"),
  plot = p_entropy,
  width = 7,
  height = 5,
  dpi = 300
)

# Predicted probability distributions
prob_long <- all_results %>%
  select(
    Cell,
    split,
    Stage_cat,
    Cell_type,
    starts_with("naive_prob_"),
    starts_with("cat_mcem_prob_")
  ) %>%
  pivot_longer(
    cols = starts_with("naive_prob_") | starts_with("cat_mcem_prob_"),
    names_to = "model_class",
    values_to = "probability"
  ) %>%
  mutate(
    model = ifelse(grepl("^naive", model_class), "Naive multinomial", "Categorical MCEM-MMIL"),
    class = gsub("^naive_prob_|^cat_mcem_prob_", "", model_class)
  )

p_prob <- ggplot(
  prob_long,
  aes(x = Stage_cat, y = probability, fill = Stage_cat)
) +
  geom_boxplot(outlier.size = 0.2) +
  facet_grid(model ~ class) +
  theme_bw() +
  labs(
    title = "Predicted latent stage-associated probabilities",
    x = "Inherited patient-level stage category",
    y = "Predicted probability"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(proj_dir, "plot_stage_categorical_MCEM_probabilities.png"),
  plot = p_prob,
  width = 12,
  height = 7,
  dpi = 300
)

# Mean MCEM probabilities by cell type
celltype_long <- all_results %>%
  group_by(Stage_cat, Cell_type) %>%
  summarise(
    n_cells = n(),
    mean_prob_I = mean(cat_mcem_prob_I, na.rm = TRUE),
    mean_prob_II_III = mean(cat_mcem_prob_II_III, na.rm = TRUE),
    mean_prob_IV = mean(cat_mcem_prob_IV, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(mean_prob_I, mean_prob_II_III, mean_prob_IV),
    names_to = "latent_class",
    values_to = "mean_probability"
  ) %>%
  mutate(
    latent_class = gsub("mean_prob_", "", latent_class)
  )

p_celltype <- ggplot(
  celltype_long,
  aes(x = reorder(Cell_type, mean_probability), y = mean_probability, fill = Stage_cat)
) +
  geom_col(position = "dodge") +
  coord_flip() +
  facet_wrap(~ latent_class) +
  theme_bw() +
  labs(
    title = "Mean categorical MCEM-MMIL latent-stage probability by cell type",
    x = "Cell type",
    y = "Mean predicted probability"
  )

ggsave(
  filename = file.path(proj_dir, "plot_stage_categorical_MCEM_by_celltype.png"),
  plot = p_celltype,
  width = 12,
  height = 7,
  dpi = 300
)

############################################################
# 13. Final console summary
############################################################

cat("\nSaved results:\n")
cat("stage_categorical_MCEM_MMIL_cell_predictions.csv\n")
cat("stage_categorical_MCEM_MMIL_summary_by_celltype.csv\n")
cat("stage_categorical_MCEM_MMIL_log.csv\n")
cat("plot_stage_categorical_MCEM_convergence.png\n")
cat("plot_stage_categorical_MCEM_entropy.png\n")
cat("plot_stage_categorical_MCEM_probabilities.png\n")
cat("plot_stage_categorical_MCEM_by_celltype.png\n")

cat("\nDone.\n")