############################################################
# Binary MMIL and MCEM-MMIL on GSE131907 stage pilot data
# Goal: Early vs Advanced stage-associated cell-state inference
############################################################

proj_dir <- "/Users/rongrong/Desktop/JHU/Academics/Course/大三/大三下/Monte Carlo Methods/Project/"
setwd(proj_dir)

library(dplyr)
library(ggplot2)
library(glmnet)
library(pROC)

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

# Binary observed weak label inherited from patient stage group
# Early = 0, Advanced = 1
df$z_obs <- ifelse(df$Stage_group == "Advanced", 1, 0)

cat("\nCell counts by observed label:\n")
print(table(df$Stage_group, df$z_obs))

cat("\nPatient counts by observed label:\n")
print(
  df %>%
    distinct(Patient, Stage_group, z_obs) %>%
    count(Stage_group, z_obs, name = "n_patients")
)

############################################################
# 2. Train/test split by patient
############################################################
# Important: split by patient, not by cell.
# Otherwise cells from the same patient could appear in both train and test.

patients <- df %>%
  distinct(Patient, Stage_group, z_obs)

early_patients <- patients$Patient[patients$z_obs == 0]
adv_patients   <- patients$Patient[patients$z_obs == 1]

test_early <- sample(early_patients, size = ceiling(0.25 * length(early_patients)))
test_adv   <- sample(adv_patients, size = ceiling(0.25 * length(adv_patients)))

test_patients <- c(test_early, test_adv)

df$split <- ifelse(df$Patient %in% test_patients, "test", "train")

cat("\nPatient split counts:\n")
print(
  df %>%
    distinct(Patient, Stage_group, split) %>%
    count(split, Stage_group, name = "n_patients")
)

cat("\nCell split counts:\n")
print(table(df$split, df$Stage_group))

train_df <- df %>% filter(split == "train")
test_df  <- df %>% filter(split == "test")

x_train <- as.matrix(train_df[, pc_cols])
x_test  <- as.matrix(test_df[, pc_cols])

z_train <- train_df$z_obs
z_test  <- test_df$z_obs

############################################################
# 3. Helper functions
############################################################

# Stable inverse-logit function
expit <- function(eta) {
  1 / (1 + exp(-eta))
}

# Fit lasso/ridge logistic model with optional soft labels.
# glmnet allows binomial response values between 0 and 1.
fit_weighted_logistic <- function(x, y_soft, alpha = 0, lambda = NULL) {
  
  # glmnet handles soft binomial labels more reliably when y is a
  # two-column matrix: column 1 = probability of class 0,
  # column 2 = probability of class 1.
  y_soft <- pmin(pmax(y_soft, 1e-6), 1 - 1e-6)
  y_mat <- cbind(class0 = 1 - y_soft, class1 = y_soft)
  
  if (is.null(lambda)) {
    cv_fit <- cv.glmnet(
      x = x,
      y = y_mat,
      family = "binomial",
      alpha = alpha,
      nfolds = 5,
      type.measure = "deviance"
    )
    
    fit <- glmnet(
      x = x,
      y = y_mat,
      family = "binomial",
      alpha = alpha,
      lambda = cv_fit$lambda.min
    )
    
    return(list(fit = fit, lambda = cv_fit$lambda.min, cv_fit = cv_fit))
    
  } else {
    fit <- glmnet(
      x = x,
      y = y_mat,
      family = "binomial",
      alpha = alpha,
      lambda = lambda
    )
    
    return(list(fit = fit, lambda = lambda, cv_fit = NULL))
  }
}
# Predict probability from glmnet fit
predict_prob <- function(fit_obj, x) {
  pred <- predict(fit_obj$fit, newx = x, type = "response", s = fit_obj$lambda)
  
  # Depending on glmnet output shape, extract class-1 probability.
  if (length(dim(pred)) == 3) {
    return(as.numeric(pred[, "class1", 1]))
  } else {
    return(as.numeric(pred))
  }
}

# Calculate log loss
log_loss <- function(y, p, eps = 1e-8) {
  p <- pmin(pmax(p, eps), 1 - eps)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

############################################################
# 4. Naive baseline
############################################################
# This model treats inherited patient-level labels as if they were true cell labels.
# Early cells = 0, Advanced cells = 1.

cat("\n==============================\n")
cat("Fitting naive logistic baseline\n")
cat("==============================\n")

naive_fit <- fit_weighted_logistic(
  x = x_train,
  y_soft = z_train,
  alpha = 0
)

naive_train_prob <- predict_prob(naive_fit, x_train)
naive_test_prob  <- predict_prob(naive_fit, x_test)

cat("\nNaive train log loss:\n")
print(log_loss(z_train, naive_train_prob))

cat("\nNaive test log loss:\n")
print(log_loss(z_test, naive_test_prob))

cat("\nNaive test AUC using inherited labels:\n")
print(as.numeric(auc(z_test, naive_test_prob)))

############################################################
# 5. Deterministic binary MMIL-style EM
############################################################
# Latent cell label y_i:
# y_i = 0 baseline-like
# y_i = 1 advanced-stage-associated
#
# Assumption:
# Early-stage cells are fixed as baseline-like, y_i = 0.
# Advanced-stage cells are a mixture.
#
# rho = assumed proportion of baseline-like cells among Advanced cells.
# So among Advanced cells, expected positive fraction is 1 - rho.

cat("\n==============================\n")
cat("Fitting deterministic EM-MMIL\n")
cat("==============================\n")

rho <- 0.70
max_iter <- 30
tol <- 1e-4

cat("\nAssumed rho, baseline-like fraction among Advanced cells:\n")
print(rho)

# Initialize soft latent labels:
# Early cells fixed at 0.
# Advanced cells initialized at 1 - rho.
y_soft <- ifelse(z_train == 1, 1 - rho, 0)

em_log <- data.frame(
  iter = integer(),
  mean_soft_advanced = numeric(),
  mean_prob_advanced = numeric(),
  delta = numeric()
)

em_fit <- NULL

for (iter in 1:max_iter) {
  
  old_y_soft <- y_soft
  
  # M-step:
  # Fit classifier using current soft labels.
  # Reuse lambda from naive model for stability and speed.
  em_fit <- fit_weighted_logistic(
    x = x_train,
    y_soft = y_soft,
    alpha = 0,
    lambda = naive_fit$lambda
  )
  
  p_hat <- predict_prob(em_fit, x_train)
  
  # E-step:
  # Early cells remain fixed as baseline-like.
  # Advanced cells get updated posterior-like probabilities.
  #
  # Simple calibrated mixture update:
  # p_hat is model-estimated P(y = 1 | x).
  # For Advanced cells, rescale so the average positive fraction is near 1-rho.
  #
  # This is a practical MMIL-style approximation for the class project.
  adv_idx <- which(z_train == 1)
  early_idx <- which(z_train == 0)
  
  y_soft[early_idx] <- 0
  
  p_adv <- p_hat[adv_idx]
  
  # Rescale advanced-cell probabilities so their mean equals 1-rho.
  # This enforces the assumed mixture proportion.
  target_mean <- 1 - rho
  current_mean <- mean(p_adv)
  
  if (current_mean > 0) {
    p_adv_scaled <- p_adv * target_mean / current_mean
  } else {
    p_adv_scaled <- rep(target_mean, length(p_adv))
  }
  
  p_adv_scaled <- pmin(pmax(p_adv_scaled, 1e-6), 1 - 1e-6)
  y_soft[adv_idx] <- p_adv_scaled
  
  delta <- mean(abs(y_soft - old_y_soft))
  
  em_log <- rbind(
    em_log,
    data.frame(
      iter = iter,
      mean_soft_advanced = mean(y_soft[adv_idx]),
      mean_prob_advanced = mean(p_hat[adv_idx]),
      delta = delta
    )
  )
  
  cat("EM iter", iter, "| delta =", round(delta, 6),
      "| mean advanced soft label =", round(mean(y_soft[adv_idx]), 4), "\n")
  
  if (delta < tol) {
    cat("EM converged.\n")
    break
  }
}

em_train_prob <- predict_prob(em_fit, x_train)
em_test_prob  <- predict_prob(em_fit, x_test)

cat("\nEM-MMIL train log loss against inherited labels:\n")
print(log_loss(z_train, em_train_prob))

cat("\nEM-MMIL test log loss against inherited labels:\n")
print(log_loss(z_test, em_test_prob))

cat("\nEM-MMIL test AUC using inherited labels:\n")
print(as.numeric(auc(z_test, em_test_prob)))

############################################################
# 6. Monte Carlo EM-MMIL
############################################################
# MCEM version:
# E-step produces probabilities for Advanced cells.
# Then sample latent labels from Bernoulli probabilities.
# M-step fits classifier using sampled latent labels.
#
# Early cells remain fixed at 0.

cat("\n==============================\n")
cat("Fitting MCEM-MMIL\n")
cat("==============================\n")

set.seed(2026)

mcem_iter <- 30
mcem_samples <- 20

# Initialize posterior probabilities
q_prob <- ifelse(z_train == 1, 1 - rho, 0)

mcem_log <- data.frame(
  iter = integer(),
  mean_q_advanced = numeric(),
  mean_sampled_y_advanced = numeric(),
  delta = numeric()
)

mcem_fit <- NULL

for (iter in 1:mcem_iter) {
  
  old_q_prob <- q_prob
  
  adv_idx <- which(z_train == 1)
  early_idx <- which(z_train == 0)
  
  sampled_y_mat <- matrix(0, nrow = length(z_train), ncol = mcem_samples)
  
  for (m in 1:mcem_samples) {
    sampled_y <- rep(0, length(z_train))
    sampled_y[early_idx] <- 0
    sampled_y[adv_idx] <- rbinom(length(adv_idx), size = 1, prob = q_prob[adv_idx])
    sampled_y_mat[, m] <- sampled_y
  }
  
  # Average sampled labels across Monte Carlo samples.
  # This approximates the E-step expectation.
  y_monte_carlo_avg <- rowMeans(sampled_y_mat)
  
  # M-step
  mcem_fit <- fit_weighted_logistic(
    x = x_train,
    y_soft = y_monte_carlo_avg,
    alpha = 0,
    lambda = naive_fit$lambda
  )
  
  p_hat <- predict_prob(mcem_fit, x_train)
  
  # Update q probabilities for next sampling step.
  p_adv <- p_hat[adv_idx]
  
  target_mean <- 1 - rho
  current_mean <- mean(p_adv)
  
  if (current_mean > 0) {
    p_adv_scaled <- p_adv * target_mean / current_mean
  } else {
    p_adv_scaled <- rep(target_mean, length(p_adv))
  }
  
  p_adv_scaled <- pmin(pmax(p_adv_scaled, 1e-6), 1 - 1e-6)
  
  q_prob[early_idx] <- 0
  q_prob[adv_idx] <- p_adv_scaled
  
  delta <- mean(abs(q_prob - old_q_prob))
  
  mcem_log <- rbind(
    mcem_log,
    data.frame(
      iter = iter,
      mean_q_advanced = mean(q_prob[adv_idx]),
      mean_sampled_y_advanced = mean(y_monte_carlo_avg[adv_idx]),
      delta = delta
    )
  )
  
  cat("MCEM iter", iter, "| delta =", round(delta, 6),
      "| mean q advanced =", round(mean(q_prob[adv_idx]), 4),
      "| mean sampled y advanced =", round(mean(y_monte_carlo_avg[adv_idx]), 4), "\n")
  
  if (delta < tol) {
    cat("MCEM converged.\n")
    break
  }
}

mcem_train_prob <- predict_prob(mcem_fit, x_train)
mcem_test_prob  <- predict_prob(mcem_fit, x_test)

cat("\nMCEM-MMIL train log loss against inherited labels:\n")
print(log_loss(z_train, mcem_train_prob))

cat("\nMCEM-MMIL test log loss against inherited labels:\n")
print(log_loss(z_test, mcem_test_prob))

cat("\nMCEM-MMIL test AUC using inherited labels:\n")
print(as.numeric(auc(z_test, mcem_test_prob)))

############################################################
# 7. Add predictions back to cell metadata
############################################################

train_results <- train_df %>%
  mutate(
    naive_prob = naive_train_prob,
    em_mmil_prob = em_train_prob,
    mcem_mmil_prob = mcem_train_prob,
    em_soft_label = y_soft,
    mcem_soft_label = q_prob,
    split = "train"
  )

test_results <- test_df %>%
  mutate(
    naive_prob = naive_test_prob,
    em_mmil_prob = em_test_prob,
    mcem_mmil_prob = mcem_test_prob,
    em_soft_label = NA_real_,
    mcem_soft_label = NA_real_,
    split = "test"
  )

all_results <- bind_rows(train_results, test_results)

############################################################
# 8. Summarize inferred advanced-associated probabilities
############################################################

cat("\nMean predicted probabilities by Stage_group and Cell_type:\n")
summary_by_celltype <- all_results %>%
  group_by(split, Stage_group, Cell_type) %>%
  summarise(
    n_cells = n(),
    mean_naive_prob = mean(naive_prob, na.rm = TRUE),
    mean_em_mmil_prob = mean(em_mmil_prob, na.rm = TRUE),
    mean_mcem_mmil_prob = mean(mcem_mmil_prob, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(split, Stage_group, desc(mean_mcem_mmil_prob))

print(summary_by_celltype)

write.csv(
  summary_by_celltype,
  file = file.path(proj_dir, "stage_binary_MMIL_summary_by_celltype.csv"),
  row.names = FALSE
)

############################################################
# 9. Save full results
############################################################

write.csv(
  all_results,
  file = file.path(proj_dir, "stage_binary_MMIL_MCEM_cell_predictions.csv"),
  row.names = FALSE
)

write.csv(
  em_log,
  file = file.path(proj_dir, "stage_binary_EM_MMIL_log.csv"),
  row.names = FALSE
)

write.csv(
  mcem_log,
  file = file.path(proj_dir, "stage_binary_MCEM_MMIL_log.csv"),
  row.names = FALSE
)

############################################################
# 10. Plots
############################################################

# Prediction distributions
plot_df <- all_results %>%
  select(Cell, split, Stage_group, Cell_type,
         naive_prob, em_mmil_prob, mcem_mmil_prob) %>%
  tidyr::pivot_longer(
    cols = c(naive_prob, em_mmil_prob, mcem_mmil_prob),
    names_to = "model",
    values_to = "probability"
  )

p1 <- ggplot(plot_df, aes(x = Stage_group, y = probability, fill = Stage_group)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~ model) +
  theme_bw() +
  labs(
    title = "Predicted advanced-associated probability by model",
    x = "Inherited patient-level stage group",
    y = "Predicted probability"
  )

ggsave(
  filename = file.path(proj_dir, "plot_stage_binary_model_probabilities.png"),
  plot = p1,
  width = 10,
  height = 5,
  dpi = 300
)

# Cell-type summary for MCEM
p2 <- all_results %>%
  group_by(Stage_group, Cell_type) %>%
  summarise(
    mean_mcem_prob = mean(mcem_mmil_prob, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = reorder(Cell_type, mean_mcem_prob), y = mean_mcem_prob, fill = Stage_group)) +
  geom_col(position = "dodge") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Mean MCEM-MMIL advanced-associated probability by cell type",
    x = "Cell type",
    y = "Mean MCEM-MMIL probability"
  )

ggsave(
  filename = file.path(proj_dir, "plot_stage_binary_mcem_by_celltype.png"),
  plot = p2,
  width = 9,
  height = 6,
  dpi = 300
)

# EM convergence
p3 <- ggplot(em_log, aes(x = iter, y = delta)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(
    title = "Deterministic EM-MMIL convergence",
    x = "Iteration",
    y = "Mean absolute change in soft labels"
  )

ggsave(
  filename = file.path(proj_dir, "plot_stage_binary_EM_convergence.png"),
  plot = p3,
  width = 7,
  height = 5,
  dpi = 300
)

p4 <- ggplot(mcem_log, aes(x = iter, y = delta)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(
    title = "MCEM-MMIL convergence",
    x = "Iteration",
    y = "Mean absolute change in posterior probabilities"
  )

ggsave(
  filename = file.path(proj_dir, "plot_stage_binary_MCEM_convergence.png"),
  plot = p4,
  width = 7,
  height = 5,
  dpi = 300
)

cat("\nSaved results:\n")
cat("stage_binary_MMIL_MCEM_cell_predictions.csv\n")
cat("stage_binary_MMIL_summary_by_celltype.csv\n")
cat("stage_binary_EM_MMIL_log.csv\n")
cat("stage_binary_MCEM_MMIL_log.csv\n")
cat("plot_stage_binary_model_probabilities.png\n")
cat("plot_stage_binary_mcem_by_celltype.png\n")
cat("plot_stage_binary_EM_convergence.png\n")
cat("plot_stage_binary_MCEM_convergence.png\n")

cat("\nDone.\n")