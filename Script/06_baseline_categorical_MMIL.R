############################################################
# Categorical deterministic MMIL-style baseline
# Dataset: GSE131907 stage pilot data
#
# Classes:
#   Stage I
#   Stage II/III
#   Stage IV
#
# No Monte Carlo is used in this script.
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

# Remove missing categorical stage labels, if any.
df <- df %>% filter(!is.na(Stage_cat))

############################################################
# 3. Train/test split by patient
############################################################
# Split at patient level to avoid cell leakage.

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

# Convert class factor into one-hot matrix.
make_onehot <- function(y, class_levels) {
  y <- factor(y, levels = class_levels)
  out <- matrix(0, nrow = length(y), ncol = length(class_levels))
  colnames(out) <- class_levels
  
  for (k in seq_along(class_levels)) {
    out[, k] <- as.numeric(y == class_levels[k])
  }
  
  return(out)
}

# Row-normalize a matrix.
row_normalize <- function(mat, eps = 1e-12) {
  mat <- pmax(mat, eps)
  mat / rowSums(mat)
}

# Multiclass log loss.
multiclass_log_loss <- function(y_true, prob_mat, class_levels, eps = 1e-8) {
  y_true <- factor(y_true, levels = class_levels)
  prob_mat <- pmin(pmax(prob_mat, eps), 1 - eps)
  
  idx <- cbind(seq_along(y_true), as.integer(y_true))
  -mean(log(prob_mat[idx]))
}

# Accuracy from predicted probability matrix.
multiclass_accuracy <- function(y_true, prob_mat, class_levels) {
  y_pred <- class_levels[max.col(prob_mat, ties.method = "first")]
  mean(y_pred == as.character(y_true))
}

# Fit multinomial glmnet with hard labels.
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
      lambda = cv_fit$lambda.1se,   # more stable than lambda.min
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
# Predict class probabilities from multinomial glmnet.
predict_multinom_prob <- function(fit_obj, x, class_levels) {
  pred <- predict(
    fit_obj$fit,
    newx = x,
    type = "response",
    s = fit_obj$lambda
  )
  
  # glmnet multinomial prediction is usually n x K x lambda.
  if (length(dim(pred)) == 3) {
    prob_mat <- pred[, , 1]
  } else {
    prob_mat <- pred
  }
  
  prob_mat <- as.matrix(prob_mat)
  
  # Make sure columns are in the desired class order.
  prob_mat <- prob_mat[, class_levels, drop = FALSE]
  prob_mat <- row_normalize(prob_mat)
  
  return(prob_mat)
}

# Fit multinomial model with soft labels by expanding data.
#
# For each cell i and each class k, create one pseudo-observation
# with label k and weight q_ik.
#
# This lets glmnet fit a deterministic soft-label multinomial model.
fit_multinom_soft <- function(x, q_mat, alpha = 0, lambda = NULL) {
  
  q_mat <- row_normalize(q_mat)
  
  n <- nrow(x)
  K <- ncol(q_mat)
  
  x_expanded <- x[rep(seq_len(n), times = K), , drop = FALSE]
  y_expanded <- rep(colnames(q_mat), each = n)
  w_expanded <- as.vector(q_mat)
  
  y_expanded <- factor(y_expanded, levels = class_levels)
  
  # Drop very tiny weights
  keep <- w_expanded > 1e-6
  
  x_expanded <- x_expanded[keep, , drop = FALSE]
  y_expanded <- y_expanded[keep]
  w_expanded <- w_expanded[keep]
  
  # Renormalize weights so glmnet is numerically happier
  w_expanded <- w_expanded / mean(w_expanded)
  
  # Use a safer lambda.
  # Reusing the naive lambda can be too weak for soft-label expanded data.
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
# 5. Naive multinomial baseline
############################################################
# This model treats inherited patient-level stage category as if it
# were the true cell-level category.

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
# 6. Deterministic categorical MMIL-style EM
############################################################
# Latent cell label:
#   Y_i in {I, II_III, IV}
#
# Observed patient-level label:
#   Stage_cat
#
# Intuition:
#   Cells from a Stage I patient are not all forced to be truly
#   Stage I-like forever. Instead, each cell gets a soft probability
#   vector over the three latent classes.
#
# But we add a weak-supervision constraint:
#   cells from class s should, on average, keep relatively high
#   probability mass on their observed patient-level class s.
#
# This is deterministic EM:
#   E-step: update q_ik using model probabilities and class constraints.
#   M-step: fit multinomial classifier using soft labels q_ik.
#
# No Monte Carlo sampling is used.
############################################################

cat("\n==============================\n")
cat("Fitting deterministic categorical MMIL baseline\n")
cat("==============================\n")

max_iter <- 40
tol <- 1e-4

# Strength of the weak patient-level label.
# Higher value means cells are more constrained to their patient's stage.
# Lower value allows more cross-stage cell-level uncertainty.
label_strength <- 0.70

cat("\nLabel strength:\n")
print(label_strength)

# Initialize q matrix using inherited patient label, but softened.
# For K = 3 and label_strength = 0.70:
#   observed class gets 0.70
#   other two classes share 0.30
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

em_log <- data.frame(
  iter = integer(),
  delta = numeric(),
  train_soft_entropy = numeric(),
  train_log_loss_inherited = numeric(),
  stringsAsFactors = FALSE
)

cat_em_fit <- NULL

for (iter in 1:max_iter) {
  
  old_q <- q_mat
  
  ##########################################################
  # M-step:
  # Fit multinomial classifier using current soft labels.
  ##########################################################
  
  cat_em_fit <- fit_multinom_soft(
    x = x_train,
    q_mat = q_mat,
    alpha = 0,
    lambda = max(naive_fit$lambda, 0.01)
  )
  
  p_hat <- predict_multinom_prob(cat_em_fit, x_train, class_levels)
  
  ##########################################################
  # E-step:
  # Update q using model probabilities + weak observed label.
  ##########################################################
  
  q_new <- p_hat
  
  # Apply a deterministic weak-label constraint:
  # each cell's posterior is multiplied by a prior that favors
  # the patient's observed class.
  label_prior <- matrix(
    (1 - label_strength) / (K - 1),
    nrow = nrow(train_df),
    ncol = K
  )
  colnames(label_prior) <- class_levels
  
  for (i in seq_len(nrow(label_prior))) {
    label_prior[i, obs_idx[i]] <- label_strength
  }
  
  q_new <- q_new * label_prior
  q_new <- row_normalize(q_new)
  
  q_mat <- q_new
  
  delta <- mean(abs(q_mat - old_q))
  
  entropy <- -mean(rowSums(q_mat * log(pmax(q_mat, 1e-12))))
  
  inherited_log_loss <- multiclass_log_loss(
    y_true = y_train_cat,
    prob_mat = p_hat,
    class_levels = class_levels
  )
  
  em_log <- rbind(
    em_log,
    data.frame(
      iter = iter,
      delta = delta,
      train_soft_entropy = entropy,
      train_log_loss_inherited = inherited_log_loss
    )
  )
  
  cat(
    "Categorical EM iter", iter,
    "| delta =", round(delta, 6),
    "| entropy =", round(entropy, 4),
    "| inherited train log loss =", round(inherited_log_loss, 4),
    "\n"
  )
  
  if (delta < tol) {
    cat("Categorical deterministic EM converged.\n")
    break
  }
}

############################################################
# 7. Evaluate categorical EM baseline
############################################################

cat_em_train_prob <- predict_multinom_prob(cat_em_fit, x_train, class_levels)
cat_em_test_prob  <- predict_multinom_prob(cat_em_fit, x_test, class_levels)

cat("\nCategorical EM train log loss against inherited labels:\n")
print(multiclass_log_loss(y_train_cat, cat_em_train_prob, class_levels))

cat("\nCategorical EM test log loss against inherited labels:\n")
print(multiclass_log_loss(y_test_cat, cat_em_test_prob, class_levels))

cat("\nCategorical EM train accuracy against inherited labels:\n")
print(multiclass_accuracy(y_train_cat, cat_em_train_prob, class_levels))

cat("\nCategorical EM test accuracy against inherited labels:\n")
print(multiclass_accuracy(y_test_cat, cat_em_test_prob, class_levels))

cat("\nCategorical EM test confusion matrix:\n")
cat_em_test_pred <- class_levels[max.col(cat_em_test_prob, ties.method = "first")]
print(table(True = y_test_cat, Predicted = cat_em_test_pred))

############################################################
# 8. Add predictions back to metadata
############################################################

# Training results
train_results <- train_df

for (k in class_levels) {
  train_results[[paste0("naive_prob_", k)]] <- naive_train_prob[, k]
  train_results[[paste0("cat_em_prob_", k)]] <- cat_em_train_prob[, k]
  train_results[[paste0("cat_em_soft_label_", k)]] <- q_mat[, k]
}

train_results$naive_pred <- class_levels[max.col(naive_train_prob, ties.method = "first")]
train_results$cat_em_pred <- class_levels[max.col(cat_em_train_prob, ties.method = "first")]
train_results$split <- "train"

# Test results
test_results <- test_df

for (k in class_levels) {
  test_results[[paste0("naive_prob_", k)]] <- naive_test_prob[, k]
  test_results[[paste0("cat_em_prob_", k)]] <- cat_em_test_prob[, k]
  test_results[[paste0("cat_em_soft_label_", k)]] <- NA_real_
}

test_results$naive_pred <- class_levels[max.col(naive_test_prob, ties.method = "first")]
test_results$cat_em_pred <- class_levels[max.col(cat_em_test_prob, ties.method = "first")]
test_results$split <- "test"

all_results <- bind_rows(train_results, test_results)

############################################################
# 9. Summarize by observed stage and cell type
############################################################

summary_by_celltype <- all_results %>%
  group_by(split, Stage_cat, Cell_type) %>%
  summarise(
    n_cells = n(),
    mean_naive_prob_I = mean(naive_prob_I, na.rm = TRUE),
    mean_naive_prob_II_III = mean(naive_prob_II_III, na.rm = TRUE),
    mean_naive_prob_IV = mean(naive_prob_IV, na.rm = TRUE),
    mean_cat_em_prob_I = mean(cat_em_prob_I, na.rm = TRUE),
    mean_cat_em_prob_II_III = mean(cat_em_prob_II_III, na.rm = TRUE),
    mean_cat_em_prob_IV = mean(cat_em_prob_IV, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(split, Stage_cat, desc(mean_cat_em_prob_IV))

cat("\nSummary by cell type:\n")
print(summary_by_celltype, n = 50)

############################################################
# 10. Save results
############################################################

write.csv(
  all_results,
  file = file.path(proj_dir, "stage_categorical_EM_MMIL_cell_predictions.csv"),
  row.names = FALSE
)

write.csv(
  summary_by_celltype,
  file = file.path(proj_dir, "stage_categorical_EM_MMIL_summary_by_celltype.csv"),
  row.names = FALSE
)

write.csv(
  em_log,
  file = file.path(proj_dir, "stage_categorical_EM_MMIL_log.csv"),
  row.names = FALSE
)

############################################################
# 11. Plots
############################################################

# Convergence plot
p_conv <- ggplot(em_log, aes(x = iter, y = delta)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(
    title = "Categorical deterministic EM-MMIL convergence",
    x = "Iteration",
    y = "Mean absolute change in soft labels"
  )

ggsave(
  filename = file.path(proj_dir, "plot_stage_categorical_EM_convergence.png"),
  plot = p_conv,
  width = 7,
  height = 5,
  dpi = 300
)

# Predicted probability distributions for each class
prob_long <- all_results %>%
  select(
    Cell,
    split,
    Stage_cat,
    Cell_type,
    starts_with("naive_prob_"),
    starts_with("cat_em_prob_")
  ) %>%
  pivot_longer(
    cols = starts_with("naive_prob_") | starts_with("cat_em_prob_"),
    names_to = "model_class",
    values_to = "probability"
  ) %>%
  mutate(
    model = ifelse(grepl("^naive", model_class), "Naive multinomial", "Categorical EM-MMIL"),
    class = gsub("^naive_prob_|^cat_em_prob_", "", model_class)
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
  filename = file.path(proj_dir, "plot_stage_categorical_probabilities.png"),
  plot = p_prob,
  width = 12,
  height = 7,
  dpi = 300
)

# Mean categorical EM probabilities by cell type
celltype_long <- all_results %>%
  group_by(Stage_cat, Cell_type) %>%
  summarise(
    n_cells = n(),
    mean_prob_I = mean(cat_em_prob_I, na.rm = TRUE),
    mean_prob_II_III = mean(cat_em_prob_II_III, na.rm = TRUE),
    mean_prob_IV = mean(cat_em_prob_IV, na.rm = TRUE),
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
    title = "Mean categorical EM-MMIL latent-stage probability by cell type",
    x = "Cell type",
    y = "Mean predicted probability"
  )

ggsave(
  filename = file.path(proj_dir, "plot_stage_categorical_EM_by_celltype.png"),
  plot = p_celltype,
  width = 12,
  height = 7,
  dpi = 300
)

############################################################
# 12. Final console summary
############################################################

cat("\nSaved results:\n")
cat("stage_categorical_EM_MMIL_cell_predictions.csv\n")
cat("stage_categorical_EM_MMIL_summary_by_celltype.csv\n")
cat("stage_categorical_EM_MMIL_log.csv\n")
cat("plot_stage_categorical_EM_convergence.png\n")
cat("plot_stage_categorical_probabilities.png\n")
cat("plot_stage_categorical_EM_by_celltype.png\n")

cat("\nDone.\n")