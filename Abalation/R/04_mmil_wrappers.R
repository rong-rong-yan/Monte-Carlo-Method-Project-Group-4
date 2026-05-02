############################################################
# 04_mmil_wrappers.R
# Binary and categorical Naive / EM-MMIL / MCEM-MMIL wrappers.
#
# Required previous scripts:
#   00_utils.R
#   03_model_backends.R
############################################################

############################################################
# Shared helpers
############################################################

sanitize_class_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

get_metadata_columns_for_output <- function(df) {
  metadata_cols <- c(
    "Cell",
    "Patient",
    "split",
    "Stage",
    "Stage_broad",
    "Stage_group",
    "z_obs",
    "Stage_cat",
    "y_cat",
    "Cell_type",
    "Cell_type.refined",
    "Cell_subtype",
    "Tissue",
    "Sample",
    "Sample_Origin"
  )
  
  intersect(metadata_cols, colnames(df))
}

############################################################
# Binary helpers
############################################################

update_binary_q <- function(
  p_hat,
  z_train,
  rho = 0.70,
  eps = 1e-6
) {
  # z_train:
  #   0 = Early
  #   1 = Advanced
  #
  # rho:
  #   assumed baseline-like fraction among Advanced cells.
  #
  # target_mean:
  #   expected advanced-associated fraction among Advanced cells.
  
  early_idx <- z_train == 0
  adv_idx <- z_train == 1
  
  target_mean <- 1 - rho
  
  q_new <- rep(0, length(z_train))
  q_new[early_idx] <- 0
  
  if (sum(adv_idx) == 0) {
    return(q_new)
  }
  
  p_adv <- p_hat[adv_idx]
  current_mean <- mean(p_adv, na.rm = TRUE)
  
  if (!is.finite(current_mean) || current_mean <= eps) {
    q_adv <- rep(target_mean, sum(adv_idx))
  } else {
    q_adv <- p_adv * target_mean / current_mean
  }
  
  q_adv <- pmin(pmax(q_adv, eps), 1 - eps)
  q_new[adv_idx] <- q_adv
  
  return(q_new)
}

monte_carlo_average_binary <- function(
  q,
  z_train,
  mcem_samples = 30
) {
  n <- length(q)
  sampled_sum <- rep(0, n)
  
  early_idx <- z_train == 0
  adv_idx <- z_train == 1
  
  for (m in seq_len(mcem_samples)) {
    sampled_y <- rep(0, n)
    sampled_y[early_idx] <- 0
    
    if (sum(adv_idx) > 0) {
      sampled_y[adv_idx] <- rbinom(
        n = sum(adv_idx),
        size = 1,
        prob = q[adv_idx]
      )
    }
    
    sampled_sum <- sampled_sum + sampled_y
  }
  
  q_mc_avg <- sampled_sum / mcem_samples
  
  return(q_mc_avg)
}

############################################################
# Binary naive baseline
############################################################

fit_binary_naive <- function(
  x_train,
  x_test,
  z_train,
  backend_config
) {
  cat("\n==============================\n")
  cat("Fitting binary naive inherited-label baseline\n")
  cat("==============================\n")
  
  naive_fit <- fit_backend_binary(
    x = x_train,
    y_soft = z_train,
    backend_config = backend_config
  )
  
  naive_train_prob <- predict_backend_binary(
    fit_obj = naive_fit,
    x = x_train
  )
  
  naive_test_prob <- predict_backend_binary(
    fit_obj = naive_fit,
    x = x_test
  )
  
  cat("\nBinary naive train log loss:\n")
  print(log_loss_binary(z_train, naive_train_prob))
  
  return(
    list(
      fit = naive_fit,
      train_prob = naive_train_prob,
      test_prob = naive_test_prob
    )
  )
}

############################################################
# Binary deterministic EM-MMIL
############################################################

fit_binary_em <- function(
  x_train,
  x_test,
  z_train,
  backend_config,
  naive_fit = NULL,
  rho = 0.70,
  max_iter = 30,
  tol = 1e-4
) {
  cat("\n==============================\n")
  cat("Fitting binary deterministic EM-MMIL\n")
  cat("==============================\n")
  
  q <- ifelse(z_train == 1, 1 - rho, 0)
  
  em_log <- data.frame(
    iter = integer(),
    delta = numeric(),
    mean_q_advanced = numeric(),
    train_logloss_inherited = numeric(),
    stringsAsFactors = FALSE
  )
  
  em_backend_config <- backend_config
  
  if (!is.null(naive_fit)) {
    em_backend_config <- reuse_glmnet_lambda_if_requested(
      backend_config = em_backend_config,
      reference_fit = naive_fit
    )
  }
  
  for (iter in seq_len(max_iter)) {
    old_q <- q
    
    em_fit_iter <- fit_backend_binary(
      x = x_train,
      y_soft = q,
      backend_config = em_backend_config
    )
    
    p_hat <- predict_backend_binary(
      fit_obj = em_fit_iter,
      x = x_train
    )
    
    q <- update_binary_q(
      p_hat = p_hat,
      z_train = z_train,
      rho = rho
    )
    
    delta <- mean(abs(q - old_q))
    
    train_logloss_inherited <- log_loss_binary(
      y = z_train,
      p = p_hat
    )
    
    mean_q_advanced <- ifelse(
      sum(z_train == 1) > 0,
      mean(q[z_train == 1]),
      NA_real_
    )
    
    em_log <- rbind(
      em_log,
      data.frame(
        iter = iter,
        delta = delta,
        mean_q_advanced = mean_q_advanced,
        train_logloss_inherited = train_logloss_inherited,
        stringsAsFactors = FALSE
      )
    )
    
    cat(
      "Binary EM iter =", iter,
      "| delta =", signif(delta, 4),
      "| mean q advanced =", signif(mean_q_advanced, 4),
      "| train logloss inherited =", signif(train_logloss_inherited, 4),
      "\n"
    )
    
    if (delta < tol) {
      cat("Binary EM converged.\n")
      break
    }
  }
  
  em_fit <- fit_backend_binary(
    x = x_train,
    y_soft = q,
    backend_config = em_backend_config
  )
  
  em_train_prob <- predict_backend_binary(
    fit_obj = em_fit,
    x = x_train
  )
  
  em_test_prob <- predict_backend_binary(
    fit_obj = em_fit,
    x = x_test
  )
  
  return(
    list(
      fit = em_fit,
      train_prob = em_train_prob,
      test_prob = em_test_prob,
      soft_label = q,
      log = em_log
    )
  )
}

############################################################
# Binary MCEM-MMIL
############################################################

fit_binary_mcem <- function(
  x_train,
  x_test,
  z_train,
  backend_config,
  naive_fit = NULL,
  rho = 0.70,
  max_iter = 30,
  tol = 1e-4,
  mcem_samples = 30,
  seed = 2026
) {
  cat("\n==============================\n")
  cat("Fitting binary MCEM-MMIL\n")
  cat("==============================\n")
  
  set.seed(seed)
  
  q <- ifelse(z_train == 1, 1 - rho, 0)
  
  mcem_log <- data.frame(
    iter = integer(),
    delta = numeric(),
    mean_q_advanced = numeric(),
    mean_q_mc_advanced = numeric(),
    train_logloss_inherited = numeric(),
    stringsAsFactors = FALSE
  )
  
  mcem_backend_config <- backend_config
  
  if (!is.null(naive_fit)) {
    mcem_backend_config <- reuse_glmnet_lambda_if_requested(
      backend_config = mcem_backend_config,
      reference_fit = naive_fit
    )
  }
  
  for (iter in seq_len(max_iter)) {
    old_q <- q
    
    # Method A:
    #   sample latent labels multiple times,
    #   average sampled labels,
    #   train one model using the averaged soft labels.
    q_mc_avg <- monte_carlo_average_binary(
      q = q,
      z_train = z_train,
      mcem_samples = mcem_samples
    )
    
    mcem_fit_iter <- fit_backend_binary(
      x = x_train,
      y_soft = q_mc_avg,
      backend_config = mcem_backend_config
    )
    
    p_hat <- predict_backend_binary(
      fit_obj = mcem_fit_iter,
      x = x_train
    )
    
    q <- update_binary_q(
      p_hat = p_hat,
      z_train = z_train,
      rho = rho
    )
    
    delta <- mean(abs(q - old_q))
    
    train_logloss_inherited <- log_loss_binary(
      y = z_train,
      p = p_hat
    )
    
    mean_q_advanced <- ifelse(
      sum(z_train == 1) > 0,
      mean(q[z_train == 1]),
      NA_real_
    )
    
    mean_q_mc_advanced <- ifelse(
      sum(z_train == 1) > 0,
      mean(q_mc_avg[z_train == 1]),
      NA_real_
    )
    
    mcem_log <- rbind(
      mcem_log,
      data.frame(
        iter = iter,
        delta = delta,
        mean_q_advanced = mean_q_advanced,
        mean_q_mc_advanced = mean_q_mc_advanced,
        train_logloss_inherited = train_logloss_inherited,
        stringsAsFactors = FALSE
      )
    )
    
    cat(
      "Binary MCEM iter =", iter,
      "| delta =", signif(delta, 4),
      "| mean q advanced =", signif(mean_q_advanced, 4),
      "| mean MC q advanced =", signif(mean_q_mc_advanced, 4),
      "| train logloss inherited =", signif(train_logloss_inherited, 4),
      "\n"
    )
    
    if (delta < tol) {
      cat("Binary MCEM converged.\n")
      break
    }
  }
  
  q_mc_final <- monte_carlo_average_binary(
    q = q,
    z_train = z_train,
    mcem_samples = mcem_samples
  )
  
  mcem_fit <- fit_backend_binary(
    x = x_train,
    y_soft = q_mc_final,
    backend_config = mcem_backend_config
  )
  
  mcem_train_prob <- predict_backend_binary(
    fit_obj = mcem_fit,
    x = x_train
  )
  
  mcem_test_prob <- predict_backend_binary(
    fit_obj = mcem_fit,
    x = x_test
  )
  
  return(
    list(
      fit = mcem_fit,
      train_prob = mcem_train_prob,
      test_prob = mcem_test_prob,
      soft_label = q,
      soft_label_mc_final = q_mc_final,
      log = mcem_log
    )
  )
}

############################################################
# Binary prediction output
############################################################

build_binary_prediction_df <- function(
  train_df,
  test_df,
  feature_cols,
  naive_result,
  em_result,
  mcem_result
) {
  keep_cols <- union(
    get_metadata_columns_for_output(train_df),
    get_metadata_columns_for_output(test_df)
  )
  
  train_out <- train_df[, keep_cols, drop = FALSE]
  test_out <- test_df[, keep_cols, drop = FALSE]
  
  for (col in setdiff(keep_cols, colnames(train_out))) {
    train_out[[col]] <- NA
  }
  
  for (col in setdiff(keep_cols, colnames(test_out))) {
    test_out[[col]] <- NA
  }
  
  train_out <- train_out[, keep_cols, drop = FALSE]
  test_out <- test_out[, keep_cols, drop = FALSE]
  
  train_out$naive_prob <- naive_result$train_prob
  train_out$em_mmil_prob <- em_result$train_prob
  train_out$mcem_mmil_prob <- mcem_result$train_prob
  
  train_out$em_soft_label <- em_result$soft_label
  train_out$mcem_soft_label <- mcem_result$soft_label
  
  test_out$naive_prob <- naive_result$test_prob
  test_out$em_mmil_prob <- em_result$test_prob
  test_out$mcem_mmil_prob <- mcem_result$test_prob
  
  test_out$em_soft_label <- NA_real_
  test_out$mcem_soft_label <- NA_real_
  
  pred_df <- rbind(train_out, test_out)
  
  if (!("Cell_type" %in% colnames(pred_df))) {
    pred_df$Cell_type <- "Unknown"
  }
  
  return(pred_df)
}

############################################################
# Main binary MMIL runner
############################################################

run_binary_mmil <- function(
  train_df,
  test_df,
  feature_cols,
  backend_config,
  rho = 0.70,
  max_iter = 30,
  tol = 1e-4,
  mcem_samples = 30,
  seed = 2026
) {
  check_required_cols(
    train_df,
    c("Stage_group", "z_obs", feature_cols),
    object_name = "train_df"
  )
  
  check_required_cols(
    test_df,
    c("Stage_group", "z_obs", feature_cols),
    object_name = "test_df"
  )
  
  x_train <- as.matrix(train_df[, feature_cols, drop = FALSE])
  x_test <- as.matrix(test_df[, feature_cols, drop = FALSE])
  
  z_train <- train_df$z_obs
  z_test <- test_df$z_obs
  
  cat("\nBinary MMIL run summary:\n")
  cat("Training cells:", nrow(x_train), "\n")
  cat("Test cells:", nrow(x_test), "\n")
  cat("Features:", length(feature_cols), "\n")
  cat("Backend:", backend_config$backend, "\n")
  
  naive_result <- fit_binary_naive(
    x_train = x_train,
    x_test = x_test,
    z_train = z_train,
    backend_config = backend_config
  )
  
  cat("\nBinary naive test log loss using inherited labels:\n")
  print(log_loss_binary(z_test, naive_result$test_prob))
  
  cat("\nBinary naive test AUC using inherited labels:\n")
  print(compute_auc_safe(z_test, naive_result$test_prob))
  
  em_result <- fit_binary_em(
    x_train = x_train,
    x_test = x_test,
    z_train = z_train,
    backend_config = backend_config,
    naive_fit = naive_result$fit,
    rho = rho,
    max_iter = max_iter,
    tol = tol
  )
  
  cat("\nBinary EM test log loss using inherited labels:\n")
  print(log_loss_binary(z_test, em_result$test_prob))
  
  cat("\nBinary EM test AUC using inherited labels:\n")
  print(compute_auc_safe(z_test, em_result$test_prob))
  
  mcem_result <- fit_binary_mcem(
    x_train = x_train,
    x_test = x_test,
    z_train = z_train,
    backend_config = backend_config,
    naive_fit = naive_result$fit,
    rho = rho,
    max_iter = max_iter,
    tol = tol,
    mcem_samples = mcem_samples,
    seed = seed
  )
  
  cat("\nBinary MCEM test log loss using inherited labels:\n")
  print(log_loss_binary(z_test, mcem_result$test_prob))
  
  cat("\nBinary MCEM test AUC using inherited labels:\n")
  print(compute_auc_safe(z_test, mcem_result$test_prob))
  
  pred_df <- build_binary_prediction_df(
    train_df = train_df,
    test_df = test_df,
    feature_cols = feature_cols,
    naive_result = naive_result,
    em_result = em_result,
    mcem_result = mcem_result
  )
  
  return(
    list(
      task_type = "binary",
      predictions = pred_df,
      naive = naive_result,
      em = em_result,
      mcem = mcem_result
    )
  )
}

############################################################
# Categorical helpers
############################################################

build_label_prior_matrix <- function(
  y_class,
  class_levels,
  label_strength = 0.70
) {
  y_class <- factor(y_class, levels = class_levels)
  
  n <- length(y_class)
  K <- length(class_levels)
  
  if (K < 2) {
    stop("Categorical MMIL requires at least two classes.")
  }
  
  off_diag <- (1 - label_strength) / (K - 1)
  
  prior <- matrix(
    off_diag,
    nrow = n,
    ncol = K
  )
  
  colnames(prior) <- class_levels
  
  for (k in seq_along(class_levels)) {
    prior[y_class == class_levels[k], k] <- label_strength
  }
  
  prior <- row_normalize(prior)
  
  return(prior)
}

sample_categorical_onehot <- function(q_mat, class_levels) {
  q_mat <- as.matrix(q_mat)
  q_mat <- row_normalize(q_mat)
  q_mat <- q_mat[, class_levels, drop = FALSE]
  
  n <- nrow(q_mat)
  K <- length(class_levels)
  
  sampled <- matrix(
    0,
    nrow = n,
    ncol = K
  )
  
  colnames(sampled) <- class_levels
  
  for (i in seq_len(n)) {
    sampled_class <- sample.int(
      n = K,
      size = 1,
      prob = q_mat[i, ]
    )
    
    sampled[i, sampled_class] <- 1
  }
  
  return(sampled)
}

monte_carlo_average_multiclass <- function(
  q_mat,
  class_levels,
  mcem_samples = 30
) {
  q_mat <- as.matrix(q_mat)
  q_mat <- row_normalize(q_mat)
  q_mat <- q_mat[, class_levels, drop = FALSE]
  
  sampled_sum <- matrix(
    0,
    nrow = nrow(q_mat),
    ncol = ncol(q_mat)
  )
  
  colnames(sampled_sum) <- class_levels
  
  for (m in seq_len(mcem_samples)) {
    sampled_sum <- sampled_sum + sample_categorical_onehot(
      q_mat = q_mat,
      class_levels = class_levels
    )
  }
  
  q_mc_avg <- sampled_sum / mcem_samples
  q_mc_avg <- row_normalize(q_mc_avg)
  
  return(q_mc_avg)
}

update_categorical_q <- function(
  p_hat,
  label_prior
) {
  p_hat <- as.matrix(p_hat)
  label_prior <- as.matrix(label_prior)
  
  if (!all(dim(p_hat) == dim(label_prior))) {
    stop("p_hat and label_prior must have the same dimensions.")
  }
  
  q_new <- p_hat * label_prior
  q_new <- row_normalize(q_new)
  
  return(q_new)
}

############################################################
# Categorical naive baseline
############################################################

fit_categorical_naive <- function(
  x_train,
  x_test,
  y_train_class,
  class_levels,
  backend_config
) {
  cat("\n==============================\n")
  cat("Fitting categorical naive inherited-label baseline\n")
  cat("==============================\n")
  
  q_naive <- make_onehot(
    y = y_train_class,
    class_levels = class_levels
  )
  
  naive_fit <- fit_backend_multiclass(
    x = x_train,
    q_mat = q_naive,
    class_levels = class_levels,
    backend_config = backend_config
  )
  
  naive_train_prob <- predict_backend_multiclass(
    fit_obj = naive_fit,
    x = x_train
  )
  
  naive_test_prob <- predict_backend_multiclass(
    fit_obj = naive_fit,
    x = x_test
  )
  
  cat("\nCategorical naive train log loss:\n")
  print(
    multiclass_log_loss(
      y_true = y_train_class,
      prob_mat = naive_train_prob,
      class_levels = class_levels
    )
  )
  
  return(
    list(
      fit = naive_fit,
      train_prob = naive_train_prob,
      test_prob = naive_test_prob,
      q_naive = q_naive
    )
  )
}

############################################################
# Categorical deterministic EM-MMIL
############################################################

fit_categorical_em <- function(
  x_train,
  x_test,
  y_train_class,
  class_levels,
  backend_config,
  naive_fit = NULL,
  label_strength = 0.70,
  max_iter = 40,
  tol = 1e-4
) {
  cat("\n==============================\n")
  cat("Fitting categorical deterministic EM-MMIL\n")
  cat("==============================\n")
  
  label_prior <- build_label_prior_matrix(
    y_class = y_train_class,
    class_levels = class_levels,
    label_strength = label_strength
  )
  
  q_mat <- label_prior
  
  em_log <- data.frame(
    iter = integer(),
    delta = numeric(),
    mean_entropy = numeric(),
    train_logloss_inherited = numeric(),
    train_accuracy_inherited = numeric(),
    stringsAsFactors = FALSE
  )
  
  em_backend_config <- backend_config
  
  if (!is.null(naive_fit)) {
    em_backend_config <- reuse_glmnet_lambda_if_requested(
      backend_config = em_backend_config,
      reference_fit = naive_fit
    )
  }
  
  for (iter in seq_len(max_iter)) {
    old_q <- q_mat
    
    em_fit_iter <- fit_backend_multiclass(
      x = x_train,
      q_mat = q_mat,
      class_levels = class_levels,
      backend_config = em_backend_config
    )
    
    p_hat <- predict_backend_multiclass(
      fit_obj = em_fit_iter,
      x = x_train
    )
    
    q_mat <- update_categorical_q(
      p_hat = p_hat,
      label_prior = label_prior
    )
    
    delta <- mean(abs(q_mat - old_q))
    
    entropy <- -rowSums(q_mat * log(pmax(q_mat, 1e-12)))
    mean_entropy <- mean(entropy)
    
    train_logloss_inherited <- multiclass_log_loss(
      y_true = y_train_class,
      prob_mat = p_hat,
      class_levels = class_levels
    )
    
    train_accuracy_inherited <- multiclass_accuracy(
      y_true = y_train_class,
      prob_mat = p_hat,
      class_levels = class_levels
    )
    
    em_log <- rbind(
      em_log,
      data.frame(
        iter = iter,
        delta = delta,
        mean_entropy = mean_entropy,
        train_logloss_inherited = train_logloss_inherited,
        train_accuracy_inherited = train_accuracy_inherited,
        stringsAsFactors = FALSE
      )
    )
    
    cat(
      "Categorical EM iter =", iter,
      "| delta =", signif(delta, 4),
      "| mean entropy =", signif(mean_entropy, 4),
      "| train logloss inherited =", signif(train_logloss_inherited, 4),
      "| train accuracy inherited =", signif(train_accuracy_inherited, 4),
      "\n"
    )
    
    if (delta < tol) {
      cat("Categorical EM converged.\n")
      break
    }
  }
  
  em_fit <- fit_backend_multiclass(
    x = x_train,
    q_mat = q_mat,
    class_levels = class_levels,
    backend_config = em_backend_config
  )
  
  em_train_prob <- predict_backend_multiclass(
    fit_obj = em_fit,
    x = x_train
  )
  
  em_test_prob <- predict_backend_multiclass(
    fit_obj = em_fit,
    x = x_test
  )
  
  return(
    list(
      fit = em_fit,
      train_prob = em_train_prob,
      test_prob = em_test_prob,
      soft_label = q_mat,
      label_prior = label_prior,
      log = em_log
    )
  )
}

############################################################
# Categorical MCEM-MMIL
############################################################

fit_categorical_mcem <- function(
  x_train,
  x_test,
  y_train_class,
  class_levels,
  backend_config,
  naive_fit = NULL,
  label_strength = 0.70,
  max_iter = 40,
  tol = 1e-4,
  mcem_samples = 30,
  seed = 2026
) {
  cat("\n==============================\n")
  cat("Fitting categorical MCEM-MMIL\n")
  cat("==============================\n")
  
  set.seed(seed)
  
  label_prior <- build_label_prior_matrix(
    y_class = y_train_class,
    class_levels = class_levels,
    label_strength = label_strength
  )
  
  q_mat <- label_prior
  
  mcem_log <- data.frame(
    iter = integer(),
    delta = numeric(),
    mean_entropy = numeric(),
    mean_mc_entropy = numeric(),
    train_logloss_inherited = numeric(),
    train_accuracy_inherited = numeric(),
    stringsAsFactors = FALSE
  )
  
  mcem_backend_config <- backend_config
  
  if (!is.null(naive_fit)) {
    mcem_backend_config <- reuse_glmnet_lambda_if_requested(
      backend_config = mcem_backend_config,
      reference_fit = naive_fit
    )
  }
  
  for (iter in seq_len(max_iter)) {
    old_q <- q_mat
    
    # Method A:
    #   sample latent categorical labels multiple times,
    #   average sampled one-hot matrices,
    #   train one model using the averaged soft labels.
    q_mc_avg <- monte_carlo_average_multiclass(
      q_mat = q_mat,
      class_levels = class_levels,
      mcem_samples = mcem_samples
    )
    
    mcem_fit_iter <- fit_backend_multiclass(
      x = x_train,
      q_mat = q_mc_avg,
      class_levels = class_levels,
      backend_config = mcem_backend_config
    )
    
    p_hat <- predict_backend_multiclass(
      fit_obj = mcem_fit_iter,
      x = x_train
    )
    
    q_mat <- update_categorical_q(
      p_hat = p_hat,
      label_prior = label_prior
    )
    
    delta <- mean(abs(q_mat - old_q))
    
    entropy <- -rowSums(q_mat * log(pmax(q_mat, 1e-12)))
    mc_entropy <- -rowSums(q_mc_avg * log(pmax(q_mc_avg, 1e-12)))
    
    mean_entropy <- mean(entropy)
    mean_mc_entropy <- mean(mc_entropy)
    
    train_logloss_inherited <- multiclass_log_loss(
      y_true = y_train_class,
      prob_mat = p_hat,
      class_levels = class_levels
    )
    
    train_accuracy_inherited <- multiclass_accuracy(
      y_true = y_train_class,
      prob_mat = p_hat,
      class_levels = class_levels
    )
    
    mcem_log <- rbind(
      mcem_log,
      data.frame(
        iter = iter,
        delta = delta,
        mean_entropy = mean_entropy,
        mean_mc_entropy = mean_mc_entropy,
        train_logloss_inherited = train_logloss_inherited,
        train_accuracy_inherited = train_accuracy_inherited,
        stringsAsFactors = FALSE
      )
    )
    
    cat(
      "Categorical MCEM iter =", iter,
      "| delta =", signif(delta, 4),
      "| mean entropy =", signif(mean_entropy, 4),
      "| mean MC entropy =", signif(mean_mc_entropy, 4),
      "| train logloss inherited =", signif(train_logloss_inherited, 4),
      "| train accuracy inherited =", signif(train_accuracy_inherited, 4),
      "\n"
    )
    
    if (delta < tol) {
      cat("Categorical MCEM converged.\n")
      break
    }
  }
  
  q_mc_final <- monte_carlo_average_multiclass(
    q_mat = q_mat,
    class_levels = class_levels,
    mcem_samples = mcem_samples
  )
  
  mcem_fit <- fit_backend_multiclass(
    x = x_train,
    q_mat = q_mc_final,
    class_levels = class_levels,
    backend_config = mcem_backend_config
  )
  
  mcem_train_prob <- predict_backend_multiclass(
    fit_obj = mcem_fit,
    x = x_train
  )
  
  mcem_test_prob <- predict_backend_multiclass(
    fit_obj = mcem_fit,
    x = x_test
  )
  
  return(
    list(
      fit = mcem_fit,
      train_prob = mcem_train_prob,
      test_prob = mcem_test_prob,
      soft_label = q_mat,
      soft_label_mc_final = q_mc_final,
      label_prior = label_prior,
      log = mcem_log
    )
  )
}

############################################################
# Categorical prediction output
############################################################

add_prob_columns <- function(df, prob_mat, prefix, class_levels) {
  for (cls in class_levels) {
    col_name <- paste0(prefix, "_", sanitize_class_name(cls))
    df[[col_name]] <- prob_mat[, cls]
  }
  
  return(df)
}

add_soft_label_columns <- function(df, q_mat, prefix, class_levels) {
  for (cls in class_levels) {
    col_name <- paste0(prefix, "_", sanitize_class_name(cls))
    df[[col_name]] <- q_mat[, cls]
  }
  
  return(df)
}

build_categorical_prediction_df <- function(
  train_df,
  test_df,
  feature_cols,
  class_levels,
  naive_result,
  em_result,
  mcem_result
) {
  keep_cols <- union(
    get_metadata_columns_for_output(train_df),
    get_metadata_columns_for_output(test_df)
  )
  
  train_out <- train_df[, keep_cols, drop = FALSE]
  test_out <- test_df[, keep_cols, drop = FALSE]
  
  for (col in setdiff(keep_cols, colnames(train_out))) {
    train_out[[col]] <- NA
  }
  
  for (col in setdiff(keep_cols, colnames(test_out))) {
    test_out[[col]] <- NA
  }
  
  train_out <- train_out[, keep_cols, drop = FALSE]
  test_out <- test_out[, keep_cols, drop = FALSE]
  
  train_out <- add_prob_columns(
    df = train_out,
    prob_mat = naive_result$train_prob,
    prefix = "naive_prob",
    class_levels = class_levels
  )
  
  train_out <- add_prob_columns(
    df = train_out,
    prob_mat = em_result$train_prob,
    prefix = "cat_em_prob",
    class_levels = class_levels
  )
  
  train_out <- add_prob_columns(
    df = train_out,
    prob_mat = mcem_result$train_prob,
    prefix = "cat_mcem_prob",
    class_levels = class_levels
  )
  
  train_out <- add_soft_label_columns(
    df = train_out,
    q_mat = em_result$soft_label,
    prefix = "cat_em_soft_label",
    class_levels = class_levels
  )
  
  train_out <- add_soft_label_columns(
    df = train_out,
    q_mat = mcem_result$soft_label,
    prefix = "cat_mcem_soft_label",
    class_levels = class_levels
  )
  
  test_out <- add_prob_columns(
    df = test_out,
    prob_mat = naive_result$test_prob,
    prefix = "naive_prob",
    class_levels = class_levels
  )
  
  test_out <- add_prob_columns(
    df = test_out,
    prob_mat = em_result$test_prob,
    prefix = "cat_em_prob",
    class_levels = class_levels
  )
  
  test_out <- add_prob_columns(
    df = test_out,
    prob_mat = mcem_result$test_prob,
    prefix = "cat_mcem_prob",
    class_levels = class_levels
  )
  
  for (cls in class_levels) {
    test_out[[paste0("cat_em_soft_label_", sanitize_class_name(cls))]] <- NA_real_
    test_out[[paste0("cat_mcem_soft_label_", sanitize_class_name(cls))]] <- NA_real_
  }
  
  pred_df <- rbind(train_out, test_out)
  
  if (!("Cell_type" %in% colnames(pred_df))) {
    pred_df$Cell_type <- "Unknown"
  }
  
  return(pred_df)
}

############################################################
# Main categorical MMIL runner
############################################################

run_categorical_mmil <- function(
  train_df,
  test_df,
  feature_cols,
  backend_config,
  class_levels,
  label_strength = 0.70,
  max_iter = 40,
  tol = 1e-4,
  mcem_samples = 30,
  seed = 2026
) {
  check_required_cols(
    train_df,
    c("Stage_cat", "y_cat", feature_cols),
    object_name = "train_df"
  )
  
  check_required_cols(
    test_df,
    c("Stage_cat", "y_cat", feature_cols),
    object_name = "test_df"
  )
  
  x_train <- as.matrix(train_df[, feature_cols, drop = FALSE])
  x_test <- as.matrix(test_df[, feature_cols, drop = FALSE])
  
  y_train_class <- as.character(train_df$Stage_cat)
  y_test_class <- as.character(test_df$Stage_cat)
  
  cat("\nCategorical MMIL run summary:\n")
  cat("Training cells:", nrow(x_train), "\n")
  cat("Test cells:", nrow(x_test), "\n")
  cat("Features:", length(feature_cols), "\n")
  cat("Backend:", backend_config$backend, "\n")
  cat("Class levels:", paste(class_levels, collapse = ", "), "\n")
  cat("Label strength:", label_strength, "\n")
  
  naive_result <- fit_categorical_naive(
    x_train = x_train,
    x_test = x_test,
    y_train_class = y_train_class,
    class_levels = class_levels,
    backend_config = backend_config
  )
  
  cat("\nCategorical naive test log loss using inherited labels:\n")
  print(
    multiclass_log_loss(
      y_true = y_test_class,
      prob_mat = naive_result$test_prob,
      class_levels = class_levels
    )
  )
  
  cat("\nCategorical naive test accuracy using inherited labels:\n")
  print(
    multiclass_accuracy(
      y_true = y_test_class,
      prob_mat = naive_result$test_prob,
      class_levels = class_levels
    )
  )
  
  em_result <- fit_categorical_em(
    x_train = x_train,
    x_test = x_test,
    y_train_class = y_train_class,
    class_levels = class_levels,
    backend_config = backend_config,
    naive_fit = naive_result$fit,
    label_strength = label_strength,
    max_iter = max_iter,
    tol = tol
  )
  
  cat("\nCategorical EM test log loss using inherited labels:\n")
  print(
    multiclass_log_loss(
      y_true = y_test_class,
      prob_mat = em_result$test_prob,
      class_levels = class_levels
    )
  )
  
  cat("\nCategorical EM test accuracy using inherited labels:\n")
  print(
    multiclass_accuracy(
      y_true = y_test_class,
      prob_mat = em_result$test_prob,
      class_levels = class_levels
    )
  )
  
  mcem_result <- fit_categorical_mcem(
    x_train = x_train,
    x_test = x_test,
    y_train_class = y_train_class,
    class_levels = class_levels,
    backend_config = backend_config,
    naive_fit = naive_result$fit,
    label_strength = label_strength,
    max_iter = max_iter,
    tol = tol,
    mcem_samples = mcem_samples,
    seed = seed
  )
  
  cat("\nCategorical MCEM test log loss using inherited labels:\n")
  print(
    multiclass_log_loss(
      y_true = y_test_class,
      prob_mat = mcem_result$test_prob,
      class_levels = class_levels
    )
  )
  
  cat("\nCategorical MCEM test accuracy using inherited labels:\n")
  print(
    multiclass_accuracy(
      y_true = y_test_class,
      prob_mat = mcem_result$test_prob,
      class_levels = class_levels
    )
  )
  
  pred_df <- build_categorical_prediction_df(
    train_df = train_df,
    test_df = test_df,
    feature_cols = feature_cols,
    class_levels = class_levels,
    naive_result = naive_result,
    em_result = em_result,
    mcem_result = mcem_result
  )
  
  return(
    list(
      task_type = "categorical",
      class_levels = class_levels,
      predictions = pred_df,
      naive = naive_result,
      em = em_result,
      mcem = mcem_result
    )
  )
}

############################################################
# Save outputs
############################################################

save_binary_mmil_outputs <- function(
  mmil_result,
  out_dir,
  pred_filename = "cell_predictions.csv"
) {
  ensure_dir(out_dir)
  
  pred_file <- file.path(out_dir, pred_filename)
  em_log_file <- file.path(out_dir, "em_mmil_log.csv")
  mcem_log_file <- file.path(out_dir, "mcem_mmil_log.csv")
  
  write.csv(
    mmil_result$predictions,
    file = pred_file,
    row.names = FALSE
  )
  
  write.csv(
    mmil_result$em$log,
    file = em_log_file,
    row.names = FALSE
  )
  
  write.csv(
    mmil_result$mcem$log,
    file = mcem_log_file,
    row.names = FALSE
  )
  
  cat("\nSaved binary MMIL outputs:\n")
  cat(pred_file, "\n")
  cat(em_log_file, "\n")
  cat(mcem_log_file, "\n")
  
  return(
    list(
      pred_file = pred_file,
      em_log_file = em_log_file,
      mcem_log_file = mcem_log_file
    )
  )
}

save_categorical_mmil_outputs <- function(
  mmil_result,
  out_dir,
  pred_filename = "cell_predictions.csv"
) {
  ensure_dir(out_dir)
  
  pred_file <- file.path(out_dir, pred_filename)
  em_log_file <- file.path(out_dir, "cat_em_mmil_log.csv")
  mcem_log_file <- file.path(out_dir, "cat_mcem_mmil_log.csv")
  
  write.csv(
    mmil_result$predictions,
    file = pred_file,
    row.names = FALSE
  )
  
  write.csv(
    mmil_result$em$log,
    file = em_log_file,
    row.names = FALSE
  )
  
  write.csv(
    mmil_result$mcem$log,
    file = mcem_log_file,
    row.names = FALSE
  )
  
  cat("\nSaved categorical MMIL outputs:\n")
  cat(pred_file, "\n")
  cat(em_log_file, "\n")
  cat(mcem_log_file, "\n")
  
  return(
    list(
      pred_file = pred_file,
      em_log_file = em_log_file,
      mcem_log_file = mcem_log_file
    )
  )
}

save_mmil_outputs <- function(
  mmil_result,
  out_dir,
  pred_filename = "cell_predictions.csv"
) {
  if (mmil_result$task_type == "binary") {
    return(
      save_binary_mmil_outputs(
        mmil_result = mmil_result,
        out_dir = out_dir,
        pred_filename = pred_filename
      )
    )
  }
  
  if (mmil_result$task_type == "categorical") {
    return(
      save_categorical_mmil_outputs(
        mmil_result = mmil_result,
        out_dir = out_dir,
        pred_filename = pred_filename
      )
    )
  }
  
  stop("Unsupported mmil_result$task_type.")
}