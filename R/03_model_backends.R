############################################################
# 03_model_backends.R
# Model backend functions for binary and categorical MMIL / MCEM.
#
# Supported backends:
#   1. glmnet binary logistic regression
#   2. glmnet multinomial logistic regression
#   3. XGBoost binary logistic model
#   4. XGBoost multiclass softprob model
#
# Binary interface:
#   fit_backend_binary()
#   predict_backend_binary()
#
# Multiclass interface:
#   fit_backend_multiclass()
#   predict_backend_multiclass()
#
# Important glmnet multinomial fix:
#   - Do NOT use cv.glmnet() for multinomial models here.
#   - Do NOT use predict() / coef() for multnet objects here.
#   - Do NOT final-fit multinomial glmnet with a single scalar lambda.
#
# Why:
#   glmnet::predict.multnet() can throw:
#     length of 'dimnames' [2] not equal to array extent
#   especially with multinomial + weighted pseudo-observations + scalar lambda.
#
# Fix:
#   - Manual CV using multnet_predict_manual().
#   - Final multinomial fit uses a lambda path, not one scalar lambda.
#   - Prediction reads glmnet internal slots directly.
#   - Reused-lambda EM / MCEM fits reuse the reference lambda path
#     when available, otherwise build a wide warm-start path.
############################################################

############################################################
# Backend configuration helper
############################################################

get_config_value <- function(config, name, default = NULL) {
  if (!is.null(config[[name]])) {
    return(config[[name]])
  }
  return(default)
}

############################################################
# Weighted soft-label expansion
############################################################

prepare_binary_weighted_data <- function(x, q, eps = 1e-8) {
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
  
  return(
    list(
      x = x_expanded,
      y = y_expanded,
      weight = w_expanded
    )
  )
}

prepare_multiclass_weighted_data <- function(
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
  K <- length(class_levels)
  
  if (nrow(q_mat) != n) {
    stop("nrow(q_mat) must match nrow(x).")
  }
  
  if (ncol(q_mat) != K) {
    stop("ncol(q_mat) must match length(class_levels).")
  }
  
  x_expanded <- x[rep(seq_len(n), times = K), , drop = FALSE]
  y_expanded <- rep(class_levels, each = n)
  w_expanded <- as.vector(q_mat)
  
  keep <- w_expanded > eps
  
  x_expanded <- x_expanded[keep, , drop = FALSE]
  y_expanded <- y_expanded[keep]
  w_expanded <- w_expanded[keep]
  
  if (length(w_expanded) == 0) {
    stop("prepare_multiclass_weighted_data: all pseudo-observation weights are zero.")
  }
  
  w_expanded <- w_expanded / mean(w_expanded)
  
  return(
    list(
      x = x_expanded,
      y = y_expanded,
      weight = w_expanded
    )
  )
}

############################################################
# Safe lambda-path helpers for glmnet multinomial
############################################################

make_safe_lambda_path <- function(lambda, relative_width = 1e-3) {
  lambda <- as.numeric(lambda)[1]
  
  if (!is.finite(lambda) || lambda <= 0) {
    stop("make_safe_lambda_path: lambda must be positive and finite.")
  }
  
  lambda_seq <- sort(
    unique(c(
      lambda * (1 + relative_width),
      lambda,
      lambda * (1 - relative_width)
    )),
    decreasing = TRUE
  )
  
  lambda_seq <- lambda_seq[is.finite(lambda_seq) & lambda_seq > 0]
  
  if (length(lambda_seq) < 2) {
    lambda_seq <- sort(
      unique(c(lambda * 1.01, lambda, lambda * 0.99)),
      decreasing = TRUE
    )
  }
  
  return(lambda_seq)
}

make_wide_lambda_path <- function(
  lambda,
  n_lambda = 60,
  upper_mult = 100,
  lower_mult = 0.01
) {
  ##########################################################
  # Used for reused-lambda multinomial glmnet fits.
  #
  # If we reuse only one selected lambda and build a tiny
  # 3-point path around it, multinomial glmnet may start at
  # an already-too-small lambda, fail to converge at the first
  # value, and throw:
  #   length of 'dimnames' [2] not equal to array extent
  #
  # This wide path starts from strong regularization and lets
  # glmnet warm-start down to the target region.
  ##########################################################
  
  lambda <- as.numeric(lambda)[1]
  
  if (!is.finite(lambda) || lambda <= 0) {
    stop("make_wide_lambda_path: lambda must be positive and finite.")
  }
  
  lambda_max <- lambda * upper_mult
  lambda_min <- lambda * lower_mult
  
  lambda_seq <- exp(
    seq(
      from = log(lambda_max),
      to = log(lambda_min),
      length.out = n_lambda
    )
  )
  
  lambda_seq <- sort(
    unique(lambda_seq[is.finite(lambda_seq) & lambda_seq > 0]),
    decreasing = TRUE
  )
  
  return(lambda_seq)
}

validate_multiclass_training_classes <- function(
  y_factor,
  class_levels,
  context = "multiclass training"
) {
  observed_classes <- unique(as.character(y_factor))
  missing_classes <- setdiff(class_levels, observed_classes)
  
  if (length(missing_classes) > 0) {
    stop(
      paste0(
        context,
        ": missing class(es) after weighted expansion: ",
        paste(missing_classes, collapse = ", "),
        ". Check patient split and q_mat columns."
      )
    )
  }
  
  invisible(TRUE)
}

############################################################
# Manual multinomial prediction (reads fit fields directly)
############################################################

multnet_predict_manual <- function(
  fit,
  x,
  class_levels,
  lambda
) {
  ##########################################################
  # Compute multinomial class probabilities by reading the
  # raw fit slots ($a0, $beta, $classnames, $lambda).
  #
  # This avoids glmnet S3 predict / coef dispatch, which can
  # trigger predict.multnet() and fail with a dimnames error.
  ##########################################################
  
  K <- length(class_levels)
  
  fit_lambdas <- fit$lambda
  
  if (is.null(fit_lambdas) || length(fit_lambdas) == 0) {
    stop("multnet_predict_manual: fit$lambda is empty.")
  }
  
  lambda <- as.numeric(lambda)[1]
  
  if (!is.finite(lambda) || lambda <= 0) {
    stop("multnet_predict_manual: lambda must be positive and finite.")
  }
  
  lambda_idx <- which.min(abs(fit_lambdas - lambda))
  
  fit_classnames <- fit$classnames
  
  if (is.null(fit_classnames)) {
    if (!is.null(names(fit$beta))) {
      fit_classnames <- names(fit$beta)
    } else {
      fit_classnames <- class_levels
    }
  }
  
  if (length(fit_classnames) != K) {
    stop(
      sprintf(
        paste0(
          "multnet_predict_manual: fit$classnames has %d entries ",
          "but class_levels has %d (%s)."
        ),
        length(fit_classnames),
        K,
        paste(class_levels, collapse = ", ")
      )
    )
  }
  
  missing_cls <- setdiff(class_levels, fit_classnames)
  
  if (length(missing_cls) > 0) {
    stop(
      sprintf(
        paste0(
          "multnet_predict_manual: class_levels (%s) missing from ",
          "fit$classnames (%s)."
        ),
        paste(missing_cls, collapse = ", "),
        paste(fit_classnames, collapse = ", ")
      )
    )
  }
  
  class_order <- match(class_levels, fit_classnames)
  
  a0_full <- fit$a0
  
  if (is.null(a0_full)) {
    stop("multnet_predict_manual: fit$a0 is missing.")
  }
  
  a0_mat <- as.matrix(a0_full)
  
  if (nrow(a0_mat) != K) {
    stop(
      sprintf(
        "multnet_predict_manual: fit$a0 has %d rows; expected %d.",
        nrow(a0_mat),
        K
      )
    )
  }
  
  if (lambda_idx > ncol(a0_mat)) {
    stop(
      sprintf(
        paste0(
          "multnet_predict_manual: lambda_idx %d exceeds number of ",
          "columns (%d) in fit$a0."
        ),
        lambda_idx,
        ncol(a0_mat)
      )
    )
  }
  
  beta_list <- fit$beta
  
  if (is.null(beta_list) || !is.list(beta_list) || length(beta_list) != K) {
    stop(
      sprintf(
        "multnet_predict_manual: fit$beta is not a list of length %d.",
        K
      )
    )
  }
  
  n <- nrow(x)
  eta <- matrix(0, nrow = n, ncol = K)
  colnames(eta) <- class_levels
  
  for (k in seq_len(K)) {
    src <- class_order[k]
    
    a0_k <- a0_mat[src, lambda_idx]
    beta_k_full <- beta_list[[src]]
    
    if (is.null(dim(beta_k_full)) || ncol(beta_k_full) < lambda_idx) {
      stop(
        sprintf(
          paste0(
            "multnet_predict_manual: fit$beta[[%d]] has unexpected ",
            "shape (lambda_idx = %d)."
          ),
          src,
          lambda_idx
        )
      )
    }
    
    beta_k <- as.numeric(beta_k_full[, lambda_idx])
    
    if (length(beta_k) != ncol(x)) {
      stop(
        sprintf(
          paste0(
            "multnet_predict_manual: beta for class '%s' has %d entries ",
            "but x has %d columns."
          ),
          class_levels[k],
          length(beta_k),
          ncol(x)
        )
      )
    }
    
    eta[, k] <- a0_k + as.numeric(x %*% beta_k)
  }
  
  row_max <- apply(eta, 1, max, na.rm = TRUE)
  eta <- eta - row_max
  
  expe <- exp(eta)
  prob_mat <- expe / rowSums(expe)
  
  colnames(prob_mat) <- class_levels
  prob_mat <- prob_mat[, class_levels, drop = FALSE]
  prob_mat <- row_normalize(prob_mat)
  
  return(prob_mat)
}

############################################################
# Manual multinomial CV (bypasses cv.glmnet)
############################################################

cv_glmnet_multiclass_manual <- function(
  x,
  y_factor,
  weights,
  class_levels,
  alpha,
  nfolds = 5,
  standardize = TRUE,
  maxit = 1e6,
  seed = NULL,
  lambda_seq = NULL,
  verbose = FALSE
) {
  ##########################################################
  # Replacement for cv.glmnet(family = "multinomial", ...).
  # Uses multnet_predict_manual() to score validation folds,
  # so it never calls predict.multnet().
  ##########################################################
  
  load_pkg("glmnet")
  
  K <- length(class_levels)
  n <- length(y_factor)
  
  validate_multiclass_training_classes(
    y_factor = y_factor,
    class_levels = class_levels,
    context = "cv_glmnet_multiclass_manual full data"
  )
  
  if (is.null(lambda_seq)) {
    full_path <- glmnet::glmnet(
      x = x,
      y = y_factor,
      weights = weights,
      family = "multinomial",
      alpha = alpha,
      standardize = standardize,
      maxit = maxit
    )
    lambda_seq <- full_path$lambda
  }
  
  lambda_seq <- as.numeric(lambda_seq)
  lambda_seq <- lambda_seq[is.finite(lambda_seq) & lambda_seq > 0]
  lambda_seq <- sort(unique(lambda_seq), decreasing = TRUE)
  
  if (length(lambda_seq) < 2) {
    lambda_seq <- make_safe_lambda_path(lambda_seq[1])
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  nfolds_eff <- min(nfolds, n)
  if (nfolds_eff < 2) {
    stop("cv_glmnet_multiclass_manual: need at least 2 observations for CV.")
  }
  
  fold_ids <- sample(rep(seq_len(nfolds_eff), length.out = n))
  
  cv_loss <- matrix(
    NA_real_,
    nrow = nfolds_eff,
    ncol = length(lambda_seq)
  )
  
  for (fold in seq_len(nfolds_eff)) {
    train_idx <- which(fold_ids != fold)
    val_idx   <- which(fold_ids == fold)
    
    if (length(train_idx) == 0 || length(val_idx) == 0) {
      next
    }
    
    if (length(unique(as.character(y_factor[train_idx]))) < K) {
      if (verbose) {
        cat("[manual cv] skipping fold ", fold, " (incomplete classes).\n", sep = "")
      }
      next
    }
    
    fold_fit <- tryCatch(
      glmnet::glmnet(
        x = x[train_idx, , drop = FALSE],
        y = y_factor[train_idx],
        weights = weights[train_idx],
        family = "multinomial",
        alpha = alpha,
        lambda = lambda_seq,
        standardize = standardize,
        maxit = maxit
      ),
      error = function(e) {
        if (verbose) {
          cat("[manual cv] glmnet failed on fold ", fold, ": ", e$message, "\n", sep = "")
        }
        NULL
      }
    )
    
    if (is.null(fold_fit)) {
      next
    }
    
    x_val      <- x[val_idx, , drop = FALSE]
    y_val_int  <- as.integer(y_factor[val_idx])
    w_val      <- weights[val_idx]
    w_val_sum  <- sum(w_val)
    
    if (w_val_sum <= 0) {
      next
    }
    
    for (l_idx in seq_along(lambda_seq)) {
      prob_mat <- tryCatch(
        multnet_predict_manual(
          fit          = fold_fit,
          x            = x_val,
          class_levels = class_levels,
          lambda       = lambda_seq[l_idx]
        ),
        error = function(e) {
          if (verbose) {
            cat(
              "[manual cv] predict failed on fold ", fold,
              ", lambda ", l_idx, ": ", e$message, "\n",
              sep = ""
            )
          }
          NULL
        }
      )
      
      if (is.null(prob_mat)) {
        next
      }
      
      idx <- cbind(seq_len(length(y_val_int)), y_val_int)
      probs_at_truth <- pmax(prob_mat[idx], 1e-12)
      
      cv_loss[fold, l_idx] <-
        -sum(w_val * log(probs_at_truth)) / w_val_sum
    }
  }
  
  mean_loss <- colMeans(cv_loss, na.rm = TRUE)
  se_loss <- apply(
    cv_loss,
    2,
    function(col) {
      col <- col[!is.na(col)]
      if (length(col) < 2) {
        return(NA_real_)
      }
      sd(col) / sqrt(length(col))
    }
  )
  
  if (all(is.na(mean_loss))) {
    stop("cv_glmnet_multiclass_manual: every CV fold failed; cannot pick lambda.")
  }
  
  best_idx <- which.min(mean_loss)
  lambda_min <- lambda_seq[best_idx]
  
  threshold <- mean_loss[best_idx] + ifelse(
    is.na(se_loss[best_idx]),
    0,
    se_loss[best_idx]
  )
  
  candidates <- which(mean_loss <= threshold)
  
  if (length(candidates) == 0) {
    lambda_1se <- lambda_min
  } else {
    # lambda_seq is decreasing. The smallest index is the strongest
    # regularization within one standard error.
    lambda_1se <- lambda_seq[min(candidates)]
  }
  
  list(
    lambda     = lambda_seq,
    cvm        = mean_loss,
    cvsd       = se_loss,
    lambda.min = lambda_min,
    lambda.1se = lambda_1se
  )
}

############################################################
# glmnet binary backend
############################################################

fit_glmnet_binary <- function(x, y_soft, backend_config = list()) {
  load_pkg("glmnet")
  
  glmnet_alpha <- get_config_value(backend_config, "glmnet_alpha", default = 0)
  lambda <- get_config_value(backend_config, "lambda", default = NULL)
  nfolds <- get_config_value(backend_config, "nfolds", default = 5)
  type_measure <- get_config_value(backend_config, "type_measure", default = "deviance")
  lambda_choice <- get_config_value(backend_config, "lambda_choice", default = "lambda.min")
  standardize <- get_config_value(backend_config, "standardize", default = TRUE)
  maxit <- get_config_value(backend_config, "maxit", default = 1000000)
  
  y_soft <- pmin(pmax(y_soft, 1e-6), 1 - 1e-6)
  y_mat <- cbind(class0 = 1 - y_soft, class1 = y_soft)
  
  if (is.null(lambda)) {
    cv_fit <- glmnet::cv.glmnet(
      x = x,
      y = y_mat,
      family = "binomial",
      alpha = glmnet_alpha,
      nfolds = nfolds,
      type.measure = type_measure,
      standardize = standardize,
      maxit = maxit
    )
    
    if (lambda_choice == "lambda.1se") {
      lambda_use <- cv_fit$lambda.1se
    } else {
      lambda_use <- cv_fit$lambda.min
    }
    
    fit <- glmnet::glmnet(
      x = x,
      y = y_mat,
      family = "binomial",
      alpha = glmnet_alpha,
      lambda = lambda_use,
      standardize = standardize,
      maxit = maxit
    )
    
    return(
      list(
        backend = "glmnet",
        task = "binary",
        fit = fit,
        cv_fit = cv_fit,
        lambda = lambda_use,
        glmnet_alpha = glmnet_alpha,
        backend_config = backend_config
      )
    )
    
  } else {
    fit <- glmnet::glmnet(
      x = x,
      y = y_mat,
      family = "binomial",
      alpha = glmnet_alpha,
      lambda = lambda,
      standardize = standardize,
      maxit = maxit
    )
    
    return(
      list(
        backend = "glmnet",
        task = "binary",
        fit = fit,
        cv_fit = NULL,
        lambda = lambda,
        glmnet_alpha = glmnet_alpha,
        backend_config = backend_config
      )
    )
  }
}

predict_glmnet_binary <- function(fit_obj, x) {
  pred <- predict(
    fit_obj$fit,
    newx = x,
    type = "response",
    s = fit_obj$lambda
  )
  
  if (length(dim(pred)) == 3) {
    prob <- as.numeric(pred[, "class1", 1])
  } else {
    prob <- as.numeric(pred)
  }
  
  prob <- pmin(pmax(prob, 1e-8), 1 - 1e-8)
  
  return(prob)
}

############################################################
# glmnet multiclass backend
############################################################

fit_glmnet_multiclass <- function(
  x,
  q_mat,
  class_levels,
  backend_config = list()
) {
  load_pkg("glmnet")
  
  glmnet_alpha <- get_config_value(backend_config, "glmnet_alpha", default = 0)
  lambda <- get_config_value(backend_config, "lambda", default = NULL)
  nfolds <- get_config_value(backend_config, "nfolds", default = 5)
  lambda_choice <- get_config_value(backend_config, "lambda_choice", default = "lambda.min")
  standardize <- get_config_value(backend_config, "standardize", default = TRUE)
  maxit <- get_config_value(backend_config, "maxit", default = 1000000)
  manual_cv_seed <- get_config_value(backend_config, "manual_cv_seed", default = 2026)
  manual_cv_verbose <- get_config_value(backend_config, "manual_cv_verbose", default = FALSE)
  
  weighted_data <- prepare_multiclass_weighted_data(
    x = x,
    q_mat = q_mat,
    class_levels = class_levels
  )
  
  y_factor <- factor(
    weighted_data$y,
    levels = class_levels
  )
  
  validate_multiclass_training_classes(
    y_factor = y_factor,
    class_levels = class_levels,
    context = "fit_glmnet_multiclass weighted data"
  )
  
  ############################################################
  # Case 1: no lambda supplied.
  # Build a full lambda path, manually CV it, and final-fit on
  # the full path. Prediction will manually select lambda_use.
  ############################################################
  
  if (is.null(lambda)) {
    full_path <- glmnet::glmnet(
      x = weighted_data$x,
      y = y_factor,
      weights = weighted_data$weight,
      family = "multinomial",
      alpha = glmnet_alpha,
      standardize = standardize,
      maxit = maxit
    )
    
    lambda_seq <- full_path$lambda
    lambda_seq <- as.numeric(lambda_seq)
    lambda_seq <- lambda_seq[is.finite(lambda_seq) & lambda_seq > 0]
    lambda_seq <- sort(unique(lambda_seq), decreasing = TRUE)
    
    if (length(lambda_seq) < 2) {
      lambda_seq <- make_safe_lambda_path(lambda_seq[1])
    }
    
    cv_fit <- cv_glmnet_multiclass_manual(
      x            = weighted_data$x,
      y_factor     = y_factor,
      weights      = weighted_data$weight,
      class_levels = class_levels,
      alpha        = glmnet_alpha,
      nfolds       = nfolds,
      standardize  = standardize,
      maxit        = maxit,
      seed         = manual_cv_seed,
      lambda_seq   = lambda_seq,
      verbose      = manual_cv_verbose
    )
    
    if (lambda_choice == "lambda.1se") {
      lambda_use <- cv_fit$lambda.1se
    } else {
      lambda_use <- cv_fit$lambda.min
    }
    
    fit <- glmnet::glmnet(
      x = weighted_data$x,
      y = y_factor,
      weights = weighted_data$weight,
      family = "multinomial",
      alpha = glmnet_alpha,
      lambda = lambda_seq,
      standardize = standardize,
      maxit = maxit
    )
    
    return(
      list(
        backend = "glmnet",
        task = "multiclass",
        fit = fit,
        cv_fit = cv_fit,
        lambda = lambda_use,
        lambda_seq = lambda_seq,
        glmnet_alpha = glmnet_alpha,
        class_levels = class_levels,
        backend_config = backend_config
      )
    )
  }
  
  ############################################################
  # Case 2: lambda supplied, usually reused from naive fit.
  # Still do NOT final-fit with scalar lambda. Prefer the full
  # reference lambda path; otherwise build a wide warm-start path.
  ############################################################
  
  lambda_use <- as.numeric(lambda)[1]
  
  if (!is.finite(lambda_use) || lambda_use <= 0) {
    stop("fit_glmnet_multiclass: supplied lambda must be positive and finite.")
  }
  
  lambda_seq_config <- get_config_value(
    backend_config,
    "lambda_seq",
    default = NULL
  )
  
  if (!is.null(lambda_seq_config)) {
    lambda_seq <- as.numeric(lambda_seq_config)
    lambda_seq <- lambda_seq[is.finite(lambda_seq) & lambda_seq > 0]
    lambda_seq <- sort(unique(lambda_seq), decreasing = TRUE)
    
    if (length(lambda_seq) < 10) {
      lambda_seq <- make_wide_lambda_path(lambda_use)
    }
  } else {
    lambda_seq <- make_wide_lambda_path(lambda_use)
  }
  
  fit <- glmnet::glmnet(
    x = weighted_data$x,
    y = y_factor,
    weights = weighted_data$weight,
    family = "multinomial",
    alpha = glmnet_alpha,
    lambda = lambda_seq,
    standardize = standardize,
    maxit = maxit
  )
  
  return(
    list(
      backend = "glmnet",
      task = "multiclass",
      fit = fit,
      cv_fit = NULL,
      lambda = lambda_use,
      lambda_seq = lambda_seq,
      glmnet_alpha = glmnet_alpha,
      class_levels = class_levels,
      backend_config = backend_config
    )
  )
}

predict_glmnet_multiclass <- function(fit_obj, x) {
  multnet_predict_manual(
    fit          = fit_obj$fit,
    x            = x,
    class_levels = fit_obj$class_levels,
    lambda       = fit_obj$lambda
  )
}

############################################################
# XGBoost binary backend
############################################################

fit_xgboost_binary <- function(x, y_soft, backend_config = list()) {
  load_pkg("xgboost")
  
  weighted_data <- prepare_binary_weighted_data(
    x = x,
    q = y_soft
  )
  
  dtrain <- xgboost::xgb.DMatrix(
    data = weighted_data$x,
    label = weighted_data$y,
    weight = weighted_data$weight
  )
  
  eta <- get_config_value(backend_config, "eta", default = 0.05)
  max_depth <- get_config_value(backend_config, "max_depth", default = 3)
  min_child_weight <- get_config_value(backend_config, "min_child_weight", default = 1)
  subsample <- get_config_value(backend_config, "subsample", default = 0.8)
  colsample_bytree <- get_config_value(backend_config, "colsample_bytree", default = 0.8)
  lambda <- get_config_value(backend_config, "xgb_lambda", default = 1)
  alpha <- get_config_value(backend_config, "xgb_alpha", default = 0)
  nthread <- get_config_value(backend_config, "nthread", default = 2)
  nrounds <- get_config_value(backend_config, "nrounds", default = 100)
  verbose <- get_config_value(backend_config, "verbose", default = 0)
  
  params <- list(
    objective = "binary:logistic",
    eval_metric = "logloss",
    eta = eta,
    max_depth = max_depth,
    min_child_weight = min_child_weight,
    subsample = subsample,
    colsample_bytree = colsample_bytree,
    lambda = lambda,
    alpha = alpha,
    nthread = nthread
  )
  
  fit <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = nrounds,
    verbose = verbose
  )
  
  return(
    list(
      backend = "xgboost",
      task = "binary",
      fit = fit,
      params = params,
      nrounds = nrounds,
      backend_config = backend_config
    )
  )
}

predict_xgboost_binary <- function(fit_obj, x) {
  dtest <- xgboost::xgb.DMatrix(data = x)
  
  prob <- predict(
    fit_obj$fit,
    newdata = dtest
  )
  
  prob <- as.numeric(prob)
  prob <- pmin(pmax(prob, 1e-8), 1 - 1e-8)
  
  return(prob)
}

############################################################
# XGBoost multiclass backend
############################################################

fit_xgboost_multiclass <- function(
  x,
  q_mat,
  class_levels,
  backend_config = list()
) {
  load_pkg("xgboost")
  
  weighted_data <- prepare_multiclass_weighted_data(
    x = x,
    q_mat = q_mat,
    class_levels = class_levels
  )
  
  y_integer <- as.integer(
    factor(weighted_data$y, levels = class_levels)
  ) - 1
  
  dtrain <- xgboost::xgb.DMatrix(
    data = weighted_data$x,
    label = y_integer,
    weight = weighted_data$weight
  )
  
  eta <- get_config_value(backend_config, "eta", default = 0.05)
  max_depth <- get_config_value(backend_config, "max_depth", default = 3)
  min_child_weight <- get_config_value(backend_config, "min_child_weight", default = 1)
  subsample <- get_config_value(backend_config, "subsample", default = 0.8)
  colsample_bytree <- get_config_value(backend_config, "colsample_bytree", default = 0.8)
  lambda <- get_config_value(backend_config, "xgb_lambda", default = 1)
  alpha <- get_config_value(backend_config, "xgb_alpha", default = 0)
  nthread <- get_config_value(backend_config, "nthread", default = 2)
  nrounds <- get_config_value(backend_config, "nrounds", default = 100)
  verbose <- get_config_value(backend_config, "verbose", default = 0)
  
  params <- list(
    objective = "multi:softprob",
    eval_metric = "mlogloss",
    num_class = length(class_levels),
    eta = eta,
    max_depth = max_depth,
    min_child_weight = min_child_weight,
    subsample = subsample,
    colsample_bytree = colsample_bytree,
    lambda = lambda,
    alpha = alpha,
    nthread = nthread
  )
  
  fit <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = nrounds,
    verbose = verbose
  )
  
  return(
    list(
      backend = "xgboost",
      task = "multiclass",
      fit = fit,
      params = params,
      nrounds = nrounds,
      class_levels = class_levels,
      backend_config = backend_config
    )
  )
}

predict_xgboost_multiclass <- function(fit_obj, x) {
  dtest <- xgboost::xgb.DMatrix(data = x)
  
  prob <- predict(
    fit_obj$fit,
    newdata = dtest
  )
  
  class_levels <- fit_obj$class_levels
  K <- length(class_levels)
  
  prob_mat <- matrix(
    prob,
    ncol = K,
    byrow = TRUE
  )
  
  colnames(prob_mat) <- class_levels
  prob_mat <- row_normalize(prob_mat)
  
  return(prob_mat)
}

############################################################
# Unified binary backend interface
############################################################

fit_backend_binary <- function(
  x,
  y_soft,
  backend_config = list()
) {
  backend <- get_config_value(backend_config, "backend", default = "glmnet")
  
  if (backend == "glmnet") {
    return(
      fit_glmnet_binary(
        x = x,
        y_soft = y_soft,
        backend_config = backend_config
      )
    )
  }
  
  if (backend == "xgboost") {
    return(
      fit_xgboost_binary(
        x = x,
        y_soft = y_soft,
        backend_config = backend_config
      )
    )
  }
  
  stop(
    paste0(
      "Unsupported backend: ",
      backend,
      ". Supported backends are 'glmnet' and 'xgboost'."
    )
  )
}

predict_backend_binary <- function(fit_obj, x) {
  if (fit_obj$backend == "glmnet") {
    return(predict_glmnet_binary(fit_obj = fit_obj, x = x))
  }
  
  if (fit_obj$backend == "xgboost") {
    return(predict_xgboost_binary(fit_obj = fit_obj, x = x))
  }
  
  stop(paste0("Unsupported fitted backend: ", fit_obj$backend))
}

############################################################
# Unified multiclass backend interface
############################################################

fit_backend_multiclass <- function(
  x,
  q_mat,
  class_levels,
  backend_config = list()
) {
  backend <- get_config_value(backend_config, "backend", default = "glmnet")
  
  if (backend == "glmnet") {
    return(
      fit_glmnet_multiclass(
        x = x,
        q_mat = q_mat,
        class_levels = class_levels,
        backend_config = backend_config
      )
    )
  }
  
  if (backend == "xgboost") {
    return(
      fit_xgboost_multiclass(
        x = x,
        q_mat = q_mat,
        class_levels = class_levels,
        backend_config = backend_config
      )
    )
  }
  
  stop(
    paste0(
      "Unsupported backend: ",
      backend,
      ". Supported backends are 'glmnet' and 'xgboost'."
    )
  )
}

predict_backend_multiclass <- function(fit_obj, x) {
  if (fit_obj$backend == "glmnet") {
    return(predict_glmnet_multiclass(fit_obj = fit_obj, x = x))
  }
  
  if (fit_obj$backend == "xgboost") {
    return(predict_xgboost_multiclass(fit_obj = fit_obj, x = x))
  }
  
  stop(paste0("Unsupported fitted backend: ", fit_obj$backend))
}

############################################################
# Backend config update helper
############################################################

reuse_glmnet_lambda_if_requested <- function(
  backend_config,
  reference_fit
) {
  reuse_lambda <- get_config_value(
    backend_config,
    "reuse_lambda",
    default = TRUE
  )
  
  if (
    reuse_lambda &&
    !is.null(reference_fit$backend) &&
    reference_fit$backend == "glmnet"
  ) {
    backend_config$lambda <- reference_fit$lambda
    
    # For multinomial glmnet, reuse the full lambda path from the
    # reference fit when available. Reusing only a scalar lambda forces
    # a tiny path in later EM/MCEM fits and can trigger glmnet's internal
    # dimnames error.
    if (!is.null(reference_fit$lambda_seq)) {
      backend_config$lambda_seq <- reference_fit$lambda_seq
    }
  }
  
  return(backend_config)
}
