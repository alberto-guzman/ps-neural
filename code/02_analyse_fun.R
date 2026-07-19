#############
## WHAT DOES THIS FUNCTION DO?
# The Analyse function is used to estimate the average treatment effect (ATE) and related metrics for a given condition.
# The function uses one of several methods, specified by the method argument, to estimate the propensity score.
# The methods used to estimate the propensity score are
# logistic regression (logit), classification and regression trees (cart), bagging (bag), random forest (forest),
# and three neural network models (nn-1, dnn-2, and dnn-3). Once the propensity score is estimated,
# the function uses survey-weighted regression to estimate the ATE, standard error of the ATE, p-value, and 95% confidence interval of the ATE.
# The function also calculates the absolute standardized average mean (ASAM) for each covariate in the data.
#
# 2026-07 revision notes (see git history for the original):
# - bag now fits CLASSIFICATION bagging (factor response); the original passed a
#   numeric 0/1 response, which silently fit regression bagging.
# - bag and forest now use out-of-bag predictions for the propensity score. The
#   original predicted in-sample (newdata = training data), which for random
#   forest returns overfit in-bag probabilities and collapses the weights.
#   CART and logit remain in-sample (no OOB analogue); NNs use a validation
#   split during training and predict in-sample, as before.
# - The three keras architectures share one helper (fit_nn_ps) instead of three
#   copy-pasted blocks. Inputs are standardized with training-set means/SDs, and
#   TF seeds are drawn from the SimDesign-managed RNG stream for reproducibility.
# - The primary ATE_se / p-value / CI now all come from the same svyglm sandwich
#   fit. The model-based weighted-lm SE is retained as ATE_se_lm for comparison.
#   A weight-truncation sensitivity (1st/99th percentiles) is saved as *_trim.
#############

# Shared keras fitting helper: trains a network with `n_hidden` hidden layers and
# returns propensity scores predicted for the full data set
fit_nn_ps <- function(dat, n_hidden) {
  # Release the previous replication's graph memory (36,000 model builds per job)
  try(keras::k_clear_session(), silent = TRUE)

  # Seed TF/Python RNGs from the SimDesign-managed R stream so runs are
  # reproducible. disable_gpu = FALSE: the default (TRUE) would silently turn
  # off the GPU this job's Slurm allocation requests; we accept the minor
  # cuDNN nondeterminism that GPU training implies.
  tensorflow::set_random_seed(sample.int(.Machine$integer.max, 1), disable_gpu = FALSE)

  # Split the data into training and validation sets (80/20)
  split <- sample(2, nrow(dat), replace = TRUE, prob = c(0.8, 0.2))
  train_data <- dat[split == 1, ]
  validation_data <- dat[split == 2, ]

  x_cols <- grep("^v", names(dat))
  x_train <- as.matrix(train_data[, x_cols])
  y_train <- as.matrix(train_data[, "T"])
  x_validation <- as.matrix(validation_data[, x_cols])
  y_validation <- as.matrix(validation_data[, "T"])

  # Standardize inputs using training-set means and SDs
  ctr <- colMeans(x_train)
  scl <- apply(x_train, 2, sd)
  scl[scl == 0] <- 1
  x_train <- scale(x_train, center = ctr, scale = scl)
  x_validation <- scale(x_validation, center = ctr, scale = scl)

  # Define model: `n_hidden` dense ReLU layers of ceiling(2p/3) units, L2-regularized
  p <- ncol(x_train)
  input_layer <- layer_input(shape = c(p))
  hidden <- input_layer
  for (i in seq_len(n_hidden)) {
    hidden <- layer_dense(
      units = ceiling(2 * p / 3), activation = "relu",
      kernel_regularizer = regularizer_l2(l = 0.01)
    )(hidden)
  }
  output_layer <- layer_dense(
    units = 1, activation = "sigmoid",
    kernel_regularizer = regularizer_l2(l = 0.01)
  )(hidden)
  model <- keras_model(inputs = input_layer, outputs = output_layer)

  model %>% compile(
    optimizer = "adam",
    loss = "binary_crossentropy",
    metrics = c("accuracy")
  )

  early_stopping <- callback_early_stopping(monitor = "val_loss", min_delta = 0.001, patience = 5)

  model %>% fit(
    x_train,
    y_train,
    epochs = 100,
    batch_size = 64,
    validation_data = list(x_validation, y_validation),
    callbacks = list(early_stopping),
    verbose = 0
  )

  # Predict propensity scores on the entire dataset (standardized the same way)
  x <- scale(as.matrix(dat[, x_cols]), center = ctr, scale = scl)
  ps <- model %>% predict(x)
  ps[, 1]
}

# function to estimate the ATE and other metrics
Analyse <- function(condition, dat, fixed_objects = NULL) {
  Attach(condition)

  if (method == "logit") {
    # estimate the propensity score using logistic regression (in-sample prediction)
    mod <- glm(T ~ . - Y - trueps, data = dat, family = binomial)
    ps <- predict(mod, newdata = dat, type = "response")
  } else if (method == "cart") {
    # estimate the propensity score using a classification tree (in-sample prediction;
    # a single tree has no out-of-bag analogue)
    mod <- rpart(T ~ . - Y - trueps, method = "class", data = dat)
    ps <- predict(mod, newdata = dat, type = "prob")[, 2]
  } else if (method == "bag") {
    # estimate the propensity score using classification bagging; omitting newdata
    # returns out-of-bag probabilities
    mod <- bagging(factor(T) ~ . - Y - trueps, data = dat, nbagg = 100)
    ps <- predict(mod, type = "prob")[, 2]
  } else if (method == "forest") {
    # estimate the propensity score using random forest; omitting newdata returns
    # out-of-bag vote proportions
    mod <- randomForest(factor(T) ~ . - Y - trueps, data = dat)
    ps <- predict(mod, type = "prob")[, 2]
  } else if (method == "gbm") {
    # estimate the propensity score using boosted trees exactly as MatchIt
    # delivers them for distance = "gbm" (replicated from MatchIt 4.7.2
    # source): bernoulli deviance, 10,000 trees, depth 3, shrinkage 0.01,
    # bag.fraction 1, iteration selected by 5-fold cross-validated deviance,
    # in-sample predictions. Same hyperparameters as WeightIt's default —
    # the packages differ only in the selection rule (CV deviance vs balance);
    # the CV rule keeps every learner in this study tuned for prediction.
    # n.cores = 1 pins gbm's internal CV parallelism (compute plumbing only —
    # SimDesign already parallelizes across replications; letting gbm detect
    # cores would oversubscribe the node)
    mod <- gbm(T ~ . - Y - trueps, data = dat, distribution = "bernoulli",
               n.trees = 10000, interaction.depth = 3, shrinkage = 0.01,
               bag.fraction = 1, cv.folds = 5, keep.data = FALSE, n.cores = 1)
    best_iter <- gbm.perf(mod, method = "cv", plot.it = FALSE)
    ps <- predict(mod, newdata = dat, n.trees = best_iter, type = "response")
  } else if (method == "bart") {
    # estimate the propensity score using Bayesian additive regression trees
    # exactly as an applied researcher gets them from WeightIt with
    # method = "bart" and all defaults (dbarts::bart2 backend, posterior-mean
    # probabilities)
    # n.threads = 1 pins dbarts' internal threading (compute plumbing only;
    # avoids oversubscription under SimDesign's worker pool)
    W <- WeightIt::weightit(
      reformulate(grep("^v", names(dat), value = TRUE), response = "T"),
      data = dat, method = "bart", estimand = "ATE", n.threads = 1L
    )
    ps <- as.numeric(W$ps)
  } else if (method == "sl") {
    # estimate the propensity score with a Super Learner stack: NNLS-weighted
    # combination of logit, lasso-logit, boosted trees (xgboost), random forest
    # (ranger), and a single-hidden-layer nnet, with 5-fold CV for the ensemble
    # weights; ranger/xgboost wrappers rather than randomForest/gbm for runtime
    x_df <- as.data.frame(dat[, grep("^v", names(dat))])
    sl <- SuperLearner(
      Y = dat$T, X = x_df, family = binomial(),
      SL.library = c("SL.glm", "SL.glmnet", "SL.xgboost", "SL.ranger", "SL.nnet"),
      cvControl = list(V = 5)
    )
    # combine the CROSS-VALIDATED library predictions (Z) with the NNLS weights
    # rather than SL.predict, which is in-sample and overfits like in-bag forest
    # predictions; this parallels the OOB choice for bag/forest
    ps <- as.numeric(sl$Z %*% sl$coef)
  } else if (method == "nn-1") {
    ps <- fit_nn_ps(dat, n_hidden = 1)
  } else if (method == "dnn-2") {
    ps <- fit_nn_ps(dat, n_hidden = 2)
  } else if (method == "dnn-3") {
    ps <- fit_nn_ps(dat, n_hidden = 3)
  }

  ##############################
  ### calculate metrics
  ##############################

  # numerical guard only: tree ensembles can emit probabilities of exactly 0 or 1
  # (undefined IPW weights); the bound is far outside any substantive range so
  # extreme-weight behavior remains visible in the diagnostics
  ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)

  dat <- dat %>%
    mutate(
      ps_pred = ps,
      ps_weights = case_when(T == 1 ~ 1 / ps, T == 0 ~ 1 / (1 - ps))
    )

  true_ATE <- 0.3

  # calculate standardized initial bias prior to weighting
  Std_In_Bias <- ((mean(dat$Y[dat$T == 1]) - mean(dat$Y[dat$T == 0])) - true_ATE) / sd(dat$Y[dat$T == 1])
  Prob_Treat <- mean(dat$T)

  # estimate the ATE with the weights; sandwich (robust) SEs from the survey package
  d.w <- svydesign(~0, weights = dat$ps_weights, data = dat)
  fit <- svyglm(Y ~ T, design = d.w)

  ATE <- unname(coef(fit)["T"])
  ATE_se <- unname(SE(fit)["T"])
  p_val <- summary(fit)$coefficients["T", "Pr(>|t|)"]

  # model-based SE from a weighted lm, retained for comparison with the sandwich SE
  modw <- lm(Y ~ T, data = dat, weights = dat$ps_weights)
  ATE_se_lm <- summary(modw)$coefficients["T", "Std. Error"]

  # calculate the 95% coverage
  conf_interval <- confint(fit, level = 0.95)["T", ]
  ci_95 <- ifelse(conf_interval[1] < true_ATE & true_ATE < conf_interval[2], 1, 0)

  # weight diagnostics
  mean_ps_weights <- mean(dat$ps_weights)
  max_ps_weight <- max(dat$ps_weights)
  ess <- sum(dat$ps_weights)^2 / sum(dat$ps_weights^2)

  # sensitivity: truncate weights at the 1st/99th percentiles and re-estimate
  qs <- quantile(dat$ps_weights, c(0.01, 0.99))
  w_trim <- pmin(pmax(dat$ps_weights, qs[1]), qs[2])
  d.w_trim <- svydesign(~0, weights = w_trim, data = dat)
  fit_trim <- svyglm(Y ~ T, design = d.w_trim)
  ATE_trim <- unname(coef(fit_trim)["T"])
  ATE_se_trim <- unname(SE(fit_trim)["T"])
  ci_trim <- confint(fit_trim, level = 0.95)["T", ]
  ci_95_trim <- ifelse(ci_trim[1] < true_ATE & true_ATE < ci_trim[2], 1, 0)

  ###############
  # calculate the ASAM for covariates
  ###############

  # subset the data into the treatment and comparison groups
  treatment_group <- dat[dat$T == 1, ]
  comparison_group <- dat[dat$T == 0, ]

  # get the names of the variables that start with "v"
  var_names <- names(dat)[grep("^v", names(dat))]

  # initialize the ASAM_list vector
  ASAM_list <- rep(NA, length(var_names))

  # loop through each covariate
  for (i in 1:length(var_names)) {
    covariate <- var_names[i]

    treatment_data <- treatment_group[[covariate]]
    comparison_data <- comparison_group[[covariate]]

    treatment_weights <- treatment_group$ps_weights
    comparison_weights <- comparison_group$ps_weights

    treatment_mean <- weighted.mean(treatment_data, treatment_weights)
    comparison_mean <- weighted.mean(comparison_data, comparison_weights)

    # standardize by the UNWEIGHTED treated-group SD so the denominator is
    # constant across methods (the cobalt convention); a weighted SD would let
    # each method's weights change its own yardstick
    treatment_sd <- sd(treatment_data)

    # absolute standardized difference of means
    ASAM_list[i] <- abs((treatment_mean - comparison_mean) / treatment_sd)
  }

  # mean and worst-covariate absolute standardized differences of means
  ASAM <- mean(ASAM_list)
  max_ASAM <- max(ASAM_list)

  ret <- c(
    Std_In_Bias = Std_In_Bias,
    Prob_Treat = Prob_Treat,
    ATE = ATE,
    ATE_se = ATE_se,
    ATE_se_lm = ATE_se_lm,
    ATE_trim = ATE_trim,
    ATE_se_trim = ATE_se_trim,
    ci_95_trim = unname(ci_95_trim),
    mean_ps_weights = mean_ps_weights,
    max_ps_weight = max_ps_weight,
    ess = ess,
    ASAM = ASAM,
    max_ASAM = max_ASAM,
    p_val = p_val,
    ci_95 = unname(ci_95)
  )
  ret
}
