#############

## WHAT DOES THIS FUNCTIOND DO?

#############

# function to estimate the ATT and other metrics
Analyse <- function(condition, dat, fixed_objects = NULL) {
  Attach(condition)
  # if the method is logit, then estimate the ATT using logistic regression
  if (method == "logit") {
    # estimate the propensity score using logistic regression
    mod <- glm(T ~ . - Y - trueps, data = dat, family = binomial)
    # save the propensity score to a vector
    ps <- mod$fitted
    # if the method is cart, then estimate the ATT using classification and regression trees
  } else if (method == "cart") {
    # estimate the propensity score using classification and regression trees
    mod <- rpart(T ~ . - Y - trueps, method = "class", data = dat)
    # save the propensity score to a vector
    ps <- predict(mod)[, 2]
    # if the method is bag, then estimate the ATT using bagging
  } else if (method == "bag") {
    # estimate the propensity score using bagging
    mod <- bagging(T ~ . - Y - trueps, data = dat, nbag = 100)
    # save the propensity score to a vector
    ps <- predict(mod, newdata = dat, type = "prob")
    # if the method is forest, then estimate the ATT using random forest
  } else if (method == "forest") {
    # estimate the propensity score using random forest
    mod <- randomForest(factor(T) ~ . - Y - trueps, ntree = 500, data = dat)
    # save the propensity score to a vector
    ps <- predict(mod, type = "prob")[, 2]
  }


  # save estimated propensity score and weights to data frame
  dat <- dat %>%
    mutate(
      ps_pred = ps,
      ps_weights = if_else(T == 0, ps / (1 - ps), 1)
    )

  # estimate the ATT with the weights
  d.w <- svydesign(~1, weights = dat$ps_weights, data = dat)
  fit <- svyglm(Y ~ T, design = d.w)

  # save the ATT and se_ATT
  ATT <- unname(coef(fit)["T"])
  vcov_matrix <- vcov(fit)
  ATT_se <- unname(sqrt(vcov_matrix["T", "T"]))

  # calculate the bias metrics
  Bias <- ATT - 0.3
  AbsBias <- abs(ATT - 0.3)

  # calculate the mean of control group weights
  dat_int <- subset(dat, T == 0)
  mean_ps_weights <- mean(dat_int$ps_weights, na.rm = TRUE)

  # calculate the 95% coverage
  lower_bound <- ATT - 1.96 * ATT_se
  upper_bound <- ATT + 1.96 * ATT_se
  ci_95 <- ifelse(lower_bound <= 0.3 && 0.3 <= upper_bound, 1, 0)

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
    # get the covariate name
    covariate <- var_names[i]

    # extract the covariate data from the treatment and comparison groups
    treatment_data <- treatment_group[[covariate]]
    comparison_data <- comparison_group[[covariate]]

    # extract the weights from the treatment and comparison groups
    treatment_weights <- treatment_group$ps_weights
    comparison_weights <- comparison_group$ps_weights

    # calculate the means of the treatment and comparison groups
    treatment_mean <- weighted.mean(treatment_data, treatment_weights)
    comparison_mean <- weighted.mean(comparison_data, comparison_weights)

    # calculate the variances of the treatment and comparison groups
    treatment_var <- wtd.var(treatment_data, treatment_weights)
    comparison_var <- wtd.var(comparison_data, comparison_weights)

    # calculate the standard deviations of the treatment and comparison groups
    treatment_sd <- sqrt(treatment_var)
    comparison_sd <- sqrt(comparison_var)

    # calculate the standardized difference of means
    sd_diff <- (treatment_mean - comparison_mean) / treatment_sd

    # take the absolute value of the standardized difference of means
    abs_sd_diff <- abs(sd_diff)

    # save the absolute standardized difference of means in the ASAM_list vector
    ASAM_list[i] <- abs_sd_diff
  }


  # calculate the mean of the absolute standardized differences of means
  ASAM <- mean(ASAM_list)
  
  ret <- c(ATT=ATT, ATT_se=ATT_se,
           Bias=Bias, AbsBias=AbsBias, mean_ps_weights=mean_ps_weights, ci_95=ci_95, ASAM=ASAM)
  ret
}



























##############
# PARKING LOT
##############
#
# ggplot(dat, aes(x = trueps, y = ps_pred)) +
#   geom_point(shape = 21, alpha = 0.2) +
#   geom_abline(slope = 1, intercept = 0) +
#   scale_x_continuous(limits = c(0, 1)) +
#   scale_y_continuous(limits = c(0, 1)) +
#   labs(x = "True PS", y = "PS predicted by main effects logistic regression")
#
#
# # Estimate propensity score nn
# neuro_n <- ceiling((2 / 3) * length(dat))
# samp <- sample(1:nrow(dat), ceiling(.70 * nrow(dat)))
# mod <- nnet(factor(T) ~ . - Y - trueps, data = dat, size = neuro_n, decay = 0.01, maxit = 200, trace = F, subset = samp)
# ps <- as.numeric(predict(mod, type = "raw"))
