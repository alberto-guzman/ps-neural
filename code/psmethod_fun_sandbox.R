#############

## WHAT DOES THIS FUNCTIOND DO?

#############



ps_methods <- function(data, method) {
  df <- data

  if (method == "logit") {
    mod <- glm(T ~ . - Y - trueps, data = df, family = binomial)
    ps <- mod$fitted
  } else if (method == "cart") {
    mod <- rpart(T ~ . - Y - trueps, method = "class", data = df)
    ps <- predict(mod)[, 2]
  } else if (method == "bag") {
    mod <- bagging(T ~ . - Y - trueps, data = df, nbag = 100)
    ps <- predict(mod, newdata = df, type = "prob")
  } else if (method == "forest") {
    mod <- randomForest(factor(T) ~ . - Y - trueps, ntree = 500, data = df)
    ps <- predict(mod, type = "prob")[, 2]
  }


  # Save estimated propensity score and weights to data frame
  df <- df %>%
    mutate(
      ps_pred = ps,
      ps_weights = if_else(T == 0, ps / (1 - ps), 1)
    )

  # estimates the ATT with the weights
  d.w <- svydesign(~1, weights = df$ps_weights, data = df)
  fit <- svyglm(Y ~ T, design = d.w)

  # saves the ATT and se_ATT
  ATT <- unname(coef(fit)["T"])
  vcov_matrix <- vcov(fit)
  ATT_se <- unname(sqrt(vcov_matrix["T", "T"]))

  # calculate bias metrics
  Bias <- ATT - 0.3
  AbsBias <- abs(ATT - 0.3)

  # calculate mean of control group weights
  df_int <- subset(df, T == 0)
  mean_ps_weights <- mean(df_int$ps_weights, na.rm = TRUE)

  # calculate 95% coverage
  lower_bound <- ATT - 1.96 * ATT_se
  upper_bound <- ATT + 1.96 * ATT_se
  ci_95 <- ifelse(lower_bound <= 0.3 && 0.3 <= upper_bound, 1, 0)

  ###############
  # calculate ASAM for covariates
  ###############

  # Subset the data into the treatment and comparison groups
  treatment_group <- df[df$T == 1, ]
  comparison_group <- df[df$T == 0, ]

  # Get the names of the variables that start with "v"
  var_names <- names(df)[grep("^v", names(df))]

  # Initialize the ASAM_list vector
  ASAM_list <- rep(NA, length(var_names))

  # Loop through each covariate
  for (i in 1:length(var_names)) {
    # Get the covariate name
    covariate <- var_names[i]

    # Extract the covariate data from the treatment and comparison groups
    treatment_data <- treatment_group[[covariate]]
    comparison_data <- comparison_group[[covariate]]

    # Extract the weights from the treatment and comparison groups
    treatment_weights <- treatment_group$ps_weights
    comparison_weights <- comparison_group$ps_weights

    # Calculate the means of the treatment and comparison groups
    treatment_mean <- weighted.mean(treatment_data, treatment_weights)
    comparison_mean <- weighted.mean(comparison_data, comparison_weights)

    # Calculate the variances of the treatment and comparison groups
    treatment_var <- wtd.var(treatment_data, treatment_weights)
    comparison_var <- wtd.var(comparison_data, comparison_weights)

    # Calculate the standard deviations of the treatment and comparison groups
    treatment_sd <- sqrt(treatment_var)
    comparison_sd <- sqrt(comparison_var)

    # Calculate the standardized difference of means
    sd_diff <- (treatment_mean - comparison_mean) / treatment_sd

    # Take the absolute value of the standardized difference of means
    abs_sd_diff <- abs(sd_diff)

    # Save the absolute standardized difference of means in the ASAM_list vector
    ASAM_list[i] <- abs_sd_diff
  }


  # Calculate the mean of the absolute standardized differences of means
  ASAM <- mean(ASAM_list)


  return(c(ATT, ATT_se, Bias, AbsBias, mean_ps_weights, ci_95, ASAM))
}




























##############
# PARKING LOT
##############
#
# ggplot(df, aes(x = trueps, y = ps_pred)) +
#   geom_point(shape = 21, alpha = 0.2) +
#   geom_abline(slope = 1, intercept = 0) +
#   scale_x_continuous(limits = c(0, 1)) +
#   scale_y_continuous(limits = c(0, 1)) +
#   labs(x = "True PS", y = "PS predicted by main effects logistic regression")
#
#
# # Estimate propensity score nn
# neuro_n <- ceiling((2 / 3) * length(df))
# samp <- sample(1:nrow(df), ceiling(.70 * nrow(df)))
# mod <- nnet(factor(T) ~ . - Y - trueps, data = df, size = neuro_n, decay = 0.01, maxit = 200, trace = F, subset = samp)
# ps <- as.numeric(predict(mod, type = "raw"))
