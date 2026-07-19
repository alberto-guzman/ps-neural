# Summarise function
# 2026-07 revision: emp_SE is the Monte Carlo SD of the ATE estimates; SE_ratio
# and SE_ratio_lm compare the average estimated SE (sandwich and weighted-lm,
# respectively) against it — values below 1 mean the estimator understates the
# true sampling variability. Trimmed-weight (*_trim) summaries mirror the primary ones.
Summarise <- function(condition, results, fixed_objects = NULL) {
  Std_In_Bias <- mean(results$Std_In_Bias)
  Prob_Treat <- mean(results$Prob_Treat)
  Bias <- bias(results$ATE, parameter = 0.3, type = "bias")
  Rel_Bias <- bias(results$ATE, parameter = 0.3, type = "relative", percent = F)
  Per_Rel_Bias <- bias(results$ATE, parameter = 0.3, type = "relative", percent = T)
  ATE <- mean(results$ATE)
  ATE_se <- mean(results$ATE_se)
  ATE_se_lm <- mean(results$ATE_se_lm)
  emp_SE <- sd(results$ATE)
  SE_ratio <- ATE_se / emp_SE
  SE_ratio_lm <- ATE_se_lm / emp_SE
  MSE <- RMSE(results$ATE, parameter = 0.3, MSE = T)
  ci_95 <- mean(results$ci_95)
  Power <- EDR(results$p_val, alpha = 0.05)
  mean_ps_weights <- mean(results$mean_ps_weights)
  max_ps_weight <- max(results$max_ps_weight)
  ess <- mean(results$ess)
  ASAM <- mean(results$ASAM)
  max_ASAM <- mean(results$max_ASAM)

  # trimmed-weight sensitivity
  Bias_trim <- bias(results$ATE_trim, parameter = 0.3, type = "bias")
  ATE_trim <- mean(results$ATE_trim)
  ATE_se_trim <- mean(results$ATE_se_trim)
  emp_SE_trim <- sd(results$ATE_trim)
  ci_95_trim <- mean(results$ci_95_trim)

  # Create a vector of the results
  ret <- c(
    Std_In_Bias = Std_In_Bias,
    Prob_Treat = Prob_Treat,
    Bias = Bias,
    Rel_Bias = Rel_Bias,
    Per_Rel_Bias = Per_Rel_Bias,
    ATE = ATE,
    ATE_se = ATE_se,
    ATE_se_lm = ATE_se_lm,
    emp_SE = emp_SE,
    SE_ratio = SE_ratio,
    SE_ratio_lm = SE_ratio_lm,
    MSE = MSE,
    ci_95 = ci_95,
    Power = Power,
    mean_ps_weights = mean_ps_weights,
    max_ps_weight = max_ps_weight,
    ess = ess,
    ASAM = ASAM,
    max_ASAM = max_ASAM,
    Bias_trim = Bias_trim,
    ATE_trim = ATE_trim,
    ATE_se_trim = ATE_se_trim,
    emp_SE_trim = emp_SE_trim,
    ci_95_trim = ci_95_trim
  )
  # Return the vector
  ret
}
