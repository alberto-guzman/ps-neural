# Summarise function
Summarise <- function(condition, results, fixed_objects = NULL) {
  Std_In_Bias <- mean(results$Std_In_Bias)
  Prob_Treat <- mean(results$Prob_Treat)
  Bias <- bias(results$ATE, parameter = 0.3, type = "bias")
  Abs_Per_Bias <- bias(results$ATE, parameter = 0.3, type = "bias", abs = T, percent = T)
  Abs_Per_Rel_Bias <- bias(results$ATE, parameter = 0.3, type = "relative", abs = T, percent = T)
  ATE_se <- mean(results$ATE_se)
  RMSE <- RMSE(results$ATE, parameter = 0.3)
  Power <- EDR(results$p_val, alpha = 0.05)
  coverage_95 <- mean(results$ci_95)
  mean_ps_weights <- mean(results$mean_ps_weights)
  ASAM <- mean(results$ASAM)
  # Create a vector of the results
  ret <- c(
    Std_In_Bias = Std_In_Bias,
    Prob_Treat = Prob_Treat,
    Bias = Bias,
    Abs_Per_Bias = Abs_Per_Bias,
    Abs_Per_Rel_Bias = Abs_Per_Rel_Bias,
    ATE_se = ATE_se,
    RMSE = RMSE,
    Power = Power,
    coverage_95 = coverage_95,
    mean_ps_weights = mean_ps_weights,
    ASAM = ASAM
  )
  # Return the vector
  ret
}
