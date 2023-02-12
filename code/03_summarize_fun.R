# Summarise function
Summarise <- function(condition, results, fixed_objects = NULL) {
  Std_In_Bias <- mean(results$Std_In_Bias)
  Prob_Treat <- mean(results$Prob_Treat)
  Bias <- bias(results$ATE, parameter = 0.3, type = "bias")
  Rel_Bias <- bias(results$ATE, parameter = 0.3, type = "relative", percent = F)
  Per_Rel_Bias <- bias(results$ATE, parameter = 0.3, type = "relative", percent = T)
  ATE <- mean(results$ATE)
  ATE_se <- mean(results$ATE_se)
  MSE <- RMSE(results$ATE, parameter = 0.3, MSE = T)
  ci_95 <- mean(results$ci_95)
  Power <- EDR(results$p_val, alpha = 0.05)
  mean_ps_weights <- mean(results$mean_ps_weights)
  ASAM <- mean(results$ASAM)

  # Create a vector of the results
  ret <- c(
    Std_In_Bias = Std_In_Bias,
    Prob_Treat = Prob_Treat,
    Bias = Bias,
    Rel_Bias = Rel_Bias,
    Per_Rel_Bias = Per_Rel_Bias,
    ATE = ATE,
    ATE_se = ATE_se,
    MSE = MSE,
    ci_95 = ci_95,
    Power = Power,
    mean_ps_weights = mean_ps_weights,
    ASAM = ASAM
  )
  # Return the vector
  ret
}
