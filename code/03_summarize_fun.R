# Summarise function
Summarise <- function(condition, results, fixed_objects = NULL) {
  # Calculate standardized initial bias
  Std_In_Bias <- mean(results$Std_In_Bias)
  # Calculate probability of treatment
  Prob_Treat <- mean(results$Prob_Treat)
  # Calculate MSE
  MSE <- RMSE(results$ATE, parameter = 0.3, MSE = T)
  # Calculate bias
  Bias <- mean(results$Bias)
  # Calculate relative bis
  RelBias <- bias(results$ATE, parameter = 0.3, type = "relative")
  # Standard error
  ATE <- mean(results$ATE)
  # Standard error
  ATE_se <- mean(results$ATE_se)
  # Calculate power
  Power <- EDR(results$p_val, alpha = 0.05)
  # Coverage
  coverage_95 <- mean(results$ci_95)
  # Mean of PS weights for control group
  ps_weight_control <- mean(results$mean_ps_weights)
  # Mean ASAM
  ASAM <- mean(results$ASAM)
  # Create a vector of the results
  ret <- c(
    Std_In_Bias = Std_In_Bias,
    Prob_Treat = Prob_Treat,
    ATE = ATE_se,
    ATE_se = ATE_se,
    Bias = Bias,
    RelBias = RelBias,
    ps_weight_control = ps_weight_control,
    ASAM = ASAM,
    Power = Power,
    coverage_95 = coverage_95
  )
  # Return the vector
  ret
}
