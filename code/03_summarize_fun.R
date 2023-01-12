# Summarise function
Summarise <- function(condition, results, fixed_objects = NULL) {
  # Calculate standardized initial bias
  Std_In_Bias <- mean(results$Std_In_Bias)
  # Calculate probability of treatment
  Prob_Treat <- mean(results$Prob_Treat)
  # Calculate the mean of the absolute bias
  AbsBias <- mean(results$AbsBias) * 100
  # Calculate abs bias from the package function
  AbsBias_fun <- bias(estimate = results$ATT, parameter = 0.3, type = "relative")
  # Standard error
  mean_ATT_se <- mean(results$ATT_se)
  # Calculate type 1 and power, detection rate
  #EDR <- mean(EDR(results$p_val, alpha = 0.05))
  # Coverage
  #ECR <- ECR(CIs = results$conf_int, parameter = 0.3)
  # Mean of PS weights for control group
  ps_weight_control <- mean(results$mean_ps_weights)
  # Mean ASAM
  ASAM <- mean(results$ASAM)
  # Create a vector of the results
  ret <- c(Std_In_Bias = Std_In_Bias, Prob_Treat = Prob_Treat, AbsBias = AbsBias, AbsBias_fun = AbsBias_fun, mean_ATT_se = mean_ATT_se, ps_weight_control = ps_weight_control, ASAM = ASAM)
  # Return the vector
  ret
}



