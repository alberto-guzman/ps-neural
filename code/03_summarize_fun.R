# Summarise function
Summarise <- function(condition, results, fixed_objects = NULL) {
  Std_In_Bias <- mean(results$Std_In_Bias)
  Prob_Treat <- mean(results$Prob_Treat)
  Bias <- bias(results$ATE, parameter = 0.3, type = "bias")
  Per_Rel_Bias <- bias(results$ATE, parameter = 0.3, type = "relative", percent = T)
  ATE <- mean(results$ATE)
  ATE_se <- mean(results$ATE_se)
  MSE <- RMSE(results$ATE, parameter = 0.3, MSE = T)
  ECR <- ECR(CIs, parameter = 0.3)
  EDR <- 1-ECR(CIs, parameter = 0.3)
  Power <- EDR(results$p_val, alpha = 0.05)
  mean_ps_weights <- mean(results$mean_ps_weights)
  ASAM <- mean(results$ASAM)
  MSRSE <- MSRSE(SE = results$ATE_se, SD = results$ATE)
  RD_SE <- RD(results$ATE_se, sd(results$ATE))
  
  
  
  # Create a vector of the results
  ret <- c(
    Std_In_Bias = Std_In_Bias,
    Prob_Treat = Prob_Treat,
    Bias = Bias,
    Per_Rel_Bias = Per_Rel_Bias,
    ATE = ATE,
    ATE_se = ATE_se,
    MSE = MSE,
    ECR = ECR,
    EDR = EDR,
    Power = Power,
    MSRSE = MSRSE,
    mean_ps_weights = mean_ps_weights,
    ASAM = ASAM,
    RD_SE = RD_SE
  )
  # Return the vector
  ret
}
