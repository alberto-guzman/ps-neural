Summarise <- function(condition, results, fixed_objects = NULL) {
  AbsBias <- mean(results$AbsBias) * 100
  Std_In_Bias <- mean(results$Std_In_Bias)
  Prob_Treat <- mean(results$Prob_Treat)
  ret <- c(AbsBias = AbsBias, Std_In_Bias = Std_In_Bias, Prob_Treat = Prob_Treat)
  ret
}
