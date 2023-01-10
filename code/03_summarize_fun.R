# Summarise function
Summarise <- function(condition, results, fixed_objects = NULL) {
  # Calculate the mean of the absolute bias
  AbsBias <- mean(results$AbsBias) * 100
  # Calculate the mean of the standard deviation of the bias
  Std_In_Bias <- mean(results$Std_In_Bias)
  # Calculate the mean of the probability of treatment
  Prob_Treat <- mean(results$Prob_Treat)
  # Create a vector of the results
  ret <- c(AbsBias = AbsBias, Std_In_Bias = Std_In_Bias, Prob_Treat = Prob_Treat)
  # Return the vector
  ret
}
