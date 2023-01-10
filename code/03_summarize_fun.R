Summarise <- function(condition, results, fixed_objects = NULL) {
  AbsBias <- mean(results$AbsBias) * 100
  ret <- c(AbsBias = AbsBias)
  ret
}
