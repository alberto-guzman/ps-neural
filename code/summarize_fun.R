Summarise <- function(condition, results, fixed_objects = NULL) {
  AbsBias <- mean(results$AbsBias)*100
  ret <- c(AbsBias=AbsBias)
  ret
}

#ret <- c(ATT=ATT, ATT_se=ATT_se,
 #        Bias=Bias, AbsBias=AbsBias, mean_ps_weights=mean_ps_weights, ci_95=ci_95, ASAM=ASAM)