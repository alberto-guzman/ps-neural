#############
## WHAT DOES THIS FUNCTION DO?
# The Summarise function (SimDesign's third stage) collapses the per-replication
# results produced by Analyse() into one row of performance metrics per design
# cell. It receives `results` — a data frame with one row per replication and
# one column per quantity returned by Analyse() — and returns a named vector.
#
# Column glossary (true ATE = 0.3 throughout; see the manuscript's Performance
# Metrics section for definitions and citations):
#   Std_In_Bias      mean standardized bias of the UNWEIGHTED treated-control
#                    outcome difference (descriptive; before any weighting)
#   Prob_Treat       mean realized treatment rate (should be ~.50 by design)
#   Bias             mean(ATE_hat) - 0.3                       [SimDesign::bias]
#   Rel_Bias         Bias / 0.3
#   Per_Rel_Bias     Rel_Bias as a percentage
#   ATE              mean point estimate
#   ATE_se           mean estimated SANDWICH standard error (svyglm; primary)
#   ATE_se_lm        mean model-based weighted-lm SE (comparison only — treats
#                    IPW weights as precision weights, a common misuse)
#   emp_SE           Monte Carlo SD of the point estimates — the "true"
#                    sampling variability under the fixed population
#   SE_ratio         ATE_se / emp_SE. Values < 1: the sandwich UNDERSTATES
#                    sampling variability (anti-conservative); > 1: conservative
#   SE_ratio_lm      same ratio for the weighted-lm SE
#   MSE              mean squared error about 0.3                [SimDesign::RMSE]
#   ci_95            95% CI coverage of 0.3 (CIs from the svyglm sandwich fit)
#   Power            share of replications rejecting H0: ATE = 0 at alpha=.05
#   mean_ps_weights  mean IPW weight (ATE weights average ~2 when healthy)
#   max_ps_weight    largest weight observed across replications (extreme-weight
#                    diagnostic; bounded above by 1e6 via the numerical guard)
#   ess              mean Kish effective sample size implied by the weights
#   ASAM / max_ASAM  mean / worst-covariate absolute standardized mean
#                    difference after weighting (unweighted-treated-SD
#                    denominator, so comparable across methods)
#   *_trim           the same estimate/SE/coverage pipeline after truncating
#                    weights at their 1st/99th percentiles (the Lee, Lessler &
#                    Stuart remedial sensitivity — NOT part of the primary
#                    practitioner workflow being evaluated)
#
# Monte Carlo standard errors for reader-facing tables can be derived from
# these: e.g., MCSE(coverage) = sqrt(ci_95 * (1 - ci_95) / R) with R = 1000.
#############

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

  # trimmed-weight sensitivity (1st/99th percentile truncation)
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
