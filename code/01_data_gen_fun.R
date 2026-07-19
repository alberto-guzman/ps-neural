#############
## WHAT DOES THIS FUNCTION DO?
# The Generate function generates simulated data based on specified conditions.
# The input is a tibble containing information about the sample size, number of covariates, and conditions for the population treatment and outcome models.
# The function first generates correlated normal variables and transforms them into normal, uniform, and Bernoulli variables (in that column order). It then selects a subset of the covariates for use in the population treatment and outcome models.
# The population treatment and outcome models are generated based on the specified conditions
# and the selected covariates, and the treatment status is also generated.
# The function returns a simulated data tibble that includes the original covariates, treatment status, and generated outcome.
#
# 2026-07 revision notes (see git history for the original):
# - The correlation matrix is now built symmetric BEFORE smoothing. The original
#   filled all p^2 cells with runif draws; psych::cor.smooth did not recognize the
#   asymmetric result as a correlation matrix and silently treated it as raw data,
#   so the correlations actually used did not match the intended U(-0.3, 0.3).
# - Complex scenarios now follow the manuscript equations (eq-psD / eq-outcome_d):
#   J quadratic terms plus J-1 adjacent-pair interactions X_j * X_{j+1} within the
#   sampled set, each reusing the variable's main-effect coefficient. The original
#   generated all pairwise combinations (2,775 interactions at p = 200), which
#   saturated the true propensity score.
# - Linear predictors are computed with matrix algebra instead of eval(parse()).
# - The POPULATION is now fixed within each design cell: the correlation matrix,
#   covariate roles, coefficients, complex-term selections, and calibrated
#   intercept are drawn under a seed derived from (n, p, scenarioT, scenarioY) —
#   deliberately NOT method, so every method faces the same population — and the
#   replication's RNG state is restored before the data draws. The original
#   redrew the population every replication, which conflates sampling variance
#   with between-population bias variation and makes CI coverage uninterpretable
#   (cf. Cannas & Arpino 2019 / Setoguchi 2008, whose DGP coefficients are fixed).
# - The treatment-model intercept is CALIBRATED per population (uniroot on a
#   seeded calibration sample of 20,000) so the marginal treatment rate is 0.50
#   in every design cell; the original fixed b0 = 0.25, which produced rates
#   from 0.22 to 0.65 across the fixed populations. Calibration centers the
#   true propensity distribution but does not change its spread — overlap in
#   the high-dimensional cells remains a (reported) design feature.
#############

Generate <- function(condition, fixed_objects = NULL) {
  # Makes the tibble of sim crossed conditions accessible to the function, from SimDesign package
  Attach(condition)

  #########################################
  # Population parameters (fixed within a design cell)
  #########################################

  # Deterministic population seed from the design cell; restore the replication
  # stream afterwards so data draws still vary across replications. The RNG
  # kind is pinned because parallel workers (L'Ecuyer-CMRG) and serial runs
  # (Mersenne-Twister) would otherwise derive DIFFERENT populations from the
  # same seed — the P and NP jobs must simulate identical populations.
  if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)
  rep_stream <- .Random.seed
  # NOTE: the additive seed formula is collision-free for this study's design
  # (n fixed at 10,000); revisit if n ever becomes a crossed factor
  set.seed(1e6 + n + 100 * p + 10 * (scenarioT == "complex_T") + (scenarioY == "complex_Y"),
           kind = "Mersenne-Twister")

  # Generate a symmetric correlation matrix with off-diagonals between -.3 and .3
  cor <- matrix(0, nrow = p, ncol = p)
  cor[upper.tri(cor)] <- runif(p * (p - 1) / 2, min = -.3, max = .3)
  cor <- cor + t(cor)
  diag(cor) <- 1

  # Smooth the correlation matrix to ensure it is positive definite
  cor <- psych::cor.smooth(cor)

  # Covariate roles: half confounders, a quarter outcome-only, a quarter treatment-only
  perm <- sample(p)
  n_confound <- floor(p / 2)
  n_rel_outcome <- floor(p / 4)

  covar_confound <- perm[1:n_confound]
  covar_rel_outcome <- perm[(n_confound + 1):(n_confound + n_rel_outcome)]
  covar_rel_treatment <- perm[(n_confound + n_rel_outcome + 1):p]

  # Column indices used by the population treatment and outcome models
  covar_for_treatment <- c(covar_confound, covar_rel_treatment)
  covar_for_outcome <- c(covar_confound, covar_rel_outcome)

  # Model coefficients (aligned to columns of X)
  beta <- runif(p, min = -0.4, max = 0.4)
  a0 <- -0.18
  g <- 0.3
  alpha <- runif(p, min = -0.2, max = 0.3)

  # Covariates carrying quadratic/interaction terms in the complex scenarios
  J_T <- floor(length(covar_for_treatment) / 2)
  cvars_T <- sample(covar_for_treatment, J_T)
  J_Y <- floor(length(covar_for_outcome) / 2)
  cvars_Y <- sample(covar_for_outcome, J_Y)

  # Shared covariate pipeline: correlated MVN draws transformed by inverse CDF
  # into half normal(0,1), a quarter uniform(0,1), a quarter Bernoulli(0.5)
  num_norm_vars <- floor(p / 2)
  num_bern_vars <- floor(p / 4)
  num_uniform_vars <- p - num_norm_vars - num_bern_vars
  make_X <- function(nn) {
    vars_unif <- pnorm(mvrnorm(nn, mu = numeric(p), Sigma = cor))
    X <- cbind(
      qnorm(vars_unif[, 1:num_norm_vars]),
      vars_unif[, (num_norm_vars + 1):(num_norm_vars + num_uniform_vars)],
      qbinom(vars_unif[, (num_norm_vars + num_uniform_vars + 1):(num_norm_vars + num_uniform_vars + num_bern_vars)], size = 1, prob = 0.5)
    )
    colnames(X) <- sprintf("v%d", 1:p)
    X
  }

  # Treatment-model linear predictor WITHOUT the intercept (base terms plus,
  # in complex cells, J quadratics and J-1 adjacent-pair interactions)
  lp_treatment <- function(X) {
    lp <- X[, covar_for_treatment, drop = FALSE] %*% beta[covar_for_treatment]
    if (scenarioT == "complex_T") {
      Xc <- X[, cvars_T, drop = FALSE]
      lp <- lp + (Xc^2) %*% beta[cvars_T]
      if (J_T >= 2) {
        lp <- lp + (Xc[, -J_T, drop = FALSE] * Xc[, -1, drop = FALSE]) %*% beta[cvars_T[-J_T]]
      }
    }
    lp
  }

  # Calibrate the intercept so the population marginal treatment rate is 0.50:
  # draw a seeded calibration sample and root-find b0 on E[plogis(b0 + lp)] = .5
  lp_cal <- lp_treatment(make_X(20000))
  b0 <- uniroot(function(b) mean(plogis(b + lp_cal)) - 0.5,
                interval = c(-50, 50), tol = 1e-6)$root

  # Population fully specified; hand the RNG back to the replication stream
  assign(".Random.seed", rep_stream, envir = .GlobalEnv)

  #########################################
  # Correlated covariates (this replication's sample)
  #########################################

  X <- make_X(n)

  #########################################
  # Population treatment model
  #########################################

  # True propensity score and binary treatment assignment
  trueps <- as.numeric(plogis(b0 + lp_treatment(X)))
  unif1 <- runif(n, 0, 1)
  T <- ifelse(unif1 < trueps, 1, 0)

  #########################################
  # Population outcome model
  #########################################

  # Error term
  e <- rnorm(n, mean = 0, sd = sqrt(0.17))

  # Base model: main effects only
  Y <- a0 + g * T + X[, covar_for_outcome, drop = FALSE] %*% alpha[covar_for_outcome] + e

  # Complex model: quadratics and adjacent-pair interactions, per eq-outcome_d
  if (scenarioY == "complex_Y") {
    Xc <- X[, cvars_Y, drop = FALSE]
    Y <- Y + (Xc^2) %*% alpha[cvars_Y]
    if (J_Y >= 2) {
      Y <- Y + (Xc[, -J_Y, drop = FALSE] * Xc[, -1, drop = FALSE]) %*% alpha[cvars_Y[-J_Y]]
    }
  }

  #########################################
  # Form simulated data tibble
  #########################################

  dat <- as_tibble(X)
  dat$T <- T
  dat$Y <- as.numeric(Y)
  dat$trueps <- trueps
  dat
}
