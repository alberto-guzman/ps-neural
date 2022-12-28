library(simhelpers)
library(dplyr)
library(tibble)
library(purrr)
library(tidyr)
library(knitr)
library(kableExtra)
library(broom)
library(ggplot2)
library(MASS)
library(Matrix)
library(psych)
library(MBESS)
library(Rlab)
library(tidyverse)





generate_dat <- function(n, p) {
  # Generate a mean vector of 0s
  mean <- rep(0, p)

  # Generate a correlation matrix with correlations between 0 and 0.6
  cor <- Matrix(runif(p^2, 0, 0.6), nrow = p)
  diag(cor) <- 1

  # Convert the correlation matrix to a numeric matrix
  cor <- as.matrix(cor)

  # Smooth the correlation matrix to ensure it is positive definite
  cor <- psych::cor.smooth(cor)

  # Generate a vector of standard deviations
  sd <- rep(1, p)

  # Convert the correlation matrix to a covariance matrix
  cov <- MBESS::cor2cov(cor, sd)

  # Generate correlated normal variables
  vars <- mvrnorm(n, mean, cov)

  # Calculate the number of normal, Bernoulli, and uniform variables to generate
  num_norm_vars <- floor(p / 2)
  num_bern_vars <- floor(p / 4)
  num_uniform_vars <- p - num_norm_vars - num_bern_vars

  # Convert all variables to uniform variables between 0 and 1
  vars_unif <- pnorm(vars)

  # Convert the first num_norm_vars variables to normal variable with mean = 1 and sd = 1
  vars_normal <- qnorm(vars_unif[, 1:num_norm_vars])

  # Convert the next num_bern_vars variables to Bernoulli variables
  vars_bern <- qbern(vars_unif[, (num_norm_vars + 1):(num_norm_vars + num_bern_vars)], 0.5)

  # The remainder are left as uniform variables
  vars_uniform <- vars_unif[, (num_norm_vars + num_bern_vars + 1):p]

  # Combine the transformed variables
  vars_transformed <- as.data.frame(cbind(vars_normal, vars_uniform, vars_bern))

  # Generate variable names and store in the master_covar list
  colnames(vars_transformed) <- paste0(rep(c("v"), p), 1:p)
  master_covar <- colnames(vars_transformed)

  # Sample half of the covariates and save to covar_confound
  covar_confound <- sample(master_covar, size = length(master_covar) / 2)

  # Sample a quarter of the covariates and save to covar_rel_outcome
  covar_rel_outcome <- sample(setdiff(master_covar, covar_confound), size = length(master_covar) / 4)

  # Save the remaining covariates to covar_rel_treatment
  covar_rel_treatment <- setdiff(master_covar, union(covar_confound, covar_rel_outcome))

  # Combine covar_confound and covar_rel_outcome, these are the covariates that will be used for the population outcome models
  covar_for_treatment <- union(covar_confound, covar_rel_treatment)

  # Combine covar_confound and covar_rel_outcome, these are the covariates that will be used for the population outcome models
  covar_for_outcome <- union(covar_confound, covar_rel_outcome)
  
  # Split the data frame into a list of numeric vectors
  list_of_vectors <- split(vars_transformed, names(vars_transformed))
  
  # Assign the vectors to variables in the current environment
  for (name in names(list_of_vectors)) {
    assign(name, list_of_vectors[[name]])
  }
  
  return(v1)
}

test <- generate_dat(n = 10000, p = 10) 
