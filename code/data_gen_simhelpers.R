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

generate_dat <- function(n, p, scenarioT) {
  # Generate a mean vector of 0s
  mean <- rep(0, p)

  # Generate a correlation matrix with correlations between -0.6 and 0.6
  cor <- Matrix(runif(p^2, -0.6, 0.6), nrow = p)
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
  vars_transformed <- cbind(vars_normal, vars_uniform, vars_bern)

  # Give the columns of vars_transformed names v1,v2,etc.
  colnames(vars_transformed) <- paste0("v", 1:p)

  # Generate variable names and store in the master_covar list
  master_covar <- dimnames(vars_transformed)[[2]]

  # Create p objects with names v1, v2, etc., and assign the corresponding column from vars_transformed to each object
  for (i in 1:p) {
    assign(paste0("v", i), vars_transformed[, i])
  }

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

  # Generate b coefficients for population treatment models
  # Initialize b0 to -1
  b0 <- -1

  # Create an empty list to store the b coefficients
  beta <- vector("list", length(master_covar))

  # Loop through all variables in the master covariate list
  for (i in seq_len(length(master_covar))) {
    # Generate a random beta value between -0.5 and 0.5
    g <- round(runif(1, min = -0.5, max = 0.5), 2)

    # Assign the value to a variable named b1, b2, etc.
    assign(paste0("b", i), g)

    # Store the variable names in the beta list
    b <- paste0("b", i)
    beta[[i]] <- b
  }

  # Extract the coefficient from the covariate name
  b <- sub(".*v", "", covar_for_treatment)

  # Create a new variable called element with the format "b * covar_for_treatment"
  element <- paste0("b", b, " * ", covar_for_treatment)

  #########################################
  # Generate base model
  #########################################
  if (scenarioT == "A") {
    # Concatenate the variables from covar_for_treatment into a single string
    equation <- paste0("(1 + exp(-(b0 + ", paste(element, collapse = " + "), ")))^-1")

    # Evaluate the equation and store the result in trueps
    trueps <- eval(parse(text = equation))
  } else

  #########################################
  # Non-linear model
  #########################################
  if (scenarioT == "B") {
    # Split covar_for_treatment into four equal-sized groups
    n <- length(covar_for_treatment)
    split_index <- round(seq(1, n, length.out = 4))
    covar_groups <- split(covar_for_treatment, cut(seq_along(covar_for_treatment), split_index))

    # Get the first group of covariates
    first_group <- covar_groups[[1]]

    # Create an empty list to store the quadratic terms
    quadratic_terms <- list()

    # Iterate over the first group of covariates and create a quadratic term for each one
    for (group in first_group) {
      # Extract the coefficient from the covariate name
      b <- sub(".*v", "", group)
      # Create the quadratic term with the format "b * covar^2"
      quadratic_terms[[group]] <- paste0("b", b, " * ", group, "^2")
    }

    # Concatenate all of the terms together and store the result in a new variable called equation
    equation <- paste0("(1 + exp(-(b0 + ", paste(c(unlist(quadratic_terms), element), collapse = " + "), ")))^-1")

    # Evaluate the equation and store the result in trueps
    trueps <- eval(parse(text = equation))
  } else

  #########################################
  # Non-additive model
  #########################################
  if (scenarioT == "C") {
    # Split covar_for_treatment into four equal-sized groups
    n <- length(covar_for_treatment)
    split_index <- round(seq(1, n, length.out = 4))
    covar_groups <- split(covar_for_treatment, cut(seq_along(covar_for_treatment), split_index))

    # Iterate over each group of covariates and create an interaction term between the two covariates
    interaction_terms <- lapply(covar_groups, function(group) {
      # Extract the coefficient from the covariate name
      b <- sub(".*v", "", group[1])
      # Create the interaction term
      paste0("b", b, " * ", group[1], " * ", group[2])
    })

    # Concatenate all of the terms together and store the result in a new variable called equation
    equation <- paste0("(1 + exp(-(b0 + ", paste(c(unlist(interaction_terms), element), collapse = " + "), ")))^-1")

    # Evaluate the equation
    trueps <- eval(parse(text = equation))
  } else

  #########################################
  # Complex model
  #########################################
  if (scenarioT == "D") {
    # Split covar_for_treatment into four equal-sized groups
    n <- length(covar_for_treatment)
    split_index <- round(seq(1, n, length.out = 4))
    covar_groups <- split(covar_for_treatment, cut(seq_along(covar_for_treatment), split_index))

    # Iterate over each group of covariates and create both a quadratic term and an interaction term for each one
    terms <- list()
    for (group in covar_groups) {
      # Extract the coefficients from the covariate names
      b <- sub(".*v", "", group[1])
      # Create the quadratic term
      terms[[group[1]]] <- paste0("b", b, " * ", group[1], "^2")
      # Create the interaction term
      terms[[paste0(group[1], " * ", group[2])]] <- paste0("b", b, " * ", group[1], " * ", group[2])
    }

    # Concatenate all of the terms together and store the result in a new variable called equation
    equation <- paste0("(1 + exp(-(b0 + ", paste(c(unlist(terms), element), collapse = " + "), ")))^-1")

    # Evaluate the equation
    trueps <- eval(parse(text = equation))
  }

  #########################################
  # Generate binary treatment T
  #########################################

  unif1 <- runif(n, 0, 1)
  T <- ifelse(trueps > unif1, 1, 0) # there is a probability of unif1 that T=1

  # return(trueps)
  return(T)
}


T <- generate_dat(n = 100000, p = 200, scenarioT = "D")
hist(T)
mean(T)
