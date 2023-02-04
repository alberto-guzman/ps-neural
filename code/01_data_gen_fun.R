#############
## WHAT DOES THIS FUNCTIOND DO?
# Based correlation, intercept, and coefficient values on CommonApp study
#############

Generate <- function(condition, fixed_objects = NULL) {
  # Makes the tibble of sim crossed conditions accessible to the function, from SimDesign package
  Attach(condition)

  # Generate a mean vector of 0s
  mean <- numeric(p)

  # Generate a correlation matrix with correlations between -.3 to .3
  cor <- matrix(runif(p^2, min = -.3, max = .3), nrow = p)
  diag(cor) <- 1

  # Smooth the correlation matrix to ensure it is positive definite
  cor <- psych::cor.smooth(cor)

  # Generate correlated normal variables
  vars <- mvrnorm(n, mean, cor)

  # Calculate the number of normal, Bernoulli, and uniform variables to generate
  num_norm_vars <- floor(p / 2)
  num_bern_vars <- floor(p / 4)
  num_uniform_vars <- p - num_norm_vars - num_bern_vars

  # Convert all variables to uniform variables between 0 and 1
  vars_unif <- pnorm(vars)

  # Convert the first num_norm_vars variables to normal variable with mean = 0 and sd = 1
  vars_normal <- qnorm(vars_unif[, 1:num_norm_vars])

  # Convert the next num_bern_vars variables to Bernoulli variables with probability of success = 0.5
  vars_bern <- qbern(vars_unif[, (num_norm_vars + 1):(num_norm_vars + num_bern_vars)], 0.5)

  # The remainder are left as uniform variables
  vars_uniform <- vars_unif[, (num_norm_vars + num_bern_vars + 1):p]

  # Combine the transformed variables
  vars_transformed <- cbind(vars_normal, vars_uniform, vars_bern)

  # Give the columns of vars names v1,v2,etc.
  colnames(vars_transformed) <- sprintf("v%d", 1:p)

  # Generate variable names and store in the master_covar list
  master_covar <- dimnames(vars_transformed)[[2]]

  # Create p objects with names v1, v2, etc. in working environment
  for (i in 1:p) {
    assign(colnames(vars_transformed)[i], vars_transformed[, i])
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

  #########################################
  #########################################
  # Population treatment models
  #########################################
  #########################################

  # Generate b coefficients for population treatment models
  # Initialize b0 to 0.25
  b0 <- 0.25

  # Create an empty list to store the b coefficients
  beta <- vector("list", length(master_covar))

  # Loop through all variables in the master covariate list
  for (i in seq_len(length(master_covar))) {
    # Generate a random number between -0.4 and 0.4
    x <- runif(1, min = -0.4, max = 0.4)

    # Assign the value to a variable named b1, b2, etc.
    assign(paste0("b", i), x)

    # Store the variable names in the beta list
    b <- paste0("b", i)
    beta[[i]] <- b
  }

  # Extract the coefficient from the covariate name
  b <- sub(".*v", "", covar_for_treatment)

  # Create a new variable called element with the format "b * covar_for_treatment"
  element <- paste0("b", b, " * ", covar_for_treatment)

  #########################################
  # Population treatment model - Generate base model
  #########################################
  if (scenarioT == "base_T") {
    # Concatenate the variables from covar_for_treatment into a single string
    equation <- paste0("(1 + exp(-(b0 + ", paste(element, collapse = " + "), ")))^-1")

    # Evaluate the equation and store the result in trueps
    trueps <- eval(parse(text = equation))
  } else

  #########################################
  # Population treatment model - Complex model
  #########################################
  if (scenarioT == "complex_T") {
    # Sample half of the variables from covar_for_treatment
    sample_vars <- sample(covar_for_treatment, length(covar_for_treatment) / 2)

    # Create a list to store the terms
    terms <- list()

    # Iterate over the sampled variables and create the quadratic terms
    for (var in sample_vars) {
      b <- sub(".*v", "", var)
      quad_term <- paste0("b", b, " * ", var, "^2")
      terms[[var]] <- quad_term
    }

    # Sample half of the variables again from covar_for_treatment
    sample_vars2 <- sample(covar_for_treatment, length(covar_for_treatment) / 2)

    # Create a list of all possible interactions between the variables
    interactions <- combn(sample_vars2, 2, paste0, collapse = "*")

    # Iterate over the interactions and create the interaction terms
    for (inter in interactions) {
      b <- sub(".*v", "", inter)
      inter_term <- paste0("b", b, " * ", inter)
      terms[[inter]] <- inter_term
    }

    # Concatenate all of the terms together and store the result in a new variable called equation
    equation <- paste0("(1 + exp(-(b0 + ", paste(c(unlist(terms), element), collapse = " + "), ")))^-1")

    # Evaluate the equation
    trueps <- eval(parse(text = equation))
  }

  #########################################
  # ~~ binary treatment T
  #########################################

  unif1 <- runif(n, 0, 1)
  T <- ifelse(unif1 < trueps, 1, 0)

  #########################################
  #########################################
  # Population outcome models
  #########################################
  #########################################

  # Generate a coefficients for population outcome models
  # Initialize a0 to -0.18
  a0 <- -0.18

  # Generate error terms for population outcome models
  e <- rnorm(n, mean = 0, sd = sqrt(0.17))
  
  alpha <- vector("list", length(master_covar))

  for (i in 1:length(master_covar)) {
    # Generate a random number between -0.2 and 0.3
    x <- runif(1, min = -0.2, max = 0.3)
    # Assign the value to a1, a2, a3, etc.
    assign(paste0("a", i), x)
    a <- paste0("a", i)
    alpha[[i + 1]] <- a
  }

  # Extract the coefficient from the covariate name
  a <- sub(".*v", "", covar_for_outcome)

  # Create a new variable called element with the format "a * covar_for_outcome"
  element <- paste0("a", a, " * ", covar_for_outcome)
  
  #########################################
  # Population outcome model - Generate base model
  #########################################
  if (scenarioY == "base_Y") {
    equation <- paste0("a0 + g * T", " + ", paste(element, collapse = " + ")," + e")
    Y <- eval(parse(text = equation))
  } else

  #########################################
  # Population outcome model - Complex model
  #########################################
  if (scenarioY == "complex_Y") {
    # Sample half of the variables from covar_for_outcome
    sample_vars <- sample(covar_for_outcome, length(covar_for_outcome) / 2)

    # Create a list to store the terms
    terms <- list()

    # Iterate over the sampled variables and create the quadratic terms
    for (var in sample_vars) {
      a <- sub(".*v", "", var)
      quad_term <- paste0("a", a, " * ", var, "^2")
      terms[[var]] <- quad_term
    }

    # Sample half of the variables again from covar_for_outcome
    sample_vars2 <- sample(covar_for_outcome, length(covar_for_outcome) / 2)

    # Create a list of all possible interactions between the variables
    interactions <- combn(sample_vars2, 2, paste0, collapse = "*")

    # Iterate over the interactions and create the interaction terms
    for (inter in interactions) {
      a <- sub(".*v", "", inter)
      inter_term <- paste0("a", a, " * ", inter)
      terms[[inter]] <- inter_term
    }

    equation <- paste0("a0 + g * T + ", paste(c(unlist(terms), element), collapse = " + "), " + e")
    Y <- eval(parse(text = equation))
    
  }

  #########################################
  # Form simulated data tibble
  #########################################

  v_list <- mget(paste0("v", 1:length(master_covar)))
  dat <- as_tibble(v_list)
  dat$T <- T
  dat$Y <- Y
  dat$trueps <- trueps
  dat
}
