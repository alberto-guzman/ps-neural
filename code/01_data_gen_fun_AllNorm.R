#############

## WHAT DOES THIS FUNCTIOND DO?

# This function generates correlated normal variables, converts them into uniform, normal, and Bernoulli variables, and then assigns each column of the resulting transformed variables to an object with a name "v1", "v2", etc.
# It also samples various subsets of these variables and assigns them to different sets (covar_confound, covar_rel_outcome, etc.).
# Finally, it generates coefficients for various population treatment and outcome models and uses them to simulate treatment and outcome data for a population, and returns this data in a tibble.


# Based correlation, intercept, and coeffecient values on CommonApp study
#############

Generate <- function(condition, fixed_objects = NULL) {
  # Makes the tibble of sim crossed conditions accessible to the function, from SimDesign package
  Attach(condition)
  # Generate a mean vector of 0s
  mean <- rep(0, p)

  # Generate a correlation matrix with correlations between -.3 to .3
  cor <- Matrix(runif(p^2, -.3, .3), nrow = p)
  diag(cor) <- 1

  # Convert the correlation matrix to a numeric matrix
  cor <- as.matrix(cor)

  # Smooth the correlation matrix to ensure it is positive definite
  cor <- psych::cor.smooth(cor)

  # Generate correlated normal variables
  vars <- mvrnorm(n, mean, cor)

  
  # Give the columns of vars_transformed names v1,v2,etc.
  colnames(vars) <- paste0("v", 1:p)

  # Generate variable names and store in the master_covar list
  master_covar <- dimnames(vars)[[2]]

  # Create p objects with names v1, v2, etc., and assign the corresponding column from vars_transformed to each object
  for (i in 1:p) {
    assign(paste0("v", i), vars[, i])
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
  # Initialize b0 to -0.25
  b0 <- 0.25

  # Create an empty list to store the b coefficients
  beta <- vector("list", length(master_covar))

  # Loop through all variables in the master covariate list
  for (i in seq_len(length(master_covar))) {
    # Generate a random number between -0.4 and 0.4
    g <- runif(1, min = -0.4, max = 0.4)

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
  # Population treatment model - Generate base model
  #########################################
  if (scenarioT == "A") {
    # Concatenate the variables from covar_for_treatment into a single string
    equation <- paste0("(1 + exp(-(b0 + ", paste(element, collapse = " + "), ")))^-1")

    # Evaluate the equation and store the result in trueps
    trueps <- eval(parse(text = equation))
  } else

  #########################################
  # Population treatment model - Non-linear model
  #########################################
  if (scenarioT == "B") {
    # Split covar_for_treatment into four equal-sized groups
    group_n <- length(covar_for_treatment)
    split_index <- round(seq(1, group_n, length.out = 4))
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
  # Population treatment model - Non-additive model
  #########################################
  if (scenarioT == "C") {
    # Split covar_for_treatment into four equal-sized groups
    group_n <- length(covar_for_treatment)
    split_index <- round(seq(1, group_n, length.out = 4))
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
  # Population treatment model - Complex model
  #########################################
  if (scenarioT == "D") {
    # Split covar_for_treatment into four equal-sized groups
    group_n <- length(covar_for_treatment)
    split_index <- round(seq(1, group_n, length.out = 4))
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
  # Initialize population treatment effect g to 0.3
  g <- 0.3

  alpha <- vector("list", length(master_covar))

  for (i in 1:length(master_covar)) {
    # Generate a random number between -0.2 and 0.3
    g <- runif(1, min = -0.2, max = 0.3)
    # Assign the value to a1, a2, a3, etc.
    assign(paste0("a", i), g)
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
  if (scenarioY == "a") {
    # Concatenate the variables from covar_for_outcome into a single string
    equation_Y <- paste0("a0 + g * T", " + ", paste(element, collapse = " + "))

    # Evaluate the equation
    Y <- eval(parse(text = equation_Y))
  } else

  #########################################
  # Population outcome model - Non-linear model
  #########################################
  if (scenarioY == "b") {
    # Split covar_for_outcome into four equal-sized groups
    group_n <- length(covar_for_outcome)
    split_index <- round(seq(1, group_n, length.out = 4))
    covar_groups <- split(covar_for_outcome, cut(seq_along(covar_for_outcome), split_index))

    # Get the first group of covariates
    first_group <- covar_groups[[1]]

    # Create an empty list to store the quadratic terms
    quadratic_terms <- list()

    # Iterate over each group of covariates and create an interaction term between the two covariates
    for (group in first_group) {
      # Extract the coefficient from the covariate name
      a <- sub(".*v", "", group)
      # Create the quadratic term with the format "a * covar^2"
      quadratic_terms[[group]] <- paste0("a", a, " * ", group, "^2")
    }

    equation_Y <- paste0("a0 + g * T + ", paste(c(unlist(quadratic_terms), element), collapse = " + "), "")

    # Evaluate the equation
    Y <- eval(parse(text = equation_Y))
  } else

  #########################################
  # Population outcome model - Non-additive model
  #########################################
  if (scenarioY == "c") {
    # Split covar_for_outcome into four equal-sized groups
    group_n <- length(covar_for_outcome)
    split_index <- round(seq(1, group_n, length.out = 4))
    covar_groups <- split(covar_for_outcome, cut(seq_along(covar_for_outcome), split_index))

    # Iterate over each group of covariates and create an interaction term between the two covariates
    interaction_terms <- lapply(covar_groups, function(group) {
      # Extract the coefficient from the covariate name
      a <- sub(".*v", "", group[1])
      # Create the interaction term
      paste0("a", a, " * ", group[1], " * ", group[2])
    })

    equation_Y <- paste0("a0 + g * T + ", paste(c(unlist(interaction_terms), element), collapse = " + "), "")

    # Evaluate the equation
    Y <- eval(parse(text = equation_Y))
  } else

  #########################################
  # Population outcome model - Complex model
  #########################################
  if (scenarioY == "d") {
    # Split covar_for_outcome into four equal-sized groups
    group_n <- length(covar_for_outcome)
    split_index <- round(seq(1, group_n, length.out = 4))
    covar_groups <- split(covar_for_outcome, cut(seq_along(covar_for_outcome), split_index))

    # Iterate over each group of covariates and create both a quadratic term and an interaction term for each one
    terms <- list()
    for (group in covar_groups) {
      # Extract the coefficients from the covariate names
      a <- sub(".*v", "", group[1])
      # Create the quadratic term
      terms[[group[1]]] <- paste0("a", a, " * ", group[1], "^2")
      # Create the interaction term
      terms[[paste0(group[1], " * ", group[2])]] <- paste0("a", a, " * ", group[1], " * ", group[2])
    }

    # Concatenate all of the terms together and store the result in a new variable called equation
    equation_Y <- paste0("a0 + g * T + ", paste(c(unlist(terms), element), collapse = " + "), "")


    # Evaluate the equation
    Y <- eval(parse(text = equation_Y))
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
