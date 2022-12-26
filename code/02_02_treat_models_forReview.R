# Generate b coefficients for population treatment models
# Initialize b0 to 0
b0 <- 0

# Create an empty list to store the b coefficients
beta <- vector("list", length(master_covar))

# Loop through all variables in the master covariate list
for (i in seq_len(length(master_covar))) {
  # Generate a random beta value
  g <- round(rbeta(1, 1, 1), 2)

  # Assign the value to a variable named b1, b2, etc.
  assign(paste0("b", i), g)

  # Store the variable names in the beta list
  b <- paste0("b", i)
  beta[[i]] <- b
}

# Add the variables from vars_transformed to the global environment
list2env(vars_transformed, .GlobalEnv)

#########################################
# Generate base model
#########################################

# Extract the coefficient from the covariate name
b <- sub(".*v", "", covar_for_treatment)

# Create a new variable called element with the format "b * covar_for_treatment"
element <- paste0("b", b, " * ", covar_for_treatment)

# Concatenate the variables from covar_for_treatment into a single string
equation <- paste0("(1 + exp(-(0 + ", paste(element, collapse = " + "), ")))^-1")

# Evaluate the equation and store the result in trueps
trueps <- eval(parse(text = equation))

#########################################
# Non-linear model
#########################################

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
equation <- paste0("(1 + exp(-(0 + ", paste(c(unlist(quadratic_terms), element), collapse = " + "), ")))^-1")

# Evaluate the equation and store the result in trueps
trueps <- eval(parse(text = equation))

#########################################
# Non-additive model
#########################################

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
equation <- paste0("(1 + exp(-(0 + ", paste(c(unlist(interaction_terms), element), collapse = " + "), ")))^-1")

# Evaluate the equation
trueps <- eval(parse(text = equation))

#########################################
# Complex model
#########################################

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
equation <- paste0("(1 + exp(-(0 + ", paste(c(unlist(terms), element), collapse = " + "), ")))^-1")

# Evaluate the equation
trueps <- eval(parse(text = equation))
