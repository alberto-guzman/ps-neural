# Generate a coefficients for population outcome models
# Initialize a0 to 1
a0 <- 1
# Initialize treatment effect g to 0.3
g <- 0.3

alpha <- vector("list", length(master_covar))

for (i in 1:length(master_covar)) {
  # Generate a random beta value
  g <- round(rbeta(1, 1, 1), 2)
  # Assign the value to a1, a2, a3, etc.
  assign(paste0("a", i), g)
  a <- paste0("a", i)
  alpha[[i + 1]] <- a
}

# Add the variables from vars_transformed to the global environment
list2env(vars_transformed, .GlobalEnv)

#########################################
# Generate base model
#########################################

# Extract the coefficient from the covariate name
a <- sub(".*v", "", covar_for_outcome)

# Create a new variable called element with the format "a * covar_for_outcome"
element <- paste0("a", a, " * ", covar_for_outcome)

# Concatenate the variables from covar_for_outcome into a single string
equation_Y <- paste0("a0 + g * T ", " + ", paste(element, collapse = " + "))

# Evaluate the equation
Y <- eval(parse(text = equation_Y))

#########################################
# Non-linear model
#########################################

# Split covar_for_outcome into four equal-sized groups
n <- length(covar_for_outcome)
split_index <- round(seq(1, n, length.out = 4))
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

#########################################
# Non-additive model
#########################################

# Split covar_for_outcome into four equal-sized groups
n <- length(covar_for_outcome)
split_index <- round(seq(1, n, length.out = 4))
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

#########################################
# Complex model
#########################################

# Split covar_for_outcome into four equal-sized groups
n <- length(covar_for_outcome)
split_index <- round(seq(1, n, length.out = 4))
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

#########################################
# Form simulated data tibble
#########################################

v_list <- mget(paste0("v", 1:length(master_covar)))
sim <- as_tibble(v_list)
sim$T <- T
sim$Y <- Y
sim$trueps <- trueps
sim$indeff <- indeff
#return(sim)
