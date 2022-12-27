# Generate a coefficients for population outcome models the intercept in Y(1) has been fixed to obtain Y(1)-Y(0)=0.3 in each scenario
# Generate coefficients for potential outcomes equation Y(0) and Y(1)

alpha <- vector("list", length(master_covar)+1)
for (i in 0:length(master_covar)) {
  # Generate a random beta value
  g <- round(rbeta(1, 1, 1), 2)
  # Assign the value to a0, a1, a2, etc.
  assign(paste0("a", i), g)
  a <- paste0("a", i)
  alpha[[i+1]] <- a
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
equation_Y0 <- paste0("a0", " + ", paste(element, collapse = " + "))
equation_Y1 <- paste0("a0 + 0.3 ", " + ", paste(element, collapse = " + "))

# Evaluate the equation
Y0 <- eval(parse(text = equation_Y0))
Y1 <- eval(parse(text = equation_Y1))

#########################################
# Non-linear model
#########################################

# Split covar_for_treatment into four equal-sized groups
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

# Concatenate the element to the equations
equation_Y0 <- paste0("a0 + ", paste(c(unlist(quadratic_terms), element), collapse = " + "), "")
equation_Y1 <- paste0("a0 + 0.3 +", paste(c(unlist(quadratic_terms), element), collapse = " + "), "")

# Evaluate the equation
Y0 <- eval(parse(text = equation_Y0))
Y1 <- eval(parse(text = equation_Y1))

# continuous outcome Y

Y <- T*Y1 + (1-T)*Y0

# true individual effect of T on Y

indeff <- Y1-Y0
