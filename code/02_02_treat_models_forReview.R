


# ~~~~~~~~~~~~~~Global Variables~~~~~~~~~~~~~~~~~~~~#

# ~~ coefficients for the population treatment and outcome models

b0 <- 0

# Generate coefficients for the population treatment and outcome models
beta <- vector("list", length(master_covar))
for (i in seq_len(length(master_covar))) {
  # Generate a random beta value
  g <- round(rbeta(1, 1, 1), 2)
  # Assign the value to a variable named b1, b2, etc.
  assign(paste0("b", i), g)
  # Store the variable names in the beta list
  b <- paste0("b", i)
  beta[[i]] <- b
}

# ~~ coefficients for potential outcomes equation Y(0) and Y(1)
# ~~ the intercept in Y(1) has been fixed to obtain Y(1)-Y(0)=0.3 in each scenario

a0 <- 3.5
a1 <- 0.5

# Generate coefficients for potential outcomes equation Y(0) and Y(1)
alpha <- vector("list", length(master_covar) - 1)
for (i in 2:length(master_covar)) {
  # Generate a random beta value
  g <- round(rbeta(1, 1, 1), 2)
  # Assign the value to a variable named a2, a3, etc.
  assign(paste0("a", i), g)
  # Store the variable names in the alpha list
  a <- paste0("a", i)
  alpha[[i - 1]] <- a
}

# Sample half of the covariates and save to covar_confound
covar_confound <- sample(master_covar, size = length(master_covar) / 2)

# Sample a quarter of the covariates and save to covar_rel_outcome
covar_rel_outcome <- sample(setdiff(master_covar, covar_confound), size = length(master_covar) / 4)

# Save the remaining covariates to covar_rel_treatment
covar_rel_treatment <- setdiff(master_covar, union(covar_confound, covar_rel_outcome))


# Combine covar_confound and covar_rel_treatment, these are the covariates that will be used for the population treatment models
covar_for_treatment <- union(covar_confound, covar_rel_treatment)

# Combine covar_confound and covar_rel_outcome, these are the covariates that will be used for the population outcome models
covar_for_outcome <- union(covar_confound, covar_rel_outcome)










#########################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~ generate scenarios for population treatment models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#########################################


#########################################
# base model
#########################################

b <- sub(".*v", "", covar_for_treatment)
element <- paste0("b", b, " * ", covar_for_treatment)

# Concatenate the variables from covar_for_treatment into a single string
equation <- paste0("(1 + exp(-(0 + ", paste(element, collapse = " + "), ")))^-1")

# Evaluate the equation
trueps <- eval(parse(text = equation))








#########################################
# non-linear model
#########################################

# Split covar_for_treatment into four equal-sized groups
n <- length(covar_for_treatment)
split_index <- round(seq(1, n, length.out = 4))
covar_groups <- split(covar_for_treatment, cut(seq_along(covar_for_treatment), split_index))

# Get the first group of covariates
first_group <- covar_groups[[1]]

# Iterate over the first group of covariates and create a quadratic term for each one
quadratic_terms <- list()
for (group in first_group) {
  # Extract the coefficient from the covariate name
  b <- sub(".*v", "", group)
  # Create the quadratic term
  quadratic_terms[[group]] <- paste0("b", b, " * ", group, "^2")
}

# Concatenate all of the terms together and store the result in a new variable called equation
equation <- paste0("(1 + exp(-(0 + ", paste(c(unlist(quadratic_terms), element), collapse = " + "), ")))^-1")



#########################################
# non-additive model
#########################################



n <- length(covar_for_treatment)
split_index <- round(seq(1, n, length.out = 4))
covar_groups <- split(covar_for_treatment, cut(seq_along(covar_for_treatment), split_index))




interaction_terms <- lapply(covar_groups, function(group) {
  # Extract the coefficients from the covariate names
  b1 <- sub(".*v", "", group[1])
  b2 <- sub(".*v", "", group[2])
  # Create the interaction term
  paste0("b", b1, " * b", b2, " * ", group[1], " * ", group[2])
})




equation <- paste0("(1 + exp(-(0 + ", paste(c(unlist(interaction_terms), element), collapse = " + "), ")))^-1")





#########################################
# complex model
#########################################




# Split covar_for_treatment into four equal-sized groups
n <- length(covar_for_treatment)
split_index <- round(seq(1, n, length.out = 4))
covar_groups <- split(covar_for_treatment, cut(seq_along(covar_for_treatment), split_index))

# Iterate over each group of covariates and create both a quadratic term and an interaction term for each one
terms <- list()
for (group in covar_groups) {
  # Extract the coefficients from the covariate names
  b1 <- sub(".*v", "", group[1])
  b2 <- sub(".*v", "", group[2])
  # Create the quadratic term
  terms[[group[1]]] <- paste0("b", b1, " * ", group[1], "^2")
  # Create the interaction term
  terms[[paste0(group[1], " * ", group[2])]] <- paste0("b", b1, " * b", b2, " * ", group[1], " * ", group[2])
}

# Concatenate all of the terms together and store the result in a new variable called equation
equation <- paste0("(1 + exp(-(0 + ", paste(c(unlist(terms), element), collapse = " + "), ")))^-1")
