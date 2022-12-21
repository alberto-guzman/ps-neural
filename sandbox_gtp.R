# Load the required packages
library(MASS)
library(Matrix)
library(psych)
library(MBESS)
library(Rlab)

# Set the number of variables (p) and the sample size (n)
p <- 10
n <- 1000

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
vars_transformed <- cbind(vars_bern, vars_uniform, vars_normal)







