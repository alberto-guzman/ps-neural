# Load the required packages
library(MASS)
library(Matrix)
library(psych)
library(MBESS)

# Set the number of variables (p) and the sample size (n)
p <- 100
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
num_bernoulli_vars <- floor(p / 4)
num_uniform_vars <- p - num_norm_vars - num_bernoulli_vars

# Convert the first num_norm_vars variables to normal variables
vars_normal <- qnorm(pnorm(vars[, 1:num_norm_vars]), mean=mean[1:num_norm_vars], sd=sd[1:num_norm_vars])

# Convert the next num_bernoulli_vars variables to Bernoulli variables
vars_bernoulli <- qbinom(pnorm(vars[, (num_norm_vars+1):(num_norm_vars+num_bernoulli_vars)]), 1, 0.5)

# Convert the remaining variables to uniform variables
vars_uniform <- qunif(pnorm(vars[, (num_norm_vars+num_bernoulli_vars+1):p]), min=min(vars[, (num_norm_vars+num_bernoulli_vars+1):p]), max=max(vars[, (num_norm_vars+num_bernoulli_vars+1):p]))

# Combine the transformed variables
vars_transformed <- cbind(vars_bernoulli, vars_uniform, vars_normal)

