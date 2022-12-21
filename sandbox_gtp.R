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
vars_transformed <- cbind(vars_normal, vars_uniform, vars_bern)

# Generate variable names and store in the master_covar list
colnames(vars_transformed) <- paste0(rep(c("v"), p), 1:p)
master_covar <- colnames(vars_transformed)


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




b <- sub(".*v","",covar_for_treatment)
element <- paste0('b',b,' * ', covar_for_treatment)
form <- paste("(1 + exp(-(0 + b1 * v1 +", paste(element, collapse = " + ")) 




#cov_20 <- sample(element,((.2)*length(element)))

if (scenarioT == "A") {
  
  #base model 
  
  f_form <- paste("(1 + exp(-(0 + b1 * v1 +", paste(element, collapse = " + "),")))^-1") 
  trueps <- eval(parse(text = f_form))
  
} else
  
  if (scenarioT == "B") {
    
    #non-linear model
    add_lin <- paste0(cov_20, '*',cov_20)
    b_2 <- sub(".*v","",cov_20)
    f_form2 <- paste(form," + ", paste0('b',b_2,'*',add_lin, collapse = " + "),")))^-1")
    trueps <- eval(parse(text = f_form2))
    
  } else
    
    if (scenarioT == "C") {
      
      #non-additive model
      n1 <- 0.5
      add_cov <- sample(cov_80,((.2)*length(cov_80)))
      b_3 <- sub(".*v","",add_cov)
      b_4 <- sample(b,length(add_cov))
      add <- paste0('b',b_3,'*',n1,'*','v',b_4,'*',add_cov)
      f_form3 <- paste(form," + ", paste(add, collapse = " + "),")))^-1") 
      trueps <- eval(parse(text = f_form3))
      
    } else
      
      if (scenarioT == "D") {
        
        #complex model
        n1 <- 0.5
        add_lin <- paste0(cov_20, '*',cov_20)
        add_cov <- sample(cov_80,((.2)*length(cov_80)))
        b_5 <- sub(".*v","",add_cov)
        b_6 <- sample(b,length(add_cov))
        add <- paste0('b',b_5,'*',n1,'*','v',b_6,'*',add_cov)
        f_form4 <- paste(form," + ", paste(add_lin, ' + ',add, collapse = " + "),")))^-1") 
        trueps <- eval(parse(text = f_form4))
      }