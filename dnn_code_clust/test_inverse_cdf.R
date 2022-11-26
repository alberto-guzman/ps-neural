library(MASS)


# Let's keep it simple, 
mu <- rep(0,100)
Sigma <- matrix(.3, nrow=100, ncol=100) + diag(100)*.3
rawvars <- mvrnorm(n=10000, mu=mu, Sigma=Sigma)


cov(rawvars); cor(rawvars)
pvars <- pnorm(rawvars)

cov(pvars); cor(pvars)



poisvars <- qpois(pvars, 5)
cor(poisvars, rawvars) 

