



# Generate a coefficients for population outcome models the intercept in Y(1) has been fixed to obtain Y(1)-Y(0)=0.3 in each scenario

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
