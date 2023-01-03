
library(rpart)
library(ipred)
library(randomForest)
library(nnet)
library(survey)

sim <- generate_dat(n = 100000, p = 20, "A","a")

mod = glm(T~ . - Y - trueps, data=sim, family=binomial)
ps = mod$fitted

# Add the predicted propensity scores and weights to the tibble
sim <- sim %>%
  mutate(ps_pred = ps,
         ps_weights = if_else(T == 0, ps / (1 - ps), 1))



sim_T1 = subset(sim, T == 0)
mean_ps_weights = mean(sim_T1$ps_weights, na.rm = TRUE)




d.w <- svydesign(~1, weights = sim$ps_weights, data = sim)
fit <- svyglm(Y ~ T, design = d.w)
summary(fit)


pred_ATT = coef(fit)["T"]
vcov_matrix = vcov(fit)
se_ATT = sqrt(vcov_matrix["T", "T"])





lower_bound = pred_ATT - 1.96 * se_ATT
upper_bound = pred_ATT + 1.96 * se_ATT
ci_95 = ifelse(lower_bound <= 0.3 && 0.3 <= upper_bound, 1, 0)


ARB = abs(pred_ATT - 0.3) / 0.3

# calculate outside the loop
MSE = mean((pred_ATT - 0.3)^2)









# Subset the data into the treatment and comparison groups
treatment_group <- sim[sim$T == 1, ]
comparison_group <- sim[sim$T == 0, ]

# Get the names of the variables that start with "v"
var_names <- names(sim)[grep("^v", names(sim))]

# Initialize the ASAM_list vector
ASAM_list <- rep(NA, length(var_names))

# Loop through each covariate
for (i in 1:length(var_names)) {
  
  # Get the covariate name
  covariate <- var_names[i]
  
  # Extract the covariate data from the treatment and comparison groups
  treatment_data <- treatment_group[[covariate]]
  comparison_data <- comparison_group[[covariate]]
  
  # Extract the weights from the treatment and comparison groups
  treatment_weights <- treatment_group$ps_weights
  comparison_weights <- comparison_group$ps_weights
  
  # Calculate the means of the treatment and comparison groups
  treatment_mean <- weighted.mean(treatment_data, treatment_weights)
  comparison_mean <- weighted.mean(comparison_data, comparison_weights)
  
  # Calculate the variances of the treatment and comparison groups
  treatment_var <- wtd.var(treatment_data, treatment_weights)
  comparison_var <- wtd.var(comparison_data, comparison_weights)
  
  # Calculate the standard deviations of the treatment and comparison groups
  treatment_sd <- sqrt(treatment_var)
  comparison_sd <- sqrt(comparison_var)
  
  # Calculate the standardized difference of means
  sd_diff <- (treatment_mean - comparison_mean) / treatment_sd
  
  # Take the absolute value of the standardized difference of means
  abs_sd_diff <- abs(sd_diff)
  
  # Save the absolute standardized difference of means in the ASAM_list vector
  ASAM_list[i] <- abs_sd_diff
}

# Calculate the mean of the absolute standardized differences of means
ASAM <- mean(ASAM_list)



# Performance Metrics I need
# ASAM DONE!
# Bias DONE!
# Weights DONE!
# Coverage DONE!
# MSE DONE!




mod = rpart(T~ . - Y - trueps,method="class",data=sim)
ps = predict(mod)[,2] 

mod = bagging(T ~ . - Y - trueps, data=sim)
ps =  predict(mod,newdata=sim,type="prob")

mod = randomForest(factor(T)~ . - Y - trueps, data=sim)
ps<-predict(mod, type="prob")[,2]






