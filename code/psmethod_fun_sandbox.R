library(rpart)
library(ipred)
library(randomForest)
library(nnet)
library(survey)
library(Hmisc)

data <- generate_data(n = 100000, p = 20, "A","a")


list <- run_methods(data,"logit")

run_methods <- function(data, method){
  
  df <- data
  
  #~~ estimate ps     
  if (method == "logit") {
    
    mod = glm(T~ . - Y - trueps, data=df, family=binomial)
    ps = mod$fitted
    
  } else
  
  if (method == "tree") {
    
    mod = rpart(T~ . - Y - trueps,method="class",data=df)
    ps = predict(mod)[,2] 
  
  }
  
  # Add the predicted propensity scores and weights to the tibble
  df <- df %>%
    mutate(ps_pred = ps,
           ps_weights = if_else(T == 0, ps / (1 - ps), 1))
  
  # estimates the ATT with the weights
  d.w <- svydesign(~1, weights = df$ps_weights, data = df)
  fit <- svyglm(Y ~ T, design = d.w)

  # saves the ATT and se_ATT
  pred_ATT = coef(fit)["T"]
  vcov_matrix = vcov(fit)
  se_ATT = sqrt(vcov_matrix["T", "T"])
  
  # calculate absolute relative bias
  ARB = abs(pred_ATT - 0.3) / 0.3
  
  # calculate mean of control gropu weights
  sim_T1 = subset(df, T == 0)
  mean_ps_weights = mean(sim_T1$ps_weights, na.rm = TRUE)
  
  # calculate 95% coverage
  lower_bound = pred_ATT - 1.96 * se_ATT
  upper_bound = pred_ATT + 1.96 * se_ATT
  ci_95 = ifelse(lower_bound <= 0.3 && 0.3 <= upper_bound, 1, 0)
  
  # calculate ASAM for covariates
  
  # Subset the data into the treatment and comparison groups
  treatment_group <- df[df$T == 1, ]
  comparison_group <- df[df$T == 0, ]
  
  # Get the names of the variables that start with "v"
  var_names <- names(df)[grep("^v", names(df))]
  
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

  return(c(pred_ATT,se_ATT,ARB,mean_ps_weights,ci_95,ASAM))
}

























mod = bagging(T ~ . - Y - trueps, data=sim)
ps =  predict(mod,newdata=sim,type="prob")

mod = randomForest(factor(T)~ . - Y - trueps, data=sim)
ps<-predict(mod, type="prob")[,2]






