#############
## WHAT DOES THIS FUNCTION DO?
# The Analyse function is used to estimate the average treatment effect (ATE) and related metrics for a given condition.
# The function uses one of several methods, specified by the method argument, to estimate the propensity score.
# The methods used to estimate the propensity score are
# logistic regression (logit), classification and regression trees (cart), bagging (bag), random forest (forest),
# and three neural network models (nn-1, dnn-2, and dnn-3). Once the propensity score is estimated,
# the function uses survey-weighted regression to estimate the ATE, standard error of the ATE, p-value, and 95% confidence interval of the ATE.
# The function also calculates the absolute standardized average mean (ASAM) for each covariate in the data.
#############

# function to estimate the ATE and other metrics
Analyse <- function(condition, dat, fixed_objects = NULL) {
  Attach(condition)

  # if the method is logit, then estimate the ATE using logistic regression
  if (method == "logit") {
    # estimate the propensity score using logistic regression
    mod <- glm(T ~ . - Y - trueps, data = dat, family = binomial)
    # predict on the entire dataframe to generate ps
    ps <- predict(mod, newdata = dat, type = "response")
    # if the method is cart, then estimate the ATE using classification and regression trees
  } else if (method == "cart") {
    # estimate the propensity score using classification and regression trees
    mod <- rpart(T ~ . - Y - trueps, method = "class", data = dat)
    # predict on the entire dataframe to generate ps
    ps <- predict(mod, newdata = dat, type = "prob")[, 2]
    # if the method is bag, then estimate the ATE using bagging
  } else if (method == "bag") {
    # estimate the propensity score using bagging
    mod <- bagging(T ~ . - Y - trueps, data = dat, nbagg = 100)
    # save the propensity score to a vector
    ps <- predict(mod, newdata = dat, type = "prob")
    # if the method is forest, then estimate the ATE using random forest
  } else if (method == "forest") {
    # estimate the propensity score using random forest
    mod <- randomForest(factor(T) ~ . - Y - trueps, data = dat)
    # save the propensity score to a vector
    ps <- predict(mod, newdata = dat, type = "prob")[, 2]
  } else if (method == "nn-1") {
    # Preprocess data
    # Split the data into training and validation sets (80/20)
    split <- sample(2, nrow(dat), replace = TRUE, prob = c(0.8, 0.2)) # random split of data
    train_data <- dat[split == 1, ]
    validation_data <- dat[split == 2, ]
    x_train <- as.matrix(train_data[, grep("^v", names(train_data))]) # select columns that start with "v" for input features
    y_train <- as.matrix(train_data[, "T"]) # select column for treatment assignment
    x_validation <- as.matrix(validation_data[, grep("^v", names(validation_data))]) # select columns that start with "v" for input features
    y_validation <- as.matrix(validation_data[, "T"]) # select column for treatment assignment

    # Define model
    p <- ncol(x_train) # number of input features
    input_layer <- layer_input(shape = c(p)) # input layer
    hidden_layer <- layer_dense(units = ceiling(2 * p / 3), activation = "relu", kernel_regularizer = regularizer_l2(l = 0.01))(input_layer)
    output_layer <- layer_dense(units = 1, activation = "sigmoid", kernel_regularizer = regularizer_l2(l = 0.01))(hidden_layer)
    model <- keras_model(inputs = input_layer, outputs = output_layer)

    # Compile model
    model %>% compile(
      optimizer = "adam",
      loss = "binary_crossentropy",
      metrics = c("accuracy")
    )

    # Define callbacks
    early_stopping <- callback_early_stopping(monitor = "val_loss", min_delta = 0.001, patience = 5)

    # Fit model
    history <- model %>% fit(
      x_train,
      y_train,
      epochs = 100,
      batch_size = 64,
      validation_data = list(x_validation, y_validation),
      callbacks = list(early_stopping),
      verbose = 0
    )

    # Preprocess data
    x <- as.matrix(dat[, grep("^v", names(dat))]) # select columns that start with "v" for input features

    # Predict propensity scores on entire dataset
    ps <- model %>% predict(x)
    ps <- ps[, 1]
  } else if (method == "dnn-2") {
    # Preprocess data
    # Split the data into training and validation sets (80/20)
    split <- sample(2, nrow(dat), replace = TRUE, prob = c(0.8, 0.2)) # random split of data
    train_data <- dat[split == 1, ]
    validation_data <- dat[split == 2, ]
    x_train <- as.matrix(train_data[, grep("^v", names(train_data))]) # select columns that start with "v" for input features
    y_train <- as.matrix(train_data[, "T"]) # select column for treatment assignment
    x_validation <- as.matrix(validation_data[, grep("^v", names(validation_data))]) # select columns that start with "v" for input features
    y_validation <- as.matrix(validation_data[, "T"]) # select column for treatment assignment

    # Define model
    p <- ncol(x_train) # number of input features
    input_layer <- layer_input(shape = c(p)) # input layer
    hidden_layer1 <- layer_dense(units = ceiling(2 * p / 3), activation = "relu", kernel_regularizer = regularizer_l2(l = 0.01))(input_layer) # first hidden layer
    hidden_layer2 <- layer_dense(units = ceiling(2 * p / 3), activation = "relu", kernel_regularizer = regularizer_l2(l = 0.01))(hidden_layer1) # second hidden layer
    output_layer <- layer_dense(units = 1, activation = "sigmoid", kernel_regularizer = regularizer_l2(l = 0.01))(hidden_layer2) # output layer
    model <- keras_model(inputs = input_layer, outputs = output_layer)

    # Compile model
    model %>% compile(
      optimizer = "adam",
      loss = "binary_crossentropy",
      metrics = c("accuracy")
    )

    # Define callbacks
    early_stopping <- callback_early_stopping(monitor = "val_loss", min_delta = 0.001, patience = 5)

    # Fit model
    history <- model %>% fit(
      x_train,
      y_train,
      epochs = 100,
      batch_size = 64,
      validation_data = list(x_validation, y_validation),
      callbacks = list(early_stopping),
      verbose = 0
    )

    # Preprocess data
    x <- as.matrix(dat[, grep("^v", names(dat))]) # select columns that start with "v" for input features

    # Predict propensity scores on entire dataset
    ps <- model %>% predict(x)
    ps <- ps[, 1]
  } else if (method == "dnn-3") {
    # Preprocess data
    # Split the data into training and validation sets (80/20)
    split <- sample(2, nrow(dat), replace = TRUE, prob = c(0.8, 0.2)) # random split of data
    train_data <- dat[split == 1, ]
    validation_data <- dat[split == 2, ]
    x_train <- as.matrix(train_data[, grep("^v", names(train_data))]) # select columns that start with "v" for input features
    y_train <- as.matrix(train_data[, "T"]) # select column for treatment assignment
    x_validation <- as.matrix(validation_data[, grep("^v", names(validation_data))]) # select columns that start with "v" for input features
    y_validation <- as.matrix(validation_data[, "T"]) # select column for treatment assignment

    # Define model
    p <- ncol(x_train) # number of input features
    input_layer <- layer_input(shape = c(p)) # input layer
    hidden_layer1 <- layer_dense(units = ceiling(2 * p / 3), activation = "relu", kernel_regularizer = regularizer_l2(l = 0.01))(input_layer) # first hidden layer
    hidden_layer2 <- layer_dense(units = ceiling(2 * p / 3), activation = "relu", kernel_regularizer = regularizer_l2(l = 0.01))(hidden_layer1) # second hidden layer
    hidden_layer3 <- layer_dense(units = ceiling(2 * p / 3), activation = "relu", kernel_regularizer = regularizer_l2(l = 0.01))(hidden_layer2) # third hidden layer
    output_layer <- layer_dense(units = 1, activation = "sigmoid", kernel_regularizer = regularizer_l2(l = 0.01))(hidden_layer3) # output layer
    model <- keras_model(inputs = input_layer, outputs = output_layer)


    # Compile model
    model %>% compile(
      optimizer = "adam",
      loss = "binary_crossentropy",
      metrics = c("accuracy")
    )

    # Define callbacks
    early_stopping <- callback_early_stopping(monitor = "val_loss", min_delta = 0.001, patience = 5)

    # Fit model
    history <- model %>% fit(
      x_train,
      y_train,
      epochs = 100,
      batch_size = 64,
      validation_data = list(x_validation, y_validation),
      callbacks = list(early_stopping),
      verbose = 0
    )

    # Preprocess data
    x <- as.matrix(dat[, grep("^v", names(dat))]) # select columns that start with "v" for input features

    # Predict propensity scores on entire dataset
    ps <- model %>% predict(x)
    ps <- ps[, 1]
  }

  ##############################
  ### calculate metrics
  ##############################

  dat <- dat %>%
    mutate(
      ps_pred = ps,
      ps_weights = case_when(T == 1 ~ 1 / ps, T == 0 ~ 1 / (1 - ps))
    )

  true_ATE <- 0.3

  # calculate standardized initial bias prior to weighting
  Std_In_Bias <- ((mean(dat$Y[dat$T == 1]) - mean(dat$Y[dat$T == 0])) - true_ATE) / sd(dat$Y[dat$T == 1])
  Prob_Treat <- mean(dat$T)

  # estimate the true_ATE with the weights
  d.w <- svydesign(~0, weights = dat$ps_weights, data = dat)
  fit <- svyglm(Y ~ T, design = d.w)

  # save the true_ATE and se_true_ATE
  ATE <- unname(coef(fit)["T"])
  # vcov_matrix <- vcov(fit)
  # ATE_se <- unname(sqrt(vcov_matrix["T", "T"]))

  modw <- lm(Y ~ T, data = dat, weights = dat$ps_weights)
  ATE_se <- summary(modw)$coefficients[c("T"), c("Std. Error")]

  # extract the p-value of T
  p_val <- summary(fit)$coefficients["T", "Pr(>|t|)"]

  # calculate the 95% coverage
  conf_interval <- confint(fit, level = 0.95)["T", ]
  lower_bound <- as.numeric(conf_interval[1])
  upper_bound <- as.numeric(conf_interval[2])
  ci_95 <- ifelse(lower_bound < true_ATE & true_ATE < upper_bound, 1, 0)

  # calculate the mean of weights
  mean_ps_weights <- mean(dat$ps_weights)

  ###############
  # calculate the ASAM for covariates
  ###############

  # subset the data into the treatment and comparison groups
  treatment_group <- dat[dat$T == 1, ]
  comparison_group <- dat[dat$T == 0, ]

  # get the names of the variables that start with "v"
  var_names <- names(dat)[grep("^v", names(dat))]

  # initialize the ASAM_list vector
  ASAM_list <- rep(NA, length(var_names))

  # loop through each covariate
  for (i in 1:length(var_names)) {
    # get the covariate name
    covariate <- var_names[i]

    # extract the covariate data from the treatment and comparison groups
    treatment_data <- treatment_group[[covariate]]
    comparison_data <- comparison_group[[covariate]]

    # extract the weights from the treatment and comparison groups
    treatment_weights <- treatment_group$ps_weights
    comparison_weights <- comparison_group$ps_weights

    # calculate the means of the treatment and comparison groups
    treatment_mean <- weighted.mean(treatment_data, treatment_weights)
    comparison_mean <- weighted.mean(comparison_data, comparison_weights)

    # calculate the variances of the treatment group
    treatment_var <- wtd.var(treatment_data, treatment_weights)

    # calculate the standard deviations of the treatment groups
    treatment_sd <- sqrt(treatment_var)

    # calculate the standardized difference of means
    sd_diff <- (treatment_mean - comparison_mean) / treatment_sd

    # take the absolute value of the standardized difference of means
    abs_sd_diff <- abs(sd_diff)

    # save the absolute standardized difference of means in the ASAM_list vector
    ASAM_list[i] <- abs_sd_diff
  }


  # calculate the mean of the absolute standardized differences of means
  ASAM <- mean(ASAM_list)

  ret <- c(
    Std_In_Bias = Std_In_Bias,
    Prob_Treat = Prob_Treat,
    ATE = ATE,
    ATE_se = ATE_se,
    mean_ps_weights = mean_ps_weights,
    ASAM = ASAM,
    p_val = p_val,
    ci_95 = ci_95
  )
  ret
}
