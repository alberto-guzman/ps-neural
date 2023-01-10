#############

## WHAT DOES THIS FUNCTIOND DO?

# this function estimates the ATT and other metrics

#############

# function to estimate the ATT and other metrics
Analyse <- function(condition, dat, fixed_objects = NULL) {
  Attach(condition)

  # if the method is logit, then estimate the ATT using logistic regression
  if (method == "logit") {
    # estimate the propensity score using logistic regression
    mod <- glm(T ~ . - Y - trueps, data = dat, family = binomial)
    # save the propensity score to a vector
    ps <- mod$fitted
    # if the method is cart, then estimate the ATT using classification and regression trees
  } else if (method == "cart") {
    # estimate the propensity score using classification and regression trees
    mod <- rpart(T ~ . - Y - trueps, method = "class", data = dat)
    # save the propensity score to a vector
    ps <- predict(mod)[, 2]
    # if the method is bag, then estimate the ATT using bagging
  } else if (method == "bag") {
    # estimate the propensity score using bagging
    mod <- bagging(T ~ . - Y - trueps, data = dat, nbag = 100)
    # save the propensity score to a vector
    ps <- predict(mod, newdata = dat, type = "prob")
    # if the method is forest, then estimate the ATT using random forest
  } else if (method == "forest") {
    # estimate the propensity score using random forest
    mod <- randomForest(factor(T) ~ . - Y - trueps, ntree = 500, data = dat)
    # save the propensity score to a vector
    ps <- predict(mod, type = "prob")[, 2]
  } else if (method == "nnet") {
    # estimate the propensity score using neural net
    neuro_n <- ceiling((2 / 3) * length(dat))
    mod <- nnet(factor(T) ~ . - Y - trueps, data = dat, size = neuro_n, decay = 0.01, maxit = 2000, trace = F, MaxNWts = 90000)
    ps <- as.numeric(predict(mod, type = "raw"))
  } else if (method == "dnn-2") {
    # Split the data into training and validation sets
    split <- sample(2, nrow(dat), replace = TRUE, prob = c(0.8, 0.2)) # random split of data
    train_data <- dat[split == 1, ]
    validation_data <- dat[split == 2, ]

    # Preprocess data
    x_train <- as.matrix(train_data[, grep("^v", names(train_data))]) # select columns that start with "v" for input features
    y_train <- as.matrix(train_data[, "T"]) # select column for treatment assignment
    x_validation <- as.matrix(validation_data[, grep("^v", names(validation_data))]) # select columns that start with "v" for input features
    y_validation <- as.matrix(validation_data[, "T"]) # select column for treatment assignment

    # Define model
    p <- ncol(x_train) # number of input features
    input_layer <- layer_input(shape = c(p)) # input layer
    hidden_layer1 <- layer_dense(units = 2 * p / 3, activation = "relu")(input_layer) # first hidden layer
    hidden_layer2 <- layer_dense(units = 2 * p / 3, activation = "relu")(hidden_layer1) # second hidden layer
    output_layer <- layer_dense(units = 1, activation = "sigmoid")(hidden_layer2) # output layer
    model <- keras_model(inputs = input_layer, outputs = output_layer)

    # Compile model
    model %>% compile(
      optimizer = "adam",
      loss = "binary_crossentropy",
      metrics = c("accuracy")
    )

    # Define callbacks
    early_stopping <- callback_early_stopping(monitor = "val_loss", min_delta = 0, patience = 5)

    # Fit model
    history <- model %>% fit(
      x_train,
      y_train,
      epochs = 100,
      batch_size = 32,
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
    # Split the data into training and validation sets
    split <- sample(2, nrow(dat), replace = TRUE, prob = c(0.8, 0.2)) # random split of data
    train_data <- dat[split == 1, ]
    validation_data <- dat[split == 2, ]

    # Preprocess data
    x_train <- as.matrix(train_data[, grep("^v", names(train_data))]) # select columns that start with "v" for input features
    y_train <- as.matrix(train_data[, "T"]) # select column for treatment assignment
    x_validation <- as.matrix(validation_data[, grep("^v", names(validation_data))]) # select columns that start with "v" for input features
    y_validation <- as.matrix(validation_data[, "T"]) # select column for treatment assignment

    # Define model
    p <- ncol(x_train) # number of input features
    input_layer <- layer_input(shape = c(p)) # input layer
    hidden_layer1 <- layer_dense(units = 2 * p / 3, activation = "relu")(input_layer) # first hidden layer
    hidden_layer2 <- layer_dense(units = 2 * p / 3, activation = "relu")(hidden_layer1) # second hidden layer
    hidden_layer3 <- layer_dense(units = 2 * p / 3, activation = "relu")(hidden_layer2) # third hidden layer
    output_layer <- layer_dense(units = 1, activation = "sigmoid")(hidden_layer3) # output layer
    model <- keras_model(inputs = input_layer, outputs = output_layer)

    # Compile model
    model %>% compile(
      optimizer = "adam",
      loss = "binary_crossentropy",
      metrics = c("accuracy")
    )

    # Define callbacks
    early_stopping <- callback_early_stopping(monitor = "val_loss", min_delta = 0, patience = 5)

    # Fit model
    history <- model %>% fit(
      x_train,
      y_train,
      epochs = 100,
      batch_size = 32,
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

  # save estimated propensity score and weights to data frame
  dat <- dat %>%
    mutate(
      ps_pred = ps,
      ps_weights = if_else(T == 0, ps / (1 - ps), 1)
    )

  # calculate standardized initial prior to weighting
  Std_In_Bias <- (mean(dat$Y[dat$T == 1]) - mean(dat$Y[dat$T == 0]) - 0.3) / sd(dat$Y[dat$T == 1])
  Prob_Treat <- mean(dat$T)

  # estimate the ATT with the weights
  d.w <- svydesign(~1, weights = dat$ps_weights, data = dat)
  fit <- svyglm(Y ~ T, design = d.w)

  # save the ATT and se_ATT
  ATT <- unname(coef(fit)["T"])
  vcov_matrix <- vcov(fit)
  ATT_se <- unname(sqrt(vcov_matrix["T", "T"]))

  # calculate the bias metrics
  Bias <- ATT - 0.3
  AbsBias <- abs(ATT - 0.3)

  # calculate the mean of control group weights
  dat_int <- subset(dat, T == 0)
  mean_ps_weights <- mean(dat_int$ps_weights, na.rm = TRUE)

  # calculate the 95% coverage
  lower_bound <- ATT - 1.96 * ATT_se
  upper_bound <- ATT + 1.96 * ATT_se
  ci_95 <- ifelse(lower_bound <= 0.3 && 0.3 <= upper_bound, 1, 0)

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

    # calculate the variances of the treatment and comparison groups
    treatment_var <- wtd.var(treatment_data, treatment_weights)
    comparison_var <- wtd.var(comparison_data, comparison_weights)

    # calculate the standard deviations of the treatment and comparison groups
    treatment_sd <- sqrt(treatment_var)
    comparison_sd <- sqrt(comparison_var)

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
    ATT = ATT,
    ATT_se = ATT_se,
    Bias = Bias,
    AbsBias = AbsBias,
    mean_ps_weights = mean_ps_weights,
    ci_95 = ci_95,
    ASAM = ASAM
  )
  ret
}
