use_condaenv("r-reticulate")
# Define weight decay values to try
weight_decay_values <- c(0.0001, 0.001, 0.01, 0.1)

# Initialize vector to store bias values
bias_values <- c()



dat_iteration <- dat


for (l in weight_decay_values) {
  
  dat_iteration <- dat
  
  
  # Split the data into training and validation sets (80/20)
  split <- sample(2, nrow(dat_iteration), replace = TRUE, prob = c(0.8, 0.2)) # random split of data
  train_data <- dat_iteration[split == 1, ]
  validation_data <- dat_iteration[split == 2, ]
  

  # Preprocess data
  x_train <- as.matrix(train_data[, grep("^v", names(train_data))]) # select columns that start with "v" for input features
  y_train <- as.matrix(train_data[, "T"]) # select column for treatment assignment
  x_validation <- as.matrix(validation_data[, grep("^v", names(validation_data))]) # select columns that start with "v" for input features
  y_validation <- as.matrix(validation_data[, "T"]) # select column for treatment assignment
  
  
  # Define model
  p <- ncol(x_train) # number of input features
  input_layer <- layer_input(shape = c(p)) # input layer
  hidden_layer1 <- layer_dense(units = ceiling(2 * p / 3), activation = "relu", kernel_regularizer = regularizer_l2(l = l))(input_layer) # first hidden layer
  hidden_layer2 <- layer_dense(units = ceiling(2 * p / 3), activation = "relu", kernel_regularizer = regularizer_l2(l = l))(hidden_layer1) # second hidden layer
  output_layer <- layer_dense(units = 1, activation = "sigmoid", kernel_regularizer = regularizer_l2(l = l))(hidden_layer2) # output layer
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
    batch_size = 64,
    validation_data = list(x_validation, y_validation),
    callbacks = list(early_stopping)
  )
  
  # Preprocess data
  x <- as.matrix(dat_iteration[, grep("^v", names(dat_iteration))]) # select columns that start with "v" for input features
  
  # Predict propensity scores on entire dataset
  ps <- model %>% predict(x)
  ps <- ps[, 1]
  
  
  dat_iteration <- dat_iteration %>%
    mutate(
      ps_pred = ps,
      ps_weights = case_when(T == 1 ~ 1 / ps, T == 0 ~ 1 / (1 - ps))
    )
  
  true_ATE <- 0.3
  
  
  # estimate the true_ATE with the weights
  d.w <- svydesign(~0, weights = dat_iteration$ps_weights, data = dat_iteration)
  fit <- svyglm(Y ~ T, design = d.w)
  
  # save the true_ATE and se_true_ATE
  ATE <- unname(coef(fit)["T"])
  bias_values <- c(bias_values, abs(ATE - true_ATE))
  
}
  