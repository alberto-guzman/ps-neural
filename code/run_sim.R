######################################################################
# Load libraries and source functions
######################################################################

# remote these before sending to cluster
library(styler)
library(grkstyle)

# library(neuralnet)
# library(keras)
# library(tensorflow)
# library(tidymodels)
# use_condaenv(condaenv = "r-reticulate", required = TRUE)

packages <- c(
  "here",
  "tidyverse",
  "MASS",
  "Rlab",
  "Matrix",
  "psych",
  "MBESS",
  "Rlab",
  "rpart",
  "ipred",
  "randomForest",
  "nnet",
  "survey",
  "Hmisc",
  "future",
  "furrr",
  "SimDesign"
)


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))


# sets working directory to root of R project
here()

######### source functions
source(here("code", "01_data_gen_fun.R"))
source(here("code", "02_analyse_fun.R"))
source(here("code", "03_summarize_fun.R"))


######################################################################
# Generate sim design dataframe
######################################################################

# fully-crossed simulation experiment
Design <- createDesign(
  n = c(100),
  p = c(20, 100),
  scenarioT = c("A", "D"),
  scenarioY = c("a", "d"),
  method = c("logit", "cart", "bag", "forest", "nnet")
)

######################################################################
# Run Simulation
######################################################################


res <- runSimulation(
  design = Design,
  replications = 1000,
  generate = Generate,
  analyse = Analyse,
  summarise = Summarise,
  parallel = T
)


res <- runSimulation(
  design = Design,
  replications = 1000,
  generate = Generate,
  analyse = Analyse,
  summarise = Summarise,
  parallel = T,
  filename = paste0("SimDesign_summary_", lubridate::today()),
  save_results = paste0("SimDesign_results_", lubridate::today())
)


######################################################################
# Testing neural network
######################################################################


df <- Design[12, ] |>
  Generate()

library(keras)
library(tensorflow)
use_condaenv("r-reticulate")

split <- sample(2, nrow(df), replace = TRUE, prob = c(0.8, 0.2)) # random split of data
train_data <- df[split == 1, ]
test_data <- df[split == 2, ]

# Preprocess data
x_train <- as.matrix(train_data[, grep("^v", names(train_data))]) # select columns that start with "v" for input features
y_train <- as.matrix(train_data[, "T"]) # select column for treatment assignment
x_test <- as.matrix(test_data[, grep("^v", names(test_data))]) # select columns that start with "v" for input features
y_test <- as.matrix(test_data[, "T"]) # select column for treatment assignment

# Define model
p <- ncol(x_train) # number of input features
input_layer <- layer_input(shape = c(p)) # input layer
hidden_layer <- layer_dense(units = 2 * p / 3, activation = "relu")(input_layer) # hidden layer
output_layer <- layer_dense(units = 1, activation = "sigmoid")(hidden_layer) # output layer
model <- keras_model(inputs = input_layer, outputs = output_layer)

# Compile model
model %>% compile(
  optimizer = "adam",
  loss = "binary_crossentropy",
  metrics = c("accuracy")
)

# Fit model
history <- model %>% fit(
  x_train,
  y_train,
  epochs = 10,
  batch_size = 32,
  validation_data = list(x_test, y_test)
)

# Preprocess data
x <- as.matrix(df[, grep("^v", names(df))]) # select columns that start with "v" for input features

# Predict propensity scores on entire dataset
ps <- model %>% predict(x)

# Add propensity scores back to original dataframe
df$ps <- ps[, 1]


#### DNN-2

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

# Fit model
history <- model %>% fit(
  x_train,
  y_train,
  epochs = 10,
  batch_size = 32,
  validation_data = list(x_test, y_test)
)

# Preprocess data
x <- as.matrix(df[, grep("^v", names(df))]) # select columns that start with "v" for input features

# Predict propensity scores on entire dataset
ps <- model %>% predict(x)

# Add propensity scores back to original dataframe
df$ps <- ps[, 1]

#### DNN-3

# Define model
p <- ncol(x_train) # number of input features
input_layer <- layer_input(shape = c(p)) # input layer
hidden_layer1 <- layer_dense(units = 2 * p / 3, activation = "relu")(input_layer) # first hidden layer
hidden_layer2 <- layer_dense(units = 2 * p / 3, activation = "relu")(hidden_layer1) # second hidden layer
hidden_layer3 <- layer_dense(units = 2 * p / 3, activation = "relu")(hidden_layer1) # third hidden layer
output_layer <- layer_dense(units = 1, activation = "sigmoid")(hidden_layer2) # output layer
model <- keras_model(inputs = input_layer, outputs = output_layer)

# Compile model
model %>% compile(
  optimizer = "adam",
  loss = "binary_crossentropy",
  metrics = c("accuracy")
)

# Fit model
history <- model %>% fit(
  x_train,
  y_train,
  epochs = 10,
  batch_size = 32,
  validation_data = list(x_test, y_test)
)

# Preprocess data
x <- as.matrix(df[, grep("^v", names(df))]) # select columns that start with "v" for input features

# Predict propensity scores on entire dataset
ps <- model %>% predict(x)

# Add propensity scores back to original dataframe
df$ps <- ps[, 1]

### early stopping





# Split the data into training and validation sets
split <- sample(2, nrow(df), replace = TRUE, prob = c(0.8, 0.2)) # random split of data
train_data <- df[split == 1, ]
validation_data <- df[split == 2, ]

# Preprocess data
x_train <- as.matrix(train_data[, grep("^v", names(train_data))]) # select columns that start with "v" for input features
y_train <- as.matrix(train_data[, "T"]) # select column for treatment assignment
x_validation <- as.matrix(validation_data[, grep("^v", names(validation_data))]) # select columns that start with "v" for input features
y_validation <- as.matrix(validation_data[, "T"]) # select column for treatment assignment

# Define model
p <- ncol(x_train) # number of input features
input_layer <- layer_input(shape = c(p)) # input layer
hidden_layer <- layer_dense(units = 2 * p / 3, activation = "relu")(input_layer) # hidden layer
output_layer <- layer_dense(units = 1, activation = "sigmoid")(hidden_layer) # output layer
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
  callbacks = list(early_stopping)
)
