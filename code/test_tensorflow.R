# open radian
update.packages(ask = FALSE)
# installed.packages()
library(tensorflow)
library(reticulate)
library(keras)
library(dplyr)

library(tfruns)
library(tfdeploy)



# force reticulate to be the correct verision
use_condaenv(condaenv = "dnn_p", conda = "auto", required = TRUE)
reticulate::use_condaenv("dnn_p", required = TRUE)


tf <- import("tensorflow")
tf$constant("Hellow Tensorflow")

# Read in MNIST data
# mnist <- dataset_mnist()
# head(mnist)
# dim(mnist)

# Read in `iris` data
iris <- read.csv(url("http://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data"), header = FALSE)
head(iris)

# convert to factor then numeric
iris[, 5] <- as.factor(iris[, 5])
iris[, 5] <- as.numeric(iris[, 5])

glimpse(iris)

names(iris) <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width", "Species")

plot(iris$Petal.Length,
     iris$Petal.Width,
     pch = 21, bg = c("red", "green3", "blue")[unclass(iris$Species)],
     xlab = "Petal Length",
     ylab = "Petal Width"
)

summary(iris)
str(iris)



# Turn `iris` into a matrix
iris <- as.matrix(iris)

# Set iris `dimnames` to `NULL`
dimnames(iris) <- NULL


# Normalize the `iris` data with keras
iris <- normalize(iris[, 1:4]
iris <- normalize(iris[, 5]-1



# Return the summary of `iris`
summary(iris)

# set seed
set.seed(5)


# Determine sample size
ind <- sample(2, nrow(iris), replace = TRUE, prob = c(0.7, 0.3))

# Split the `iris` data
iris.training <- iris[ind == 1, 1:4]
iris.test <- iris[ind == 2, 1:4]

# Split the class attribute
iris.trainingtarget <- iris[ind == 1, 5]
iris.testtarget <- iris[ind == 2, 5]



# One hot encode training target values, wont need this for my code since outcome is binary
iris.trainLabels <- to_categorical(iris.trainingtarget)

# One hot encode test target values
iris.testLabels <- to_categorical(iris.testtarget)

# Print out the iris.testLabels to double check the result
print(iris.testLabels)



########## KERAS MODELING



# Initialize a sequential model
model <- keras_model_sequential()

# Add layers to the model
model %>%
     layer_dense(units = 8, activation = "relu", input_shape = c(4)) %>%
     layer_dense(units = 2, activation = "softmax")

# Print a summary of a model
summary(model)

# Get model configuration
get_config(model)

# Get layer configuration
get_layer(model, index = 1)

# List the model's layers
model$layers

# List the input tensors
model$inputs

# List the output tensors
model$outputs



# Compile the model
model %>% compile(
     loss = "binary_crossentropy",
     optimizer = "adam",
     metrics = "accuracy"
)



# Fit the model
model %>% fit(
     iris.training,
     iris.trainLabels,
     #how many times your run epoch
     epochs = 200,
     batch_size = 5,
     validation_split = 0.2
)


# Store the fitting history in `history`
history <- model %>% fit(
     iris.training,
     iris.trainLabels,
     epochs = 200,
     batch_size = 5,
     validation_split = 0.2
)

# Plot the history
plot(history)

# Prediction
prob <- model %>%
     predict_proba(test)



# Prediction
prob <- model %>%
     predict_classes(test)
table(Predicted = pred, Actual = testtarget)     

# need to put model for entire sample 


cbind(prob, pred, testtarget)


