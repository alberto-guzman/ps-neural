# Read in `data` data
head(data)

# convert to factor then numeric
# data[, 5] <- as.factor(data[, 5])
# data[, 11] <- as.numeric(data[, 11])

str(data)



# Turn `data` into a matrix
data <- as.matrix(data)

# Set data `dimnames` to `NULL`
dimnames(data) <- NULL


# Normalize the `data` data with keras
# data <- normalize(data[, 1:10]
# data <- normalize(data[, 5]-1

# drop the inf var
data <- data[, -13]
# drop outcome
data <- data[, -12]

# Return the summary of `data`
summary(data)
head(data)
# set seed
set.seed(5)


# Determine sample size
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3))

# Split the `data` data
data.training <- data[ind == 1, 1:10]
data.test <- data[ind == 2, 1:10]

# Split the class attribute
data.trainLabels <- data[ind == 1, 11]
data.testLabels <- data[ind == 2, 11]



# code broke here....

# One hot encode training target values, wont need this for my code since outcome is binary
# data.trainLabels <- to_categorical(data.trainingtarget)

# One hot encode test target values
# data.testLabels <- to_categorical(data.testtarget)

# Print out the data.testLabels to double check the result
print(data.testLabels)



########## KERAS MODELING



# Initialize a sequential model
model <- keras_model_sequential()

# Add layers to the model
model %>%
     layer_dense(units = 8, activation = "relu", input_shape = c(11)) %>%
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
     data.training,
     data.trainLabels,
     # how many times your run epoch
     epochs = 200,
     batch_size = 5,
     validation_split = 0.2
)


# Store the fitting history in `history`
history <- model %>% fit(
     data.training,
     data.trainLabels,
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