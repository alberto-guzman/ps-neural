########################
## DNN Paper
## Alberto Guzman-Alvarez
## Last Worked on: 090120
########################


list.of.packages <- c(
    "Matching", "rpart", "randomForest", "gbm", "twang", "ipred", "neuralnet",
    "nnet", "e1071", "klaR", "xtable", "flexmix", "AUC", "Hmisc", "Kendall", "lattice",
    "keras", "clusterGeneration", "PoisBinOrdNor", "devtools", "matrixcalc", "tensorflow"
) # replace xx and yy with package names
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages) > 0) {
    install.packages(new.packages)
}
lapply(list.of.packages, require, character.only = T)

library(tensorflow)
library(reticulate)
library(keras)
library(dplyr)
