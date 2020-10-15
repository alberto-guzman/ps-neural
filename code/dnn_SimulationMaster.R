########################
## DNN Paper
## Alberto Guzman-Alvarez
########################

#devtools::install_github("rstudio/reticulate")

list.of.packages <- c(
  "Matching", "rpart", "randomForest", "gbm", "twang", "ipred", "neuralnet",
  "nnet", "e1071", "klaR", "xtable", "flexmix", "AUC", "Hmisc", "Kendall",
  "lattice", "keras", "clusterGeneration", "PoisBinOrdNor", "devtools",
  "matrixcalc", "tensorflow", "dplyr", "reticulate"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages) > 0) {
  install.packages(new.packages)
}

lapply(list.of.packages, require, character.only = T)

update.packages(ask = FALSE)


tf <- import("tensorflow")
tf$constant("Hellow Tensorflow")

os <- import("os")
os$listdir(".")