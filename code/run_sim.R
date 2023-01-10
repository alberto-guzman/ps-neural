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
source(here("code", "data_gen_fun.R"))
source(here("code", "psmethod_fun_sandbox.R"))
source(here("code", "summarize_fun.R"))


######################################################################
# Generate sim design dataframe
######################################################################

# fully-crossed simulation experiment
Design <- createDesign(
  n = c(10000),
  p = c(20, 100),
  scenarioT = c("A", "D"),
  scenarioY = c("a", "d"),
  method = c("logit")
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
  parallel = T,
  filename = paste0("SimDesign_summary_", lubridate::today()),
  save_results = paste0("SimDesign_results_", lubridate::today())
)
