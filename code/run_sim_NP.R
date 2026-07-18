######################################################################
# Canonical runner: keras-based methods (nn-1, dnn-2, dnn-3), serial
# (TensorFlow does not play well with SimDesign's parallel workers)
# Full crossed design; tree/logit methods run via run_sim_P.R
######################################################################

packages <- c(
  "here",
  "tidyverse",
  "MASS",
  "Matrix",
  "psych",
  "rpart",
  "ipred",
  "randomForest",
  "survey",
  "Hmisc",
  "SimDesign",
  "keras",
  "tensorflow",
  "reticulate"
)

lapply(packages, library, character.only = TRUE)

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
  n = c(10000),
  p = c(20, 100, 200),
  scenarioT = c("base_T", "complex_T"),
  scenarioY = c("base_Y", "complex_Y"),
  method = c("nn-1", "dnn-2", "dnn-3")
)

######################################################################
# Run Simulation
######################################################################

use_virtualenv("/ihome/xqin/alg223/.virtualenvs/r-reticulate")
# use_condaenv("r-reticulate")

res <- runSimulation(
  design = Design,
  replications = 1000,
  generate = Generate,
  analyse = Analyse,
  summarise = Summarise,
  parallel = FALSE,
  seed = 43000 + seq_len(nrow(Design)),
  filename = "sim_results_v2_NP.rds",
  save_results = TRUE
)
