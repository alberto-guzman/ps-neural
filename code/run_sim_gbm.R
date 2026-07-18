######################################################################
# Canonical runner: GBM only — split into its own Slurm job because 5-fold
# CV of 10,000 trees makes it the longest block (~30 min/rep at p=200)
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
  "gbm",
  "glmnet",
  "nnet",
  "ranger",
  "xgboost",
  "SuperLearner",
  "WeightIt",
  "dbarts",
  "survey",
  "Hmisc",
  "SimDesign"
)

lapply(packages, library, character.only = TRUE)

# sets working directory to root of R project (save_results paths are cwd-relative)
setwd(here())

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
  method = c("gbm")
)

######################################################################
# Run Simulation
######################################################################

# SimDesign does not create parent directories for save_results
dir.create(here("data"), showWarnings = FALSE)

res <- runSimulation(
  design = Design,
  replications = 1000,
  generate = Generate,
  analyse = Analyse,
  summarise = Summarise,
  parallel = TRUE,
  packages = packages,
  seed = 44000 + seq_len(nrow(Design)),
  filename = "sim_results_v2_gbm.rds",
  save_results = TRUE,
  save_details = list(save_results_dirname = "data/sim_results_v2_gbm")
)
