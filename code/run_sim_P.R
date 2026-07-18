######################################################################
# Canonical runner: parallel-safe methods (logit, cart, bag, forest)
# Full crossed design; keras-based methods run via run_sim_NP.R
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
  "future",
  "furrr",
  "SimDesign"
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
  method = c("logit", "cart", "bag", "forest", "gbm", "bart", "sl")
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
  seed = 42000 + seq_len(nrow(Design)),
  filename = "sim_results_v2_P.rds",
  save_results = TRUE,
  save_details = list(save_results_dirname = here("data", "sim_results_v2_P"))
)
