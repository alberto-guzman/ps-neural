######################################################################
# PRODUCTION RUNNER 2 of 3: GBM only.
#
# GBM lives in its own Slurm job because even under validation-split
# selection (see 02_analyse_fun.R) its 10,000-tree fits are the study's
# largest single block: ~5.2 min/rep at p = 200 (measured 2026-07-19),
# ~550 core-hours for its 12 cells x 1,000 replications (~9h on a
# 64-core node).
#
# Populations are identical to the other jobs (cell-derived seeds in
# Generate()). Submit with submit_gbm.slurm from the repo root.
#
# Outputs: data/sim_results_v2_gbm/results-row-*  and sim_results_v2_gbm_val_20.rds
# Combine with code/combine_res_fun.R after all three jobs finish.
######################################################################

packages <- c(
  "here",
  "tidyverse",
  "MASS",
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
  p = c(20),
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
  replications = 2,
  generate = Generate,
  analyse = Analyse,
  summarise = Summarise,
  parallel = TRUE,
  ncores = 2,
  packages = packages,
  seed = 44000 + seq_len(nrow(Design)),
  filename = "sim_results_v2_gbm_val_20.rds",
  save_results = TRUE,
  save_details = list(save_results_dirname = "data/sim_results_v2_gbm_val_20")
)
