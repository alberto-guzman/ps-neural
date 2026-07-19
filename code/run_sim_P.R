######################################################################
# PRODUCTION RUNNER 1 of 3: the "P" (parallel) job.
#
# Runs the six non-keras, non-GBM methods — logit, cart, bag, forest,
# bart, sl — over the full crossed design (p = 20/100/200 x base/complex
# treatment x base/complex outcome = 12 cells per method, 72 design rows),
# 1,000 replications each, parallelized across replications by SimDesign.
#
# Companions: run_sim_gbm.R (GBM alone — its 10,000-tree fits remain the
# study's largest block, so it keeps its own job) and
# run_sim_NP.R (the three keras neural networks). All three jobs simulate
# IDENTICAL populations because the population seed inside Generate()
# depends only on the design cell, never on the method or the job.
#
# Expected cost: ~600 core-hours (~10h on a 64-core node; timings measured
# in the 2026-07 AWS validation, see infra/aws/validation_2026-07-18/).
# Submit on Slurm with submit_dnn_R_P.slurm. Run from the repo root (the
# script setwd()s to the project root via here()).
#
# Outputs: data/sim_results_v2_P/results-row-*  (per-cell replication data)
#          sim_results_v2_P.rds                 (summarised, one row per cell)
# Combine all three jobs' outputs with code/combine_res_fun.R.
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
  p = c(20, 100, 200),
  scenarioT = c("base_T", "complex_T"),
  scenarioY = c("base_Y", "complex_Y"),
  method = c("logit", "cart", "bag", "forest", "bart", "sl")
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
  seed = 42000 + seq_len(nrow(Design)),
  filename = "sim_results_v2_P.rds",
  save_results = TRUE,
  save_details = list(save_results_dirname = "data/sim_results_v2_P")
)
