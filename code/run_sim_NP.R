######################################################################
# PRODUCTION RUNNER 3 of 3: the "NP" (non-parallel) keras job.
#
# Runs the three neural network methods (nn-1, dnn-2, dnn-3) over the
# full crossed design (36 design rows x 1,000 replications). Historically
# serial because TensorFlow cannot survive FORKED workers; note the 2026-07
# AWS pilot demonstrated TF runs fine in SimDesign's fresh-process PSOCK
# workers (parallel = TRUE, each worker owns its own TF session), which cuts
# wall time from ~87 serial hours to a few — validate with a 2-replication
# cluster test before switching parallel to TRUE here.
#
# Environment: needs a python with TensorFlow reachable by reticulate.
# On Pitt CRC the virtualenv below is used if present; elsewhere, build the
# environment with infra/aws/bootstrap.sh (which encodes the hard-won
# version pins: python 3.11, tensorflow-cpu 2.16.2, tf-keras 2.16.0,
# protobuf < 5) and export TF_USE_LEGACY_KERAS=1. GPU is NOT required —
# these are small MLPs (~7-14s per replication on one CPU core).
#
# Outputs: data/sim_results_v2_NP/results-row-*  and sim_results_v2_NP.rds
# Combine with code/combine_res_fun.R after all three jobs finish.
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
  method = c("nn-1", "dnn-2", "dnn-3")
)

######################################################################
# Run Simulation
######################################################################

# Pitt CRC virtualenv if present; elsewhere (e.g., AWS) reticulate discovers
# the environment set up by keras::install_keras() on its own
crc_venv <- "/ihome/xqin/alg223/.virtualenvs/r-reticulate"
if (dir.exists(crc_venv)) use_virtualenv(crc_venv)
# use_condaenv("r-reticulate")

# SimDesign does not create parent directories for save_results
dir.create(here("data"), showWarnings = FALSE)

res <- runSimulation(
  design = Design,
  replications = 1000,
  generate = Generate,
  analyse = Analyse,
  summarise = Summarise,
  parallel = FALSE,
  packages = packages,
  seed = 43000 + seq_len(nrow(Design)),
  filename = "sim_results_v2_NP.rds",
  save_results = TRUE,
  save_details = list(save_results_dirname = "data/sim_results_v2_NP")
)
