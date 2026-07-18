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
