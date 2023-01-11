######################################################################
# Load libraries and source functions
######################################################################

packages <- c(
  "here",
  "tidyverse",
  "MASS",
  "Rlab",
  "Matrix",
  "psych",
  "Rlab",
  "rpart",
  "ipred",
  "randomForest",
  "nnet",
  "survey",
  "Hmisc",
  "future",
  "furrr",
  "SimDesign",
  "keras",
  "tensorflow",
  "reticulate"
)



invisible(lapply(packages, library, character.only = TRUE))


# sets working directory to root of R project
here()

######### source functions
source(here("code","01_data_gen_fun.R"))
source(here("code","02_analyse_fun.R"))
source(here("code","03_summarize_fun.R"))

######################################################################
# Generate sim design dataframe
######################################################################

# fully-crossed simulation experiment
Design <- createDesign(
  n = c(10000),
  p = c(20, 100),
  scenarioT = c("A", "B", "C", "D"),
  scenarioY = c("a", "b", "c", "d"),
  method = c("logit", "cart", "bag", "forest", "nnet", "dnn-2", "dnn-3")
)


######################################################################
# Run Simulation
######################################################################

use_virtualenv("/ihome/xqin/alg223/.virtualenvs/r-reticulate/bin/python")


res <- runSimulation(
  design = Design,
  replications = 100,
  generate = Generate,
  analyse = Analyse,
  summarise = Summarise,
  parallel = F
)

#
# res <- runSimulation(
#   design = Design,
#   replications = 1000,
#   generate = Generate,
#   analyse = Analyse,
#   summarise = Summarise,
#   parallel = T,
#   filename = paste0("SimDesign_summary_", lubridate::today()),
#   save_results = paste0("SimDesign_results_", lubridate::today())
# )
#
