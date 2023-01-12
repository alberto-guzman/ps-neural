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



lapply(packages, library, character.only = TRUE)

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
  n = c(100000),
  p = c(20, 100, 200),
  scenarioT = c("A", "B", "C", "D"),
  scenarioY = c("a", "b", "c", "d"),
  method = c("logit","cart","bag","forest")
)


######################################################################
# Run Simulation
######################################################################

#use_virtualenv("/ihome/xqin/alg223/.virtualenvs/r-reticulate")


res <- runSimulation(
  design = Design,
  replications = 1000,
  generate = Generate,
  analyse = Analyse,
  summarise = Summarise,
  parallel = T,
  filename = "sim_results_2.rds"
)
