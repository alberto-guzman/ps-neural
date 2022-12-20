######################################################################
# This script:
# - 
######################################################################

# Load packages, managed by renv  --------------------------------------------------------------------------------------------------------------------------------------------------------------------


packages <- c(
  "here","tidyverse"
)

invisible(lapply(packages, library, character.only = TRUE))

# sets working directory to root of R project 
here()

######### load functions
source(here("dnn_code_clust","dnn_datagen.R"))
source(here("dnn_code_clust","dnn_psmethod_tidy.R"))



######### SIMULATION STUDY  ############

######### set simulation conditions


size <- c(3000)
num_variables <- c(10)
nrep <- 1:3
scenarioT <- c("A", "D")
scenarioY <- c("a", "d")
method <- c("logit", "dnn")

conditions <- crossing(
  size,
  num_variables,
  nrep,
  scenarioT,
  scenarioY,
  method
)

conditions

######### run simulation
ncores <- parallelly::availableCores() - 1
plan(multisession, workers = ncores, gc = T)
ncores


tic()

results <-
  conditions %>%
  group_by(size, num_variables, nrep, scenarioT, scenarioY) %>%
  mutate(datasets = pmap(list(size, num_variables, nrep, scenarioT, scenarioY), possibly(funcov, NA)))

results$values <- future_map2(results$datasets, results$method, possibly(funsim, NA), .options = furrr_options(seed = 123))


results_tidy <- results %>%
  dplyr::select(-(datasets)) %>%
  separate(values, c("rmse", "hatg", "absrbias", "varhatg", "hatgsew", "covw"), sep = ",")

results_tidy$rmse <- substring(results_tidy$rmse, 3)
results_tidy$covw <- substring(results_tidy$covw, 1, nchar(results_tidy$covw) - 1)

results_tidy <- results_tidy %>%
  mutate_at(c("rmse", "hatg", "absrbias", "varhatg", "hatgsew", "covw"), as.numeric)



 results_summary <- results_tidy %>%
   group_by(num_variables,scenarioT, scenarioY, method) %>%
   summarise(mean_rmse = mean(rmse, na.rm = T))

 # results_summary %>%
 #   ggplot (aes (x = as.factor(num_variables), y = mean_mse, group = method, color = method, shape = method)) +
 #   geom_line () + geom_point() + facet_grid (scenarioT ~  scenarioY)

 # New facet label names for dose variable
 t.labs <- c("Base", "Interactions", "Quad Terms","Iteractions + Quad")
 names(t.labs) <- c("A", "B", "C","D")

 # New facet label names for supp variable
 y.labs <- c("Base", "Interactions", "Quad Terms","Iteractions + Quad")
 names(y.labs) <- c("a", "b", "c","d")

 toc()

 results_summary %>%
   ggplot (aes (x = method, y = mean_rmse, fill = as.factor(num_variables))) +
   xlab('Method') + ylab('RMSE') +
   geom_bar (position = 'dodge', stat ='identity') + facet_grid (scenarioT ~  scenarioY,labeller = labeller(scenarioT = t.labs, scenarioY = y.labs), scales = 'fixed') +
   scale_fill_discrete(name = "Number of Covars") +
   theme(legend.position="top") +
   theme_minimal()
toc()
