######################################################################

######################################################################

# remote these before sending to cluster
library(styler)
library(grkstyle)

# library(neuralnet)
# library(keras)
# library(tensorflow)
# library(tidymodels)
# use_condaenv(condaenv = "r-reticulate", required = TRUE)

packages <- c(
  "here",
  "tidyverse",
  "MASS",
  "Rlab",
  "Matrix",
  "psych",
  "MBESS",
  "Rlab",
  "rpart",
  "ipred",
  "randomForest",
  "nnet",
  "survey",
  "Hmisc",
  "future",
  "furrr"
)


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE))


# sets working directory to root of R project
here()

######### load functions
source(here("code", "data_gen_fun.R"))
source(here("code", "psmethod_fun_sandbox.R"))


######### SIMULATION STUDY  ############


n <- c(1000)
p <- c(20, 100)
nrep <- 1:20
scenarioT <- c("A", "B", "C", "D")
scenarioY <- c("a", "b", "c", "d")
method <- c("logit", "cart", "bag", "forest")

conditions <- crossing(
  n,
  p,
  nrep,
  scenarioT,
  scenarioY,
  method
)

######### run simulation
ncores <- parallelly::availableCores() - 1
plan(multisession, workers = ncores, gc = T)

results <-
  conditions %>%
  group_by(n, p, nrep, scenarioT, scenarioY) %>%
  mutate(datasets = pmap(list(n, p, nrep, scenarioT, scenarioY), possibly(generate_data, NA)))

results$values <- future_map2(results$datasets, results$method, possibly(ps_methods, NA), .options = furrr_options(seed = 123))








results_tidy <- results %>%
  dplyr::select(-(datasets)) %>%
  separate(values, c("ATT", "ATT_se", "Bias", "AbsBias", "mean_ps_weights", "ci_95", "ASAM"), sep = ",")

results_tidy$ATT <- substring(results_tidy$ATT, 3)
results_tidy$ASAM <- substring(results_tidy$ASAM, 1, nchar(results_tidy$ASAM) - 1)

results_tidy <- results_tidy %>%
  mutate_at(c("ATT", "ATT_se", "Bias", "AbsBias", "mean_ps_weights", "ci_95", "ASAM"), as.numeric)









results_summary <- results_tidy %>%
  group_by(method, p, scenarioT, scenarioY) %>%
  summarise(AbsBias = mean(AbsBias, na.rm = T))


# New facet label names for dose variable
t.labs <- c("Base", "Interactions", "Quad Terms", "Complex")
names(t.labs) <- c("A", "B", "C", "D")

# New facet label names for supp variable
y.labs <- c("Base", "Interactions", "Quad Terms", "Complex")
names(y.labs) <- c("a", "b", "c", "d")


results_summary %>%
  ggplot(aes(x = method, y = AbsBias,  fill = as.factor(p))) +
  xlab("Method") +
  ylab("AbsBias") +
  geom_bar(position = "dodge", stat = "identity") +
  facet_grid(scenarioT ~ scenarioY, labeller = labeller(scenarioT = t.labs, scenarioY = y.labs), scales = "fixed") +
  scale_fill_discrete(name = "# Covars") +
  theme(legend.position = "top") +
  theme_minimal()
