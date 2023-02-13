





library(tibble)
library(dplyr)
library(simhelpers)
library(rsimsum)
library(ggplot2)

combine_results <- function(dir_path) {
  # Get a list of all the "results-row-*.rds" files in the directory
  files <- list.files(dir_path, pattern = "results-row-\\d+\\.rds", full.names = TRUE)
  
  # Read in the "results" and "condition" dataframes from each file and combine them into a single tibble
  results_list <- lapply(files, function(file) {
    res <- readRDS(file)[["results"]] %>% as_tibble()
    cond <- readRDS(file)[["condition"]] %>% as_tibble()
    bind_cols(res, cond)
  })
  
  # Combine all the resulting tibbles into a single tibble
  combined_results <- bind_rows(results_list)
  
  return(combined_results)
}


res1 <- combine_results("~/Projects/inProgress/2018_propensity_neuralnet_paper/data/sim_results_n10000_r1000_P_e_all_old/")
res2 <- combine_results("~/Projects/inProgress/2018_propensity_neuralnet_paper/data/dnn_n10000_r1000_P_e_0212_all/")
res <- bind_rows(res1,res2)



result1 <- combine_results("~/Projects/inProgress/2018_propensity_neuralnet_paper/data/sim_results_n10000_r1000_P_e_all/")
result2 <- combine_results("~/Projects/inProgress/2018_propensity_neuralnet_paper/data/sim_results_n10000_r1000_NP_e_all/")


res <- bind_rows(result1,result2)
res <- res %>% mutate(true_para = 0.3)

s1 <- simsum(data = res, estvarname = "ATE", true = 0.3, se = "ATE_se", methodvar = "method", by = c("p", "scenarioT", "scenarioY"), ref = "logit", x = TRUE)
summary(s1)
autoplot(summary(s1), type = "forest", stats = "bias")
autoplot(summary(s1), type = "lolly", stats = "bias")

autoplot(s1, type = "zip")
autoplot(s1, type = "est_density")
autoplot(s1, type = "est_ridge") +
  xlim(0,.4)

tidy <- tidy(summary(s1))

tidy <- tidy |> 
  filter(stat == "bias", )


res %>%
  ggplot(aes(x = method, y = ATE, fill = as.factor(p))) + 
  geom_violin(trim = T) +
  ylim(-1,1) +
  stat_summary(fun.data=ATE, mult=1, 
                 geom="pointrange", color="red")
  