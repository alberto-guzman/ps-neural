





library(tibble)
library(dplyr)
library(simhelpers)
library(rsimsum)

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


result1 <- combine_results("~/Projects/inProgress/2018_propensity_neuralnet_paper/data/sim_results_n10000_r1000_P_e_all/")
result2 <- combine_results("~/Projects/inProgress/2018_propensity_neuralnet_paper/data/sim_results_n10000_r1000_NP_e_all/")


res <- bind_rows(result1,result2)
res <- res %>% mutate(true_para = 0.3)


test <- res %>%
  group_by(p,method) %>% # grouping 
  do(calc_absolute(., estimates = ATE, true_param = true_para)) # run the function

s1 <- simsum(data = res, estvarname = "ATE", true = 0.3, se = "ATE_se", methodvar = "method", by = c("p"), ref = "logit", dropbig = F)
summary(s1)
