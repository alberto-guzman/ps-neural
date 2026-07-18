library(tidyverse)
library(here)

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
  bind_rows(results_list)
}

# SimDesign's save_results writes per-condition files into these directories
# (created by run_sim_P.R, run_sim_gbm.R, and run_sim_NP.R respectively)
result_p <- combine_results(here("data", "sim_results_v2_P"))
result_gbm <- combine_results(here("data", "sim_results_v2_gbm"))
result_np <- combine_results(here("data", "sim_results_v2_NP"))

res <- bind_rows(result_p, result_gbm, result_np)
saveRDS(res, file = here("data", "res_all_v2.rds"))
