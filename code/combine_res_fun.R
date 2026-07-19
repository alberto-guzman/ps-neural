######################################################################
# COMBINE STEP: merges the per-replication results of the three production
# jobs (run_sim_P.R, run_sim_gbm.R, run_sim_NP.R) into one analysis-ready
# tibble, data/res_all_v2.rds — one row per replication, with the design-cell
# identifiers (p, scenarioT, scenarioY, method) bound to each row. This is
# the file the manuscript (propensity_over.qmd) reads.
#
# Run AFTER all three jobs finish: Rscript code/combine_res_fun.R
# Note: SimDesign writes per-cell files named results-row-<k> — classic RDS
# in older SimDesign versions, qs2 "qdata" format (no extension) in newer
# ones; read_result() handles both.
######################################################################

library(tidyverse)
library(here)

# SimDesign's per-row files are classic RDS in older versions and qs2 qdata
# format in newer ones (which also drop the .rds extension) — read either
read_result <- function(file) {
  tryCatch(readRDS(file), error = function(e) qs2::qd_read(file))
}

combine_results <- function(dir_path) {
  # Get a list of all the "results-row-*" files in the directory
  files <- list.files(dir_path, pattern = "results-row-\\d+", full.names = TRUE)
  stopifnot("no results-row files found in dir_path" = length(files) > 0)

  # Read in the "results" and "condition" dataframes from each file and combine them into a single tibble
  results_list <- lapply(files, function(file) {
    x <- read_result(file)
    bind_cols(as_tibble(x[["results"]]), as_tibble(x[["condition"]]))
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
