#!/usr/bin/env bash
# PILOT run: 100 reps/cell for ALL ten methods — enough for usable coverage
# (±2.2pp MCSE) and SE-calibration reads to draft Results/Discussion; the
# 1000-rep Pitt production run then tightens the numbers. GBM's
# validation-split selection (~5x faster than the July 18 pilot's CV) makes
# 100 GBM reps affordable. NP runs parallel (PSOCK workers tolerate TF —
# proven in the July 18 pilot). Sized for a 64-core box by measured
# core-hours (P 56, gbm 41, NP 5): P 34 workers, gbm 24, NP 6 => ~2h wall.
# Outputs use *_pilot names.
set -euo pipefail
cd "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export TF_USE_LEGACY_KERAS=1

mkdir -p aws_pilot data
V=aws_pilot

sed 's/replications = 1000/replications = 100/;
     s/parallel = TRUE,/parallel = TRUE,\n  ncores = 34,/;
     s/sim_results_v2_P.rds/sim_results_v2_P_pilot.rds/;
     s|"data/sim_results_v2_P"|"data/sim_results_v2_P_pilot"|' \
  code/run_sim_P.R > $V/run_P_pilot.R

sed 's/replications = 1000/replications = 100/;
     s/parallel = TRUE,/parallel = TRUE,\n  ncores = 24,/;
     s/sim_results_v2_gbm.rds/sim_results_v2_gbm_pilot.rds/;
     s|"data/sim_results_v2_gbm"|"data/sim_results_v2_gbm_pilot"|' \
  code/run_sim_gbm.R > $V/run_gbm_pilot.R

# NP: parallel PSOCK (validated in the July 18 pilot)
sed 's/replications = 1000/replications = 100/;
     s/parallel = FALSE,/parallel = TRUE,\n  ncores = 6,/;
     s/sim_results_v2_NP.rds/sim_results_v2_NP_pilot.rds/;
     s|"data/sim_results_v2_NP"|"data/sim_results_v2_NP_pilot"|' \
  code/run_sim_NP.R > $V/run_NP_pilot.R

echo "launching 3 pilot processes: $(date)"
declare -A PIDS
Rscript $V/run_P_pilot.R   > $V/P.log   2>&1 & PIDS[P]=$!
Rscript $V/run_gbm_pilot.R > $V/gbm.log 2>&1 & PIDS[gbm]=$!
Rscript $V/run_NP_pilot.R  > $V/NP.log  2>&1 & PIDS[NP]=$!

FAIL=0
for k in "${!PIDS[@]}"; do
  if ! wait "${PIDS[$k]}"; then echo "PROCESS FAILED: $k (see $V/$k.log)"; FAIL=1; fi
done
echo "all processes finished: $(date)"
[ "$FAIL" -eq 0 ] || echo "=== AT LEAST ONE PROCESS FAILED — report covers survivors ==="

Rscript - > $V/pilot_report.txt 2>&1 <<'EOF'
suppressMessages({library(tidyverse); library(here)})
setwd(here::here())
read_val <- function(f) tryCatch(readRDS(f), error = function(e) qs2::qd_read(f))
files <- Filter(file.exists, c("sim_results_v2_P_pilot.rds",
                               "sim_results_v2_gbm_pilot.rds",
                               "sim_results_v2_NP_pilot.rds"))
res <- bind_rows(lapply(files, function(f) as_tibble(read_val(f))))
cat("== PILOT REPORT (100 reps, all methods) ==\n")
cat("rows:", nrow(res), "of 120 | methods:", paste(sort(unique(res$method)), collapse=", "), "\n")
cat("any NA ATE:", any(is.na(res$ATE)), "\n\n")
res |> mutate(mins_per_rep = round(as.numeric(SIM_TIME)/60/REPLICATIONS, 2)) |>
  group_by(method, p) |>
  summarise(bias = round(mean(ATE - 0.3), 3), emp_SE = round(mean(emp_SE), 3),
            SE_ratio = round(mean(SE_ratio), 2), cover = round(mean(ci_95), 2),
            cover_trim = round(mean(ci_95_trim), 2), ess = round(mean(ess)),
            mins_rep = round(mean(mins_per_rep), 2), .groups = "drop") |>
  arrange(method, p) |> print(n = 40, width = 200)
EOF
cat $V/pilot_report.txt
echo "=== pilot complete ==="
