#!/usr/bin/env bash
# Full-scale validation run: 2 replications of every design cell at the real
# n = 10,000 across all p, exercising the ACTUAL runner scripts (sed-derived
# copies that change only replications, worker count, and output names).
# GBM is sharded by p into three concurrent processes; P and NP run alongside
# (5 processes total; each SimDesign cluster capped at 2 workers to match the
# 2 replications). All validation output uses *_val names so nothing collides
# with a real run. Expected wall: ~2.5-4h, dominated by the gbm p=200 shard.
set -euo pipefail
cd "$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export TF_USE_LEGACY_KERAS=1

mkdir -p aws_validation data
V=aws_validation

# --- derive validation copies of the real runners ---
sed 's/replications = 1000/replications = 2/;
     s/parallel = TRUE,/parallel = TRUE,\n  ncores = 2,/;
     s/sim_results_v2_P.rds/sim_results_v2_P_val.rds/;
     s|"data/sim_results_v2_P"|"data/sim_results_v2_P_val"|' \
  code/run_sim_P.R > $V/run_P_val.R

for PP in 20 100 200; do
  sed "s/replications = 1000/replications = 2/;
       s/parallel = TRUE,/parallel = TRUE,\n  ncores = 2,/;
       s/p = c(20, 100, 200)/p = c(${PP})/;
       s/sim_results_v2_gbm.rds/sim_results_v2_gbm_val_${PP}.rds/;
       s|\"data/sim_results_v2_gbm\"|\"data/sim_results_v2_gbm_val_${PP}\"|" \
    code/run_sim_gbm.R > $V/run_gbm_val_${PP}.R
done

sed 's/replications = 1000/replications = 2/;
     s/sim_results_v2_NP.rds/sim_results_v2_NP_val.rds/;
     s|"data/sim_results_v2_NP"|"data/sim_results_v2_NP_val"|' \
  code/run_sim_NP.R > $V/run_NP_val.R

# --- run: P + three gbm shards + NP concurrently; fail if any process fails ---
echo "launching 5 validation processes: $(date)"
declare -A PIDS
Rscript $V/run_P_val.R       > $V/P.log      2>&1 & PIDS[P]=$!
Rscript $V/run_gbm_val_20.R  > $V/gbm20.log  2>&1 & PIDS[gbm20]=$!
Rscript $V/run_gbm_val_100.R > $V/gbm100.log 2>&1 & PIDS[gbm100]=$!
Rscript $V/run_gbm_val_200.R > $V/gbm200.log 2>&1 & PIDS[gbm200]=$!
Rscript $V/run_NP_val.R      > $V/NP.log     2>&1 & PIDS[NP]=$!

FAIL=0
for k in "${!PIDS[@]}"; do
  if ! wait "${PIDS[$k]}"; then echo "PROCESS FAILED: $k (see $V/$k.log)"; FAIL=1; fi
done
echo "all processes finished: $(date)"
[ "$FAIL" -eq 0 ] || { echo "=== VALIDATION FAILED — inspect logs ==="; exit 1; }

# --- diagnostic report from the five summarised outputs (includes SIM_TIME) ---
Rscript - > $V/validation_report.txt 2>&1 <<'EOF'
suppressMessages({library(tidyverse); library(here)})
setwd(here::here())
read_result <- function(f) tryCatch(readRDS(f), error = function(e) qs2::qd_read(f))
files <- c("sim_results_v2_P_val.rds",
           sprintf("sim_results_v2_gbm_val_%d.rds", c(20, 100, 200)),
           "sim_results_v2_NP_val.rds")
res <- bind_rows(lapply(files, function(f) as_tibble(read_result(f))))
cat("== AWS VALIDATION REPORT ==\n")
cat("design rows completed:", nrow(res), "of 132\n")
cat("methods:", paste(sort(unique(res$method)), collapse = ", "), "\n")
cat("any NA ATE:", any(is.na(res$ATE)), " | any NA ATE_se:", any(is.na(res$ATE_se)), "\n\n")
res |>
  mutate(mins_per_rep = round(as.numeric(SIM_TIME) / 60 / 2, 2)) |>
  group_by(method, p) |>
  summarise(cells = n(), mean_mins_per_rep = round(mean(mins_per_rep), 2),
            max_mins_per_rep = max(mins_per_rep),
            mean_abs_bias = round(mean(abs(ATE - 0.3)), 3),
            max_w = signif(max(max_ps_weight), 3), .groups = "drop") |>
  arrange(method, p) |> print(n = 40)
cat("\nProjected full-run core-hours by method (1000 reps):\n")
res |> mutate(mins = as.numeric(SIM_TIME) / 60 / 2) |>
  group_by(method) |>
  summarise(core_hours_1000reps = round(sum(mins) * 1000 / 60), .groups = "drop") |>
  print(n = 12)
EOF
cat $V/validation_report.txt
echo "=== validation complete; report at $V/validation_report.txt ==="
