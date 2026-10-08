# Claude Code instructions: ML propensity score simulation (JEBS)

README.md has the design, the repo map, and the reproduction steps. This file
holds the standing constraints and the gotchas that cost time. The global
rules apply (R scripts, writing voice); where this file is narrower, it wins.

## Constraints that do not move

- No Common App data in this paper. The author has no permission to use it
  for an empirical application. The DGP calibration may cite summary values
  already public in the dissertation. The co-author owns the real-data
  application section; a placeholder HTML comment sits before Discussion in
  the qmd.
- Target journal is JEBS. The paper is co-authored; names are pending in the
  YAML author block. Voice is we/our throughout.
- Design philosophy: evaluate methods with out-of-the-box practitioner
  configurations. Logit and GBM are literal software defaults (Stata
  `teffects ipw`, WeightIt); the trees are literature configurations with
  out-of-bag prediction; Super Learner and the neural networks are the
  informed end. The one stated deviation is GBM iteration selection on a
  single 20% validation split instead of 5-fold CV, and the design-philosophy
  paragraph names it.
- Fixed populations: `Generate()` draws population parameters once per cell
  under a deterministic seed with the RNG kind pinned to Mersenne-Twister,
  then calibrates the intercept to a 0.50 treatment rate. Do not revert to a
  per-replication redraw.
- Primary SEs, p-values, and CIs come from the svyglm sandwich; the
  weighted-lm SE is kept as a comparison column. Trimmed weights (1st/99th
  percentile) are a sensitivity, not the primary.
- ASAM uses the unweighted treated-SD denominator (cobalt convention), and
  max_ASAM is reported alongside it.
- Declined by the author and not to be re-proposed: a standalone LASSO-logit
  arm, an AIPW doubly-robust arm, xgboost in place of GBM, per-replication
  population redraw, balance-aware NN training (dropout was tested and
  rejected in 2023, commit e8048e5), and a population-seed robustness
  appendix until a reviewer asks for one.

## Running the simulation

- The three jobs (P, gbm, NP) share no state and run from the repo root.
  `Rscript code/audit_dgp_properties.R` should print PASS on any new machine.
- One thread per worker: gbm `n.cores = 1`, dbarts `n.threads = 1`, OMP and
  OPENBLAS threads set to 1 in the Slurm scripts.
- SimDesign leaves `SIMDESIGN-TEMPFILE_*` files that abort a relaunch; delete
  them first. `save_results_dirname` must be relative. Newer SimDesign writes
  extension-less qs2 row files; `combine_res_fun.R` handles both formats.
- Keras: `TF_USE_LEGACY_KERAS=1`, tensorflow-cpu 2.16.2, tf-keras 2.16.0,
  protobuf below 5, Python 3.11. `infra/aws/bootstrap.sh` encodes all of it.
  The NN arm runs in PSOCK workers without a GPU (7 to 14 seconds per fit at
  n = 10,000).
- AWS sizing: 64 workers on a 128 GB instance OOM-thrashed on the pilot; the
  512 GB r7a.16xlarge ran clean. GBM p = 200 cells run about three times the
  2-rep extrapolation. Keep the dead-man shutdown and the S3 heartbeat in
  `run_pilot.sh`. Pilot SIM_TIME is wall time; multiply by workers when
  budgeting core-hours.

## Manuscript

- Every citation is DOI-verified before it enters `references.bib`; notes go
  in `research/citations/`.
- Results and Discussion are written from v2 evidence. The dissertation's
  logit-bias claim was a DGP artifact and the corrected implementation
  supersedes it.
- Keep method rankings for the corner cell (complex treatment, complex
  outcome, p = 200) out of the text until the population-seed check has run.
