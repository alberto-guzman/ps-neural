# Evaluating Modern Propensity Score Estimation Methods with High-Dimensional Data

Monte Carlo simulation comparing ten propensity score estimation methods —
logistic regression, CART, bagged CART, random forest, gradient boosted trees,
BART, Super Learner, and neural networks of one to three hidden layers — for
inverse-probability-of-treatment weighting (IPW) estimation of the ATE with
high-dimensional data calibrated to education administrative records
(p = 20/100/200 covariates, n = 10,000, 120 design cells x 1,000 replications).

Manuscript: `propensity_over.qmd` (Quarto; renders with `references.bib` +
`apa.csl`). Author: Alberto Guzman-Alvarez. Derived from the author's
dissertation (University of Pittsburgh); the implementation here supersedes
the dissertation's (see the revision notes in `code/01_data_gen_fun.R` and the
git history for a complete record of corrections).

## Repository map

```
propensity_over.qmd        the manuscript (reads data/res_all_v2.rds)
references.bib, apa.csl    bibliography (every entry DOI-verified; notes in research/citations/)
code/
  01_data_gen_fun.R        Generate(): the data-generating process (see header)
  02_analyse_fun.R         Analyse(): the ten PS methods + IPW estimation,
                           sandwich SEs, trimming sensitivity, diagnostics
  03_summarize_fun.R       Summarise(): per-cell performance metrics (glossary in header)
  run_sim_P.R              production job 1: logit/cart/bag/forest/bart/sl  (~10h on 64 cores)
  run_sim_gbm.R            production job 2: GBM alone                      (~9h on 64 cores)
  run_sim_NP.R             production job 3: keras neural networks          (see header)
  combine_res_fun.R        merges the three jobs into data/res_all_v2.rds
  audit_dgp_properties.R   rerunnable empirical correctness tests for Generate()
submit_dnn_R_P.slurm       Slurm submissions for the three jobs
submit_gbm.slurm             (update the module loads for your cluster)
submit_dnn_R_NP.slurm
infra/aws/
  bootstrap.sh             full environment build for Ubuntu 24.04 (r2u + pinned
                           TensorFlow stack; encodes every version pitfall we hit)
  run_validation.sh        2-replication full-scale shakedown of the real runners
  run_pilot.sh             tiered pilot (100 reps; 50 for GBM)
  validation_2026-07-18/   archived validation report and logs (measured timings)
research/                  literature-refresh log, per-citation verification
                           notes, v1-vs-v2 comparison figures
example_sim_code/          reference implementations from prior literature
                           (incl. Cannas & Arpino 2019 replication code)
```

## Reproducing the simulation

1. **Environment.** R >= 4.5 with the packages listed at the top of
   `code/run_sim_P.R`, plus a TensorFlow-backed keras for the NP job.
   `infra/aws/bootstrap.sh` builds everything from a bare Ubuntu 24.04 machine
   and documents the required version pins (Python 3.11 — R's embedded
   interpreter lacks `distutils` shims on 3.12; `tensorflow-cpu==2.16.2`;
   `tf-keras==2.16.0`; `protobuf<5`; `TF_USE_LEGACY_KERAS=1` at runtime).
2. **Verify the DGP** (recommended on any new machine):
   `Rscript code/audit_dgp_properties.R` — every check should print PASS.
3. **Run the three jobs** from the repository root (they share no state and
   may run concurrently; identical populations are guaranteed by cell-derived
   seeds inside `Generate()`):
   `Rscript code/run_sim_P.R`, `Rscript code/run_sim_gbm.R`,
   `Rscript code/run_sim_NP.R` — or submit the corresponding Slurm scripts.
   Total compute ~1,300 core-hours, ~40% of it GBM.
4. **Combine:** `Rscript code/combine_res_fun.R` → `data/res_all_v2.rds`.
5. **Render:** `quarto render propensity_over.qmd`.

## Design decisions worth knowing before modifying anything

- **Fixed populations.** Each of the 12 design cells (p x treatment-model x
  outcome-model) defines ONE population: correlation matrix, covariate roles,
  coefficients, complexity terms, and calibrated intercept are drawn once
  under a cell-derived seed (RNG kind pinned so parallel and serial jobs
  agree), then only data are redrawn per replication. Coverage is therefore a
  true repeated-sampling property.
- **50% treatment rate by calibration.** The treatment-model intercept is
  root-found per population so E[e(X)] = 0.50 in every cell.
- **Methods as practitioners encounter them.** Logit and BART are the literal
  software defaults (Stata `teffects`, WeightIt, MatchIt — sources audited in
  the manuscript); GBM uses the MatchIt/WeightIt shared default
  hyperparameters with ONE stated deviation (iteration selection on a single
  20% validation split instead of MatchIt's 5-fold CV — a ~5x feasibility
  saving mirroring the networks' early stopping); CART and random forest
  match MatchIt's implementations (including out-of-bag prediction for
  forest); bagged CART follows Lee et al. (2010); Super Learner and the
  networks are the deliberately-configured tier. No learner is tuned on the
  balance metrics used for evaluation.
- **Weights are used untrimmed as the primary analysis** (a numerical bound at
  1e-6 only); 1st/99th-percentile truncation is reported as a labeled
  remedial sensitivity. All inference treats estimated weights as fixed
  (the field's default), which the manuscript evaluates rather than assumes.
