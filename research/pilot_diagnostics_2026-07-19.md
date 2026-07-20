# Pilot diagnostics — 2026-07-19 (100 reps/cell, all ten methods)

Post-pilot checks run on `infra/aws/pilot_2026-07-19/`. Per-cell table with
MCSE columns: `infra/aws/pilot_2026-07-19/pilot_percell_diagnostics.csv`.

## 1. Coverage failures are under-adjustment (balance-mediated)

Across all non-bag cells: cor(coverage, max_ASAM) = -0.44;
cor(|bias|, max_ASAM) = +0.54.

Forest at p=100 leaves residual max_ASAM 0.24-0.27 (logit 0.04, SL 0.05)
and the per-cell biases range -0.118 to +0.233 with coverage 0.00 in 3 of
4 cells (the aggregate 0.02 masked sign-varying per-population bias).
CART p=100 similar (max_ASAM 0.32). Mechanism for the Discussion: smooth
high-dimensional PS -> residual confounding -> biased with tight SEs.

## 2. Bagging arm: clip-floor mechanism confirmed

All 12 bag cells have max weight EXACTLY 1e6 (= the 1e-6 PS clip guard)
and WARNINGS = 100/100 reps; ESS as low as 8. Granular OOB probabilities
(k/25 with nbagg = 25) emit exact 0s; the clip converts them to 1e6
weights and a handful of units hijack the estimate. Per-cell bias ranges
-1.69 to +0.12. DECISION NEEDED: report as cautionary defaults finding
(recommended; on-message for the defaults-anchored framing) vs. raise the
clip floor as a labeled sensitivity. Diagnose-before-production item.

## 3. Production budget (worker-corrected from pilot SIM_TIME)

Pilot SIM_TIME is wall time under parallel workers; core time = wall x
workers (P methods x34, gbm x24, NP x6). CAVEAT: pilot ran 64 workers on
64 cores, so these include contention inflation (upper bound-ish).

| method | core-h @1000 reps |
|---|---:|
| gbm | 1,224 |
| sl | 822 |
| bag | 371 |
| forest | 150 |
| bart | 82 |
| NNs (3) | 81 |
| cart+logit | 26 |
| **total** | **~2,760** |

vs the shakedown-derived ~1,300 estimate: the 100-rep cells run ~2x the
2-rep extrapolation (fixed overhead amortizes differently + contention).
Pitt planning: gbm-only job ~= 1,224/64 ~= 19h wall on 64 cores; P job
(sl+bag+forest+bart+cart+logit ~= 1,470) ~= 23h; NP ~= 81/6 ~= 13h serial
lanes or trivially parallel. All fit 3-day windows comfortably.

## 4. Queued AWS runs (specs, ready to launch)

**(a) Population-seed sensitivity (un-parked, sentinel cells):** 5
alternative population seeds x 100 reps for: all four p=200 scenario
cells + forest's p=100 cells, methods {logit, sl, forest, dnn-2}. Q: does
the corner's shared -0.08 bias and forest's p=100 bias replicate across
draws, or is it population idiosyncrasy? Drives how production claims are
phrased + whether the appendix enters the paper. Est. ~2-4h on a
c7a/r7a.8xlarge, < $10.

**(b) Bag clip autopsy:** 1 cell (p=20 base x base, worst bias), 50 reps,
storing per-rep clip counts (n units at floor/ceiling, share of weight
mass on clipped units) + a variant with clip floor at 1e-3. Confirms the
mechanism quantitatively for the writeup. Est. < 1h, ~$2. Can share the
instance with (a).

## 5. Other decisions closed this session

- Logit M-estimation SE column: RECOMMENDED add (analysis-side, ties to
  kostouraki2024).
- g-computation arm: RECOMMENDED against (scope paragraph already
  excludes estimator families).
- Bootstrap sub-study: RECOMMENDED as targeted sub-experiment (logit,
  nn-1, dnn-3, sl at p=100/200; ~200 reps x B=200) rather than the
  draft's dangling promise. Decision pending.
- Cross-platform seed check (Pitt vs AWS, 2-3 cells) before production
  venue commit. Decision pending.
