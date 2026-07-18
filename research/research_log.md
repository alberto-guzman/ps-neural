# Research Log — PS Neural Net Manuscript (JEBS)

Citation protocol (inherited from the CommonApp/JREE workflow): every new
citation must be DOI/landing-page verified before entering references.bib.
Per-citation verification notes live in research/citations/.

## Run 1 — Post-2023 refresh sweep (2026-07-18, in progress)

Three parallel search agents launched:

1. **ML/NN for PS estimation lane** — neural network and deep learning PS
   estimation, tree ensembles / boosting / super learner for PS, systematic
   reviews and best-practice guidance (2023–2026). Also tasked with checking
   whether anything scoops the high-dimensional-DNN-in-education contribution.
2. **High-dimensional / DML lane** — hdPS and ML extensions, regularized PS,
   double/debiased ML, AIPW/TMLE, PS calibration, high-dimensional CBPS.
3. **IPW variance / weights lane** — SE estimation for IPW with estimated PS
   (sandwich vs bootstrap vs influence function), weight trimming/truncation
   guidance, overlap weights, CI coverage in weighting simulations.

Findings independently re-verified before entering references.bib.

### Verified earlier today (pre-agent, during journal targeting)

- **Leite et al. (2025), Psychological Methods** — "Machine Learning for
  Propensity Score Estimation: A Systematic Review and Reporting Guidelines."
  Online ahead of print Oct 2025. 179 applications reviewed across 40 years;
  GBM most used, then RF; 48% fail to report balance; reporting guidelines
  provided. Verified via PubMed (PMID 41100266). Authors: Leite, Zhang,
  Collier, Chawla, Kong, Lee, Quan, Soyoye.
- **Autenrieth, Levine, Fan & Guarcello, Journal of Educational Data Mining**
  — "Stacked Ensemble Learning for Propensity Score Methods in Observational
  Studies" (2021). GBM-meta-learner stack beat individual learners in
  education-calibrated Monte Carlo scenarios; balance measures can mislead
  model selection. Landing page verified (JEDM article 525); full biblio
  details to be confirmed in Run 1.

Awaiting agent reports.
