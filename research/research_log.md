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

### Run 1 results (2026-07-18, all three lanes reported)

Agent-verified via DOI landing pages / Crossref metadata; re-verify each item
once more at the moment it enters references.bib (protocol).

#### Novelty assessment — the headline

**No 2023–2026 study does a high-dimensional, education-calibrated, multi-method
(logit/CART/ensembles/SL/multi-layer NN) IPW-ATE simulation. No 2023–2026 JEBS
article on ML/NN propensity scores was found — good for fit and novelty.**
Two close competitors MUST be cited and differentiated:

1. **kim2026** — Kim, Lee & Jung (2026), "On the classification–causal tradeoff
   in neural network propensity score estimation," *Stats* 9(2):37,
   DOI 10.3390/stats9020037. DNN/CNN vs logit Monte Carlo, 36 conditions,
   high-dimensional. Nearest neighbor. Differentiators: our broad comparator
   set (theirs is logit only), education calibration, coverage/SE focus. Their
   "aggressive classification optimization destabilizes weighting" framing
   SUPPORTS our discussion of why deeper NNs can hurt IPW.
2. **liJSSM2025** — Li et al. (2025), "Comparative effectiveness of PS
   estimation methods for IPTW with complex survey data," *J Survey Stat &
   Methodology*, DOI 10.1093/jssam/smaf003. Same learner menu (logit, CBPS,
   CART, RF, GBM, SL) + IPW, but survey-weights niche, moderate p, no deep
   NNs. RF best under poor overlap. Cite and position on high-p + DNN +
   education.

Also near-neighbor (preprint): **karim2026pre** — Karim & Hu (2026), arXiv
2607.07065, plasmode NHANES p≈167 regularized-PS vs TMLE. Regularization-
focused, health data; mild threat only; flag as preprint if cited.

#### Anchor citations for the revision

- **leite2025** — Leite et al. (2025), *Psych Methods*, DOI 10.1037/met0000789.
  Systematic review, 179 ML-PS applications; GBM most used; severe
  underreporting (48% no balance report, ~23% sensitivity analysis). THE
  framing citation; the revision should visibly follow its reporting
  guidelines (a differentiator JEBS reviewers will value).
- **naimi2023 + balzer2023** — Naimi, Mishler & Kennedy (2023) *AJE* 192(9):
  1536–44, DOI 10.1093/aje/kwab201 (ML plug-ins in singly robust estimators →
  bias + invalid inference; DR + sample splitting fixes) and the Balzer &
  Westling invited commentary, DOI 10.1093/aje/kwab200 (context-specific
  simulations legitimized — directly supports our design). Canonical citations
  for our poor-coverage finding.
- **kostouraki2024** — *Stat Med* 43(13):2672–94, DOI 10.1002/sim.10078.
  M-estimation tutorial for IPTW variance across weight types; explains why
  ignoring PS estimation is conservative for ATE (matches our logit base/base
  over-coverage) but can be anti-conservative for ATT.
- **austin2022** — *Stat Med* 41(22):4426–43, DOI 10.1002/sim.9519. Bootstrap
  vs asymptotic SEs for PS weighting; benchmark for our planned bootstrap arm.

#### Supporting citations by manuscript section

Lit review (high-dimensional framing): karim2025 (*Am Stat* 79(1):72–90, DOI
10.1080/00031305.2024.2368794, hdPS + ML extensions); balde2023 (*Biometrics*
79(1):514–20, GOAL outcome-adaptive lasso; represents the penalized-selection
branch we don't test); NOTE: no 2023–26 high-dimensional CBPS extension exists
(gap worth a sentence; Ning & Imai 2020 Biometrika is still the reference).

Methods justification: phillips2023 (*IJE* 52(4):1276–85, DOI
10.1093/ije/dyad023, SL specification guidance — supports our library/folds);
collier2022 (*JXE* 90(4):1003–20, ANN-PS tutorial — correct year, online
2021); autenrieth2021 (*JEDM* 13(1), stacked ensembles in education — YEAR IS
2021, fix if cited as recent).

Discussion (why ML-IPW misbehaves + remedies): gutman2024 (*Epidemiology*
35(4):473–80, DOI 10.1097/EDE.0000000000001733 — post-calibrating PS before
IPW, biggest gains for tree learners; most on-point calibration cite);
vanderlaan2025 (CLeaR/PMLR 275, arXiv 2411.06342 — isotonic calibration of
weights, principled alternative to trimming); ballinari2026 (*J Applied
Econometrics* early view, DOI 10.1002/jae.70057 — DML + PS calibration; the
arXiv 2409.04874 item, now published); shang2025 (*Stat Methods Med Res*
34(3):457–72 — balance-penalized loss for NN PS); peng2024 (arXiv 2404.04794,
still preprint as of v2 Apr 2026, title changed to "Local balance calibration
for nonparametric PS estimation"; LBCNet).

Discussion (DML as future direction): ahrens2026 (*J Econ Literature* in
press, arXiv 2504.08324 — THE accessible DML intro); moccia2024 (*Eur J
Epidemiol* 39(10):1097–1108 — TMLE/AIPW/DML triad for applied readers);
ellul2025 (*Stat Med* 44:e70272 — AIPW/TMLE + SL in high-dimensional
confounding; cross-fitting mainly improves coverage); bach2024 (*JSS* 108(3),
DoubleML R package); bach2024tuning (CLeaR/PMLR 236 — tuning propagates into
causal estimates; predictive performance imperfect guide).

Discussion (weights/positivity): matsouaka2024 (*Biometrical J* 66(4):2300156
— positivity violations, overlap weights); liu2025 (*Stat Med* 44:e70329 —
tutorial on weighting under positivity violations, ChiPS package; recommends
overlap-type weights over IPW+trimming); rizk2025 (*J Clin Epi* 187:111942 —
overlap weighting changes the estimand to ATO; caution against using it to
rescue ATE); lewis2024 (*Biometrika* 111(3):989–1011 — separation breaks
logistic MLE in high dimensions; statistical grounding for our logit
weight-explosion mechanism); li2025sb (*Statistics in Biosciences* online
first, DOI 10.1007/s12561-025-09503-7 — bootstrap vs sandwich for weighted
ATEs, bootstrap failure modes under poor overlap); kosko2024 (*Stat Med*
43(15):2894–2927 — bag-of-little-bootstraps for IPW at scale; makes our
bootstrap arm computationally feasible); orihara2023 (arXiv 2310.20182,
preprint — when sandwich is conservative vs anti-conservative);
reifeis2022 (*AJE* 191(6):1092–97 — ATT sandwich can be anti-conservative;
stacked estimating equations correction).

Secondary/optional: lourenco2024 (IJERPH 21(11):1484, health-policy scoping
review — redundant with leite2025); wang2026 (*Eur J Epidemiol* 41(3):317–28 —
trial-emulation benchmark where logit beat GBM; "prediction ≠ causal accuracy"
counterweight); yang2023 (BMC Med Res Methodol 23:41 — multi-task NN for PS
with missing data; drop if space-constrained).

#### Excluded / unverified (do not cite without further verification)

Quantum NN PS (arXiv 2506.19973); preprints.org AI+PSM review; Rabenseifner
calibration strategies (arXiv 2503.17290, authors unverified); data-adaptive
covariate balancing (arXiv 2512.18069); AIPW + outcome-adaptive lasso
comparison (arXiv 2405.11522); Venkatesan PDS 2025 dimensionality-reduction
(unfetched); DML bootstrap consistency (arXiv 2604.17239); survival-specific
overlap-weight applications.

#### Corrections to prior beliefs

- Autenrieth JEDM stacked-ensemble paper is **2021**, Collier & Leite JXE
  tutorial is **2022 (online 2021)** — both predate the manuscript's 2023
  cutoff and may already be cited; check before treating as "new."
- vanderLaan2007 Super Learner DOI is 10.2202/1544-6115.1309 (10.1515 is wrong).
