# JEBS Structure Survey (2026-07-19)

Seven verified articles (2020-2025), full texts archived in session scratchpad;
survey run for structural alignment of the ps-neural manuscript. Closest
analogs: Keller 2020 (estimator lineup: IPW via GBM/twang, CBPS, BART, TMLE)
and Suk 2024 (structural template: sim designs + ECLS-K application).

## Modal JEBS structure for a simulation-based methods paper
1. Introduction (prose aims, ends with roadmap; no formal RQ blocks)
2. Setup/Notation, Estimand, and Assumptions
3. Existing methods / estimators to be compared
4. Proposed or focal method(s) [3-4 may merge for pure comparisons]
5. Simulation Study — canonical subsections: Data Generation Process
   (often ENUMERATED STEPS), Conditions/Experimental Design, Estimators,
   Evaluation Measures (formulas displayed), Results by design
   (design-first, then estimator, then metric; headline metric in main
   text, secondary metrics in supplement)
6. Empirical/Real Data Application (Data + Results subsections)
7. Discussion or Conclusions (combined OK)
Back matter: Authors' Note -> Acknowledgments -> Conflicting Interests ->
Funding -> ORCID -> Notes -> References -> bios -> received/accepted dates.

## Near-universal conventions (7/7 or 6/7)
1. UNSTRUCTURED abstract ~120-260 words (mode 150-200) + "Keywords:" line,
   4-9 lowercase keywords.
2. Numbered decimal headings in the published version.
3. NO separate Limitations subsection — signposted paragraphs in
   Discussion ("There are some limitations... First,... Despite these...").
4. Online supplement expected: secondary metrics, full per-condition
   results, proofs, code. Figures outnumber tables for performance results.
5. SAGE declarations block; replication materials as a subsection (Vembye),
   Data Availability section (Gilbert 2025), or GitHub footnote (Suk).

## Real-data illustration: 6/7 include one (ECLS-K x2, NLSY, TIMSS, reading
RCT, Swedish equating data), own top-level section, usually AFTER the
simulation; the exception used a hypothetical computational example.

## Per-article structures: see agent report archived in this file's git
history / session log; key exemplars:
- Vembye et al. 2023: Simulation Study subsections = Data Generation
  Process / Estimators / Experimental Design / Performance Assessment /
  Assumptions / REPLICATION MATERIALS / Results; three enumerated aims.
- Suk 2024: numbered sections; DGP as steps 1-4 per design; |Bias| in main
  text, SD/RMSE in supplemental appendices; GitHub footnote; ECLS-K
  application with sensitivity analyses; ends "Conclusions".
- Wallin & Wiberg 2023: empirical example BEFORE simulation as motivation;
  13 PS specifications; full results in supplement.
- Keller 2020: numbered a-priori expectations before results; case study
  (ECLS-K); appendices with variable list + pedagogy example.
- Pashley & Miratrix 2022 (version of record): back-matter order above;
  159-word abstract, 6 keywords.
- Nguyen & Stuart 2020: metric-first results (balance -> bias -> MSE) — the
  one metric-first organizer; 257-word abstract, 9 keywords.
- Gilbert et al. 2025: compact sim section; Data Availability + R code
  appendix.
