# kim2026tradeoff — Classification-causal tradeoff in NN PS estimation

Kim, S., Lee, J., & Jung, K. (2026). On the classification-causal tradeoff in
neural network propensity score estimation. *Stats*, 9(2), Article 37.
DOI: 10.3390/stats9020037

**Verified:** Run 1 agent (Crossref + Semantic Scholar, 2026-07-18);
re-verified via Crossref API at bib entry (2026-07-19). MDPI landing page
blocks automated fetch — full-text details from indexed summaries.

**Key facts:** Monte Carlo, 36 conditions, DNN/CNN vs logistic regression
(only comparator), high-dimensional emphasis; NNs cut bias/imbalance where
logit worsened them; thesis: aggressive classification optimization pushes PS
toward 0/1, destabilizing IPW ("stability-aware training", e.g., restricted
epochs; PropensityNet = 5 dense layers, 30% dropout, TF/Keras).

**Used for:** Intro/Lit Review (nearest DNN-for-PS neighbor — cite &
differentiate on comparator breadth, education calibration, coverage/SE
analysis); Discussion (their tradeoff framing supports our depth gradient).
