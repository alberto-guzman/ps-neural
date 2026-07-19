# liJSSM2025 — PS estimation methods for IPTW with complex survey data

Li, L., Yang, C., Hu, L., Zhang, W., Aldridge, M., Liu, B., & Mazumdar, M.
(2025). Comparative effectiveness of propensity score estimation methods for
IPTW analysis with complex survey data: A simulation study. *JSSM*.
DOI: 10.1093/jssam/smaf003

**Verified:** DOI via Crossref (Run 1 agent, 2026-07-18); full text read via
PMC (PMC12721855, 2026-07-19).

**Key facts:** compares logit, CBPS, CART, RF, GBM, SL for IPTW across 12
scenarios; survey-weights niche, no deep NNs, moderate p. GBM fit with gbm
package "suggested default values" — exact hyperparameters NOT reported
(their reproducibility gap; cited in our Method as the contrast to our
fully-reported configuration). RF best under poor overlap in their study.

**Used for:** Method (design-philosophy paragraph — reporting-gap contrast);
Lit Review (nearest methodological neighbor — cite & differentiate on
high-p + DNNs + education calibration + coverage/SE analysis).
