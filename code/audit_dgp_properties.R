suppressMessages({library(tidyverse); library(MASS); library(psych); library(SimDesign)})
setwd("~/Projects/inProgress/2018_propensity_neuralnet_paper")
source("code/01_data_gen_fun.R")
ok <- function(name, cond) cat(sprintf("%-62s %s\n", name, ifelse(isTRUE(cond), "PASS", "** FAIL **")))

n <- 50000; p <- 200
set.seed(42)
t0 <- Sys.time()
d <- suppressWarnings(Generate(data.frame(n=n, p=p, scenarioT="base_T", scenarioY="base_Y")))
gen_secs <- as.numeric(difftime(Sys.time(), t0, units="secs"))
X <- as.matrix(d[, 1:p])

## replay the population block exactly (same draw order as Generate)
set.seed(1e6 + n + 100*p + 0 + 0, kind="Mersenne-Twister")
corm <- matrix(0,p,p); corm[upper.tri(corm)] <- runif(p*(p-1)/2,-.3,.3)
corm <- corm + t(corm); diag(corm) <- 1
sm <- suppressWarnings(psych::cor.smooth(corm))          # deterministic, no RNG
perm <- sample(p)
n_conf <- floor(p/2); n_out <- floor(p/4)
conf <- perm[1:n_conf]; rel_out <- perm[(n_conf+1):(n_conf+n_out)]; rel_trt <- perm[(n_conf+n_out+1):p]
beta <- runif(p,-.4,.4); alpha <- runif(p,-.2,.3)
cft <- c(conf, rel_trt); cfo <- c(conf, rel_out)
J_T <- floor(length(cft)/2); cvars_T <- sample(cft, J_T)
J_Y <- floor(length(cfo)/2); cvars_Y <- sample(cfo, J_Y)

## A. marginals and type counts
ok("A1 normal block ~N(0,1)", max(abs(colMeans(X[,1:100]))) < .02 && max(abs(apply(X[,1:100],2,sd)-1)) < .02)
ok("A2 uniform block in [0,1], mean .5", min(X[,101:150]) >= 0 && max(X[,101:150]) <= 1 && max(abs(colMeans(X[,101:150])-.5)) < .012)
ok("A3 bernoulli block {0,1}, p=.5", all(X[,151:200] %in% c(0,1)) && max(abs(colMeans(X[,151:200])-.5)) < .012)

## B. correlation fidelity on the latent-faithful normal block
emp <- cor(X[,1:100]); tgt <- sm[1:100,1:100]
ok("B1 correlations track smoothed target (r > .97)", cor(emp[upper.tri(emp)], tgt[upper.tri(tgt)]) > .97)

## C. role exclusions with correct specifications
sfT <- summary(glm(d$T ~ X[, c(cft, rel_out)], family=binomial))$coefficients
zT <- abs(sfT[-1, "z value"]); names(zT) <- c(cft, rel_out)
ok("C1 outcome-only covariates ~null in treatment model (mean |z| < 1.3)",
   mean(zT[as.character(rel_out)]) < 1.3 && mean(zT[as.character(cft)]) > 2)
sfY <- summary(lm(d$Y ~ d$T + X[, c(cfo, rel_trt)]))$coefficients
tY <- abs(sfY[-(1:2), "t value"]); names(tY) <- c(cfo, rel_trt)
ok("C2 treatment-only covariates ~null in outcome model (mean |t| < 1.3)",
   mean(tY[as.character(rel_trt)]) < 1.3 && mean(tY[as.character(cfo)]) > 3)

## D. assignment validity: T | trueps is Bernoulli(trueps)
dec <- cut(d$trueps, quantile(d$trueps, 0:10/10), include.lowest=TRUE)
cal <- tibble(ps=d$trueps, T=d$T, dec=dec) |> group_by(dec) |> summarise(m_ps=mean(ps), m_T=mean(T))
ok("D1 realized T matches trueps within deciles (max gap < .015)", max(abs(cal$m_ps - cal$m_T)) < .015)

## E. constant effect: correctly specified regression recovers 0.3
ok("E1 coef on T in correct outcome spec = .3 (+/- .01)", abs(sfY["d$T","Estimate"] - 0.3) < .01)

## F. beta/alpha recovery: fitted coefs match population values
fitcoef <- sfY[-(1:2), "Estimate"]; names(fitcoef) <- c(cfo, rel_trt)
ok("F1 outcome-model coefficients recover alpha (r > .995)",
   cor(fitcoef[as.character(cfo)], alpha[cfo]) > .995)

## G. calibration determinism + rate
d2 <- suppressWarnings(Generate(data.frame(n=20000, p=p, scenarioT="base_T", scenarioY="base_Y")))
ok("G1 marginal treatment rate = .50 (+/- .01, two draws)", abs(mean(d$T)-.5) < .01 && abs(mean(d2$T)-.5) < .01)

## H. generation cost at worst case
cat(sprintf("%-62s %.1fs (n=50k, p=200; production n=10k will be ~5x cheaper)\n", "H1 Generate() runtime", gen_secs))

## C1 diagnosis
cat("\nC1 decomposition:\n")
cat("  rows in coef table:", nrow(sfT)-1, "of", length(c(cft, rel_out)), "(any dropped => misalignment)\n")
cat("  mean |z| outcome-only:", round(mean(zT[as.character(rel_out)]), 2),
    "| mean |z| treatment-set:", round(mean(zT[as.character(cft)]), 2), "\n")
cat("  share of outcome-only with |z|>2:", round(mean(zT[as.character(rel_out)] > 2), 2), "\n")
cat("  glm converged:", glm(d$T ~ X[, c(cft, rel_out)], family=binomial)$converged, "\n")
