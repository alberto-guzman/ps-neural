suppressMessages({library(tidyverse); library(rsimsum); library(patchwork); library(qs2)})
source("~/Projects/inProgress/dissertation/ggplot_theme_pub.R")
outdir <- "~/Projects/inProgress/2018_propensity_neuralnet_paper/research/v1_v2_figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

meth_lvls <- c("logit","nn-1","dnn-2","dnn-3","cart","bag","forest")

v1 <- readRDS("~/Projects/inProgress/2018_propensity_neuralnet_paper/data/res_all_nnmod.rds") |>
  mutate(method = fct_relevel(as.factor(method), meth_lvls))

read_row <- function(f) tryCatch(readRDS(f), error = function(e) qs2::qd_read(f))
rows <- list()
for (d in list.dirs("pilot_salvage/data", recursive = FALSE)) {
  for (f in list.files(d, pattern = "results-row", full.names = TRUE)) {
    x <- read_row(f)
    rows[[length(rows)+1]] <- bind_cols(as_tibble(x$results), as_tibble(x$condition))
  }
}
v2 <- bind_rows(rows) |> filter(method %in% meth_lvls) |>
  mutate(method = fct_relevel(as.factor(method), meth_lvls))

simstat <- function(res, stat_, by_ = "p") {
  res |> simsum(estvarname="ATE", true=0.3, se="ATE_se", methodvar="method", ref="logit", by=by_) |>
    summary() |> tidy() |> filter(stat == stat_) |>
    mutate(method = fct_relevel(as.factor(method), meth_lvls))
}
bar_fig <- function(df, yl, ylim = NULL, hlines = list(), pr = TRUE) {
  g <- ggplot(df, aes(x = method, y = est, fill = as.factor(p))) +
    ylab(yl) + geom_bar(position="dodge", stat="identity") + labs(fill="Covariates") +
    theme_Publication() + scale_fill_Publication()
  if (!is.null(ylim)) g <- g + coord_cartesian(ylim = ylim)
  for (h in hlines) g <- g + geom_hline(yintercept=h[[1]], linetype=h[[2]], color=h[[3]], size=h[[4]])
  if (pr) g <- g + geom_pointrange(aes(ymin=lower, ymax=upper), position=position_dodge(0.9), size=0.1, alpha=0.5)
  g
}
save_pair <- function(g1, g2, name, w = 13, h = 5) {
  gg <- (g1 + ggtitle("v1 (dissertation)")) + (g2 + ggtitle("v2 (corrected pilot)")) +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  ggsave(file.path(outdir, paste0(name, ".png")), gg, width = w, height = h, dpi = 150)
  cat("saved:", name, "\n")
}

## ASAM
asam <- function(res) res |> group_by(p, method) |> summarise(est = mean(ASAM), .groups="drop")
save_pair(bar_fig(asam(v1), "ASAM", c(0,.5), list(list(.1,2,"black",.25), list(.2,2,"black",.5)), pr=FALSE),
          bar_fig(asam(v2), "ASAM", c(0,.5), list(list(.1,2,"black",.25), list(.2,2,"black",.5)), pr=FALSE), "fig-asam")

## Bias
save_pair(bar_fig(simstat(v1,"bias"), "Bias"), bar_fig(simstat(v2,"bias"), "Bias"), "fig-bias")

## Relative bias (MCSE via n())
relb <- function(res) res |> mutate(rel_bias=(ATE-0.3)/0.3) |> group_by(p, method) |>
  summarise(est=mean(rel_bias), mcse=sd(rel_bias)/sqrt(n()), .groups="drop") |>
  mutate(lower=est-1.96*mcse, upper=est+1.96*mcse)
save_pair(bar_fig(relb(v1), "Relative Bias", c(-.6,.6), list(list(-.1,"dashed","grey",.5), list(.1,"dashed","grey",.5))),
          bar_fig(relb(v2), "Relative Bias", c(-.6,.6), list(list(-.1,"dashed","grey",.5), list(.1,"dashed","grey",.5))), "fig-rel-bias")

## Relative bias faceted by scenario
relbc <- function(res) res |> mutate(rel_bias=(ATE-0.3)/0.3) |> group_by(p, method, scenarioT, scenarioY) |>
  summarise(est=mean(rel_bias), mcse=sd(rel_bias)/sqrt(n()), .groups="drop") |>
  mutate(lower=est-1.96*mcse, upper=est+1.96*mcse)
f1 <- bar_fig(relbc(v1), "Relative Bias", c(-.6,.6), list(list(-.1,"dashed","grey",.5), list(.1,"dashed","grey",.5))) +
  facet_grid(scenarioT ~ scenarioY) + theme(panel.spacing = unit(.5,"cm"))
f2 <- bar_fig(relbc(v2), "Relative Bias", c(-.6,.6), list(list(-.1,"dashed","grey",.5), list(.1,"dashed","grey",.5))) +
  facet_grid(scenarioT ~ scenarioY) + theme(panel.spacing = unit(.5,"cm"))
save_pair(f1, f2, "fig-rel-bias-cond", w = 15, h = 7)

## SE, MSE, coverage, relerror, power (simsum stats)
save_pair(bar_fig(simstat(v1,"modelse"), "Standard Error (SE)"), bar_fig(simstat(v2,"modelse"), "Standard Error (SE)"), "fig-se")
save_pair(bar_fig(simstat(v1,"mse"), "Mean Squared Error (MSE)", c(0,10)), bar_fig(simstat(v2,"mse"), "Mean Squared Error (MSE)", c(0,10)), "fig-mse")
save_pair(bar_fig(simstat(v1,"cover"), "95% C.I. Coverage", c(0,1), list(list(.95,2,"black",.5))),
          bar_fig(simstat(v2,"cover"), "95% C.I. Coverage", c(0,1), list(list(.95,2,"black",.5))), "fig-cover")
save_pair(bar_fig(simstat(v1,"relerror"), "Relative % error in standard error", NULL, list(list(-5,2,"black",.5))),
          bar_fig(simstat(v2,"relerror"), "Relative % error in standard error", NULL, list(list(-5,2,"black",.5))), "fig-ate-relerr")
save_pair(bar_fig(simstat(v1,"power"), "Power", c(0,1)), bar_fig(simstat(v2,"power"), "Power", c(0,1)), "fig-power")

## Weights
wts <- function(res) res |> group_by(p, method) |> summarise(est = mean(mean_ps_weights), .groups="drop")
save_pair(bar_fig(wts(v1), "IPTW Weights", c(0,5), pr=FALSE), bar_fig(wts(v2), "IPTW Weights", c(0,5), pr=FALSE), "fig-weights")

cat("ALL FIGURES DONE\n")
