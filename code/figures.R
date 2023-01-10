
# New facet label names for dose variable
t.labs <- c("Base", "Interactions", "Quad Terms", "Complex")
names(t.labs) <- c("A", "B", "C", "D")

# New facet label names for supp variable
y.labs <- c("Base", "Interactions", "Quad Terms", "Complex")
names(y.labs) <- c("a", "b", "c", "d")


results_summary %>%
  ggplot(aes(x = method, y = AbsBias, fill = as.factor(p))) +
  xlab("Method") +
  ylab("AbsBias") +
  geom_bar(position = "dodge", stat = "identity") +
  facet_grid(scenarioT ~ scenarioY, labeller = labeller(scenarioT = t.labs, scenarioY = y.labs), scales = "fixed") +
  scale_fill_discrete(name = "# Covars") +
  theme(legend.position = "top") +
  theme_minimal()



# ggplot(dat, aes(x = trueps, y = ps_pred)) +
#   geom_point(shape = 21, alpha = 0.2) +
#   geom_abline(slope = 1, intercept = 0) +
#   scale_x_continuous(limits = c(0, 1)) +
#   scale_y_continuous(limits = c(0, 1)) +
#   labs(x = "True PS", y = "PS predicted by main effects logistic regression")
#
