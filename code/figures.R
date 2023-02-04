######################################################################
# Load libraries and source functions
######################################################################

packages <- c(
  "tidyverse",
  "gt"
)

lapply(packages, library, character.only = TRUE)


# append simulation for rep=100
df1 <- readRDS("~/Projects/inProgress/2018_propensity_neuralnet_paper/data/sim_results_n10000_r1000_NP.rds")
df2 <- readRDS("~/Projects/inProgress/2018_propensity_neuralnet_paper/data/sim_results_n10000_r1000_P.rds")
res <- bind_rows(df1, df2)



# generate overall table
over_df <- res |> 
  dplyr::select(p, scenarioT, scenarioY, method, Bias, Abs_Per_Bias, Abs_Per_Rel_Bias, ATE_se, MSE, Power, coverage_95, ASAM)  

gt(over_df) |>
  tab_spanner(
    label = "Condition",
    columns = c(p, scenarioT, scenarioY)
  ) |>
  tab_spanner(
    label = "Performance Measure",
    columns = c(Bias, Abs_Per_Bias, Abs_Per_Rel_Bias, ATE_se, MSE, Power, coverage_95, ASAM)
  ) |>
  fmt_number(columns = c(Bias, Abs_Per_Bias, Abs_Per_Rel_Bias, ATE_se, MSE, Power, coverage_95, ASAM), decimals = 3) |>
  cols_align_decimal()

# table for weights




# table for standardized initial bias and probability of treatment
init_df <-
  res %>%
  group_by(scenarioT, scenarioY) %>%
  summarise(
    Std_In_Bias = mean(Std_In_Bias),
    Prob_Treat = mean(Prob_Treat),
    .groups = "rowwise"
  ) |>
  ungroup()


gt(init_df) |>
  tab_spanner(
    label = "Condition",
    columns = c(scenarioT, scenarioY)
  ) |>
  fmt_number(columns = c(Std_In_Bias, Prob_Treat), decimals = 2) |>
  cols_label(
    scenarioT = md("Treatment<br />Model"),
    scenarioY = md("Outcome<br />Model"),
    Std_In_Bias = md("Standardized<br />Initial Bias"),
    Prob_Treat = "P(Z=1)"
  ) |>
  cols_align_decimal()



# overall figure

over_df <-
  res %>%
  group_by(p, method) %>%
  summarise(
    Bias = mean(Bias),
    Abs_Per_Bias = mean(Abs_Per_Bias),
    Abs_Per_Rel_Bias = mean(Abs_Per_Rel_Bias),
    ATE_se = mean(ATE_se),
    MSE = mean(MSE),
    Power = mean(Power),
    PS_Weights = mean(mean_ps_weights),
    ASAM = mean(ASAM),
    coverage_95 = mean(coverage_95),
    .groups = "rowwise"
  ) %>%
  ungroup() %>%
  pivot_longer(-c(p, method), names_to = "metric", values_to = "value") %>%
  mutate(
    p = as.factor(p),
    method = as.factor(method),
    metric = as.factor(metric)
  )


ggplot(data = over_df, aes(x = method, y = value, fill = as.factor(p))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~metric, scales = "free") +
  theme_minimal() +
  guides(fill = guide_legend(title = "p")) +
  theme(
    legend.position = "right",
    panel.spacing = unit(1, "cm")
  ) +
  ylab(element_blank()) +
  xlab(element_blank()) 

# figures by condition



res_sum_df <-
  res %>%
  group_by(p, method, scenarioT, scenarioY) %>%
  summarise(
    Bias = mean(Bias),
    Abs_Per_Bias = mean(Abs_Per_Bias),
    Abs_Per_Rel_Bias = mean(Abs_Per_Rel_Bias),
    ATE_se = mean(ATE_se),
    MSE = mean(MSE),
    Power = mean(Power),
    PS_Weights = mean(mean_ps_weights),
    ASAM = mean(ASAM),
    coverage_95 = mean(coverage_95),
    .groups = "rowwise"
  ) %>%
  ungroup()

res_sum_df %>%
  ggplot(aes(x = method, y = Bias, fill = as.factor(p))) +
  ylab("Bias") +
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_grid(scenarioT ~ scenarioY, scales = "fixed") +
  scale_fill_discrete(name = "# Covars") +
  theme(legend.position = "top") +
  theme_minimal()

res_sum_df %>%
  ggplot(aes(x = method, y = Abs_Per_Bias, fill = as.factor(p))) +
  ylab("Abs_Per_Bias") +
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_grid(scenarioT ~ scenarioY, scales = "fixed") +
  scale_fill_discrete(name = "# Covars") +
  theme(legend.position = "top") +
  theme_minimal() +
  coord_cartesian(ylim=c(0,10))

res_sum_df %>%
  ggplot(aes(x = method, y = ASAM, fill = as.factor(p))) +
  ylab("ASAM") +
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "black") +
  facet_grid(scenarioT ~ scenarioY, scales = "fixed") +
  scale_fill_discrete(name = "# Covars") +
  theme(legend.position = "top") +
  theme_minimal()


res_sum_df %>%
  ggplot(aes(x = method, y = MSE, fill = as.factor(p))) +
  ylab("MSE") +
  geom_bar(position = "dodge", stat = "identity") +
  facet_grid(scenarioT ~ scenarioY, scales = "fixed") +
  scale_fill_discrete(name = "# Covars") +
  theme(legend.position = "top") +
  theme_minimal()

res_sum_df %>%
  ggplot(aes(x = method, y = ATE_se, fill = as.factor(p))) +
  ylab("SE") +
  geom_bar(position = "dodge", stat = "identity") +
  facet_grid(scenarioT ~ scenarioY, scales = "fixed") +
  scale_fill_discrete(name = "# Covars") +
  theme(legend.position = "top") +
  theme_minimal()

res_sum_df %>%
  ggplot(aes(x = method, y = PS_Weights, fill = as.factor(p))) +
  ylab("PS_Weights") +
  geom_bar(position = "dodge", stat = "identity") +
  facet_grid(scenarioT ~ scenarioY, scales = "fixed") +
  scale_fill_discrete(name = "# Covars") +
  theme(legend.position = "top") +
  theme_minimal()

ggplot(dat, aes(x = trueps, y = ps_pred)) +
  geom_point(shape = 21, alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "True PS", y = "PS predicted by main effects logistic regression")