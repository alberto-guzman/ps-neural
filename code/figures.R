
library(gt)
library(tidyverse)

#####
# table for standardized intial bias and probabilyt of treatment
#####

init_df <-
  res %>% 
  group_by(scenarioT, scenarioY) %>% 
  summarise(Std_In_Bias = mean(Std_In_Bias),
            Prob_Treat = mean(Prob_Treat), .groups = "rowwise") |> 
  ungroup() 


gt(init_df) |> 
  tab_spanner(label = "Condition",
              columns = c(scenarioT, scenarioY)) |> 
  fmt_number(columns = c(Std_In_Bias, Prob_Treat), decimals = 2) |> 
  cols_label(
    scenarioT = md("Treatment<br />Model"),
    scenarioY = md("Outcome<br />Model"),
    Std_In_Bias = md("Standardized<br />Initial Bias"),
    Prob_Treat = "P(Z=1)") |> 
  cols_align_decimal()




over_df <-
  res %>% 
  group_by(p, method) %>% 
  summarise(Abs_Per_Bias = mean(Abs_Per_Bias),
            RMSE = mean(RMSE),
            .groups = "rowwise") %>%
  ungroup() %>%
  pivot_longer(-c(p, method), names_to = "metric", values_to = "value") %>%
  mutate(p = as.factor(p),
         method = as.factor(method),
         metric = as.factor(metric))


ggplot(data = over_df, aes(x = method, y = value, fill = as.factor(p))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ metric, scales = "free") +
  theme_minimal() +
  guides(fill = guide_legend(title = "p")) +
  theme(legend.position = "right",
        panel.spacing = unit(1, "cm")) +
  ylab(element_blank()) +
  xlab(element_blank()) 








#####
# overall figures
#####

over_df <-
  res %>% 
  group_by(p, method) %>% 
  summarise(Bias = mean(Bias),
            RelBias = mean(RelBias),
            SE = mean(ATE_se),
            ATT = mean(ATE),
            ASAM = mean(ASAM),
            coverage_95 = mean(coverage_95),
            MSE = mean(MSE),
            .groups = "rowwise") %>%
  ungroup() %>%
  pivot_longer(-c(p, method), names_to = "metric", values_to = "value") %>%
  mutate(p = as.factor(p),
         method = as.factor(method),
         metric = as.factor(metric))


ggplot(data = over_df, aes(x = method, y = value, fill = as.factor(p))) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  facet_wrap(~ metric, scales = "free") +
  theme_minimal() +
  guides(fill = guide_legend(title = "p")) +
  theme(legend.position = "right",
        panel.spacing = unit(1, "cm")) +
  ylab(element_blank()) +
  xlab(element_blank()) 



#####
# overall figures
#####


res_sum_df <-
  res %>% 
  group_by(p, method, scenarioT, scenarioY) %>% 
  summarise(Bias = mean(Bias),
            RelBias = mean(RelBias),
            ATE = mean(ATE),
            ATE_se = mean(ATE_se),
            ASAM = mean(ASAM),
            coverage_95 = mean(coverage_95),
            MSE = mean(MSE),
            .groups = "rowwise") %>%
  ungroup() 

# New facet label names for dose variable
t.labs <- c("Base", "Interactions", "Quad Terms", "Complex")
names(t.labs) <- c("A", "B", "C", "D")

# New facet label names for supp variable
y.labs <- c("Base", "Interactions", "Quad Terms", "Complex")
names(y.labs) <- c("a", "b", "c", "d")


res_sum_df %>%
  ggplot(aes(x = method, y = ATE, fill = as.factor(p))) +
  ylab("ATE") +
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "black") +
  facet_grid(scenarioT ~ scenarioY, labeller = labeller(scenarioT = t.labs, scenarioY = y.labs), scales = "fixed") +
  scale_fill_discrete(name = "# Covars") +
  theme(legend.position = "top") +
  theme_minimal()


res_sum_df %>%
  ggplot(aes(x = method, y = Bias, fill = as.factor(p))) +
  ylab("Bias") +
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_grid(scenarioT ~ scenarioY, labeller = labeller(scenarioT = t.labs, scenarioY = y.labs), scales = "fixed") +
  scale_fill_discrete(name = "# Covars") +
  theme(legend.position = "top") +
  theme_minimal()

res_sum_df %>%
  ggplot(aes(x = method, y = RelBias, fill = as.factor(p))) +
  ylab("Relative Bias") +
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_grid(scenarioT ~ scenarioY, labeller = labeller(scenarioT = t.labs, scenarioY = y.labs), scales = "fixed") +
  scale_fill_discrete(name = "# Covars") +
  theme(legend.position = "top") +
  theme_minimal()

res_sum_df %>%
  ggplot(aes(x = method, y = ATE_se, fill = as.factor(p))) +
  ylab("Standard Error") +
  geom_bar(position = "dodge", stat = "identity") +
  facet_grid(scenarioT ~ scenarioY, labeller = labeller(scenarioT = t.labs, scenarioY = y.labs), scales = "fixed") +
  scale_fill_discrete(name = "# Covars") +
  theme(legend.position = "top") +
  theme_minimal()


res_sum_df %>%
  ggplot(aes(x = method, y = ASAM, fill = as.factor(p))) +
  ylab("ASAM") +
  geom_bar(position = "dodge", stat = "identity") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "black") +
  facet_grid(scenarioT ~ scenarioY, labeller = labeller(scenarioT = t.labs, scenarioY = y.labs), scales = "fixed") +
  scale_fill_discrete(name = "# Covars") +
  theme(legend.position = "top") +
  theme_minimal()


res_sum_df %>%
  ggplot(aes(x = method, y = MSE, fill = as.factor(p))) +
  ylab("MSE") +
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
