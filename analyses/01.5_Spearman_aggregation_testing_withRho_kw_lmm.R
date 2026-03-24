# ==============================
# 1. Setup
# ==============================
library(dplyr)
library(rstatix)
library(psych)   # for correlation comparison
library(tidyverse)

# Function to compute Spearman rho per treatment
spearman_by_group <- function(df, response) {
  df %>%
    group_by(rank, context, treatment) %>%
    cor_test(complexity, !!sym(response), method = "spearman") %>%
    mutate(metric = response) %>%
    ungroup()
}

# Collect correlations for all metrics
spearman_results <- bind_rows(
  spearman_by_group(ppv_df, "precision"),
  spearman_by_group(ppv_df, "sensitivity"),
  spearman_by_group(ppv_df, "F1")
)



library(dplyr)
library(rstatix)
library(lme4)
library(lmerTest)

ppv_long <- ppv_df %>%
  pivot_longer(cols = c(precision, sensitivity, F1),
               names_to = "metric", values_to = "value")

ppv_long$rank <- factor(ppv_long$rank)
ppv_long$treatment <- factor(ppv_long$treatment, levels = c("no_damage", "mudline_damage",
                                                  "middle_damage", "bottom_damage", "all_damage"))
ppv_long$context <- factor(ppv_long$context, levels = c("simple", "combined"))

# ---------------------------
# 1. Kruskal-Wallis per group
# ---------------------------
agg_test <- ppv_long %>%
  group_by(rank, context, metric) %>%
  nest() %>%
  mutate(
    test = map(data, ~ kruskal_test(.x, value ~ treatment))
  ) %>%
  unnest(test) %>%
  ungroup() %>%
  select(rank, context, metric, n, statistic, df, p, method) %>%
  mutate(kw_decision = ifelse(p > 0.05, "OK to aggregate", "Keep separate"))
write.csv(agg_test, "kw_aggregate.csv")

# ---------------------------
# 2. Global Linear Mixed Model
# ---------------------------
fit <- lmer(value ~ treatment + (1|rank) + (1|context) + (1|metric),
            data = ppv_long)

lmm_anova <- anova(fit) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Effect") %>%
  select(Effect, NumDF, DenDF, `F value`, `Pr(>F)`)

lmm_decision <- ifelse(lmm_anova$`Pr(>F)`[lmm_anova$Effect == "treatment"] > 0.05,
                       "OK to aggregate", "Keep separate")
write.csv(lmm_anova, "lmm_anova.csv")

# ---------------------------
# 3. Final summary
# ---------------------------
final_summary <- agg_test %>%
  mutate(global_lmm_decision = lmm_decision)

final_summary

write.csv(final_summary, "final_aggregation_decision.csv")

# ---------------------------
# 1. Kruskal-Wallis per group
# ---------------------------
agg_test <- ppv_long %>%
  group_by(rank, metric) %>%
  nest() %>%
  mutate(
    test = map(data, ~ kruskal_test(.x, value ~ context))
  ) %>%
  unnest(test) %>%
  ungroup() %>%
  select(rank, metric, n, statistic, df, p, method) %>%
  mutate(kw_decision = ifelse(p > 0.05, "OK to aggregate", "Keep separate"))

# ---------------------------
# 2. Global Linear Mixed Model
# ---------------------------
fit <- lmer(value ~ context + (1|rank) + (1|metric),
            data = ppv_long)

lmm_anova <- anova(fit) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Effect") %>%
  select(Effect, NumDF, DenDF, `F value`, `Pr(>F)`)

lmm_decision <- ifelse(lmm_anova$`Pr(>F)`[lmm_anova$Effect == "treatment"] > 0.05,
                       "OK to aggregate", "Keep separate")


# ---------------------------
# 3. Final summary
# ---------------------------
final_summary <- agg_test %>%
  mutate(global_decision = lmm_decision)

final_summary

# ---------------------------
  # 1. Kruskal-Wallis per group
  # ---------------------------
agg_test <- ppv_long %>%
  group_by(context, treatment, metric) %>%
  nest() %>%
  mutate(
    test = map(data, ~ kruskal_test(.x, value ~ rank))
  ) %>%
  unnest(test) %>%
  ungroup() %>%
  select(treatment, context, metric, n, statistic, df, p, method) %>%
  mutate(kw_decision = ifelse(p > 0.05, "OK to aggregate", "Keep separate"))

# ---------------------------
# 2. Global Linear Mixed Model
# ---------------------------
fit <- lmer(value ~ rank + (1|treatment) + (1|context) + (1|metric),
            data = ppv_long)

lmm_anova <- anova(fit) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Effect") %>%
  select(Effect, NumDF, DenDF, `F value`, `Pr(>F)`)

lmm_decision <- ifelse(lmm_anova$`Pr(>F)`[lmm_anova$Effect == "rank"] > 0.05,
                       "OK to aggregate", "Keep separate")


# ---------------------------
# 3. Final summary
# ---------------------------
final_summary <- agg_test %>%
  mutate(global_decision = lmm_decision)

final_summary
