# ==============================
# 1. Setup
# ==============================
library(dplyr)
library(ARTool)
library(emmeans)
library(rstatix)  # for Dunn test, tidy KW
library(tidyr)

#########load in data #################### 
metrics_df <- read_xlsx("Precision_sensitivity_calculations.xlsx")

metrics_long <- metrics_df %>% 
  pivot_longer(
    cols = c(no_damage, all_damage, mudline_damage, 
             middle_damage, bottom_damage,
             spiked_mudline_damage, spiked_middle_damage, spiked_bottom_damage),
    names_to = "damage",
    values_to = "value"
  )

metrics_f1 <- metrics_long %>%
  pivot_wider(names_from = stat_type, values_from = value) %>%
  mutate(F1 = 2 * (precision * sensitivity) / (precision + sensitivity))

ppv_df <- metrics_f1 %>%
  rename(
    rank = Level,
    complexity = input_fragment
  ) %>%
  select(rank, complexity, damage, precision, sensitivity, F1)


#because I also have two different contexts, need to take that into account (e.g. simple vs. combined IODP)
ppv_df <- ppv_df %>%
  mutate(
    context = ifelse(grepl("^spiked", damage), "combined", "simple"),
    treatment = case_when(
      damage %in% c("no_damage", "all_damage") ~ damage,
      TRUE ~ gsub("spiked_", "", damage)  # strips "spiked_" prefix
    )
  )
ppv_df <- ppv_df %>%
  select(-damage)
# ==============================
# 1b. Factor fixed and random effects safely
# ==============================
ppv_df$complexity <- factor(ppv_df$complexity,
                            levels = c(10, 25, 50, 100, 250, 500, 1000, 2500, 5000))
ppv_df$context <- factor(ppv_df$context, levels = c("simple", "combined"))
ppv_df$treatment <- factor(ppv_df$treatment,
                           levels = c("no_damage", "mudline_damage", "middle_damage", "bottom_damage", "all_damage"))
ppv_df$rank <- factor(ppv_df$rank)

metric <- c("precision", "sensitivity", "F1")


# ------------------------------
# Helper: ART wrappers
# ------------------------------
run_art_tests <- function(df, response, ref_complexity = "5000", adjust_method = "BH") {
  formula <- as.formula(paste(response, "~ complexity * treatment"))
  
  art_model <- art(formula, data = df)
  
  # ANOVA
  art_anova <- anova(art_model)
  
  # Estimated marginal means (complexity pooled across factors)
  em <- emmeans(artlm(art_model, "complexity"),
                pairwise ~ complexity,
                adjust = adjust_method,
                ref = ref_complexity, 
                infer = c(TRUE, TRUE))  # adds conf.low and conf.high)
  
  contrasts <- as.data.frame(em$contrasts) %>%
    mutate(sig = ifelse(p.value < 0.05, "*", "ns"))
  
  return(list(
    anova = art_anova,
    posthoc = contrasts,
    emmeans = as.data.frame(em$emmeans),
    ref_complexity = ref_complexity
  ))
}

run_art_tests_paired <- function(df, response, ref_complexity = "5000", adjust_method = "BH") {
  formula <- as.formula(paste(response, "~ complexity * treatment * context"))
  
  art_model <- art(formula, data = df)
  
  # ANOVA
  art_anova <- anova(art_model)
  
  # Estimated marginal means (complexity pooled across factors)
  em <- emmeans(artlm(art_model, "complexity"),
                pairwise ~ complexity,
                adjust = adjust_method,
                ref = ref_complexity, 
                infer = c(TRUE, TRUE)  # adds conf.low and conf.high
                )
  
  contrasts <- as.data.frame(em$contrasts) %>%
    mutate(sig = ifelse(p.value < 0.05, "*", "ns"))
  
  return(list(
    anova = art_anova,
    posthoc = contrasts,
    emmeans = as.data.frame(em$emmeans),
    ref_complexity = ref_complexity
  ))
}

# ------------------------------
# Helper: Format ART results
# ------------------------------

art_to_tibble <- function(posthoc, metric, model_type) {
  posthoc %>%
    as_tibble() %>%
    select(any_of(c("contrast", "estimate", "SE", "df", "p.value", "conf.low", "conf.high"))) %>%
    mutate(metric = metric, model = model_type)
}
# ------------------------------
# Helper: KW + Dunn wrapper
# ------------------------------
kw_dunn_tibble <- function(df, response, group_var = "complexity") {
  # KW
  kw <- df %>%
    group_by(rank) %>%
    kruskal_test(as.formula(paste(response, "~ complexity"))) %>%
    mutate(metric = response, test = "Kruskal-Wallis") %>%
    ungroup()
  
  # Dunn
  dunn <- df %>%
    group_by(rank) %>%
    dunn_test(as.formula(paste(response, "~ complexity")),
              p.adjust.method = "BH") %>%
    mutate(metric = response, test = "Dunn") %>%
    ungroup()
  
  # Tidy
  tibble(
    metric = response,
    test = "Kruskal-Wallis",
    statistic = kw$statistic,
    df = kw$df,
    p.value = kw$p
  ) %>%
    bind_rows(
      dunn %>%
        select(group1, group2, p.adj, p = p) %>%
        mutate(
          contrast = paste(group1, "-", group2),
          metric = response,
          test = "Dunn"
        ) %>%
        select(metric, test, contrast, p, p.adj)
    )
}


# ==============================
# 2. Run tests
# ==============================
# KW + Dunn
kw_dunn_all <- bind_rows(
  kw_dunn_tibble(ppv_df, "precision"),
  kw_dunn_tibble(ppv_df, "sensitivity"),
  kw_dunn_tibble(ppv_df, "F1")
)

# Pooled ART
ppv_results   <- run_art_tests(ppv_df, "precision")
recall_results <- run_art_tests(ppv_df, "sensitivity")
f1_results    <- run_art_tests(ppv_df, "F1")

# Paired ART (subset treatment)
paired_df <- ppv_df %>%
  filter(treatment %in% c("mudline_damage", "middle_damage", "bottom_damage"))

ppv_results_paired   <- run_art_tests_paired(paired_df, "precision")
recall_results_paired <- run_art_tests_paired(paired_df, "sensitivity")
f1_results_paired    <- run_art_tests_paired(paired_df, "F1")

# ART pooled + paired in one tibble
art_all <- bind_rows(
  art_to_tibble(ppv_results$posthoc, "precision", "pooled"),
  art_to_tibble(recall_results$posthoc, "sensitivity", "pooled"),
  art_to_tibble(f1_results$posthoc, "F1", "pooled"),
  art_to_tibble(ppv_results_paired$posthoc, "precision", "paired"),
  art_to_tibble(recall_results_paired$posthoc, "sensitivity", "paired"),
  art_to_tibble(f1_results_paired$posthoc, "F1", "paired")
)

# ==============================
# 3. Export
# ==============================
write.csv(kw_dunn_all, "kw_dunn_results.csv", row.names = FALSE)
write.csv(art_all, "art_results.csv", row.names = FALSE)
