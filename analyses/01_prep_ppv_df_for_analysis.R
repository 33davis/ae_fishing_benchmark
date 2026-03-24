# ==============================
# 1. Setup
# ==============================
library(dplyr)
library(readxl)
library(tidyr)
library(rstatix)   # has correlation functions
library(ggplot2)

# --------- Load in data ----------
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

# Add context and treatment
ppv_df <- ppv_df %>%
  mutate(
    context = ifelse(grepl("^spiked", damage), "combined", "simple"),
    treatment = case_when(
      damage %in% c("no_damage", "all_damage") ~ damage,
      TRUE ~ gsub("spiked_", "", damage)  # strips "spiked_" prefix
    )
  ) %>%
  select(-damage)

# Factor categorical variables
ppv_df$context   <- factor(ppv_df$context, levels = c("simple", "combined"))
ppv_df$treatment <- factor(ppv_df$treatment,
                           levels = c("no_damage", "mudline_damage",
                                      "middle_damage", "bottom_damage",
                                      "all_damage"))
ppv_df$rank <- factor(ppv_df$rank)

# Ensure complexity is numeric (not factor!)
ppv_df$complexity <- as.numeric(ppv_df$complexity)

# --------- Save for analysis ----------
write.csv(ppv_df, "ppv_df.csv")
