# Load necessary libraries
library(ggplot2)
library(dplyr)
library(betareg)  # For beta regression
library(lme4)     # For linear mixed models
library(tidyverse)
library(readxl)
library (RColorBrewer)
#navigate to F1 stats file and upload 
setwd("C:/Users/davisee/OneDrive - University of Tasmania/Documents/Bioinformatics_Documentation/Post-MEGAN_stats")

f1_table <- read_xlsx(path ="F1_statistics_clean_table.xlsx")

f1_long <- f1_table %>% 
  pivot_longer(
    cols = -c(Level, input_fragment), 
    names_to = "damage_treatment", 
    values_to = "f1_score"
  )




# # Adjust p-values
# kruskal_results <- kruskal_results %>%
#   adjust.pvalue(method = "BH") %>%
#   arrange(p.adj)
# 
# 
# glm_result <- glm(cbind(successes, failures) ~ Treatment * as.factor(n) * Taxonomic_Level, 
#                   family = binomial(link = "logit"), data = data)
# summary(glm_result)
# 
# # 3. Kruskal-Wallis Test (Non-Parametric Alternative)
# kruskal_treatment <- kruskal.test(Proportion ~ Treatment, data = data)
# kruskal_n <- kruskal.test(Proportion ~ as.factor(n), data = data)
# kruskal_taxonomic <- kruskal.test(Proportion ~ Taxonomic_Level, data = data)
# 
# # --- Post-hoc Tests --- #
# # Pairwise comparisons if significant differences found
# pairwise_treatment <- pairwise.wilcox.test(data$Proportion, data$Treatment, p.adjust.method = "bonferroni")
# pairwise_taxonomic <- pairwise.wilcox.test(data$Proportion, data$Taxonomic_Level, p.adjust.method = "bonferroni")
# 
# # GLM analysis
# glm_results <- glm(f1_score ~ Level * damage_type, data = df_f1, family = gaussian())
# summary(glm_results)
# 
# dunn_test_BH <- dunnTest(Indentity ~ Gene_mapped, data = spec_iden, method = "bh")
# 
# # Boxplot with statistical comparisons
# ggplot(df_f1, aes(x = Level, y = f1_score, fill = Level)) +
#   geom_boxplot() +
#   facet_wrap(~damage_type, scales = "free_y") +
#   stat_compare_means(method = "kruskal.test", label = "p.signif") +
#   theme_minimal() +
#   labs(title = "F1 Score Comparisons Across Levels and Damage Types",
#        x = "Taxonomic Level", y = "F1 Score") +
#   scale_fill_brewer(palette = "Set2")
# 

#bind these together to do stat testing? 

# --- Data Formatting Instructions --- #
# Ensure both datasets have the same column structure:
# - Treatment: Categorical variable representing different experimental groups
# - n: Numeric variable representing sample sizes (e.g., 1000, 100, etc.)
# - Taxonomic_Level: Categorical variable (Species, Genus, Family)
# - Sensitivity, Precision: Numeric values between 0 and 1
# - F1: Harmonic mean of Sensitivity & Precision
# - Source: Factor variable indicating dataset origin ("Dataset1" or "Dataset2")

# # Simulated dataset example (Replace with actual data)
# data1 <- data.frame(
#   Treatment = rep(c("A", "B", "C", "D", "E"), each = 9),
#   n = rep(c(1000, 100, 10, 2500, 250, 25, 5000, 500, 50), times = 5),
#   Taxonomic_Level = rep(c("Species", "Genus", "Family"), length.out = 45),
#   Sensitivity = runif(45, 0.5, 1),
#   Precision = runif(45, 0.5, 1)
# )
# data1$F1 <- 2 * (data1$Sensitivity * data1$Precision) / (data1$Sensitivity + data1$Precision)
# data1$Source <- "Dataset1"
# 
# # Simulated second dataset
# set.seed(123)  # Ensures reproducibility
# data2 <- data1
# # Modify values slightly to simulate a different dataset
# data2$Sensitivity <- runif(45, 0.4, 0.95)
# data2$Precision <- runif(45, 0.4, 0.95)
# data2$F1 <- 2 * (data2$Sensitivity * data2$Precision) / (data2$Sensitivity + data2$Precision)
# data2$Source <- "Dataset2"
# 
# # Combine datasets
# combined_data <- rbind(data1, data2)

# --- Statistical Tests --- #
# 
# # 1. Beta Regression for F1 Score (Recommended for Proportions between 0 and 1)
# glm_f1 <- betareg(F1 ~ Source * Treatment * as.factor(n) * Taxonomic_Level, data = combined_data)
# summary(glm_f1)
# 
# # 2. Linear Mixed Model (Alternative for Repeated Measures)
# lmm_f1 <- lmer(F1 ~ Source * Treatment * as.factor(n) * Taxonomic_Level + (1 | Treatment), data = combined_data)
# summary(lmm_f1)
# 
# # 3. Kruskal-Wallis Test (Non-Parametric Alternative)
# kruskal_test_f1 <- kruskal.test(F1 ~ Source, data = combined_data)
# 
# # 4. Pairwise Wilcoxon Test (if significant differences are found)
# pairwise_wilcox_f1 <- pairwise.wilcox.test(combined_data$F1, combined_data$Source, p.adjust.method = "bonferroni")
# 
# # --- Visualizations --- #
# 
# # Boxplot of F1 Scores by Dataset
# ggplot(combined_data, aes(x = Source, y = F1, fill = Source)) +
#   geom_boxplot() +
#   theme_minimal() +
#   labs(x = "Dataset", y = "F1 Score", fill = "Dataset")
# 
# # Interaction Plot of F1 Score by n Value
# ggplot(combined_data, aes(x = as.factor(n), y = F1, color = Source, group = Source)) +
#   stat_summary(fun = mean, geom = "line", size = 1) +
#   stat_summary(fun = mean, geom = "point", size = 2) +
#   theme_minimal() +
#   labs(x = "n Value", y = "Mean F1 Score", color = "Dataset") +
#   scale_x_discrete(limits = as.character(c(10, 25, 50, 100, 250, 500, 1000, 2500, 5000)))
# 
# # Boxplot of F1 Scores by Taxonomic Level
# ggplot(combined_data, aes(x = Taxonomic_Level, y = F1, fill = Source)) +
#   geom_boxplot() +
#   theme_minimal() +
#   labs(x = "Taxonomic Level", y = "F1 Score", fill = "Dataset")
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(lme4)
library(emmeans)
library(multcompView)

# Step 1: Reshape your original wide table into long format (if not already done)
# f1_long <- f1_table %>%
#   pivot_longer(
#     cols = -c(Level, input_fragment),
#     names_to = "damage_treatment",
#     values_to = "f1_score"
#   )

# # Load necessary libraries
# library(dplyr)
# library(ggplot2)
# library(tidyr)
# library(lme4)
# library(emmeans)
# library(multcompView)

# Step 1: Reshape your original wide table into long format (if not already done)
# f1_long <- f1_table %>%
#   pivot_longer(
#     cols = -c(Level, input_fragment),
#     names_to = "damage_treatment",
#     values_to = "f1_score"
#   )

# Step 2: Descriptive Statistics
summary_stats <- f1_long %>%
  group_by(Level, damage_treatment, input_fragment) %>%
  summarise(
    mean_f1 = mean(f1_score, na.rm = TRUE),
    sd_f1 = sd(f1_score, na.rm = TRUE),
    median_f1 = median(f1_score, na.rm = TRUE),
    iqr_f1 = IQR(f1_score, na.rm = TRUE),
    min_f1 = min(f1_score, na.rm = TRUE),
    max_f1 = max(f1_score, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

# Step 3: ANOVA - test for effects of input_fragment and damage_treatment
anova_model <- aov(f1_score ~ input_fragment * damage_treatment, data = f1_long)
summary(anova_model)

# Step 4: Post-hoc Tests (Tukey's HSD for damage_treatment)
tukey_results <- TukeyHSD(anova_model, which = "damage_treatment")
print(tukey_results)

# Step 5: Visualizations

## Boxplot - f1 score by input fragment across treatments
ggplot(f1_long, aes(x = factor(input_fragment), y = f1_score, fill = damage_treatment)) +
  geom_boxplot() +
  facet_wrap(~ Level) +
  labs(title = "F1 Scores by Fragment Amount and Damage Treatment",
       x = "Input Fragment Amount",
       y = "F1 Score") +
  theme_minimal()

## Bar plot with error bars
ggplot(summary_stats, aes(x = factor(input_fragment), y = mean_f1, fill = damage_treatment)) +
  geom_col(position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_f1 - sd_f1, ymax = mean_f1 + sd_f1), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~ Level) +
  labs(title = "Mean F1 Scores by Damage Treatment and Fragment Amount",
       x = "Input Fragment Amount",
       y = "Mean F1 Score") +
  theme_minimal()

## Interaction plot
ggplot(f1_long, aes(x = factor(input_fragment), y = f1_score, color = damage_treatment, group = damage_treatment)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun = mean, geom = "point") +
  facet_wrap(~ Level) +
  labs(title = "Interaction Between Fragment Amount and Damage Treatment on F1 Score",
       x = "Input Fragment Amount",
       y = "F1 Score") +
  theme_minimal()

## Heatmap
ggplot(f1_long, aes(x = factor(input_fragment), y = damage_treatment, fill = f1_score)) +
  geom_tile() +
  facet_wrap(~ Level) +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Heatmap of F1 Scores by Fragment Amount and Damage Treatment",
       x = "Input Fragment Amount",
       y = "Damage Treatment") +
  theme_minimal()

# Step 6: Mixed Effects Model (if appropriate)
# Including random intercept for Level if multiple levels exist and are not main focus
mixed_model <- lmer(f1_score ~ input_fragment * damage_treatment + (1 | Level), data = f1_long)
summary(mixed_model)

# Optional: Post-hoc tests from mixed model
emm <- emmeans(mixed_model, ~ damage_treatment | input_fragment)
pairs(emm, adjust = "tukey")




library(glmmTMB)
library(purrr)
library(broom)

# Ensure F1 is within (0,1)
f1_long <- f1_long %>%
  mutate(f1_adj = ifelse(f1_score == 0, 0.0001,
                         ifelse(f1_score == 1, 0.9999, f1_score)))


model_beta <- glmmTMB(f1_adj ~ damage_treatment + input_fragment + (1 | Level),
                      family = beta_family(), data = f1_long_test)

summary(model_beta)

#check model fit 
library(DHARMa)

sim_res <- simulateResiduals(model_beta)
plot(sim_res)  # Check for uniformity, overdispersion, outliers


### the KS and combined outlier test were significant 
## try again? 
model_beta_interaction <- glmmTMB(f1_adj ~ damage_treatment * input_fragment + (1 | Level),
                                  family = beta_family(), data = f1_long_test)
summary(model_beta_interaction)

sim_res <- simulateResiduals(model_beta_interaction)
plot(sim_res) 


library(broom.mixed)

lrt_results <- f1_long_test %>%
  group_by(Level) %>%
  group_split() %>%
  map_dfr(~ {
    full_model <- glmmTMB(f1_adj ~ damage_treatment + input_fragment,
                          family = beta_family(), data = .x)
    null_model <- glmmTMB(f1_adj ~ input_fragment, family = beta_family(), data = .x)
    
    test <- anova(null_model, full_model, test = "LRT")
    tidy(test) %>% mutate(Level = unique(.x$Level))
  })

lrt_results



library(glmmTMB)
library(dplyr)
library(broom)

# Function to run beta regression and check significance of damage_treatment per Level
test_damage_treatment_per_level <- function(data) {
  # Fit the model for each Level
  model <- glmmTMB(f1_adj ~ damage_treatment + input_fragment + (1 | Level),
                   family = beta_family(), data = data)
  
  # Perform ANOVA to test significance of damage_treatment
  anova_results <- anova(model, test = "Chisq") # Likelihood ratio test
  tidy(anova_results) %>%
    mutate(Level = unique(data$Level))
}

# Apply the function to each Level group
lrt_results <- f1_long_test %>%
  group_by(Level) %>%
  group_split() %>%
  map_dfr(test_damage_treatment_per_level)

# View the results of the likelihood ratio tests
lrt_results


library(glmmTMB)
library(dplyr)
library(broom)

library(glmmTMB)
library(dplyr)
library(broom)

# Function to run beta regression and test significance of damage_treatment per Level
test_damage_treatment_per_level <- function(data) {
  # Fit the full model with damage_treatment
  full_model <- glmmTMB(f1_adj ~ damage_treatment + input_fragment + (1 | Level),
                        family = beta_family(), data = data)
  
  # Fit the null model without damage_treatment
  null_model <- glmmTMB(f1_adj ~ input_fragment + (1 | Level),
                        family = beta_family(), data = data)
  
  # Calculate log-likelihood values for both models
  logLik_full <- logLik(full_model)
  logLik_null <- logLik(null_model)
  
  # Compute the likelihood ratio statistic
  lr_stat <- 2 * (logLik_full - logLik_null)
  
  # Get the degrees of freedom for the chi-squared distribution
  df <- length(coef(full_model)) - length(coef(null_model))
  
  # Compute the p-value from the chi-squared distribution
  p_value <- pchisq(lr_stat, df, lower.tail = FALSE)
  
  # Return the results
  tibble(
    Level = unique(data$Level),
    lr_stat = lr_stat,
    df = df,
    p_value = p_value
  )
}

# Apply the function to each Level group
lrt_results <- f1_long_test %>%
  group_by(Level) %>%
  group_split() %>%
  map_dfr(test_damage_treatment_per_level)

# View the results of the likelihood ratio tests
lrt_results

library(emmeans)
library(ggpubr)

# Fit model for each Level
model_per_level <- glmmTMB(f1_adj ~ damage_treatment + input_fragment + (1 | Level),
                           family = beta_family(), data = f1_long_test)

# Pairwise comparisons for damage_treatment
emmeans_results <- emmeans(model_per_level, pairwise ~ damage_treatment | Level, type = "response")
summary(emmeans_results)



kruskal_results <- f1_long %>%
  group_by(damage_treatment) %>%
  kruskal_test(f1_score ~ Level)

# GLM analysis
glm_results <- glm(f1_score ~ Level * damage_treatment, data = f1_long, family = gaussian())
summary(glm_results)

# Boxplot with statistical comparisons
ggplot(f1_long, aes(x = Level, y = f1_score, fill = Level)) +
  geom_boxplot() +
  facet_wrap(~damage_treatment, scales = "free_y") +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +
  theme_minimal() +
  labs(title = "F1 Score Comparisons Across Levels and Damage Types",
       x = "Taxonomic Level", y = "F1 Score") +
  scale_fill_brewer(palette = "Set2")


################################
# Load necessary packages
library(tidyverse)
library(rstatix)     # For Kruskal-Wallis and Dunn's test with adjustments
library(ggpubr)       # For visualization and stat_pvalue_manual
library(ggsignif)     # Optional for adding significance brackets

# ==== 1. Load Data ====
# Replace this with your actual data load
# df <- read_csv("your_data.csv")

# ==== 2. FILTER DATA: Compare spiked damage vs no_damage ====
spiked_df <- df %>%
  filter(str_detect(damage_treatment, "spiked") | damage_treatment == "no_damage")

# ==== 3. Non-parametric test (Kruskal-Wallis) ====
spiked_kruskal <- spiked_df %>%
  kruskal_test(f1_score ~ damage_treatment)

# ==== 4. Post-hoc pairwise test (Dunn’s Test with BH correction) ====
spiked_dunn <- spiked_df %>%
  dunn_test(f1_score ~ damage_treatment, p.adjust.method = "BH")

# Add significance letters
spiked_letters <- spiked_dunn %>%
  add_xy_position(x = "damage_treatment")

# ==== 5. PLOT: F1 Score by Damage Treatment ====
ggplot(spiked_df, aes(x = damage_treatment, y = f1_score, fill = damage_treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  stat_pvalue_manual(spiked_letters, label = "p.adj.signif", tip.length = 0.01) +
  theme_minimal(base_size = 14) +
  labs(title = "F1 Score by Spiked vs No Damage Treatment",
       y = "F1 Score", x = "Damage Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")

# ==== 6. TEST: Does F1 differ by Taxonomic Level? ====
level_kruskal <- df %>%
  kruskal_test(f1_score ~ Level)

level_dunn <- df %>%
  dunn_test(f1_score ~ Level, p.adjust.method = "BH")

level_letters <- level_dunn %>%
  add_xy_position(x = "Level")

# ==== 7. PLOT: F1 Score by Taxonomic Level ====
ggplot(df, aes(x = Level, y = f1_score, fill = Level)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.6) +
  stat_pvalue_manual(level_letters, label = "p.adj.signif", tip.length = 0.01) +
  theme_minimal(base_size = 14) +
  labs(title = "F1 Score by Taxonomic Level",
       y = "F1 Score", x = "Taxonomic Level") +
  scale_fill_brewer(palette = "Pastel1")
#################################
# Load required libraries
library(tidyverse)
library(rstatix)     # Kruskal-Wallis, Dunn’s Test
library(ggpubr)       # For stat_pvalue_manual
library(ggsignif)     # For significance brackets if needed

# === 1. Load your data ===
# df <- read_csv("your_data.csv")

# === 2. NON-PARAMETRIC TEST: F1 Score by Damage Treatment within each Level ===
damage_within_level_tests <- df %>%
  group_by(Level) %>%
  kruskal_test(f1_score ~ damage_treatment)

damage_within_level_posthoc <- df %>%
  group_by(Level) %>%
  dunn_test(f1_score ~ damage_treatment, p.adjust.method = "BH")

# Add y-position for plotting
damage_within_level_letters <- damage_within_level_posthoc %>%
  add_xy_position(x = "damage_treatment")

# === 3. Plot: F1 Score by Damage Treatment, Faceted by Taxonomic Level ===
ggplot(df, aes(x = damage_treatment, y = f1_score, fill = damage_treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_pvalue_manual(damage_within_level_letters, label = "p.adj.signif", 
                     tip.length = 0.01, hide.ns = TRUE) +
  facet_wrap(~ Level, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(title = "F1 Score by Damage Treatment within Taxonomic Level",
       y = "F1 Score", x = "Damage Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")

# === 4. TEST EFFECT OF INPUT AMOUNT ===
# Kruskal test for input_fragment across full dataset
kruskal_input <- df %>%
  kruskal_test(f1_score ~ as.factor(input_fragment))

# Optional: test by Level
kruskal_input_by_level <- df %>%
  group_by(Level) %>%
  kruskal_test(f1_score ~ as.factor(input_fragment))

# === 5. Interaction test using linear model (non-parametric workaround limited) ===
# Use logit_f1 to allow parametric model for interaction detection
interaction_model <- lm(logit_f1 ~ damage_treatment * input_fragment * Level, data = df)
summary(interaction_model)

# === 6. Optional: Post-hoc from interaction model ===
library(emmeans)
emm_interaction <- emmeans(interaction_model, ~ damage_treatment | Level * input_fragment)
pairs(emm_interaction, adjust = "BH")

