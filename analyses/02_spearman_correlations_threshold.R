# =========================================
# Title: Spearman correlation thresholds (ranked)
# Purpose: Identify global (pooled across treatment) thresholds
# Input: ppv_df.csv (wide format)
# Output: Thresholds (global) + optional metric-strict variants
# =========================================

# ------------------------------
# 0. Load packages
# ------------------------------
library(dplyr)
library(rstatix)
library(ggplot2)
library(readr)
library(purrr)
library(tidyr)
library(zoo)
# ------------------------------
# 1. Load & prepare data
# ------------------------------
ppv_df <- read_csv("ppv_df.csv")

ppv_long <- ppv_df %>%
  pivot_longer(cols = c(precision, sensitivity, F1),
               names_to = "metric", values_to = "value") %>%
  mutate(
    rank = factor(rank, 
                  levels = c("Family", "Genus", "Species")),
    treatment = factor(treatment,
                       levels = c("no_damage", "mudline_damage",
                                  "middle_damage", "bottom_damage", "all_damage")),
    context = factor(context, levels = c("simple", "combined"))
  )

# =========================================
# 2. Global Spearman correlations (aggregated across treatment)
# =========================================
# Based on 1.5 aggregation testing: aggregate across treatment,
# keep rank × metric × context structure.

spearman_results <- ppv_long %>%
  group_by(rank, metric, context) %>%
  reframe(
    cor_test = list(cor.test(complexity, value, method = "spearman")),
    mean_metric_value = mean(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    rho = map_dbl(cor_test, ~ .x$estimate),
    p.value = map_dbl(cor_test, ~ .x$p.value),
    p.adj = p.adjust(p.value, method = "BH"),
    threshold_type = "Global"
  ) %>%
  select(rank, metric, context, rho, p.value, p.adj, threshold_type, mean_metric_value)

write_csv(spearman_results, "spearman_results_ranked.csv")

# ------------------------------
# 3. Significant correlations only
# ------------------------------
sig_rho <- spearman_results %>%
  filter(p.adj < 0.05, rho > 0)  # only positive correlations

# sig_rho_ne <- spearman_results %>%
#   filter(p.adj < 0.05, rho < 0)  # only negative correlations
# 
# sig_rho_all <- spearman_results %>%
#   filter(p.adj < 0.05)  # all significant correlations

write_csv(sig_rho, "spearman_significant_rho_ranked.csv")

# =========================================
# 4. Metric-strict Global thresholds 
# =========================================
# Define all strict cut-offs
strict_cutoffs <- c(0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5)

# Prepare an empty list to store results
thresholds_list <- list()
previous_cutoff <- 1.0  # initialize previous cutoff to 1.0
# Loop through each cut-off
for (cutoff in strict_cutoffs) {
  sig_strict <- sig_rho %>%
    left_join(ppv_long, by = c("rank", "context", "metric")) %>%
    filter(value > cutoff, 
           value <= previous_cutoff, 
           metric %in% c("precision", "sensitivity"))  # ensure metric is valid
  
  thresholds_strict <- sig_strict %>%
    group_by(rank, metric, context) %>%
    summarise(
      threshold = min(complexity, na.rm = TRUE),
      rho_at_threshold = rho[which.min(complexity)],
      pvalue_adjusted_BH = p.adj[which.min(complexity)],
      metric_value = value[which.min(complexity)],
      threshold_type = paste0("Global metric (value > ", cutoff, ")"),
      expected_error = paste0("Error rate (", (1-cutoff)*100, "%)"),
      .groups = "drop"
    ) %>%
    arrange(rank, metric, context)
  
  thresholds_list[[as.character(cutoff)]] <- thresholds_strict
  previous_cutoff <- cutoff  # update upper bound cutoff
}

# Combine all thresholds into one data frame
thresholds_strict_all <- bind_rows(thresholds_list)

print(thresholds_strict_all)

# Write to CSV
write_csv(thresholds_strict_all, "thresholds_all_metricstrict_combined.csv")

# =========================================
# 5. Local (Stepwise) Spearman corelations
# identify threshold/ "transitional areas" 
# where rho changes sign or significance 
# =========================================
# ---------------------------------
# Function: compute local Spearman correlations within a moving window
# ---------------------------------
local_spearman <- function(df, window = 3) {
  complexities <- sort(unique(df$complexity), decreasing = TRUE)
  
  results <- map_df(seq_along(complexities), function(i) {
    # define window range centered at current complexity
    start <- max(1, i - floor(window / 2))
    end <- min(length(complexities), i + floor(window / 2))
    subset <- df %>% filter(complexity %in% complexities[start:end])
    
    if (nrow(subset) > 2) {
      ct <- suppressWarnings(cor.test(subset$complexity, subset$value, method = "spearman"))
      tibble(
        local_center = complexities[i],
        rho_local = ct$estimate,
        p_local = ct$p.value,
        p.adj = p.adjust(ct$p.value, method = "BH"),
        mean_metric_value = mean(subset$value, na.rm = TRUE)  # include metric value
      )
    } else {
      tibble(local_center = complexities[i], rho_local = NA, p_local = NA, mean_value = NA)
    }
  })
  
  return(results)
}

# ---------------------------------
# Apply by rank × metric × context
# ---------------------------------
local_results <- ppv_long %>%
  group_by(rank, metric, context) %>%
  group_modify(~ local_spearman(.x, window = 3)) %>%
  ungroup()

# ---------------------------------
# Identify transitional areas: sign changes or local significance
# ---------------------------------
local_trends <- local_results %>%
  group_by(rank, metric, context) %>%
  arrange(desc(local_center), rank) %>%
  mutate(
    sign_change = ifelse(
      is.na(rho_local) | is.na(lag(rho_local)),
      FALSE,  # treat missing as no change
      sign(rho_local) != sign(lag(rho_local))
    ),
    significant = p.adj < 0.05, 
    threshold_type = "Local"
  ) %>%
  ungroup() %>%
 filter(sign_change & significant) %>%
  # filter(significant) %>%
  #filter(sign_change) %>%
  select(rank, metric, context, local_center, rho_local, p_local, p.adj, mean_metric_value, threshold_type)


print(local_trends)

local_strictcutoff <- 0.6

sig_strict_local <- local_trends %>%
  filter(mean_metric_value > local_strictcutoff) %>%
  group_by(context, metric) %>%
  arrange(rank) 

View(sig_strict_local)

# ========================================= 
# Testing  multiple strict cutoffs
# =========================================
# Define all strict cut-offs
strict_cutoffs <- c(0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5)

# Prepare an empty list to store results
thresholds_list_local <- list()
previous_cutoff <- 1.0  # initialize previous cutoff to 1.0
# Loop through each cut-off
for (cutoff in strict_cutoffs) {
  sig_strict_loca <- local_trends %>%
    filter(mean_metric_value > cutoff, 
           mean_metric_value <= previous_cutoff) %>%
    group_by(rank, context, metric) %>%
    mutate(
      threshold = local_center,
      rho_at_threshold = rho_local,
      pvalue_adjusted_BH = p.adj,
      metric_value = mean_metric_value,
      threshold_type = paste0("Local metric (value > ", cutoff, ")"),
      expected_error = paste0("Error rate (", (1-cutoff)*100, "%)")) %>%
    select(rank, metric, context, threshold, rho_at_threshold,
           pvalue_adjusted_BH, metric_value, threshold_type, expected_error) %>%
    arrange(rank, desc(threshold))
            
  thresholds_list_local[[as.character(cutoff)]] <- sig_strict_loca
  
  #update upper bound cutoff
  previous_cutoff <- cutoff
}

# Combine all cutoffs into one table
thresholds_strict_local_all <- bind_rows(thresholds_list_local)


print(thresholds_strict_local_all)

# ---------------------------
# merge with global trends?
# ---------------------------
# Add a 'source' column to each dataset
thresholds_strict_all <- thresholds_strict_all %>%
  mutate(source = "Global")

thresholds_strict_local_all <- thresholds_strict_local_all %>%
  mutate(source = "Local")

# Combine both tables
thresholds_all_combined <- bind_rows(thresholds_strict_all, thresholds_strict_local_all)

# Arrange for readability
thresholds_all_combined <- thresholds_all_combined %>%
  arrange(rank, metric, context, threshold_type)


thresholds_all_combined <- thresholds_all_combined %>%
  select(source, rank, metric, context, threshold, rho_at_threshold,
         pvalue_adjusted_BH, metric_value, threshold_type, expected_error)

thresholds_all_combined <- thresholds_all_combined %>%
  mutate(
    expected_error = factor(expected_error,
                            levels = c("Error rate (1%)", "Error rate (5%)",
                                       "Error rate (10%)", "Error rate (20%)",
                                       "Error rate (30%)", "Error rate (40%)",
                                       "Error rate (50%)"))
  )

# View merged table
print(thresholds_all_combined)


write_csv(thresholds_all_combined, "thresholds_all_combined_ranked.csv")

# =========================================
# 6. Plot combined thresholds (optional)
# =========================================
library(ggplot2)

mytheme <- theme_minimal() +
  theme(
    # Legend
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.width = unit(0.6, "cm"),
    legend.spacing.x = unit(0.1, 'cm'),
    
    # Axis text and titles
    axis.text.x = element_text(size = 12, angle = 90, hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12, vjust = -0.5),
    axis.title.y = element_text(size = 12, vjust = 1.5),
    
    # Facet strip
    strip.text = element_text(size = 10),
    strip.background = element_blank(),
    
    # Panel spacing (between facets)
    panel.spacing = unit(1, "lines"),
    # Background with light grid
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", size = 0.3),
    panel.grid.minor = element_line(color = "grey95", size = 0.2)
  )




ggplot(
  thresholds_all_combined %>% 
    filter(
      threshold != 5000,
      metric_value >= 0.9),
  aes(x = threshold, y = rho_at_threshold, color = rank, shape = expected_error)
) +
  geom_point(size = 3) +
  facet_grid(metric~context) +
  mytheme+
  labs(
    title = "Comparison of Global (combined) vs Local (simple) Thresholds",
    x = "Complexity Threshold",
    y = "Spearman rho", 
    color = "Rank",              
    shape = "Expected Error"     
  )
# =======================
# End of script
# =======================
