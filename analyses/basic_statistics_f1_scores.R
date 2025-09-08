# #############################################
# # Non-Parametric Analysis of Precision, Recall, and F1
# # Author: E.E. Davis
# # Goal: Test if metrics differ by Rank (taxonomic level),
# #       Complexity (input fragments), or Damage treatment
# #############################################
# ========================================
# Load libraries
# ========================================
library(tidyverse)
library(readxl)
library(rstatix)
library(ggpubr)
library(ggsignif)
library(purrr)
library(viridis)
library(emmeans)
library(Hmisc)

# ========================================
# 1. Load and format data
# ========================================
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

# Factor ordering
ppv_df$complexity <- factor(ppv_df$complexity,
                            levels = c(10, 25, 50, 100, 250, 500, 1000, 2500, 5000))
ppv_df$damage <- factor(ppv_df$damage,
                        levels = c("no_damage", "mudline_damage","spiked_mudline_damage", 
                                   "middle_damage", "spiked_middle_damage", 
                                   "bottom_damage", "spiked_bottom_damage", "all_damage"))

# ========================================
# 2. Kruskal-Wallis + Dunn testing
#    - Per rank (taxonomic level)
# ========================================

metrics <- c("precision", "sensitivity", "F1")

# ==========================================================
# Kruskal-Wallis tests: complexity ~ metric, grouped by rank
# ==========================================================
run_kw_per_rank <- function(df, metric_name) {
  df %>%
    group_by(rank) %>%
    kruskal_test(as.formula(paste(metric_name, "~ complexity"))) %>%
    mutate(metric = metric_name)
}

kw_results <- map_dfr(metrics, ~run_kw_per_rank(ppv_df, .x))
write_csv(kw_results, "KW_results_by_complexity.csv")

# ==========================================================
# Dunn’s post-hoc tests: complexity pairwise within rank × metric
# ==========================================================
run_dunn_per_rank <- function(df, metric_name) {
  df %>%
    group_by(rank) %>%
    dunn_test(as.formula(paste(metric_name, "~ complexity")), p.adjust.method = "BH") %>%
    mutate(metric = metric_name) %>%
    add_xy_position(x = "complexity")   # for plotting
}

dunn_results <- map_dfr(metrics, ~run_dunn_per_rank(ppv_df, .x))
write_csv(dunn_results, "Dunn_results_by_complexity.csv")

# ========================================
# Custom theme for plots
# ========================================
mytheme_bigheatmap_facetgrid_light <- theme_minimal() +
  theme(
    # Legend
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.width = unit(0.7, "cm"),
    legend.spacing.x = unit(0.2, 'cm'),
    
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
# ==========================================================
# Plot Precision
# ==========================================================
dunn_sig_precision <- dunn_results %>%
  filter(metric == "precision", p.adj <= 0.05)

plot_precision <- ggplot(ppv_df, aes(x = complexity, y = precision, 
                                     color = damage, group = damage)) +
  stat_summary(fun = median, geom = "line", linewidth = 1) +
  stat_summary(fun.data = median_hilow, geom = "errorbar", width = 0.2) +
  facet_wrap(~rank, scales = "free_y") +
  scale_colour_viridis_d(option = "C") +
  labs(y = "Precision (PPV)", x = "Input Fragment Amount") +
  stat_pvalue_manual(
    dunn_sig_precision,
    label = "p.adj.signif",
    hide.ns = TRUE
  )+ 
  mytheme_bigheatmap_facetgrid_light +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom", 
    strip.text = element_text(size = 9)) 
plot_precision # review plot before saving

ggsave("Precision_by_complexity.png", plot_precision, width = 10, height = 6, dpi = 300)

# ==========================================================
# Plot Sensitivity
# ==========================================================
dunn_sig_sensitivity <- dunn_results %>%
  filter(metric == "sensitivity", p.adj <= 0.05)

plot_sensitivity <- ggplot(ppv_df, aes(x = complexity, y = sensitivity, 
                                       color = damage, group = damage)) +
  stat_summary(fun = median, geom = "line", size = 1) +
  stat_summary(fun.data = median_hilow, geom = "errorbar", width = 0.2) +
  facet_wrap(~rank, scales = "free_y") +
  scale_color_viridis_d(option = "C") +
  labs(y = "Sensitivity (Recall)", x = "Input Fragment Amount") +
  stat_pvalue_manual(
    dunn_sig_sensitivity,
    label = "p.adj.signif",
    hide.ns = TRUE
  )+  
  mytheme_bigheatmap_facetgrid_light +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom", 
    strip.text = element_text(size = 9))
plot_sensitivity # review plot before saving

ggsave("Sensitivity_by_complexity.png", plot_sensitivity, width = 10, height = 6, dpi = 300)

# ==========================================================
# Plot F1
# ==========================================================
dunn_sig_f1 <- dunn_results %>%
  filter(metric == "F1", p.adj <= 0.05)

plot_f1 <- ggplot(ppv_df, aes(x = complexity, y = F1, 
                              color = damage, group = damage)) +
  stat_summary(fun = median, geom = "line", size = 1) +
  stat_summary(fun.data = median_hilow, geom = "errorbar", width = 0.2) +
  facet_wrap(~rank, scales = "free_y") +
  mytheme_bigheatmap_facetgrid_light +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom", 
    strip.text = element_text(size = 9)) +
  labs(y = "F1 Score", x = "Input Fragment Amount") +
  stat_pvalue_manual(
    dunn_sig_f1,
    label = "p.adj.signif",
    hide.ns = TRUE
  )
plot_f1 # review plot before saving
ggsave("F1_by_complexity.png", plot_f1, width = 10, height = 6, dpi = 300)


# ========================================
# 4. Optional: ART interaction test (Complexity × Damage × Rank)
# ========================================
 library(ARTool)
art_model <- art(F1 ~ complexity * damage + (1|rank), data = ppv_df)
anova(art_model)
#Optional post-hoc: 
  emmeans(art_model, pairwise ~ complexity | damage)
  