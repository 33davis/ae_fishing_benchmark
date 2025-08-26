#this is a draft script for Precision - Recall 

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(readxl)
library(viridis)
library(metR)

#goal make versions of plots that mimic the ones in manuscript
#comparing all synthetic data
#comparing all spiked data
#direct synthetic vs spiked at same damage type

# Replace with your actual dataframe name
# Expected columns: tool, rank, complexity (numeric), damage (character),
#                   precision (PPV), sensitivity (Recall)
#excel sheet is sheet 2 of supplementary tables, but a copy within R project folder
metrics_df <- read_xlsx("Precision_sensitivity_calculations.xlsx")

#mutate long first 
metrics_long <- metrics_df %>% 
  pivot_longer(
    cols = c(no_damage, all_damage, mudline_damage, 
             middle_damage, bottom_damage,
             spiked_mudline_damage, spiked_middle_damage, spiked_bottom_damage),
    names_to = "damage",
    values_to = "value"
  )
#mutate wide from long to calculate f1 stats
metrics_f1 <- metrics_long %>%
  pivot_wider(names_from = stat_type, values_from = value) %>%
  mutate(F1 = 2 * (precision * sensitivity) / (precision + sensitivity))

### need to get data in format that will work for script
ppv_df <- metrics_f1 %>%
  rename(
    rank = Level,
    complexity = input_fragment
  ) %>%
  select(rank, complexity, damage, precision, sensitivity, F1)

# Convert complexity to factor for proper facet ordering
ppv_df$complexity <- factor(ppv_df$complexity,
                            levels = c(10, 25, 50, 100, 250, 500, 1000, 2500, 5000))
# Convert damage to factor with ordered levels
ppv_df$damage <- factor(ppv_df$damage,
                        levels = c("no damage", "mudline_damage","spiked_mudline_damage", "middle_damage", "spiked_middle_damage", "bottom_damage", "all_damage",
                                  "spiked_bottom_damage"))
# -------------------------------
# 6. Calculate AUC for each combination
# -------------------------------
ppv_df <- ppv_df %>%
  arrange(complexity, damage, sensitivity) %>%
  group_by(complexity, damage) %>%
  mutate(
    auc = sum(diff(sensitivity) * head(precision, -1), na.rm = TRUE))%>%
  ungroup()

# # -------------------------------
# # 7. Generate F1 score contour lines
# # -------------------------------
# get_f1_contours <- function(f1_scores = c(0.2, 0.4, 0.6, 0.8, 0.9)) {
#   recall_vals <- seq(0.01, 1, by = 0.01)
#   contour_data <- do.call(rbind, lapply(f1_scores, function(f1) {
#     precision_vals <- (f1 * recall_vals) / (2 * recall_vals - f1)
#     df <- data.frame(sensitivity = recall_vals,
#                      precision = precision_vals,
#                      F1 = f1)
#     df <- df[precision_vals <= 1 & precision_vals >= 0, ]  # keep valid values
#     return(df)
#   }))
#   return(contour_data)
# }
# 
# f1_df <- get_f1_contours()


# # Set target rank for plotting
# rank_to_plot <- "Species"
# 
# # Filter dataset for selected rank
# df_rank <- ppv_df %>%
#   filter(rank == rank_to_plot)
# 
# # Compute AUC (trapezoidal) and prepare tool label
# pr_combined <- df_rank %>%
#   arrange(tool, complexity, damage, sensitivity) %>%
#   group_by(tool, complexity, damage) %>%
#   mutate(
#     auc = sum(diff(sensitivity) * head(precision, -1), na.rm = TRUE),
#     tool_label = paste0(tool, " (AUC=", round(unique(auc), 3), ")")
#   ) %>%
#   ungroup()
# 
# # Define a function to generate F1 score contour data
# get_f1_contours <- function(f1_scores = c(0.2, 0.4, 0.6, 0.8, 0.9)) {
#   recall_vals <- seq(0.01, 1, by = 0.01)
#   contour_data <- do.call(rbind, lapply(f1_scores, function(f1) {
#     precision_vals <- (f1 * recall_vals) / (2 * recall_vals - f1)
#     df <- data.frame(sensitivity = recall_vals,
#                      precision = precision_vals,
#                      F1 = f1)
#     df <- df[precision_vals <= 1 & precision_vals >= 0, ]  # Keep valid values only
#     return(df)
#   }))
#   return(contour_data)
# }
# 
# # Generate contour lines for F1 scores
# f1_df <- get_f1_contours()
# 
# # Convert complexity to factor with ordered levels for proper facet ordering
# pr_combined$complexity <- factor(pr_combined$complexity,
#                                  levels = c(10, 25, 50, 100, 250, 500, 1000, 2500, 5000))
# 
# # Plot with PR curves + F1 contours
# ggplot(pr_combined, aes(x = sensitivity, y = precision, color = tool_label, group = tool_label)) +
#   # PR curves and points
#   geom_line(size = 1.2) +
#   geom_point(size = 2, alpha = 0.8) +
#   
#   # F1 contours
#   geom_line(data = f1_df, aes(x = sensitivity, y = precision, group = F1),
#             inherit.aes = FALSE, linetype = "dashed", color = "gray50") +
#   
#   # Labels for F1 contours
#   geom_text(data = f1_df %>% group_by(F1) %>% filter(sensitivity == max(sensitivity)),
#             aes(x = sensitivity, y = precision, label = paste0("F1=", F1)),
#             inherit.aes = FALSE, hjust = 1.1, vjust = 0.5, size = 3, color = "gray30") +
#   
#   # Facet by numeric complexity (rows) and damage (columns)
#   facet_grid(rows = vars(complexity), cols = vars(damage)) +
#   
#   # Labels and theme
#   labs(title = paste("Precision-Recall Curves with F1 Score Contours at", rank_to_plot, "Level"),
#        subtitle = "Faceted by Species Complexity (numeric) and Deamination Damage",
#        x = "Sensitivity (Recall)",
#        y = "Precision (PPV)",
#        color = "Tool (AUC)") +
#   theme_minimal(base_size = 14) +
#   theme(strip.text = element_text(size = 12),
#         axis.title = element_text(size = 13),
#         legend.position = "bottom",
#         legend.title = element_text(size = 12))


# Function to generate F1 contour data
make_f1_contours <- function(n_points = 200, f1_values = c(0.2, 0.4, 0.6, 0.8)) {
  sens <- seq(0.001, 1, length.out = n_points)
  purrr::map_dfr(f1_values, function(f1) {
    prec <- (f1 * sens) / (2 * sens - f1)
    tibble(
      sensitivity = sens,
      precision = prec,
      F1 = f1
    ) %>%
      filter(!is.nan(precision), precision >= 0, precision <= 1)
  })
}

f1_contours <- make_f1_contours()

# Plot: Precision vs Sensitivity with F1 contours and AUC curves
ggplot() +
  # F1 score contours
  geom_line(data = f1_contours,
            aes(x = sensitivity, y = precision, group = F1),
            linetype = "dashed", color = "grey50") +
  # geom_text(data = f1_contours %>% group_by(F1) %>% slice_tail(n=1),
  #           aes(x = sensitivity, y = precision, label = paste0("F1=", F1)),
  #           hjust = -0.1, vjust = 0.5, size = 3, color = "grey40") +
  geom_text(
    data = f1_contours %>% group_by(F1) %>% slice_tail(n = 1),
    aes(x = sensitivity, y = precision, label = paste0("F1=", F1)),
    hjust = -0.1, vjust = 0.5, size = 3, color = "grey40"
  )+
  
  
  # Data curves: AUC curves colored by complexity
  geom_line(data = ppv_df,
            aes(x = sensitivity, y = precision, color = complexity, group = rank),
            size = 1) +
  
  # Points at observed values
  geom_point(data = ppv_df,
             aes(x = sensitivity, y = precision, color = complexity, shape = rank),
             size = 2) +
  
  # Facet by rank and damage
  facet_grid(damage ~ rank) +
  # facet_grid(damage ~ rank, labeller = labeller(
  #   damage = c("no damage" = "Control",
  #              "light damage" = "Mild",
  #              "heavy damage" = "Severe")
  # ))
  # 
  
  scale_color_viridis_d(option = "C") +
  labs(
    title = "Precision-Recall with F1 Score Contours",
    x = "Sensitivity (Recall)",
    y = "Precision (PPV)",
    color = "Complexity",
    shape = "Rank"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  )




# ROC-style PPV/AUC plot
p1 <- ggplot(ppv_df, aes(x = sensitivity, y = precision, 
                         color = complexity, group = complexity)) +
  geom_line(size = 1) +
  facet_wrap(~ damage) +
  geom_contour(aes(z = F1), color = "grey50", bins = 5, alpha = 0.5) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,1)) +
  theme_bw() +
  labs(
    title = "Precision vs Sensitivity (ROC-style)",
    x = "Sensitivity (Recall)",
    y = "Precision (PPV)"
  )
p1

# AUC summary plot
p2 <- ggplot(ppv_df, aes(x = damage, y = auc, 
                         color = complexity, group = complexity)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ rank) +
  theme_bw() +
  labs(
    title = "AUC by Damage level and Complexity",
    x = "Damamge level",
    y = "AUC"
  )
p2


####################
# Option 1: AUC curves with F1 contours (labelled)
ggplot(ppv_df, aes(x = sensitivity, y = precision, color = complexity, group = complexity)) +
  geom_path(size = 1) +
  geom_contour_filled(aes(z = F1), bins = 6, alpha = 0.3) +   # contour shading
  geom_text_contour(aes(z = F1), bins = 6, stroke = 0.2, size = 3) + # F1 labels
  facet_wrap(~damage) +
  scale_color_viridis_d(option = "C") + # <-- change complexity colours here
  labs(x = "Sensitivity (Recall)",
       y = "Precision (PPV)",
       color = "Complexity",
       fill = "F1 Score") +
  theme_minimal(base_size = 14)

# Option 2: PPV vs AUC curves with F1 contours (labelled)
ggplot(ppv_df, aes(x = auc, y = precision, color = complexity, group = complexity)) +
  geom_path(size = 1) +
  geom_contour_filled(aes(z = F1), bins = 6, alpha = 0.3) +
  geom_text_contour(aes(z = F1), bins = 6, stroke = 0.2, size = 3) +
  facet_wrap(~rank) +
  scale_color_viridis_d(option = "C") + # <-- change complexity colours here
  labs(x = "AUC",
       y = "Precision (PPV)",
       color = "Complexity",
       fill = "F1 Score") +
  theme_minimal(base_size = 14)

########## so clean up, different ideas for structuring, and information about diffiernet stuff
# ================================================================
# Precision-Recall Analysis Script
# Bundled plotting functions for ppv_df
# ================================================================

# ------------------------------------------------
# Example: relabel facets
# Adjust this as needed for your dataset
damage_labels <- c("no damage" = "0%", 
  "mudline_damage" = "3% (0 mbsf - 0 kya)", 
  "middle_damage" = "7% (14.55 mbsf - 14.5 kya)", 
  "bottom_damage" = "18% (61.5 mbsf - 214 kya)", 
  "all_damage" = "100%",
  "spiked_mudline_damage" = "combined 3% (0 mbsf - 0 kya)",
  "spiked_middle_damage" = "combined 7% (14.55 mbsf - 14.5 kya)",
  "spiked_bottom_damage" = "combined 18% (61.5 mbsf - 214 kya)")

# ------------------------------------------------
# Example: preferred colour mapping
# Replace with your desired categories
preferred_colors <- c(
  "shark"   = "darkgreen",
  "lobster" = "blue"
)

# ================================================================
# 1. Precision-Recall with F1 contours
# ------------------------------------------------
# Interpretation:
# - Shows trade-off between precision (PPV) and sensitivity (recall).
# - Dashed lines are F1 score contours (balance of precision and recall).
# - Curves represent model performance by complexity and rank.
# ================================================================

# Generate F1 contour data
f1_scores <- seq(0.1, 0.9, by = 0.1)
f1_contours <- expand.grid(
  sensitivity = seq(0.01, 1, by = 0.01),
  F1 = f1_scores
) %>%
  mutate(precision = (F1 * sensitivity) / (2 * sensitivity - F1)) %>%
  filter(precision >= 0 & precision <= 1)

p1 <- ggplot() +
  geom_line(data = f1_contours,
            aes(x = sensitivity, y = precision, group = F1),
            linetype = "dashed", color = "grey50") +
  geom_text(data = f1_contours %>% group_by(F1) %>% slice_tail(n = 1),
            aes(x = sensitivity, y = precision, label = paste0("F1=", F1)),
            hjust = -0.1, vjust = 0.5, size = 3, color = "grey40") +
  geom_line(data = ppv_df,
            aes(x = sensitivity, y = precision, color = complexity,
                group = interaction(complexity, rank)),
            size = 1) +
  geom_point(data = ppv_df,
             aes(x = sensitivity, y = precision, color = complexity, shape = rank),
             size = 2) +
  facet_grid(damage ~ rank, labeller = labeller(damage = damage_labels)) +
  scale_color_manual(values = preferred_colors) +
  labs(
    title = "Precision-Recall with F1 Score Contours",
    x = "Sensitivity (Recall)",
    y = "Precision (PPV)",
    color = "Complexity",
    shape = "Rank"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

# ================================================================
# 2. Simple Precision-Recall curves
# ------------------------------------------------
# Interpretation:
# - Similar to ROC, but better for imbalanced datasets.
# - Curves closer to the top-right indicate better trade-offs.
# ================================================================

p2 <- ggplot(ppv_df, aes(x = sensitivity, y = precision,
                         color = complexity,
                         group = interaction(complexity, rank))) +
  geom_line(size = 1) +
  geom_point(aes(shape = rank), size = 2) +
  facet_grid(damage ~ rank, labeller = labeller(damage = damage_labels)) +
  scale_color_manual(values = preferred_colors) +
  labs(
    title = "Precision-Recall Curves by Rank and Damage",
    x = "Sensitivity (Recall)",
    y = "Precision (PPV)",
    color = "Complexity",
    shape = "Rank"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

# ================================================================
# 3. Aggregate comparison plot (by complexity)
# ------------------------------------------------
# Interpretation:
# - Collapses across ranks to highlight complexity differences.
# - Useful when rank is secondary to complexity comparison.
# ================================================================

p3 <- ggplot(ppv_df, aes(x = sensitivity, y = precision,
                         color = complexity)) +
  geom_line(size = 1, alpha = 0.8) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = preferred_colors) +
  labs(
    title = "Aggregate Precision-Recall by Complexity",
    x = "Sensitivity (Recall)",
    y = "Precision (PPV)",
    color = "Complexity"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

# ================================================================
# Output: arrange plots or save individually
# ================================================================

# Example to display them together if you have patchwork:
# library(patchwork)
# p1 / p2 / p3
