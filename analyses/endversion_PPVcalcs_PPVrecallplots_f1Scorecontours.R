#this is a draft script for PPV plotting and visualization

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(readxl)
library(viridis)
library(metR)
library(grid)
library(patchwork)
library(ggrepel)

################# publication themes for plots 
mytheme <-   theme(legend.position = "bottom", 
                   axis.text.x = element_text(size = 12, angle = 45), 
                   legend.text = element_text(size = 12), 
                   legend.title = element_text(size = 15), 
                   axis.text.y = element_text(size = 12), 
                   axis.title = element_text(size =12), 
                   strip.text.x.top = element_text(size = 12), 
                   legend.key.width = unit(0.5, "cm"),
                   legend.spacing.x = unit(0.2, 'cm'))

mytheme_test <- theme(
  legend.position = "bottom",
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 12),
  legend.key.width = unit(0.7, "cm"),
  legend.spacing.x = unit(0.2, 'cm'),
  
  axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
  axis.text.y = element_text(size = 10),
  axis.title.x = element_text(size = 12, vjust = -0.5),
  axis.title.y = element_text(size = 12, vjust = 1.5),
  
  strip.text = element_text(size = 11) #, face = "bold"),
  # strip.background = element_rect(fill = "gray90", color = NA)
)

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

mytheme_nofacet <-   theme(legend.position = "bottom", 
                           axis.text.x = element_text(size = 12, angle = 45), 
                           legend.text = element_text(size = 12), 
                           legend.title = element_text(size = 14), 
                           axis.text.y = element_text(size = 12), 
                           axis.title = element_text(size =12), 
                           legend.key.width = unit(2.5, "cm"),
                           legend.spacing.x = unit(0.3, 'cm'))

mytheme_test_nofacet <- theme(
  legend.position = "bottom",
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 12),
  legend.key.width = unit(1.5, "cm"),
  legend.spacing.x = unit(0.2, 'cm'),
  
  axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
  axis.text.y = element_text(size = 12),
  axis.title.x = element_text(size = 12, vjust = -0.5),
  axis.title.y = element_text(size = 12, vjust = 1.5))
#goal make versions of plots that mimic the ones in manuscript
#comparing all synthetic data
#comparing all spiked data
#direct synthetic vs spiked at same damage type


mytheme <-   theme(legend.position = "bottom", 
                   axis.text.x = element_text(size = 12, angle = 45), 
                   legend.text = element_text(size = 12), 
                   legend.title = element_text(size = 15), 
                   axis.text.y = element_text(size = 12), 
                   axis.title = element_text(size =12), 
                   strip.text.x.top = element_text(size = 12), 
                   legend.key.width = unit(0.5, "cm"),
                   legend.spacing.x = unit(0.2, 'cm'))
# Replace with your actual dataframe name
# Expected columns: tool, rank, complexity (numeric), damage (character),
#                   precision (PPV), sensitivity (Recall)
##########excel sheet is sheet 2 of supplementary tables, but a copy within R project folder #########
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
                        levels = c("no_damage", "mudline_damage","spiked_mudline_damage", "middle_damage", "spiked_middle_damage", "bottom_damage",
                                   "spiked_bottom_damage", "all_damage"))
# # -------------------------------
# # 6. Calculate AUC for each combination (this isn't actually helpful at all, because there is not PR curve to get the area under, just summarised values)
# # -------------------------------
# ppv_df <- ppv_df %>%
#   arrange(complexity, damage, sensitivity) %>%
#   group_by(complexity, damage) %>%
#   mutate(
#     auc = sum(diff(sensitivity) * head(precision, -1), na.rm = TRUE))%>%
#   ungroup()
# 

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


# ------------------------------------------------
# Example: relabel facets
# Adjust this as needed for your dataset
damage_labels <- c("no_damage" = "0%", 
                   "mudline_damage" = "3% (0 mbsf - 0 kya)", 
                   "middle_damage" = "7% (14.55 mbsf - 14.5 kya)", 
                   "bottom_damage" = "18% (61.5 mbsf - 214 kya)", 
                   "all_damage" = "100%",
                   "spiked_mudline_damage" = "combined 3%",
                   "spiked_middle_damage" = "combined 7% ",
                   "spiked_bottom_damage" = "combined 18%")


########## so clean up, different ideas for structuring, and information about diffiernet stuff

# ================================================================
# 1. Precision-Recall with F1 contours
# ------------------------------------------------
# Interpretation:
# - Shows trade-off between precision (PPV) and sensitivity (recall).
# - Dashed lines are F1 score contours (balance of precision and recall).
# - Curves represent model performance by complexity and rank.
# ================================================================

# Plot: Precision vs Sensitivity with F1 contours and AUC curves (currently commented out)
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
    hjust = 0.9, vjust = 0.9, size = 2, color = "grey40")+
  # Data curves: AUC curves colored by complexity
  # geom_line(data = ppv_df,
  #           aes(x = sensitivity, y = precision, color = complexity),
  #           size = 1) +
  
  # Points at observed values
  geom_point(data = ppv_df,
             aes(x = sensitivity, y = precision, color = complexity),
             size = 2) +
  
  # Facet by rank and damage
  facet_grid(rank ~ damage, labeller = labeller(damage = damage_labels)) +
  
  scale_color_viridis_d(option = "C") +
  labs(
    #title = "Precision-Recall with F1 Score Contours",
    x = "Sensitivity (Recall)",
    y = "Precision (PPV)",
    color = "Complexity"
  ) +
  mytheme_bigheatmap_facetgrid_light +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom", 
    strip.text = element_text(size = 9))


# # AUC summary plot ### fix naming of x axis 
# ggplot(ppv_df, aes(x = damage, y = auc, 
#                    color = complexity, group = complexity)) +
#   geom_line(size = 1) +
#   geom_point(size = 2) +
#   facet_wrap(~ rank, labeller = labeller(damage = damage_labels)) +
#   mytheme_bigheatmap_facetgrid +
#   labs(
#     title = "AUC by Damage level and Complexity",
#     x = "Damamge level",
#     y = "AUC"
#   )



# ================================================================
# 2. Simple Precision-Recall curves -----> not possible with current summarised data counts, would need to go back to raw data
# ------------------------------------------------
# Interpretation:
# - Similar to ROC, but better for imbalanced datasets.
# - Curves closer to the top-right indicate better trade-offs.
# ================================================================

# ggplot(ppv_df, aes(x = sensitivity, y = precision,
#                         color = complexity,
#                         group = interaction(complexity, rank))) +
#  geom_line(size = 1) +
#  geom_point(aes(shape = rank), size = 2) +
#  facet_grid(damage ~ rank, labeller = labeller(damage = damage_labels)) +
#  scale_color_viridis_d(option = "C") +
#  labs(
#    title = "Precision-Recall Curves by Rank and Damage",
#    x = "Sensitivity (Recall)",
#    y = "Precision (PPV)",
#    color = "Complexity",
#    shape = "Rank"
#  ) +
#  theme_bw() +
#  theme(panel.grid = element_blank(),
#        legend.position = "bottom")

# ============================================================
# PPV / Recall Plotting Script with F1 Contours
# ============================================================

# -------------------------
# Subset Data
# -------------------------

# 1. Synthetic-only damage types
synthetic_df <- ppv_df %>%
  filter(damage %in% c("no_damage", "all_damage", 
                       "mudline_damage", "middle_damage", "bottom_damage"))
# Labels for synthetic plot
synthetic_damage_labels <- c("no_damage" = "0%", 
                             "mudline_damage" = "3% (0 mbsf - 0 kya)", 
                             "middle_damage" = "7% (14.55 mbsf - 14.5 kya)", 
                             "bottom_damage" = "18% (61.5 mbsf - 214 kya)", 
                             "all_damage" = "100%")

# 2. Spiked-only damage types
spiked_df <- ppv_df %>%
  filter(damage %in% c("spiked_mudline_damage", 
                       "spiked_middle_damage", "spiked_bottom_damage"))
# Labels for spiked plot
spiked_damage_labels <- c("spiked_mudline_damage" = "combined 3% (0 mbsf - 0 kya",
                          "spiked_middle_damage" = "combined 7% (14.55 mbsf - 14.5 kya",
                          "spiked_bottom_damage" = "combined 18% (61.5 mbsf - 214 kya)")

# 3. Paired synthetic vs spiked damage
paired_df <- ppv_df %>%
  filter(damage %in% c("mudline_damage", "spiked_mudline_damage",
                       "middle_damage", "spiked_middle_damage",
                       "bottom_damage", "spiked_bottom_damage")) %>%
  mutate(pair = case_when(
    damage %in% c("mudline_damage", "spiked_mudline_damage") ~ "mudline",
    damage %in% c("middle_damage", "spiked_middle_damage") ~ "middle",
    damage %in% c("bottom_damage", "spiked_bottom_damage") ~ "bottom"
  ))

# Labels for paired plot
paired_damage_labels <- c(
  "mudline_damage"        = "3% (0 mbsf – 0 kya)",
  "spiked_mudline_damage" = "Combined 3%",
  "middle_damage"         = "7% (14.55 mbsf – 14.5 kya)",
  "spiked_middle_damage"  = "Combined 7%",
  "bottom_damage"         = "18% (61.5 mbsf – 214 kya)",
  "spiked_bottom_damage"  = "Combined 18%")
# ============================================================
# PLOTS ;)
# ============================================================

############## 1. Synthetic-only ##################
p_synthetic <- 
  ggplot() +
  # F1 score contours
  geom_line(data = f1_contours,
            aes(x = sensitivity, y = precision, group = F1),
            linetype = "dashed", color = "black") +
  # geom_text(data = f1_contours %>% group_by(F1) %>% slice_tail(n=1),
  #           aes(x = sensitivity, y = precision, label = paste0("F1=", F1)),
  #           hjust = -0.1, vjust = 0.5, size = 3, color = "grey40") +
  geom_text(
    data = f1_contours %>% group_by(F1) %>% slice_tail(n = 1),
    aes(x = sensitivity, y = precision, label = paste0("F1=", F1)),
    hjust = 0.8, vjust = -0.6, size = 3, color = "black")+
  # Data curves: AUC curves colored by complexity
  # geom_line(data = ppv_df,
  #           aes(x = sensitivity, y = precision, color = complexity),
  #           size = 1) +
  
  # Points at observed values
  geom_point(data = synthetic_df,
             aes(x = sensitivity, y = precision, color = complexity),
             size = 4, position = "dodge") +
  
  # Facet by rank and damage
  facet_grid(rank~ damage, labeller = labeller(damage = synthetic_damage_labels)) +
  
  scale_color_viridis_d(option = "C") +
  labs(
    #title = "Precision-Recall with F1 Score Contours",
    x = "Sensitivity (Recall)",
    y = "Precision (PPV)",
    color = "Complexity") +
  mytheme_bigheatmap_facetgrid_light +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom", 
    strip.text = element_text(size = 12))

print(p_synthetic)
###################### 2. Spiked-only ################
p_spiked <- 
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
    hjust = 0.9, vjust = 0.9, size = 2, color = "grey40")+
  # Data curves: AUC curves colored by complexity
  # geom_line(data = ppv_df,
  #           aes(x = sensitivity, y = precision, color = complexity),
  #           size = 1) +
  
  # Points at observed values
  geom_point(data = spiked_df,
             aes(x = sensitivity, y = precision, color = complexity),
             size = 2) +
  
  # Facet by rank and damage
  facet_grid(rank ~ damage, labeller = labeller(damage = spiked_damage_labels)) +
  
  scale_color_viridis_d(option = "C") +
  labs(
    #title = "Precision-Recall with F1 Score Contours",
    x = "Sensitivity (Recall)",
    y = "Precision (PPV)",
    color = "Complexity"
  ) +
  mytheme_bigheatmap_facetgrid +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom", 
    strip.text = element_text(size = 9))

################## 3. Paired comparisons #####################

p_paired <-
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
    hjust = 0.9, vjust = 0.9, size = 2, color = "grey40")+
  # Data curves: AUC curves colored by complexity
  # geom_line(data = ppv_df,
  #           aes(x = sensitivity, y = precision, color = complexity),
  #           size = 1) +
  
  # Points at observed values
  geom_point(data = paired_df,
             aes(x = sensitivity, y = precision, color = complexity),
             size = 2) +
  
  # Facet by rank and damage
  facet_grid(rank ~ damage, labeller = labeller(damage = paired_damage_labels)) +
  
  scale_color_viridis_d(option = "C") +
  labs(
    #title = "Precision-Recall with F1 Score Contours",
    x = "Sensitivity (Recall)",
    y = "Precision (PPV)",
    color = "Complexity"
  ) +
  mytheme_bigheatmap_facetgrid +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom", 
    strip.text = element_text(size = 9))
# ============================================================
# INTERPRETATION NOTES (COMMENTED)
# ============================================================
# - Synthetic plots: compare performance of classifiers under different 
#   artificial ("clean") damage scenarios. Helps establish baseline.
# - Spiked plots: show effect of additional spiked complexity.
# - Paired plots: directly contrast synthetic vs spiked for mudline, middle, 
#   and bottom. Useful for identifying performance degradation or improvement.
# - F1 contours: lines of equal balance between PPV and recall (via F1).
#   Higher F1 contour levels indicate better overall performance.
# - Color (complexity): complexity level of analysis (e.g., species/genus/family).
#   Comparing colors within a facet shows how complexity influences outcomes.

# ============================================================
# Example: Print plots
# ============================================================
print(p_synthetic)
print(p_spiked)
print(p_paired)

# ================================================================
# Output: arrange plots or save individually
# ================================================================

# Example to display them together if you have patchwork:
p_synthetic / p_spiked / p_paired

###################### extra version ##############
#this is a draft script for PPV plotting and visualization

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(readxl)
library(viridis)
library(metR)
library(grid)
library(patchwork)
library(ggrepel)

################# publication themes for plots 
mytheme_bigheatmap_facetgrid <- theme_minimal() +
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
    plot.background  = element_rect(fill = "white", color = NA)
  )

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

#goal make versions of plots that mimic the ones in manuscript
#comparing all synthetic data
#comparing all spiked data
#direct synthetic vs spiked at same damage type


mytheme <-   theme(legend.position = "bottom", 
                   axis.text.x = element_text(size = 12, angle = 45), 
                   legend.text = element_text(size = 12), 
                   legend.title = element_text(size = 15), 
                   axis.text.y = element_text(size = 12), 
                   axis.title = element_text(size =12), 
                   strip.text.x.top = element_text(size = 12), 
                   legend.key.width = unit(0.5, "cm"),
                   legend.spacing.x = unit(0.2, 'cm'))
# Replace with your actual dataframe name
# Expected columns: tool, rank, complexity (numeric), damage (character),
#                   precision (PPV), sensitivity (Recall)
##########excel sheet is sheet 2 of supplementary tables, but a copy within R project folder #########
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
                        levels = c("no_damage", "mudline_damage","spiked_mudline_damage", "middle_damage", "spiked_middle_damage", "bottom_damage",
                                   "spiked_bottom_damage", "all_damage"))
# -------------------------------
# 6. Calculate AUC for each combination (this isn't actually helpful at all, because there is not PR curve to get the area under, just summarised values)
# -------------------------------
ppv_df <- ppv_df %>%
  arrange(complexity, damage, sensitivity) %>%
  group_by(complexity, damage) %>%
  mutate(
    auc = sum(diff(sensitivity) * head(precision, -1), na.rm = TRUE))%>%
  ungroup()


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


# ------------------------------------------------
# Example: relabel facets
# Adjust this as needed for your dataset
damage_labels <- c("no_damage" = "0%", 
                   "mudline_damage" = "3% (0 mbsf - 0 kya)", 
                   "middle_damage" = "7% (14.55 mbsf - 14.5 kya)", 
                   "bottom_damage" = "18% (61.5 mbsf - 214 kya)", 
                   "all_damage" = "100%",
                   "spiked_mudline_damage" = "combined 3%",
                   "spiked_middle_damage" = "combined 7% ",
                   "spiked_bottom_damage" = "combined 18%")


########## so clean up, different ideas for structuring, and information about diffiernet stuff

# ================================================================
# 1. Precision-Recall with F1 contours
# ------------------------------------------------
# Interpretation:
# - Shows trade-off between precision (PPV) and sensitivity (recall).
# - Dashed lines are F1 score contours (balance of precision and recall).
# - Curves represent model performance by complexity and rank.
# ================================================================

# Plot: Precision vs Sensitivity with F1 contours and AUC curves (currently commented out)
ggplot() +
  # F1 score contours
  geom_line(data = f1_contours,
            aes(x = sensitivity, y = precision, group = F1),
            linetype = "dashed", color = "black") +
  # geom_text(data = f1_contours %>% group_by(F1) %>% slice_tail(n=1),
  #           aes(x = sensitivity, y = precision, label = paste0("F1=", F1)),
  #           hjust = -0.1, vjust = 0.5, size = 3, color = "grey40") +
  geom_text(
    data = f1_contours %>% group_by(F1) %>% slice_tail(n = 1),
    aes(x = sensitivity, y = precision, label = paste0("F1=", F1)),
    hjust = 0.8, vjust = -0.6, size = 2, color = "black")+
  # Data curves: AUC curves colored by complexity
  # geom_line(data = ppv_df,
  #           aes(x = sensitivity, y = precision, color = complexity),
  #           size = 1) +
  
  # Points at observed values
  geom_point(data = ppv_df,
             aes(x = sensitivity, y = precision, color = complexity),
             size = 3) +
  
  # Facet by rank and damage
  facet_grid(rank ~ damage, labeller = labeller(damage = damage_labels)) +
  
  scale_color_viridis_d(option = "C") +
  labs(
    #title = "Precision-Recall with F1 Score Contours",
    x = "Sensitivity (Recall)",
    y = "Precision (PPV)",
    color = "Complexity"
  ) +
  mytheme_bigheatmap_facetgrid_light +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom", 
    strip.text = element_text(size = 9))



# ================================================================
# 2. Simple Precision-Recall curves -----> not possible with current summarised data counts, would need to go back to raw data
# ------------------------------------------------
# Interpretation:
# - Similar to ROC, but better for imbalanced datasets.
# - Curves closer to the top-right indicate better trade-offs.
# ================================================================
# AUC summary plot ### fix naming of x axis 
# ggplot(ppv_df, aes(x = damage, y = auc, 
#                    color = complexity, group = complexity)) +
#   geom_line(size = 1) +
#   geom_point(size = 2) +
#   facet_wrap(~ rank, labeller = labeller(damage = damage_labels)) +
#   mytheme_bigheatmap_facetgrid +
#   labs(
#     title = "AUC by Damage level and Complexity",
#     x = "Damamge level",
#     y = "AUC"
#   )
# 

# ggplot(ppv_df, aes(x = sensitivity, y = precision,
#                         color = complexity,
#                         group = interaction(complexity, rank))) +
#  geom_line(size = 1) +
#  geom_point(aes(shape = rank), size = 2) +
#  facet_grid(damage ~ rank, labeller = labeller(damage = damage_labels)) +
#  scale_color_viridis_d(option = "C") +
#  labs(
#    title = "Precision-Recall Curves by Rank and Damage",
#    x = "Sensitivity (Recall)",
#    y = "Precision (PPV)",
#    color = "Complexity",
#    shape = "Rank"
#  ) +
#  theme_bw() +
#  theme(panel.grid = element_blank(),
#        legend.position = "bottom")

# ============================================================
# PPV / Recall Plotting Script with F1 Contours
# ============================================================

# -------------------------
# Subset Data
# -------------------------

# 1. Synthetic-only damage types
synthetic_df <- ppv_df %>%
  filter(damage %in% c("no_damage", "all_damage", 
                       "mudline_damage", "middle_damage", "bottom_damage"))
# Labels for synthetic plot
synthetic_damage_labels <- c("no_damage" = "0%", 
                             "mudline_damage" = "3% (0 mbsf - 0 kya)", 
                             "middle_damage" = "7% (14.55 mbsf - 14.5 kya)", 
                             "bottom_damage" = "18% (61.5 mbsf - 214 kya)", 
                             "all_damage" = "100%")

# 2. Spiked-only damage types
spiked_df <- ppv_df %>%
  filter(damage %in% c("spiked_mudline_damage", 
                       "spiked_middle_damage", "spiked_bottom_damage"))
# Labels for spiked plot
spiked_damage_labels <- c("spiked_mudline_damage" = "combined 3% (0 mbsf - 0 kya",
                          "spiked_middle_damage" = "combined 7% (14.55 mbsf - 14.5 kya",
                          "spiked_bottom_damage" = "combined 18% (61.5 mbsf - 214 kya)")

# 3. Paired synthetic vs spiked damage
paired_df <- ppv_df %>%
  filter(damage %in% c("mudline_damage", "spiked_mudline_damage",
                       "middle_damage", "spiked_middle_damage",
                       "bottom_damage", "spiked_bottom_damage")) %>%
  mutate(pair = case_when(
    damage %in% c("mudline_damage", "spiked_mudline_damage") ~ "mudline",
    damage %in% c("middle_damage", "spiked_middle_damage") ~ "middle",
    damage %in% c("bottom_damage", "spiked_bottom_damage") ~ "bottom"
  ))

# Labels for paired plot
paired_damage_labels <- c(
  "mudline_damage"        = "Synthetic 3% (0 mbsf – 0 kya)",
  "spiked_mudline_damage" = "Combined 3%",
  "middle_damage"         = "Synthetic 7% (14.55 mbsf – 14.5 kya)",
  "spiked_middle_damage"  = "Combined 7% ",
  "bottom_damage"         = "Synthetic 18% (61.5 mbsf – 214 kya)",
  "spiked_bottom_damage"  = "Combined 18%")
# ============================================================
# PLOTS ;)
# ============================================================

############## 1. Synthetic-only ##################
p_synthetic <- 
  ggplot() +
  # F1 score contours
  geom_line(data = f1_contours,
            aes(x = sensitivity, y = precision, group = F1),
            linetype = "dashed", color = "black") +
  # geom_text(data = f1_contours %>% group_by(F1) %>% slice_tail(n=1),
  #           aes(x = sensitivity, y = precision, label = paste0("F1=", F1)),
  #           hjust = -0.1, vjust = 0.5, size = 3, color = "grey40") +
  geom_text(
    data = f1_contours %>% group_by(F1) %>% slice_tail(n = 1),
    aes(x = sensitivity, y = precision, label = paste0("F1=", F1)),
    hjust = 0.8, vjust = -0.6, size = 3, color = "black")+
  # Data curves: AUC curves colored by complexity
  # geom_line(data = ppv_df,
  #           aes(x = sensitivity, y = precision, color = complexity),
  #           size = 1) +
  
  # Points at observed values
  geom_point(data = synthetic_df,
             aes(x = sensitivity, y = precision, color = complexity),
             size = 4, position = "dodge") +
  
  # Facet by rank and damage
  facet_grid(rank~ damage, labeller = labeller(damage = synthetic_damage_labels)) +
  
  scale_color_viridis_d(option = "C") +
  labs(
    #title = "Precision-Recall with F1 Score Contours",
    x = "Sensitivity (Recall)",
    y = "Precision (PPV)",
    color = "Complexity") +
  mytheme_bigheatmap_facetgrid_light +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom", 
    strip.text = element_text(size = 12))

print(p_synthetic)
###################### 2. Spiked-only ################
p_spiked <- 
  ggplot() +
  # F1 score contours
  geom_line(data = f1_contours,
            aes(x = sensitivity, y = precision, group = F1),
            linetype = "dashed", color = "black") +
  # geom_text(data = f1_contours %>% group_by(F1) %>% slice_tail(n=1),
  #           aes(x = sensitivity, y = precision, label = paste0("F1=", F1)),
  #           hjust = -0.1, vjust = 0.5, size = 3, color = "grey40") +
  geom_text(
    data = f1_contours %>% group_by(F1) %>% slice_tail(n = 1),
    aes(x = sensitivity, y = precision, label = paste0("F1=", F1)),
    hjust = 0.8, vjust = -0.6, size = 3, color = "black")+
  # Data curves: AUC curves colored by complexity
  # geom_line(data = ppv_df,
  #           aes(x = sensitivity, y = precision, color = complexity),
  #           size = 1) +
  
  # Points at observed values
  geom_point(data = spiked_df,
             aes(x = sensitivity, y = precision, color = complexity),
             size = 4) +
  
  # Facet by rank and damage
  facet_grid(rank ~ damage, labeller = labeller(damage = spiked_damage_labels)) +
  
  scale_color_viridis_d(option = "C") +
  labs(
    #title = "Precision-Recall with F1 Score Contours",
    x = "Sensitivity (Recall)",
    y = "Precision (PPV)",
    color = "Complexity"
  ) +
  mytheme_bigheatmap_facetgrid_light +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom", 
    strip.text = element_text(size = 12))
print(p_spiked)
################## 3. Paired comparisons #####################

p_paired <-
  ggplot() +
  # F1 score contours
  geom_line(data = f1_contours,
            aes(x = sensitivity, y = precision, group = F1),
            linetype = "dashed", color = "black") +
  # geom_text(data = f1_contours %>% group_by(F1) %>% slice_tail(n=1),
  #           aes(x = sensitivity, y = precision, label = paste0("F1=", F1)),
  #           hjust = -0.1, vjust = 0.5, size = 3, color = "grey40") +
  geom_text(
    data = f1_contours %>% group_by(F1) %>% slice_tail(n = 1),
    aes(x = sensitivity, y = precision, label = paste0("F1=", F1)),
    hjust = 0.8, vjust = -0.6, size = 3, color = "black")+
  # Data curves: AUC curves colored by complexity
  # geom_line(data = ppv_df,
  #           aes(x = sensitivity, y = precision, color = complexity),
  #           size = 1) +
  
  # Points at observed values
  geom_point(data = paired_df,
             aes(x = sensitivity, y = precision, color = complexity),
             size = 3) +
  
  # Facet by rank and damage
  facet_grid(rank ~ damage, labeller = labeller(damage = paired_damage_labels)) +
  
  scale_color_viridis_d(option = "C") +
  labs(
    #title = "Precision-Recall with F1 Score Contours",
    x = "Sensitivity (Recall)",
    y = "Precision (PPV)",
    color = "Complexity"
  ) +
  mytheme_bigheatmap_facetgrid_light +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom", 
    strip.text = element_text(size = 9))
print(p_paired)
# ============================================================
# INTERPRETATION NOTES (COMMENTED)
# ============================================================
# - Synthetic plots: compare performance of classifiers under different 
#   artificial ("clean") damage scenarios. Helps establish baseline.
# - Spiked plots: show effect of additional spiked complexity.
# - Paired plots: directly contrast synthetic vs spiked for mudline, middle, 
#   and bottom. Useful for identifying performance degradation or improvement.
# - F1 contours: lines of equal balance between PPV and recall (via F1).
#   Higher F1 contour levels indicate better overall performance.
# - Color (complexity): complexity level of analysis (e.g., species/genus/family).
#   Comparing colors within a facet shows how complexity influences outcomes.

# ============================================================
# Example: Print plots
# ============================================================
print(p_synthetic)
print(p_spiked)
print(p_paired)

# ================================================================
# Output: arrange plots or save individually
# ================================================================

# Example to display them together if you have patchwork:
p_synthetic / p_spiked / p_paired