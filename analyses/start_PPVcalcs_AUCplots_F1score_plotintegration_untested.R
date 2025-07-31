#this is a draft script for Precision - Recall 

# Load required libraries
library(ggplot2)
library(dplyr)

# Replace with your actual dataframe name
# Expected columns: tool, rank, complexity (numeric), damage (character),
#                   precision (PPV), sensitivity (Recall)

# Set target rank for plotting
rank_to_plot <- "Species"

# Filter dataset for selected rank
df_rank <- metrics_df %>%
  filter(rank == rank_to_plot)

# Compute AUC (trapezoidal) and prepare tool label
pr_combined <- df_rank %>%
  arrange(tool, complexity, damage, sensitivity) %>%
  group_by(tool, complexity, damage) %>%
  mutate(
    auc = sum(diff(sensitivity) * head(precision, -1), na.rm = TRUE),
    tool_label = paste0(tool, " (AUC=", round(unique(auc), 3), ")")
  ) %>%
  ungroup()

# Define a function to generate F1 score contour data
get_f1_contours <- function(f1_scores = c(0.2, 0.4, 0.6, 0.8, 0.9)) {
  recall_vals <- seq(0.01, 1, by = 0.01)
  contour_data <- do.call(rbind, lapply(f1_scores, function(f1) {
    precision_vals <- (f1 * recall_vals) / (2 * recall_vals - f1)
    df <- data.frame(sensitivity = recall_vals,
                     precision = precision_vals,
                     F1 = f1)
    df <- df[precision_vals <= 1 & precision_vals >= 0, ]  # Keep valid values only
    return(df)
  }))
  return(contour_data)
}

# Generate contour lines for F1 scores
f1_df <- get_f1_contours()

# Convert complexity to factor with ordered levels for proper facet ordering
pr_combined$complexity <- factor(pr_combined$complexity,
                                 levels = c(10, 25, 50, 100, 250, 500, 1000, 2500, 5000))

# Plot with PR curves + F1 contours
ggplot(pr_combined, aes(x = sensitivity, y = precision, color = tool_label, group = tool_label)) +
  # PR curves and points
  geom_line(size = 1.2) +
  geom_point(size = 2, alpha = 0.8) +
  
  # F1 contours
  geom_line(data = f1_df, aes(x = sensitivity, y = precision, group = F1),
            inherit.aes = FALSE, linetype = "dashed", color = "gray50") +
  
  # Labels for F1 contours
  geom_text(data = f1_df %>% group_by(F1) %>% filter(sensitivity == max(sensitivity)),
            aes(x = sensitivity, y = precision, label = paste0("F1=", F1)),
            inherit.aes = FALSE, hjust = 1.1, vjust = 0.5, size = 3, color = "gray30") +
  
  # Facet by numeric complexity (rows) and damage (columns)
  facet_grid(rows = vars(complexity), cols = vars(damage)) +
  
  # Labels and theme
  labs(title = paste("Precision-Recall Curves with F1 Score Contours at", rank_to_plot, "Level"),
       subtitle = "Faceted by Species Complexity (numeric) and Deamination Damage",
       x = "Sensitivity (Recall)",
       y = "Precision (PPV)",
       color = "Tool (AUC)") +
  theme_minimal(base_size = 14) +
  theme(strip.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.position = "bottom",
        legend.title = element_text(size = 12))
