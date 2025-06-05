# Load libraries
library(tidyverse)
library(viridis)

# Customize these paths and theme functions as needed
plots_dir <- "your/plots/directory/" # e.g., "plots/"
# source("path/to/your/custom_theme.R") # if using custom themes

# Replace with your own dataset
# Input dataframe must contain: read_name, file_name, assignment_level, read_count
input_data <- your_dataframe  # e.g., ha_reads

# ---------- FULL HEATMAP PLOT ----------
ggplot(input_data, aes(x = assignment_level, y = file_name, fill = read_count)) + 
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  facet_grid(~read_name, scales = "free", space = "free") +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Taxonomic Assignment Level", 
    y = "Input File",
    title = "Taxonomic Assignment Heatmap by Marker Gene",
    fill = "Read Count"
  )

ggsave(filename = paste0(plots_dir, "full_taxonomic_assignment_heatmap.png"), width = 12, height = 7, units = "in", dpi = 300)


# ---------- FILTERED HEATMAP (Shared Assignment Levels) ----------
# Filter to include only assignment levels present across all file_names for each read_name
filtered_data <- input_data %>%
  group_by(read_name) %>%
  mutate(total_files = n_distinct(file_name)) %>%
  group_by(read_name, assignment_level) %>%
  mutate(n_files = n_distinct(file_name)) %>%
  filter(n_files == total_files) %>%
  ungroup()

ggplot(filtered_data, aes(x = assignment_level, y = file_name, fill = read_count)) + 
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  facet_wrap(~read_name, scales = "free") +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Taxonomic Assignment Level", 
    y = "Input File",
    title = "Filtered Heatmap (Shared Assignments Only)",
    fill = "Read Count"
  )

ggsave(filename = paste0(plots_dir, "filtered_taxonomic_assignment_heatmap.png"), width = 5.5, height = 9, units = "in", dpi = 300)
