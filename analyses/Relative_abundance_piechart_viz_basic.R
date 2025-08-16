#visualisation script for Linda presentation, to be in the style of 2021 pie chart 
#relative abunbances of taxonomic assignments across all samples

#this script was created with the help of chatgpt
#required packages
library(tidyverse)
library(ggrepel)
#load in summarised .txt file from megan @desired taxonomic level 
phy_sum <- read_delim("C:/Users/davisee/OneDrive - University of Tasmania/Documents/Bioinformatics_Documentation/28_08_2024-MARES_BAR_PR2v5_20240216/MALT/Comparison_IODP_CLEAN_phylumlevel_counts.txt", delim = "\t", col_names = TRUE)
head(phy_sum)

#sum taxa counts across all datasets
phy_sum_totals <- phy_sum %>%
  mutate(Total = rowSums(across(-`#Datasets`))) %>%
  select(Taxa = `#Datasets`, Total) %>%
  arrange(desc(Total))

#convert to relative abundances (%)
phy_sum_rel <- phy_sum_totals %>%
  mutate(Rel_abun = Total/ sum(Total) * 100) 

ggplot(phy_sum_rel, aes(x = "", y = Rel_abun, fill = Taxa)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  labs(title = "Relative Abundance of Phyla") +
  theme_void() +
  theme(legend.position = "none")


# Identify top 10 taxa + "Chordata"
top_taxa <- phy_sum_rel %>%
  arrange(desc(Rel_abun)) %>%
  slice_head(n = 5) %>%
  pull(Taxa)

taxa_to_label <- union(top_taxa, "Chordata")

# Create label column: show only for selected taxa
phy_sum_rel <- phy_sum_rel %>%
  mutate(label = ifelse(Taxa %in% taxa_to_label,
                        paste0(Taxa, " (", round(Rel_abun, 1), "%)"),
                        ""))

# Pie chart with top ten samples and "Chordata" labeled
ggplot(phy_sum_rel, aes(x = "", y = Rel_abun, fill = Taxa)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  # geom_text(aes(label = label),
  #           position = position_stack(vjust = 0.5), size = 3) +
  theme_void() +
  theme(legend.position = "none")

  )
#############################################
library(tidyverse)
library(ggrepel)

# Start from your phy_sum table
phy_sum_totals <- phy_sum %>%
  mutate(Total = rowSums(across(-`#Datasets`))) %>%
  select(Taxa = `#Datasets`, Total) %>%
  arrange(desc(Total)) %>%
  mutate(Rel_abun = Total / sum(Total) * 100)

# Top 10 + Chordata
top_taxa <- phy_sum_totals %>%
  slice_max(order_by = Rel_abun, n = 5) %>%
  pull(Taxa)

taxa_to_label <- union(top_taxa, "Chordata")

# Compute slice positions
plot_data <- phy_sum_totals %>%
  arrange(Rel_abun) %>%
  mutate(
    ymax = cumsum(Rel_abun),
    ymin = lag(ymax, default = 0),
    label_pos = (ymax + ymin) / 2,
    label = ifelse(Taxa %in% taxa_to_label,
                   paste0(Taxa, " (", round(Rel_abun, 1), "%)"),
                   NA)
  )

# Pie chart with outside labels
ggplot(plot_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = Taxa)) +
  geom_rect(color = "white") +
  coord_polar(theta = "y") +
  xlim(c(2, 6)) +  # gives more space for labels
  # geom_text_repel(
  #   aes(
  #     x = 4.5,
  #     y = label_pos,
  #     label = label
  #   ),
  #   segment.size = 0.5,
  #   na.rm = TRUE,
  #   direction = "y",
  #   hjust = 0
  # ) +
  theme_void() +
  theme(legend.position = "none")

