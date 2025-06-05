library(vegan)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(readr)
library(stringr)

# 1. Read the data from CSV file
setwd("C:/Users/davisee/OneDrive - University of Tasmania/Documents/Bioinformatics_Documentation/upscale_runs/newparams/single_runs")
data <- read_tsv("Comparison_allnonspiked_genus_summarised.txt")
data_fam <- read_tsv("Comparison_allnonspiked_family_summarised.txt")
data_spec <- read_tsv("")

# 2. Remove the Genus column 
data_no_genus <- data[, -which(names(data) == "#Datasets")]

# 3. Apply Shannon Diversity Index calculation for each column
shannon_indices <- apply(data_no_genus, 2, function(x) diversity(x, index = "shannon"))

# 4. Arrange the indices from lowest to highest 
df_plot <- data.frame(Method = names(shannon_indices), ShannonIndex = shannon_indices) %>%
  arrange(ShannonIndex)

#4.5 clean up and rename columns for better readership from pipeline testing
df_plot <- df_plot %>%
      mutate(
    Method = str_replace(Method, "newparams", "n10"), 
    Method = str_replace(Method, "n1012092024", "n10"))

df_plot <- df_plot %>%
  mutate(# Replace "newparams" with "n10"
    n_value = str_extract(Method, "n\\d+"),                          # Extract the "n" value (e.g., "n500", "n1000")
    sample_type = str_extract(Method, "(deam_)?(iterative|\\d+h\\d)"), # Extract "deam_iterative", "deam_1h1", etc.
    value = str_extract(n_value,"\\d+" )# extract just integer value with "n"
      )

row.names(df_plot) <- NULL

df_plot$sample_type <- replace(df_plot$sample_type, df_plot$sample_type == "iterative", "no_damage")
df_plot[df_plot$Method == "deam_n10", "sample_type"] <- "deam"
df_plot[df_plot$Method == "n10", "sample_type"] <- "no_damage"


df_plot <- df_plot %>% 
  mutate(sample_type = factor(sample_type, 
                            levels = c("no_damage", 
                                       "1h1", 
                                       "2h6", 
                                       "8h2",
                                       "deam")))
df_plot$value <- as.integer(df_plot$value)

#5. Plot a connected scatter plot 
library(ggplot2)
 
mytheme_nofacet <-   theme(legend.position = "bottom", 
                           axis.text.x = element_text(size = 12, angle = 45), 
                           legend.text = element_text(size = 12), 
                           legend.title = element_text(size = 14), 
                           axis.text.y = element_text(size = 12), 
                           axis.title = element_text(size =12), 
                           legend.key.width = unit(0.5, "cm"),
                           legend.spacing.x = unit(0.2, 'cm'))

 ggplot(df_plot, aes(x = reorder(n_value, value), y = ShannonIndex, group = sample_type)) +
   geom_point(aes(color = sample_type), size = 2) +
   scale_color_brewer(name = "% Deaminated fragments", 
                     labels = c("0%", "3% - mudline", "7% - middle", "18% - bottom", "100%"),
                     palette = "Set2") + 
   geom_line(aes(color = sample_type)) + # custom colors + 
   labs(x= "Input Fragment size per species", y = "Shannon Index", title = "Shannon Diversity Index of Synthetic Data - Genus level") +
   theme_minimal()+ 
    mytheme_nofacet
 
 ##### family plot ###

 # 2. Remove the Genus column 
 data_no_family <- data_fam[, -which(names(data_fam) == "#Datasets")]
 
 # 3. Apply Shannon Diversity Index calculation for each column
 shannon_indices_fam <- apply(data_no_family, 2, function(x) diversity(x, index = "shannon"))
 
 # 4. Arrange the indices from lowest to highest 
 df_plot_fam <- data.frame(Method = names(shannon_indices_fam), ShannonIndex = shannon_indices_fam) %>%
   arrange(ShannonIndex)
 
 #4.5 clean up and rename columns for better readership from pipeline testing
 df_plot_fam <- df_plot_fam %>%
   mutate(
     Method = str_replace(Method, "newparams", "n10"), 
     Method = str_replace(Method, "n1012092024", "n10"))
 
 df_plot_fam <- df_plot_fam %>%
   mutate(# Replace "newparams" with "n10"
     n_value = str_extract(Method, "n\\d+"),                          # Extract the "n" value (e.g., "n500", "n1000")
     sample_type = str_extract(Method, "(deam_)?(iterative|\\d+h\\d)"), # Extract "deam_iterative", "deam_1h1", etc.
     value = str_extract(n_value,"\\d+" )# extract just integer value with "n"
   )
 
 row.names(df_plot_fam) <- NULL
 
 df_plot_fam$sample_type <- replace(df_plot_fam$sample_type, df_plot_fam$sample_type == "iterative", "no_damage")
 df_plot_fam[df_plot_fam$Method == "deam_n10", "sample_type"] <- "deam"
 df_plot_fam[df_plot_fam$Method == "n10", "sample_type"] <- "no_damage"
 
 
 df_plot_fam <- df_plot_fam %>% 
   mutate(sample_type = factor(sample_type, 
                               levels = c("no_damage", 
                                          "1h1", 
                                          "2h6", 
                                          "8h2",
                                          "deam")))
 df_plot_fam$value <- as.integer(df_plot_fam$value)
 
 #5. Plot a connected scatter plot 
 ggplot(df_plot_fam, aes(x = reorder(n_value, value), y = ShannonIndex, group = sample_type)) +
   geom_point(aes(color = sample_type), size = 2) +
   scale_color_brewer(name = "% Deaminated fragments", 
                      labels = c("0%", "3% - mudline", "7% - middle", "18% - bottom", "100%"),
                      palette = "Set2") + 
   geom_line(aes(color = sample_type)) + # custom colors + 
   labs(x= "Input Fragment size per species", y = "Shannon Index", title = "Shannon Diversity Index of Synthetic Data - Family level") +
   theme_minimal()+ 
   mytheme_nofacet
 
 
 ####create combined dataset for facet plotting
 df_plot$TaxLevel <- "Genus"
 df_plot_fam$TaxLevel <- "Family"
 
 df_combined <- rbind(df_plot, df_plot_fam)
 
 df_combined$sample_type <- factor(df_combined$sample_type, levels = c("no_damage", 
                                                                                  "1h1", 
                                                                                  "2h6", 
                                                                                  "8h2",
                                                                                  "deam"))
 df_combined$TaxLevel <- factor(df_combined$TaxLevel, levels = c("Genus", "Family"))

 ### plots with paiwise comparisons 
 
 ggplot(df_combined, aes(x = reorder(n_value, ShannonIndex), y = ShannonIndex, group = sample_type)) +
   geom_point(aes(color = sample_type), size = 2) +
   geom_line(aes(color = sample_type)) +
   scale_color_brewer(
     name = "% Deaminated fragments",
     labels = c("0%", "3% - mudline", "7% - middle", "18% - bottom", "100%"),
     palette = "Set2"
   ) +
   facet_wrap(~TaxLevel, scales = "free_y") +
   labs(
     x = "Input Fragment size per species",
     y = "Shannon Index",
     title = "Shannon Diversity Index of Synthetic Data"
   ) +
   theme_minimal() +
   mytheme_nofacet
 
