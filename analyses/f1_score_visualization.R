library(tidyverse)
library(ggplot2)
library(readxl)
library (RColorBrewer)
#navigate to F1 stats file and upload 
setwd("C:/Users/davisee/OneDrive - University of Tasmania/Documents/Bioinformatics_Documentation/Post-MEGAN_stats")

f1_table <- read_xlsx(path = "F1_statistics_clean_table.xlsx")



#create a custom theme
mytheme <-   theme(legend.position = "bottom", 
                   axis.text.x = element_text(size = 12, angle = 45), 
                   legend.text = element_text(size = 12), 
                   legend.title = element_text(size = 15), 
                   axis.text.y = element_text(size = 12), 
                   axis.title = element_text(size =12), 
                   strip.text.x.top = element_text(size = 12), 
                   legend.key.width = unit(0.5, "cm"),
                   legend.spacing.x = unit(0.2, 'cm'))

mytheme_nofacet <-   theme(legend.position = "bottom", 
                           axis.text.x = element_text(size = 12, angle = 45), 
                           legend.text = element_text(size = 12), 
                           legend.title = element_text(size = 14), 
                           axis.text.y = element_text(size = 12), 
                           axis.title = element_text(size =12), 
                           legend.key.width = unit(0.5, "cm"),
                           legend.spacing.x = unit(0.2, 'cm'))

############### look at only synthetic data

f1_syn <- f1_table[,1:7]
View(f1_syn)

## reshape data into a long format, to make it easier to plot
f1_long_syn <- f1_syn %>%
  pivot_longer(cols = c(no_damage, mudline_damage, middle_damage, bottom_damage, all_damage),  # Add other damage columns here if necessary
               names_to = "Damage_Type", 
               values_to = "F1_Score") 

#add levels to make it a bit easier to understand proportional damage from 0 - 100% 
f1_long_syn <- f1_long_syn %>% mutate(Damage_Type = factor(Damage_Type, 
                                                   levels = c("no_damage", "mudline_damage", 
                                                              "middle_damage", "bottom_damage", "all_damage")))  # Add other damage types as needed
# 
# #clean up level names to look more cohesion with stupid english 
# f1_long_syn$Level <- replace(f1_long_syn$Level, f1_long_syn$Level == "family", "Family")
# f1_long_syn$Level <- replace(f1_long_syn$Level, f1_long_syn$Level == "species", "Species")

View(f1_long_syn)

#create a barplot that compares f1 results for only synthetic runs
f1_long_syn %>% 
  ggplot(aes(x = as.factor(input_fragment), y = F1_Score, fill = Damage_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Level, scales = "fixed") +
  labs(x = "Input Fragment size per species", y = "F1 Score", title = "Performance assessment at Across Input Fragment Sizes of Synthetic only Runs") +
  scale_fill_brewer(name = "% Deaminated fragments", 
                    labels = c("0%", "3% - mudline", "7% - middle", "18% - bottom", "100%"),
                    palette = "Set2") +
  theme_minimal() +
  mytheme


f1_long_syn %>%
  ggplot(aes(x = as.factor(input_fragment), y = F1_Score, group = Damage_Type, colour = as.factor(Damage_Type))) + 
  geom_line() + 
  geom_point(size = 3)+ 
  facet_wrap(~Level, scales = "fixed") + 
  labs(x = "Input Fragment size per species", y = "F1 Score", title = "Performance assessment at Across Input Fragment Sizes of Synthetic only Runs") +
  theme_minimal() +
  mytheme + 
  scale_color_discrete(name = "% Deaminated Fragments", 
                       labels=c("0%","3% (0 mbsf - 0 kya)", "7% (14.55 mbsf - 14.5 kya)", "18% (61.35 mbsf - 214 kya)", "100%"), 
                       type =c("aquamarine3","darkorange1", "blueviolet", "deeppink1","darkolivegreen3")) 

#####################333for initial analyses, not going to include the differential damaged, but non spiked data
f1_small <- f1_table[,-(5:7)]
View(f1_small)

#create stacked barplot to show f1 scores
no_damage_smallf1_plot <- f1_small %>% 
ggplot(aes(x = as.factor(input_fragment), y = no_damage, fill = Level)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Input Fragment", y = "F1 Score", title = "Undamaged profiling performance assessment") +
  theme_minimal()+ 
  mytheme

no_damage_smallf1_plot
#reshape data to long format
f1_long <- f1_small %>%
  pivot_longer(cols = c(no_damage, spiked_mudline_damage, spiked_middle_damage, spiked_bottom_damage, all_damage),  # Add other damage columns here if necessary
               names_to = "Damage_Type", 
               values_to = "F1_Score") 
  
  f1_long <- f1_long %>% mutate(Damage_Type = factor(Damage_Type, 
                              levels = c("no_damage", "spiked_mudline_damage", 
                                         "spiked_middle_damage", "spiked_bottom_damage", "all_damage")))  # Add other damage types as needed

# Create the grouped barplot

 f1_long %>% 
  ggplot(aes(x = as.factor(input_fragment), y = F1_Score, fill = Damage_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Level, scales = "fixed") +
  labs(x = "Input Fragment size per species", y = "F1 Score", title = "Performance assessment at Across Input Fragment Sizes") +
 scale_fill_brewer(name = "% Deaminated fragments", 
                   labels = c("0%", "combined 3% (0 mbsf - 0 kya)", "combined 7% (14.55 - 14.5 kya)", "combined 18% (61.35 mbsf - 214 kya)", "100%"),
                   palette = "Set2") +
   theme_minimal() +
  mytheme
 

####################################
 
f1_comparespiked <- f1_table[,-(3:4)]

# Transforming from wide to long format
f1_long_compare_spiked <- f1_comparespiked %>%
  pivot_longer(
    cols = c(mudline_damage, spiked_mudline_damage, middle_damage, spiked_middle_damage, bottom_damage, spiked_bottom_damage),
    names_to = "Damage_Type",
    values_to = "F1_Score"
  ) 

# Factor the Damage_Type with levels specified as strings
f1_long_compare_spiked <- f1_long_compare_spiked %>% 
  mutate(Damage_Type = factor(Damage_Type, 
                              levels = c("mudline_damage", "spiked_mudline_damage", 
                                         "middle_damage", "spiked_middle_damage", 
                                         "bottom_damage", "spiked_bottom_damage")))
# f1_long_compare_spiked$Level <- replace(f1_long_compare_spiked$Level, f1_long_compare_spiked$Level == "family", "Family")
# f1_long_compare_spiked$Level <- replace(f1_long_compare_spiked$Level, f1_long_compare_spiked$Level == "species", "Species")


f1_plot3 <-
  f1_long_compare_spiked %>% 
  ggplot(aes(x = as.factor(input_fragment), y = F1_Score, fill = Damage_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Level, scales = "fixed") +
  labs(x = "Input Fragments per species", y = "F1 Score", title = "Performance assessment of synthetic only vs. combined IDOP + synthetic data") +
  scale_fill_brewer(name = "% Deaminated fragments", palette = "Paired", 
                    labels = c("3% (0 mbsf - 0 kya)", "combined 3% (0 mbsf - 0 kya)", 
                               "7% (14.55 mbsf - 14.5 kya)", "combined 7% (14.55 - 14.5 kya)", 
                               "18% (61.35 mbsf - 214 kya)", "combined 18% (61.35 mbsf - 214 kya)")) +
theme_minimal() +
  mytheme

# Exporting with high resolution
ggsave("3_f1_comp_plot_high_res_APRIL.png", plot = f1_plot3, width = 9, height = 4, dpi = 600, bg = "white")




###########################################
#only non_spiked_data
####################################

f1_comparesyn <- f1_table[,-(8:10)]

# Transforming from wide to long format
f1_long_compare_small_synthetic <- f1_comparesyn %>%
  pivot_longer(
    cols = c(no_damage, all_damage, mudline_damage, middle_damage, bottom_damage),
    names_to = "Damage_Type",
    values_to = "F1_Score"
  ) 

# Factor the Damage_Type with levels specified as strings
f1_long_compare_small_synthetic <- f1_long_compare_small_synthetic %>% 
  mutate(Damage_Type = factor(Damage_Type, 
                              levels = c("no_damage", 
                                         "mudline_damage", 
                                         "middle_damage", 
                                         "bottom_damage",
                                         "all_damage")))


f1_plot2 <- f1_long_compare_small_synthetic %>% 
  ggplot(aes(x = as.factor(input_fragment), y = F1_Score, fill = Damage_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Level, scales = "fixed") +
  labs(x = "Input Fragments per species", y = "F1 Score", title = "Performance assessment of synthetic taxonomic assignment") +
  scale_fill_brewer(name = "% Deaminated fragments", 
                    labels = c("0%", "3% (0 mbsf - 0 kya)", "7% (14.55 mbsf - 14.5 kya)", "18% (61.35 mbsf - 214 kya)", "100%"),
                    palette = "Set2") +
  # theme_minimal() +
  # theme(legend.position = "bottom", 
  #       axis.text.x = element_text(size = 6, angle = 45), 
  #       legend.text = element_text(size = 6),
  #       legend.title = element_text(size = 8)) 
theme_minimal() +
 mytheme

# Exporting with high resolution
ggsave("2_f1_comp_non_syn_high_res_APRIL.png", plot = f1_plot2, width = 10, height = 5, dpi = 600, bg = "white")

#################### just look at all or nothing
f1_compare_allno <- f1_table[,1:4]
# Transforming from wide to long format
f1_long_compare_allno <- f1_compare_allno %>%
  pivot_longer(
    cols = c(no_damage, all_damage),
    names_to = "Damage_Type",
    values_to = "F1_Score"
  ) 

# Factor the Damage_Type with levels specified as strings
f1_long_compare_allno <- f1_long_compare_allno %>% 
  mutate(Damage_Type = factor(Damage_Type, 
                              levels = c("no_damage", 
                                         "all_damage")))


f1_long_compare_allno %>% 
  ggplot(aes(x = as.factor(input_fragment), y = F1_Score, fill = Damage_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Level, scales = "fixed") +
  labs(x = "Input Fragment size per species", y = "F1 Score", title = "Performance assessment of no deaminated vs all deaminated") +
  theme_minimal() +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(size = 12, angle = 45)) + 
  scale_fill_discrete(name = "% Fragment Deamination", labels=c("0", "100"), type =c("cadetblue1", "cornsilk4")) 



f1_long_compare_allno %>%
  ggplot(aes(x = as.factor(input_fragment), y = F1_Score, group = Damage_Type, colour = as.factor(Damage_Type))) + 
  geom_line() + 
  geom_point(size = 3)+ 
  facet_wrap(~Level, scales = "fixed") + 
  labs(x = "Input Fragment size per species", y = "F1 Score", title = "Performance assessment of Damage Treatments") +
  theme_minimal() +
  mytheme + 
  scale_color_discrete(name = "% Deaminated Fragments", labels=c("0%", "100%"), type =c("aquamarine3", "darkolivegreen2")) 

  



##### Create a section that only looks at the spiked samples, perhaps only a line graph? 

only_spiked <- f1_table[,-(3:7)] 
# only_spiked$Level <- replace(only_spiked$Level, only_spiked$Level == "family", "Family")
# only_spiked$Level <- replace(only_spiked$Level, only_spiked$Level == "species", "Species")
# 

spiked_long <- only_spiked %>% 
  pivot_longer(cols = c(spiked_mudline_damage, spiked_middle_damage, spiked_bottom_damage), 
               names_to = "Damage_Type", 
               values_to = "F1_Score")

spiked_long %>%
  mutate(Damage_Type = factor(Damage_Type, 
                              levels = c("spiked_mudline_damage", 
                                         "spiked_middle_damage", 
                                         "spiked_bottom_damage")))

spiked_long %>% 
  ggplot(aes(x = as.factor(input_fragment), y = F1_Score, fill = Damage_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Level, scales = "fixed") + 
  labs(x = "Input Fragment size per species", y = "F1 Score", title = "Peformance assessment of Spiked Damage Treatments") +
  theme_minimal() +
  mytheme + 
  scale_fill_discrete(name = "% Deaminated Fragments", 
                       labels=c("3% (0 mbsf - 0 kya)", "7% (14.55 mbsf - 14.5 kya)", "18% (61.35 mbsf - 214 kya)"), 
                       type =c("darkorange", "blueviolet", "deeppink")) 

spiked_long %>%
  ggplot(aes(x = as.factor(input_fragment), y = F1_Score, group = Damage_Type, colour = as.factor(Damage_Type))) + 
  geom_line() + 
  geom_point(size = 3)+ 
  facet_wrap(~Level, scales = "fixed") + 
  labs(x = "Input Fragment size per species", y = "F1 Score", title = "Peformance assessment of combined IODP + synthetic Damage Treatments") +
  theme_minimal() +
  mytheme + 
  scale_color_discrete(name = "% Deaminated Fragments", 
                       labels=c("3% (0 mbsf - 0 kya)", "7% (14.55 mbsf - 14.5 kya)", "18% (61.35 mbsf - 214 kya)"), 
                       type =c("darkorange", "blueviolet", "deeppink")) 

