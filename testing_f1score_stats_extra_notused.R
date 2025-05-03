# Load necessary libraries
library(ggplot2)
library(dplyr)
library(betareg)  # For beta regression
library(lme4)     # For linear mixed models
library(tidyverse)
library(readxl)
library (RColorBrewer)
library(emmeans)
library(multcompView)
library(car)

mytheme <-   theme(legend.position = "bottom", 
                   axis.text.x = element_text(size = 12, angle = 45), 
                   legend.text = element_text(size = 12), 
                   legend.title = element_text(size = 15), 
                   axis.text.y = element_text(size = 12), 
                   axis.title = element_text(size =12), 
                   strip.text.x.top = element_text(size = 12), 
                   legend.key.width = unit(0.5, "cm"),
                   legend.spacing.x = unit(0.2, 'cm'))


library(grid)  # for unit()

mytheme_test <- theme(
  legend.position = "bottom",
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 12),
  legend.key.width = unit(0.5, "cm"),
  legend.spacing.x = unit(0.2, 'cm'),
  
  axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
  axis.text.y = element_text(size = 10),
  axis.title.x = element_text(size = 12, vjust = -0.5),
  axis.title.y = element_text(size = 12, vjust = 1.5),
  
  strip.text = element_text(size = 11) #, face = "bold"),
  # strip.background = element_rect(fill = "gray90", color = NA)
  
)

mytheme_test_nofacet <- theme(
  legend.position = "bottom",
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 14),
  legend.key.width = unit(0.5, "cm"),
  legend.spacing.x = unit(0.2, 'cm'),
  
  axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
  axis.text.y = element_text(size = 12),
  axis.title.x = element_text(size = 12, vjust = -0.5),
  axis.title.y = element_text(size = 12, vjust = 1.5))
#navigate to F1 stats file and upload 
setwd("C:/Users/davisee/OneDrive - University of Tasmania/Documents/Bioinformatics_Documentation/Post-MEGAN_stats")

f1_table <- read_xlsx(path ="F1_statistics_clean_table.xlsx")


f1_long <- f1_table %>% 
  pivot_longer(
    cols = -c(Level, input_fragment), 
    names_to = "damage_treatment", 
    values_to = "f1_score"
  )
    
    
    # Main ANOVA models
    anova_model <- aov(f1_score ~ Level, data = f1_long)
    anova_frag <- aov(f1_score ~ Level * input_fragment, data = f1_long)
    anova_damage <- aov(f1_score ~ Level * damage_treatment, data = f1_long)
    
    # Summary outputs
    summary(anova_model)
    summary(anova_frag)
    summary(anova_damage)
    
    # Check assumptions
    # Homogeneity of variance
    leveneTest(f1_score ~ Level, data = f1_long)
    
    # Normality of residuals
    shapiro.test(residuals(anova_model))
    
    # # Effect sizes
    # etaSquared(anova_model, type = 2)
    # etaSquared(anova_frag, type = 2)
    # etaSquared(anova_damage, type = 2)
    # 
    # # Post-hoc test with Tukey HSD
    # tukey_result <- TukeyHSD(anova_model)
    # print(tukey_result)
    # 
    # # Extract compact letter display for Tukey groups
    # tukey_letters <- multcompView::multcompLetters4(anova_model, tukey_result)
    # 
    # # Estimated marginal means with 95% CIs
    # emm <- emmeans(anova_model, pairwise ~ Level)
    # summary(emm, infer = TRUE)
    # 
    # # Save compact summary table
    # library(knitr)
    # anova_table <- tidy(anova_model)
    # kable(anova_table, caption = "One-way ANOVA results for F1 Score by Level")
    
    library(tidyverse)
    library(rstatix)
    library(FSA)        # for dunnTest
    library(ggpubr)
    
    # Kruskal-Wallis tests
    kruskal_level <- f1_long %>% kruskal_test(f1_score ~ Level)
    kruskal_frag <- f1_long %>% kruskal_test(f1_score ~ input_fragment)
    kruskal_damage <- f1_long %>% kruskal_test(f1_score ~ damage_treatment)
    
    # Pairwise post hoc (Dunn's test)
    dunn_level <- dunnTest(f1_score ~ Level, data = f1_long, method = "bh")
    
    f1_long$input_fragment <- as.factor(f1_long$input_fragment)
    
    dunn_frag <- dunnTest(f1_score ~ input_fragment, data = f1_long, method = "bh")
    dunn_damage <- dunnTest(f1_score ~ damage_treatment, data = f1_long, method = "bh")
    
    # Create a simple result summary
    list(
      kruskal_level = kruskal_level,
      kruskal_frag = kruskal_frag,
      kruskal_damage = kruskal_damage,
      dunn_level = dunn_level$res,
      dunn_frag = dunn_frag$res,
      dunn_damage = dunn_damage$res
    )
    
    # Prepare p-value annotation data
    pvals <- dunn_level %>%
      filter(P.adj < 0.05) %>%
      mutate(
        x = gsub(" - .*", "", Comparison),
        xend = gsub(".*- ", "", Comparison),
        y = seq(0.95, 1.0, 1.05)  # Adjust as needed
      )
    
    # Ensure levels are ordered for plotting
    f1_long$Level <- factor(f1_long$Level, levels = c("Family", "Genus", "Species"))
    
    # Final styled plot
    ggplot(f1_long, aes(x = Level, y = f1_score, fill = Level)) +
      geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
      geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.9) +
      geom_segment(data = pvals, 
                   aes(x = as.numeric(factor(x, levels = c("Family", "Genus", "Species"))),
                       xend = as.numeric(factor(xend, levels = c("Family", "Genus", "Species"))),
                       y = y - 0.04, yend = y - 0.04),
                   inherit.aes = FALSE) +
      geom_text(data = pvals,
                aes(x = (as.numeric(factor(x)) + as.numeric(factor(xend))) / 2,
                    y = y - 0.03,
                    label = paste("p =", format.pval(P.adj, digits = 3))),
                inherit.aes = FALSE, size = 5) +
      scale_fill_manual(values = c("Family" = "#0d0887", 
                                   "Genus" = "#cc4778", 
                                   "Species" = "#f0f921")) +
      labs(
        title = "F1 Score Differences by Taxonomic Level",
        x = "Taxonomic Level",
        y = "F1 Score"
      ) +
      theme_minimal() +
     mytheme_test_nofacet 
 
    
    # Combine all pairwise results into a single data frame with a group column
    library(dplyr)
    library(readr)
    
    # Add group info and combine
    significant_dunn <- bind_rows(
      dunn_level %>% filter(P.adj < 0.05) %>% mutate(group = "Taxonomic Level"),
      dunn_frag %>% filter(P.adj < 0.05) %>% mutate(group = "Fragment Length"),
      dunn_damage %>% filter(P.adj < 0.05) %>% mutate(group = "Damage Category")
    )
    
    # Export to CSV
    write_csv(significant_dunn, "significant_pairwise_results.csv")
    