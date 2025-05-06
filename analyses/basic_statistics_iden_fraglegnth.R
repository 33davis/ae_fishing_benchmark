#### Install necessary packages if not already installed ####
if (!require(ggplot2)) install.packages("ggplot2", dependencies = TRUE)
if (!require(nortest)) install.packages("nortest", dependencies = TRUE)
if (!require(moments)) install.packages("moments", dependencies = TRUE)
if (!require(FSA)) install.packages("FSA", dependencies = TRUE)
if (!require(rstatix)) install.packages("rstatix", dependencies = TRUE)
if (!require(ggsignif)) install.packages("ggsignif", dependencies = TRUE)
# Load packages
library(ggplot2)
library(nortest)
library(moments)
library(readxl)
library(FSA)
library(dplyr)
library(rstatix)
library(ggsignif)
#### set working directoy and load in desired table #### 
setwd("C:/Users/davisee/OneDrive - University of Tasmania/Documents/Bioinformatics_Documentation/Post-MEGAN_stats")
spec_iden <- read_excel("C:/Users/davisee/OneDrive - University of Tasmania/Documents/Bioinformatics_Documentation/Post-MEGAN_stats/Draft_stats_table_updated.xlsx", sheet = "Sheet2")
#### clean up datasheet for testing #### 
#check structure of download
str(spec_iden)
#check gene names that you want to test by
unique(spec_iden$Gene_mapped)

#### checking for normality visually ####
# Histogram
ggplot(spec_iden, aes(x = Indentity)) +
  geom_histogram(binwidth = 1, fill = "blue", alpha = 0.5, color = "black") +
  labs(title = "Histogram of Percent Identity Scores", x = "Percent Identity", y = "Frequency") +
  theme_minimal()

# Q-Q Plot
ggplot(spec_iden, aes(sample = Indentity)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(title = "Q-Q Plot of Percent Identity Scores") +
  theme_minimal()

# Density Plot
ggplot(spec_iden, aes(x = Indentity)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density Plot of Percent Identity Scores", x = "Percent Identity", y = "Density") +
  theme_minimal()

#### statisical tests for normality ####
# Shapiro-Wilk Test (best for small to moderate sample sizes n < 5000)
shapiro.test(spec_iden$Indentity)
#if p-value > 0.05 -> fail to reject null hypothesis (data is likely normal)
#if p-value < 0.05 -> reject null hypothesis (data is not normal)

# Skewness & Kurtosis Test
skewness(spec_iden$Indentity)  # Close to 0 if normal
kurtosis(spec_iden$Indentity)  # Close to 3 if normal

#### comparing percent identitiy scores across Gene_mappeds #### 
# One-Way ANOVA (if normality assumption holds)
# anova_result <- aov(Indentity ~ Gene_mapped, data = spec_iden)
# summary(anova_result)
# 
# # Post-hoc test: Tukey’s HSD for pairwise comparisons
# TukeyHSD(anova_result)

# Kruskal-Wallis Test (if data is non-normal)
KW_test <- kruskal.test(Indentity ~ Gene_mapped, data = spec_iden)
print(KW_test)

dunn_test_bonf <- dunnTest(Indentity ~ Gene_mapped, data = spec_iden, method = "bonferroni")# Adjusted p-values
dunn_test_holm <- dunnTest(Indentity ~ Gene_mapped, data = spec_iden, method = "holm")
dunn_test_BH <- dunnTest(Indentity ~ Gene_mapped, data = spec_iden, method = "bh")
print(dunn_test_bonf)
print(dunn_test_holm)
print(dunn_test_BH)


## visualise this
# Create a boxplot
ggplot(spec_iden, aes(x = Gene_mapped, y = Indentity, fill = Gene_mapped)) +
  geom_boxplot() +
  #stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white") +  # Mean point
  theme_minimal() +
  labs(title = "Percent Identity Score by Gene",
       y = "Percent Identity",
       x = "Gene") +
  scale_fill_manual(values = c("16S" = "lightblue", "18S" = "red", "COI" = "green")) +
  geom_signif(
    comparisons = list(c("16S", "18S"), c("16S", "COI"), c("18S", "COI")),
    map_signif_level = TRUE,
    step_increase = 0.1
  )

#### compare each gene's mean to the global mean #### 
# Compute global mean
global_mean <- mean(spec_iden$Indentity)

# # One-sample t-test for each Gene_mapped (if normality assumption holds)
# t.test(spec_iden$Indentity[spec_iden$Gene_mapped == "16S"], mu = global_mean)
# t.test(spec_iden$Indentity[spec_iden$Gene_mapped == "18S"], mu = global_mean)
# t.test(spec_iden$Indentity[spec_iden$Gene_mapped == "COI"], mu = global_mean)

# Wilcoxon Signed-Rank Test for each Gene_mapped (if data is non-normal)
wilcox.test(spec_iden$Indentity[spec_iden$Gene_mapped == "16S"], mu = global_mean)
wilcox.test(spec_iden$Indentity[spec_iden$Gene_mapped == "18S"], mu = global_mean)
wilcox.test(spec_iden$Indentity[spec_iden$Gene_mapped == "COI"], mu = global_mean)



######### testing the same but with fragment length ####### 
ggplot(spec_iden, aes(x = length)) +
  geom_histogram(binwidth = 1, fill = "blue", alpha = 0.5, color = "black") +
  labs(title = "Histogram of Percent Identity Scores", x = "Percent Identity", y = "Frequency") +
  theme_minimal()

# Q-Q Plot
ggplot(spec_iden, aes(sample = length)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(title = "Q-Q Plot of Percent Identity Scores") +
  theme_minimal()

# Density Plot
ggplot(spec_iden, aes(x = length)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density Plot of Percent Identity Scores", x = "Percent Identity", y = "Density") +
  theme_minimal()

#### statisical tests for normality ####
# Shapiro-Wilk Test (best for small to moderate sample sizes n < 5000)
shapiro.test(spec_iden$length)
#if p-value > 0.05 -> fail to reject null hypothesis (data is likely normal)
#if p-value < 0.05 -> reject null hypothesis (data is not normal)

# Skewness & Kurtosis Test
skewness(spec_iden$length)  # Close to 0 if normal
kurtosis(spec_iden$length)  # Close to 3 if normal

#### comparing percent identitiy scores across Gene_mappeds #### 
# One-Way ANOVA (if normality assumption holds)
# anova_result <- aov(length ~ Gene_mapped, data = spec_iden)
# summary(anova_result)
# 
# # Post-hoc test: Tukey’s HSD for pairwise comparisons
# TukeyHSD(anova_result)

# Kruskal-Wallis Test (if data is non-normal)
KW_test_length <- kruskal.test(length ~ Gene_mapped, data = spec_iden)
print(KW_test_length)

# dunn_test_bonf_L <- dunnTest(length ~ Gene_mapped, data = spec_iden, method = "bonferroni")# Adjusted p-values
# dunn_test_holm_L <- dunnTest(length ~ Gene_mapped, data = spec_iden, method = "holm")
# dunn_test_BH_L <- dunnTest(length ~ Gene_mapped, data = spec_iden, method = "bh")
# print(dunn_test_bonf_L)
# print(dunn_test_holm_L)
# print(dunn_test_BH_L)


## visualise this
ggplot(spec_iden, aes(x = Gene_mapped, y = length, fill = Gene_mapped)) +
  geom_boxplot() +
  #stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white") +  # Mean point
  theme_minimal() +
  labs(title = "Fragment Length by assigned gene",
       y = "Fragment Length",
       x = "Gene") +
  scale_fill_manual(values = c("16S" = "lightblue", "18S" = "red", "COI" = "green")) +
  geom_signif(
    comparisons = list(c("16S", "18S"), c("16S", "COI"), c("18S", "COI")),
    map_signif_level = TRUE,
    step_increase = 0.1
  )

#### compare each gene's mean to the global mean #### 
# Compute global mean
global_mean <- mean(spec_iden$length)

# # One-sample t-test for each Gene_mapped (if normality assumption holds)
# t.test(spec_iden$length[spec_iden$Gene_mapped == "16S"], mu = global_mean)
# t.test(spec_iden$length[spec_iden$Gene_mapped == "18S"], mu = global_mean)
# t.test(spec_iden$length[spec_iden$Gene_mapped == "COI"], mu = global_mean)

# Wilcoxon Signed-Rank Test for each Gene_mapped (if data is non-normal)
wilcox.test(spec_iden$length[spec_iden$Gene_mapped == "16S"], mu = global_mean)
wilcox.test(spec_iden$length[spec_iden$Gene_mapped == "18S"], mu = global_mean)
wilcox.test(spec_iden$length[spec_iden$Gene_mapped == "COI"], mu = global_mean)
