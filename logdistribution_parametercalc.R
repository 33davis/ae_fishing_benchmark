#load data working packages
library(ggplot2)
library(tidyr)
library(tidyselect)
library(readr)
library(dplyr)
library(ggpmisc)
library(broom)
library(DescTools)
library(fitdistrplus)
library(readr)
#set working directory 

#load in U1538 collapsed filted tsv data from multiqc report

reads_table <- read_tsv("./fastqc_sequence_length_distribution_plot.tsv")
clean_read_table <- na.omit(reads_table)


#check that tsv download was all right
View(reads_table)


############ create a visualization of all of the U1538 samples in graph (this is available from MUtliQC, no real reason to do this) #########
#########isolate names of samples 
## tidying the data taken from the MultiQC 
sample_names <- colnames(reads_table) 
samples <- sample_names[2:66] 

## this pivots the table to a more readable format
reads_long <- reads_table %>%
  pivot_longer(c(samples)) 


################################shray's orginal script for simulating a probability distribution ##############
# Set the parameters for the log-normal distribution
# location <- 4.04772894  # Location parameter (μ)
# scale <- 0.588874723     # Scale parameter (σ)


#good aproximation of one of the samples 1H by eye!!!! important to know, just playing around with paramters  
location <-  4.03772894
scale <- 0.408874723


location <- 10.713426
scale <- 2.322461 

# alternative values during testing phase ######
# #good approximation of the other 
# location <-  4.16772894
# scale <- 0.578874723 
# 
# ### meanlog for all samples, unweighted
# location <- 4.43218
# scale <- 0.5311455
# 
# #u1538  weighted mean, grouped by depth
# location <- 3.91973
# scale <- 0.03565758
# 
# #site weighted mean, sum total reads
# location <- 3.987025
# scale <- 0.0351081
# 
# 
# 
# location <- 4.624973
# scale <- 0.9740629

#Generate a range of values for the x-axis -- proof of concept for probability density function #############
x <- seq(0.01, 170, length.out = 100000) ## added min/max range from actual data? 

# Calculate the probability density function (PDF) for each x
pdf_values <- dlnorm(x, meanlog = location, sdlog = scale)


# Plot the log-normal distribution
plot(x, pdf_values, type = "l", col = "blue", lty = 1,
     main = "Log-normal Distribution",
     xlab = "Fragment Size", ylab = "Probability Density Function (PDF)")
#abline(v = c(42, 47, 160), col = c("black", "red","black"), lty = 3)
#add labels to ablines?
# text(20, 0.003, "25")
# text(42, 0.003, "47")
# text(75, 0.003, "70")

legend("topright", legend = paste("Location =", location, "\nScale =", scale),
       col = "blue", lty = 1)

grid()


############## create subset list for further analysis and visualization ##############
spiked_punch <- reads_long %>% ### reads_long contains all reported data Exp382_Deep-collapsed (three sites)
  filter(name %in% c("IODP.23137_382_U1538C_1H_1_0_5cm.collapsed",
                     "IODP.23153_382_U1538C_5H_4_145_150cm.collapsed",
                     "IODP.23169_382_U1538D_14H_1_145_150cm.collapsed")) %>%
  mutate(mbsf = case_when( ### this adds extra data for visualization
    grepl("^IODP.23137", name)~"0", 
    grepl("^IODP.23153", name)~"40.55", 
    grepl("^IODP.23169", name)~"118.35")) %>%
  mutate(kya = case_when( ### add more useful metadata information for isolated samples 
    grepl("^IODP.23137", name)~"0", 
    grepl("^IODP.23153", name)~"118", 
    grepl("^IODP.23169", name)~"392")) %>%
  rename(seqlength = `Sequence Length (bp)`) %>% #better naming convention than multiQC
  na.omit(value) %>% #imputation to take out values where there were no recorded sequences by bp 
  mutate(log_length = log(seqlength)) #this shows one meanlog and one meansd for all values
 #        ,mean_loglength = mean(log_length), #calculate "location" parameter to use in pdf
#          sd_loglength = sd(log_length)) %>% #caculate "scale" parameter for pdf
#   group_by(name) %>%
#   mutate(percent_reads = value/sum(value, na.rm = TRUE)) ## create a proportion of reads/total reads to act as weight
# ungroup() #good practice for data
# 
##################### check that percent reads is caculated by group, sum_reads == 1 for each group######################
# v_spiked <- spiked_punch %>%
#   group_by(name) %>%
#   mutate(
#     sum_reads = sum(percent_reads),
#     sum_value = sum(value),
#     mean_value = mean(value),
#     max_value = max(value),
#     min_value = min(value),
#     q25 = quantile(value, 0.25),
#     q75 = quantile(value, 0.75),
#     iqr_value = q75 - q25
#     ) %>%
#   filter(value >= q25,
#          value <= q75)

### add a weighted mean and weighted sd by group
# wspiked_punch <- spiked_punch %>%
#  group_by(name) %>%
#   summarize(
#     weighted_mean = weighted.mean(log_length, percent_reads),
#     sum_weights = sum(percent_reads),
#     n = n(),
#     sd_weighted_mean = sqrt(sum(percent_reads * (log_length - weighted_mean)^2) / ((n - 1) * sum_weights)), 
#     mean_readlength = exp(weighted_mean),
#     q10 = exp(quantile(log_length, probs = 0.10, weights = percent_reads)),
#     q90 = exp(quantile(log_length, probs = 0.90, weights = percent_reads)),
#     iqr_weighted_mean = q90 - q10
#     ) %>%
#   ungroup()
# View(wspiked_punch)

# 
# wspiked_punch <- spiked_punch %>%
#   group_by(name) %>%
#   summarize(
#     weighted_mean = weighted.mean(log_length, percent_reads),
#     sum_weights = sum(percent_reads),
#     n = n(),
#     sd_weighted_mean = sqrt(sum(percent_reads * (log_length - weighted_mean)^2) / ((n - 1) * sum_weights)), 
#     mean_readlength = exp(weighted_mean),
#     q10 = exp(quantile(log_length, probs = 0.10, weights = percent_reads)),
#     q90 = exp(quantile(log_length, probs = 0.90, weights = percent_reads)),
#     iqr_weighted_mean = q90 - q10,
#     # Calculate the interquartile range (IQR) of the y values and find corresponding x values
#     iqr_y = quantile(value, probs = c(0.25, 0.75)),
#     iqr_x_max = log_length[which(value == iqr_y[1])], 
#     iqr_x_min = log_length[which(pdf_values == iqr_y[2])]
#   ) %>%
#   ungroup()
# # running this option instead of group_by will get the weighted mean for all samples rather than by depth 
# ### sums all of the reads together 
# wsum_spiked_1 <- spiked_punch[,1:8] %>% #leaves out per sample calculated percent reads
#   ungroup() %>%
#   mutate(percent_reads = value/sum(value, na.rm = TRUE)) %>%
#   summarize(
#     weighted_mean = weighted.mean(log_length, percent_reads),
#     sum_weights = sum(percent_reads),
#     n = n(),
#     sd_weighted_mean = sqrt(sum(percent_reads * (log_length - weighted_mean)^2) / ((n - 1) * sum_weights)), 
#     mean_readlength = exp(weighted_mean),
#     q10 = exp(quantile(log_length, probs = 0.10, weights = percent_reads)),
#     q90 = exp(quantile(log_length, probs = 0.90, weights = percent_reads)),
#     iqr_weighted_mean = q90 - q10
#   ) 
# View(wsum_spiked_1)
################# summary tables with geometric mean // sd ###############
geom_spiked <- spiked_punch %>% 
  group_by(name) %>% 
  summarise(
    mean_med = median(log_length),
   arith_mean = mean(log_length),
   arith_sd = sd(log_length),
    mean_geo = Gmean(log_length), 
    sd_geo = Gsd(log_length), 
   sd_theory = arith_mean / mean_med # use the calculated value rather that DescTools::package 
         )


####Theory --- > 
# E(x) --> arithmetic mean of data
#GM[x] --> geometric mean // median of data 
#Gsd[x] --> geometric standard deviation 

# E(x) = G[x] * Gsd[x]
location <- 4.624973 # geom_spiked median log data transform
scale <- 0.9740629 #calculated geom_sd with median 

################## iqr ############
#Define the mean, standard deviation, and number of data points
mean_value <- location
sd_value <- scale
n <- 800000

# Generate random normal data with the specified mean and standard deviation
set.seed(123)  # for reproducibility
data <- rnorm(n, mean = mean_value, sd = sd_value)

# Calculate the first and third quartiles of the normal distribution
q1 <- qnorm(0.05, mean = mean_value, sd = sd_value)
q3 <- qnorm(0.95, mean = mean_value, sd = sd_value)

# Calculate the interquartile range (IQR)
iqr <- q3 - q1

# Print the IQR
print(iqr)

exp#

########### Generate a range of values for the x-axis #############
x <- seq(0.01, 170, length.out = 100000) ## added min/max range from actual data? 

# Calculate the probability density function (PDF) for each x
pdf_values <- dlnorm(x, meanlog = location, sdlog = scale)


# Plot the log-normal distribution
plot(x, pdf_values, type = "l", col = "blue", lty = 1,
     main = "Log-normal Distribution",
     xlab = "Fragment Size", ylab = "Probability Density Function (PDF)")
#abline(v = c(42, 47, 160), col = c("black", "red","black"), lty = 3)
#add labels to ablines?
# text(20, 0.003, "25")
# text(42, 0.003, "47")
# text(75, 0.003, "70")

legend("topright", legend = paste("Location =", location, "\nScale =", scale),
       col = "blue", lty = 1)

grid()



################ EX: create pdf distributions stacked in one plot ############### 
#good approximation of the other 
meanlog <- c(4.624973, 4.472671)
sdlog <- c(1.132280, 1.132280)


# Generate a range of values for the x-axis
x <- seq(0.01, 160, length.out = 1000)

# Calculate the probability density function (PDF) for each x
pdf_values <- dlnorm(x, meanlog = location, sdlog = scale)

# Plot the log-normal distribution
# Plot the log-normal distributions
plot(x, dlnorm(x, meanlog = meanlog[1], sdlog = sdlog[1]), type = "l", col = "blue", lty = 1,
     main = "Log-normal Distribution",
     xlab = "Fragment Size", ylab = "Probability Density Function (PDF)")
lines(x, dlnorm(x, meanlog = meanlog[2], sdlog = sdlog[2]), col = "red", lty = 1)
#lines(x, dlnorm(x, meanlog = meanlog[3], sdlog = sdlog[3]), col = "black", lty = 1)

# Add vertical lines at specific points
abline(v = c(25, 47, 70), col = c("black", "red","black"), lty = 3)
text(20, 0.003, "25")
text(75, 0.003, "70")

# Add a legend
legend("topright", legend = paste("Meanlog =", meanlog, "\nSdlog =", sdlog),
       col = c("blue", "red", "black"), lty = 1)

# Add a grid
grid()



