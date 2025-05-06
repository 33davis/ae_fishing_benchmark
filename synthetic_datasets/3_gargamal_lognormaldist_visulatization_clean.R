#set working directory 
setwd('~/Bioinformatics_Documentation/Dataset_Setup/datafiles/U1538/')

#load data working packages
library(ggplot2)
library(tidyr)
library(tidyselect)
library(readr)
library(dplyr)
library(broom)


#load in U1538 collapsed filted tsv data from multiqc report

reads_table <- read_tsv("./fastqc_sequence_length_distribution_plot.tsv")


#check that tsv download was all right
View(reads_table)


############ create a visualization of all of the U1538 samples in graph (this is available from MUtliQC, no real reason to do this) #########
#########isolate names of samples 
sample_names <- colnames(reads_table) 
samples <- sample_names[2:66] 

## this pivots the table to a more readable format
reads_long <- reads_table %>%
  pivot_longer(c(samples)) 

## after the pivot, I need to get rid of all the other samples, I just want U1538 
reads_long_U1538 <- reads_long %>%
  filter(grepl("U1538", name))
View(reads_long_U1538)


# Overview of reads
all <- ggplot(reads_long_U1538, aes(x = `Sequence Length (bp)`, y = value, color = name)) + 
  geom_line( show.legend = FALSE) + 
  ggtitle(label = "U1538 Samples: Collapsed Read Length Distribution") + 
  theme_bw()

#################### Spiked sampled version #####################
### now lets do the same thing but only looking at the samples we are spiking#

spiked_punch <- reads_long %>% ### reads_long contains all reported data Exp382_Deep-collapsed (three sites)
  filter(name %in% c("IODP.23137_382_U1538C_1H_1_0_5cm.collapsed",
                     "IODP.23146_382_U1538C_2H_6_95_100cm.collapsed",
                     "IODP.23164_382_U1538C_8H_2_145_150cm.collapsed")) %>%
  mutate(mbsf = case_when( ### this adds extra data for visualization
    grepl("^IODP.23137", name)~"0", 
    grepl("^IODP.23146", name)~"14.55", 
    grepl("^IODP.23164", name)~"61.35")) %>%
  mutate(kya = case_when( ### add more useful metadata information for isolated samples 
    grepl("^IODP.23137", name)~"0", 
    grepl("^IODP.23146", name)~"14.5", 
    grepl("^IODP.23164", name)~"214")) %>%
  rename(sum_seqlength = `Sequence Length (bp)`) %>% #better naming convention than multiQC
  na.omit(value) %>% #imputation to take out values where there were no recorded sequences by bp 
   mutate(log_length = log(sum_seqlength)) %>% #this shows one meanlog and one meansd for all values
  mutate(percent_reads = value/sum(value, na.rm = TRUE)) %>% ## create a proportion of reads/total reads to act as weight
  mutate(short_name = case_when( ### shorten the sample names for visualization purposes
    grepl("^IODP.23137", name) ~ "1H_1", 
    grepl("^IODP.23146", name) ~ "2H_6", 
    grepl("^IODP.23164", name) ~ "8H_2"))%>% 
  ungroup() #good practice for data 

View(spiked_punch)

########## create a stats summary table ######

spiked_punch %>% 
  group_by(name) %>%
  slice_max(value)

spiked_punch %>%
  group_by(name) %>%
  slice_min(value)


### creates a plot with all three samples, with annotation #################
over <- spiked_punch %>% 
  ggplot(aes(x = sum_seqlength, y = value, color = short_name)) +
  geom_line() +
  geom_point() +
  ggtitle(label = "U1538 Selected Samples: Read Length Distribution") + 
  ylab(label = "Number of Fragments") +
  xlab (label = "Fragment length")+
  labs(color = "Chosen Samples") +
  annotate("text", x = 110, y = 2.0e+05, label = "61.35 mbsf - 214 kya", color = "blue", size = 4) +  #consider changing the names to match mbsf or estimated age range
 annotate("text", x = 35, y = 7.9e+05, label = "0 mbsf - 0 kya", color = "red", size = 4) +
  annotate("text", x = 78, y = 1.1e+05, label = "14.55 mbsf - 14.5 kya", color = "darkgreen", size = 4) +
 geom_vline(xintercept = 47, linetype = "dashed", color = "red") + 
 annotate("text", x = 42,y = 6.1e+05, label = "47 bp", color = "red", size = 3)+
geom_vline(xintercept = 62, linetype = "dashed", color = "blue") + 
annotate("text", x = 69,y = 6.5e+05, label = "62 bp", color = "blue", size = 3) + 
  geom_vline(xintercept = 52, linetype = "dashed", color = "darkgreen") + 
  annotate("text", x = 56,y = 5.2e+05, label = "52 bp", color = "darkgreen", size = 3)
#over + theme_classic()
over + theme_bw()


#### visualization of raw values to see if there is log normal distribution for reads length########

spiked_punch %>%
  filter(name == "IODP.23137_382_U1538C_1H_1_0_5cm.collapsed") %>% 
  ggplot(aes(x = sum_seqlength, y = value)) + 
  geom_line(linewidth = 0.6)+
  geom_vline(xintercept = 36, linetype = "dashed", color = "blue") + 
  annotate("text", x = 31, y = 1.1e+05, label = "36", color = "blue", size = 5) + 
  geom_vline(xintercept = 65, linetype = "dashed", color = "blue") + 
  annotate("text", x = 70, y = 1.1e+05, label = "65", color = "blue", size = 5) +
  labs(title = "U1538C_1H_1_0_5cm read length distribution")

spiked_punch %>%
  filter(name == "IODP.23153_382_U1538C_5H_4_145_150cm.collapsed") %>% 
  ggplot(aes(x =sum_seqlength, y = value)) + 
  geom_line(linwidth = 0.6) + 
  geom_vline(xintercept = 33, linetype = "dashed", color = "blue") + 
  annotate("text", x = 30, y = 1.5e+05, label = "33", color = "blue", size = 5) + 
  geom_vline(xintercept = 60, linetype = "dashed", color = "blue") + 
  annotate("text", x = 65, y = 1.5e+05, label = "65", color = "blue", size = 5) +
  labs(title = "U1538C_5H_2_145_150cm read length distribution")

spiked_punch %>%
  filter(name == "IODP.23166_382_U1538C_10H_5_145_150cm.collapsed") %>% 
  ggplot(aes(x =sum_seqlength, y = value)) + 
  geom_line(linewidth = 0.7)+
 # geom_vline(xintercept = 29, linetype = "dashed", color = "blue") + 
  #annotate("text", x = 24, y = 2500, label = "29", color = "blue", size = 5) + 
 # geom_vline(xintercept = 53, linetype = "dashed", color = "blue") + 
  #annotate("text", x = 58, y = 2500, label = "53", color = "blue", size = 5) +
  labs(title = "U1538C_10H_5_145_150cm read length distribution" ) 



