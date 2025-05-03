setwd("./")
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

counts <- read_tsv('./Comparison-2_eukfamily.txt.txt')
View(counts)
counts <- counts %>% rename(classification = `#Datasets`)


#########isolate names of samples 
sample_names <- colnames(counts) 

counts <- counts %>% colnames(c(sample_names))
remove_these <- c("IODP\\.","text")
sample_names <- str_remove_all(sample_names, paste(remove_these, collapse = "|"))

newcounts <- counts %>% setNames(sample_names)
samples <- sample_names[2:34] 

## this pivots the table to a more readable format
reads_long <- newcounts %>%
  pivot_longer(c(samples)) 

  reads_long %>% 
    ggplot(aes(name, value, fill = classification)) + 
    geom_bar(position= "stack", stat = "identity")
  