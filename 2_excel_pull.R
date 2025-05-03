# load packages for dataset upload, clean up, and reading excel 
library(tidyverse)
library(tidyselect)
library(tidyr)
library(broom)
library(tidyverse)
library(tidyr)
library(readxl)

#Set working directory for uploading excel data
setwd("~/Bioinformatics_Documentation/Bioinformatic_Opti_1")

##upload metadata file 
#metadata <- read_excel("~/RefSeq_Genomes/notho_metadata.xlsx")
#View(metadata)

#upload excel file, sheet, and define any NA values to be defined as new sheet 
clean <- read_excel("~/RefSeq_Genomes/notho_metadata.xlsx", sheet = 5, na = " - ") 
view(clean)
# double check structure of new dataframe
str(clean)

#isolate accession values columns into new text vectors 
accession_18S <- clean %>% 
  filter(gene_name_short == "18S") %>%
  pull(genbank_accession)

acession_COI <- clean %>%
  filter(gene_name_short == "COX1") %>%
  pull(genbank_accession)