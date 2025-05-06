library(tidyverse)
library(readr)
library(stringr)

synthetic_file <- read.csv("./deam_syn_names.txt", header = FALSE)
synthetic_file <- (lapply())
str(synthetic_file)
synthetic_file <- synthetic_file[,2:3]
colnames(synthetic_file) <- c("taxon", "gene_types")


# Split taxon column into separate columns and create a new taxa column
df_split <- synthetic_file %>%
  mutate(
    Accession = str_split_fixed(taxon, " ", 4)[,1],
    Genus = str_split_fixed(taxon, " ", 4)[,2],
    Species = str_split_fixed(taxon, " ", 4)[,3],
    Gene = str_split_fixed(taxon, " ", 4)[,4]
  ) %>%
  unite(taxa, Genus, Species, sep = " ") %>%
  select(taxa, Gene, Accession) # Adjust 'sum' or 'V4' to the correct column name


# Define the interested values and mapping
interested_values <- c("18S", "COI", "COX1")
gene_mapping <- c(
  "COI" = "COI",
  "COX1" = "COI",
  "18S" = "18S"
)

# Simplify the code to create the unified_gene_factor column directly
df_split <- df_split %>%
  mutate(
    gene_name = Gene %>%
      str_extract(paste(interested_values, collapse = "|")) %>%
      recode(!!!gene_mapping) %>%
      factor(levels = unique(gene_mapping))
  )

#summarise fragment counts by taxa and gene type
begin_sum <- df_split %>%
  group_by(taxa, gene_name) %>%
  summarise(input_fragments = n(), .groups = 'drop')


# Read and prepare other datasets
#change directory 
all_damage_run <- read_tsv("./deam_synthetic-ex.txt", col_names = FALSE)
colnames(all_damage_run) <- c("taxa", "alldamage_sum")

oneh1_damage_run <- read_tsv("./IODP_1h1_combo-ex.txt", col_names = FALSE)
colnames(oneh1_damage_run) <- c("taxa", "1h1_sum")

twoh6_damage_run <- read_tsv("./IODP_2h6_combo-ex.txt", col_names = FALSE)
colnames(twoh6_damage_run) <- c("taxa", "2h6_sum")

eighth2_damage_run <- read_tsv("./IODP_8h2_combo-ex.txt", col_names = FALSE)
colnames(eighth2_damage_run) <- c("taxa", "8h2_sum")

# Left join the datasets based on taxa
full_left_syn_comp <- begin_sum %>%
  left_join(oneh1_damage_run, by = 'taxa') %>%
  left_join(twoh6_damage_run, by = 'taxa') %>%
  left_join(eighth2_damage_run, by = 'taxa') %>%
  left_join(all_damage_run, by = 'taxa') 


# Ensure full_left_syn_comp includes gene_factor from df_split
full_left_syn_comp <- full_left_syn_comp %>%
  left_join(df_split %>% select(taxa, gene_factor), by = 'taxa')

#create gpplot visualization 
library(ggplot2)
library(tidyr)


# Gather the data for plotting
plot_data <- full_left_syn_comp %>%
  pivot_longer(cols = -taxa, names_to = "Sample", values_to = "Sum") %>%
  filter(!is.na(Sum)) # Remove missing values, will cause bar graph to be kinda funky

sp <- c("Artedidraco lonnbergi", "Champsocephalus gunnari", "Dissostichus eleginoides","Dissostichus mawsoni", "Harpagifer antarcticus",  "Harpagifer bispinis", 
        "Harpagifer georgianus", "Harpagifer kerguelensis", "Trematomus loennbergii", "Trematomus scotti","Aptenodytes forsteri", "Aptenodytes patagonicus", 
        "Eudyptes chrysocome", "Eudyptes chrysolophus", "Eudyptes filholi", "Eudyptes schlegeli", "Eudyptula minor", "Megadyptes antipodes", 
        "Pygoscelis adeliae", "Pygoscelis antarcticus")
plot_data <- plot_data %>% 
  mutate(taxa = factor(taxa, levels = sp)) %>%
  arrange(taxa)

# Create the barplot
ggplot(plot_data, aes(x = taxa, y = Sum, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF", "#F564E3", "#00C19F")) + # Adjust palette for 5 colors
  theme_gray() +
  labs(x = "Taxa",
       y = "Sum",
       fill = "Sample",
       title = "Comparison of synthetic taxa assignment with MARES_PR2 database") +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  coord_flip() %>%
  annotate()
