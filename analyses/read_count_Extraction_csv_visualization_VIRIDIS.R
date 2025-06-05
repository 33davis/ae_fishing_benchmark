### this is for extracted read count data exported from MEGAN CE per sample file, collating
#together into merged csv files
#data cleanup, visualization 
library(dplyr)
library(data.table)
library(ggplot2)
library(grid) # for unit()
library(stringr)
library(tidyr)
############# create a custom theme for chapter viz #############
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
  legend.key.width = unit(0.7, "cm"),
  legend.spacing.x = unit(0.2, 'cm'),
  
  axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
  axis.text.y = element_text(size = 10),
  axis.title.x = element_text(size = 12, vjust = -0.5),
  axis.title.y = element_text(size = 12, vjust = 1.5),
  
  strip.text = element_text(size = 11) #, face = "bold"),
  # strip.background = element_rect(fill = "gray90", color = NA)
)

mytheme_bigheatmap_facetgrid <- theme_minimal() +
  theme(
    # Legend
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.key.width = unit(0.7, "cm"),
    legend.spacing.x = unit(0.2, 'cm'),
    
    # Axis text and titles
    axis.text.x = element_text(size = 7, angle = 90, hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 12, vjust = -0.5),
    axis.title.y = element_text(size = 12, vjust = 1.5),
    
    # Facet strip
    strip.text = element_text(size = 10),
    strip.background = element_blank(),
    
    # Panel spacing (between facets)
    panel.spacing = unit(1, "lines"),
    # 
    # # Plot margins (optional: tighter layout)
    # plot.margin = margin(5, 5, 5, 5), 
      # add a background so it isn't transparent
      panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

mytheme_nofacet <-   theme(legend.position = "bottom", 
                           axis.text.x = element_text(size = 12, angle = 45), 
                           legend.text = element_text(size = 12), 
                           legend.title = element_text(size = 14), 
                           axis.text.y = element_text(size = 12), 
                           axis.title = element_text(size =12), 
                           legend.key.width = unit(2.5, "cm"),
                           legend.spacing.x = unit(0.3, 'cm'))

mytheme_test_nofacet <- theme(
  legend.position = "bottom",
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 12),
  legend.key.width = unit(1.5, "cm"),
  legend.spacing.x = unit(0.2, 'cm'),
  
  axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
  axis.text.y = element_text(size = 12),
  axis.title.x = element_text(size = 12, vjust = -0.5),
  axis.title.y = element_text(size = 12, vjust = 1.5))


######### process extracted fasta files from MEGAN CE reads-extractedtaxa.fasta ########
process_fasta_files <- function(file_directory, output_directory = NULL) {  
  # List all FASTA files in the specified folder
  fasta_files <- list.files(path = file_directory, pattern = "\\.fasta$", full.names = TRUE)
  
  # Initialize a list to store data for each file
  fasta_data <- list()
  
  # Loop through each FASTA file and extract data
  for (file in fasta_files) {
    # Read the FASTA file
    fasta_content <- readLines(file)
    
    # Check if the file contains any FASTA headers (lines starting with ">")
    if (length(fasta_content) == 0 || !any(grepl("^>", fasta_content))) {
      message(paste("Skipping file:", file, "- No FASTA content found"))
      next
    }
    
    # Extract read names (lines starting with ">"), removing ">"
    read_names <- gsub("^>", "", fasta_content[grepl("^>", fasta_content)])
    
    # Create a data.table for the current file
    temp_dt <- data.table(
      read_name = read_names,
      file_name = basename(file)  # Extract filename only
    )
    
    # Add the data.table to the list
    fasta_data[[file]] <- temp_dt
  }
  
  # Combine all data.tables into one
  if (length(fasta_data) > 0) {
    combined_data <- rbindlist(fasta_data)
    
    # Remove everything after the first `.` in `read_name`, if applicable
    combined_data[, read_name := gsub("\\..*", "", read_name)]
    
    # Extract assignment_read_name from file_name (update regex as needed)
    combined_data[, assignment_read_name := sub(".*reads-([^.]+)\\.fasta.*", "\\1", file_name)]
    
    # Create a summary table with counts
    summary_table <- combined_data %>%
      group_by(assignment_read_name, read_name) %>%
      summarize(read_count = n(), .groups = "drop") %>%
      mutate(file_name = file_directory)
    
    # Handle output directory
    if (is.null(output_directory)) {
      output_directory <- file_directory
    }
    if (!dir.exists(output_directory)) {
      dir.create(output_directory, recursive = TRUE)
    }
    
    # Define output file name
    output_file <- file.path(output_directory, paste0(basename(file_directory), "_summary_table.csv"))
    
    # Write the summary table to a CSV file
    write.csv(summary_table, output_file, row.names = FALSE)
    
    message(paste("Summary table saved to:", output_file))
    return(summary_table)
  } else {
    message("No valid FASTA files found in the directory.")
    return(NULL)
  }
}

# Example usage:
# summary_table <- process_fasta_files(
#   file_directory = "C:/path/to/input/files",
#   output_directory = "C:/path/to/output/directory"
# )


### currently is not set up to deal with multiple subdirectories, such as the deamination treatment folders within each subdirectory
### double check that working directory is set, otherwise you WILL get an error
main_directory <-"C:/Users/davisee/OneDrive - University of Tasmania/Documents/Bioinformatics_Documentation/upscale_runs/newparams/all_spiked_extracted/n5000"
setwd("C:/Users/davisee/OneDrive - University of Tasmania/Documents/Bioinformatics_Documentation/upscale_runs/newparams/all_spiked_extracted/n5000")
folder_names <- list.dirs(path = main_directory, full.names = FALSE, recursive = FALSE)
folder_names

for (folder in folder_names) {
  ## message which folder is being processed: 
  message(paste("This folder is being processed:", folder))
  
process_fasta_files(
  file_directory = folder,
  output_directory = "n5000")
}
message("completed processing folders")


############ putting all of the new csv files together to view the data #########
#establish parent csv directory, output directory created above
# Load required package
library(dplyr) # for bind_rows (optional but robust for varied structures)

# Define the function
merge_csv_files <- function(csv_directory, pattern = "\\.csv$", output_file = NULL) {
  # Step 1: List all CSV files in the directory
  csv_files <- list.files(path = csv_directory, pattern = pattern, full.names = TRUE)
  
  # Check if any files are found
  if (length(csv_files) == 0) {
    stop("No CSV files found in the directory with the specified pattern.")
  }
  
  # Step 2: Read each file and store in a list
  csv_data <- lapply(csv_files, read.csv)
  
  # Step 3: Combine all data frames into one
  merged_data <- bind_rows(csv_data) # dplyr's bind_rows handles inconsistent column orders
  
  # Step 4: Save to output file if specified
  if (!is.null(output_file)) {
    write.csv(merged_data, file = output_file, row.names = FALSE)
    message(paste("Merged data saved to:", output_file))
  }
  
  # Step 5: Return the merged data
  return(merged_data)
}

csv_dir <- "C:/Users/davisee/OneDrive - University of Tasmania/Documents/Bioinformatics_Documentation/upscale_runs/newparams/all_spiked_extracted/n1000/n1000"
output_file <- "C:/Users/davisee/OneDrive - University of Tasmania/Documents/Bioinformatics_Documentation/upscale_runs/newparams/all_spiked_extracted/n1000/summary_n1000.csv"

merged_n1000_data <- 
  merge_csv_files(
  csv_directory = csv_dir, 
  pattern = "\\.csv$", 
  output_file = output_file
)
##################### function to merge all created  summary csv files from parent directory ##############################

library(dplyr)
library(data.table)

merge_csv_files <- function(parent_directory, output_file = "merged_all.csv") {
  # Define allowed folder names
  valid_folders <- c("n10", "n25", "n50", "n100", "n250", "n500", "n1000", "n2500", "n5000")
  
  # Get list of subdirectories
  subdirs <- list.dirs(parent_directory, recursive = FALSE)
  subdirs <- subdirs[basename(subdirs) %in% valid_folders]  # Keep only valid folders
  
  # Initialize a list to store data
  csv_data <- list()
  
  # Loop through each valid folder and check for merged_n$.csv
  for (subdir in subdirs) {
    # Construct the expected file name based on folder name
    folder_name <- basename(subdir)
    csv_file <- file.path(subdir, paste0("summary_", folder_name, ".csv"))
    
    # Read and store the file if it exists
    if (file.exists(csv_file)) {
      temp_data <- fread(csv_file)  # Read CSV efficiently
      temp_data[, folder_name := folder_name]  # Add folder info
      csv_data[[folder_name]] <- temp_data
    } else {
      message(paste("Skipping:", csv_file, "- File not found"))
    }
  }
  
  # Combine all data.tables if any files were found
  if (length(csv_data) > 0) {
    merged_data <- rbindlist(csv_data, fill = TRUE)
    
    # Save the merged file
    fwrite(merged_data, file.path(parent_directory, output_file))
    
    message("Merged CSV saved to:", file.path(parent_directory, output_file))
    return(merged_data)
  } else {
    message("No matching CSV files found in the specified folders.")
    return(NULL)
  }
}

# Example usage
merge_csv_files("C:/Users/davisee/OneDrive - University of Tasmania/Documents/Bioinformatics_Documentation/upscale_runs/newparams/all_spiked_extracted")  # Change this to your actual directory

######################## load in new summary csv files to do visualizations #########################################
setwd("C:/Users/davisee/OneDrive - University of Tasmania/Documents/Bioinformatics_Documentation/upscale_runs/newparams/single_runs/reads_extracted/summary_csv")
merged_data <- read.csv("summary_allsynthetic.csv") # used for the synthetic version 

setwd("C:/Users/davisee/OneDrive - University of Tasmania/Documents/Bioinformatics_Documentation/upscale_runs/newparams/all_spiked_extracted")
merged_combo_data <- read.csv("merged_all.csv")
merged_combo_data <- merged_combo_data %>% 
  rename("assignment_level" = "assignment_read_name") ## will make consistent with previous code

plots_dir <- "C:/Users/davisee/OneDrive - University of Tasmania/Documents/chapter1/figures"

################ Do some preliminary synthetic data factoring ###################
merged_data <- merged_data %>% 
  mutate(file_name = factor(file_name, 
      levels = c("n10_nodeam", "n10_lowdeam", 
                "n10_middledeam", "n10_highdeam", "n10_deam", 
                "n25_nodeam", "n25_lowdeam", 
                "n25_middledeam", "n25_highdeam", "n25_deam", 
                "n50_nodeam", "n50_lowdeam",
                "n50_middledeam", "n50_highdeam", "n50_deam", 
                "n100_nodeam", "n100_lowdeam", 
                "n100_middledeam", "n100_highdeam", "n100_deam",
                "n250_nodeam", "n250_lowdeam", 
                "n250_middledeam", "n250_highdeam", "n250_deam", 
                "n500_nodeam", "n500_lowdeam", 
                "n500_middledeam", "n500_highdeam", "n500_deam", 
                "n1000_nodeam", "n1000_lowdeam", 
                "n1000_middledeam", "n1000_highdeam", "n1000_deam", 
                "n2500_nodeam", "n2500_lowdeam", 
                "n2500_middledeam", "n2500_highdeam", "n2500_deam", 
                "n5000_nodeam", "n5000_lowdeam", 
                "n5000_middledeam", "n5000_highdeam", "n5000_deam"))) 
merged_data$read_name <- factor(merged_data$read_name, unique(merged_data$read_name))

# all_assignment <- merged_data %>% 
#   group_by(assignment_level, read_name, file_name) %>%
#   mutate(read_count = n())
# 
# all_assignment <- merged_data %>%
#   group_by(assignment_level, read_name, file_name) %>%
#   summarise(read_count = sum(read_count), .groups = "drop")  

#Group and count properly
all_assignment <- merged_data %>%
  group_by(assignment_level, read_name, file_name) %>%
  summarise(read_count = n(), .groups = "drop")  # Count number of reads per group# Sum read counts instead of recounting
View(all_assignment)

S_reads <- all_assignment %>%
  filter(read_name %in% c("AF518190","KP875235"))

coi_reads <- all_assignment %>%
  filter(!read_name %in% c("AF518190","KP875235"))


### combined IODP + synthetic data version preliminary factoring #####
merged_combo_data <- merged_combo_data %>% 
  mutate(file_name = factor(file_name, 
                            levels = c("n10_lowdeam", "n10_middledeam", "n10_highdeam",  
                                      "n25_lowdeam", "n25_middledeam", "n25_highdeam",
                                      "n50_lowdeam","n50_middledeam", "n50_highdeam",  
                                      "n100_lowdeam", "n100_middledeam", "n100_highdeam",
                                      "n250_lowdeam", "n250_middledeam", "n250_highdeam",
                                      "n500_lowdeam", "n500_middledeam", "n500_highdeam",
                                      "n1000_lowdeam", "n1000_middledeam", "n1000_highdeam",
                                      "n2500_lowdeam", "n2500_middledeam", "n2500_highdeam",
                                      "n5000_lowdeam", "n5000_middledeam", "n5000_highdeam")))
merged_combo_data$read_name <- factor(merged_combo_data$read_name, unique((merged_combo_data$read_name)))

############ exploratory analysis to understand outputs 
S_reads <- merged_combo_data %>%
  filter(read_name %in% c("AF518190","KP875235"))


coi_reads <- merged_combo_data %>%
  filter(!read_name %in% c("AF518190","KP875235"))


########################## rework for combined IODP + synthetic analysis ##################
all_assignment <- merged_combo_data %>% 
  group_by(assignment_level, read_name, file_name) %>%
  mutate(read_count = n())

all_assignment <- merged_combo_data %>%
  group_by(assignment_level, read_name, file_name) %>%
  summarise(read_count = sum(read_count), .groups = "drop")  # Sum read counts instead of recounting

View(all_assignment)

##################### harpagifer antarcticus 18s v. COI #######################
ha_reads <- all_assignment %>%
  filter(read_name %in% c("AF518190", "ON891147")) %>%
  mutate(read_name = case_when(
    read_name == "AF518190" ~ "18S - AF518190",
    read_name == "ON891147" ~ "COI - ON891147",
    TRUE ~ read_name # In case of any other values, retain them unchanged
  ))

## review plot of 
ha_reads %>%
  #filter(read_name == "18S - AF518190") %>%
  #filter(read_name == "COI - ON891147") %>% 
ggplot(aes(assignment_level, file_name, fill = read_count)) + 
  geom_tile()+
  scale_fill_viridis_c(option = "plasma", direction = -1)+ #more colourblind accessible
  #facet_wrap(~read_name, scales = "free")+ 
  scale_x_discrete(expand = c(0, 0)) +  # remove x padding around tiles
  facet_grid(~read_name, scales = "free", space = "free") + # enables varying widths
  mytheme_bigheatmap_facetgrid+ 
  labs(
    x = "MALT Taxonomic Assignment", 
    y = "Damage and Fragment Input Run", 
    #title = "Harpagifer antarcticus - 18S vs. COI taxonomic assignment sytnthic datasets", #title for synthetic_only data
    title = "Harpagifer antarcticus - 18S vs. COI taxonomic assignment combined IODP + synthetic datasets", #title for merged_combo data
    fill = "Total count\nfragment assignment"
   # title = "Harpagifer antarcticus - 18S taxonomic assignment"
    #title = "Harpagifer antarcticus - COI taxonomic assignment"
  ) 

  ggsave(filename = paste0(plots_dir, "HA_18SCOI_synthetic_fullheatmap.png"), width = 12, height = 7, units = "in", dpi = 300) #for synthetic plot
  ggsave(filename = paste0(plots_dir, "HA_18SCOI_combinediodpsynthetic_fullheatmap.png"), width = 12, height = 7, units = "in", dpi = 300) #for combined plot
#filter data for paper to have ONLY instances were it showed up for all file_names (still faceted)
# Filter to include only assignment levels present in all file_names per read_name
filtered_ha_reads <- ha_reads %>%
  group_by(read_name) %>%
  mutate(total_files = n_distinct(file_name)) %>%
  group_by(read_name, assignment_level) %>%
  mutate(n_files = n_distinct(file_name)) %>%
  ungroup() %>%
  filter(n_files == total_files)

# Plot
filtered_ha_reads %>%
  ggplot(aes(assignment_level, file_name, fill = read_count)) + 
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", direction = -1)+ #more colourblind accessible
  # scale_fill_gradient(low = "purple", high = "yellow")+
  #scale_fill_gradient(low = "#440154", high = "#FDE725")+
  facet_wrap(~read_name, scales = "free") + 
  theme_minimal() +
   mytheme_test + 
  theme(panel.background = element_rect(fill = "white", color = NA),
         plot.background = element_rect(fill = "white", color = NA)) + 
  labs(
    x = "MALT Taxonomic Assignment", 
    y = "Damage and Fragment Input Run", 
    #title = "Harpagifer antarcticus - 18S vs. COI taxonomic assignment synthetic data",
    fill = "Total count\nfragment assignment"
  )
ggsave( "HA_18SCOI_synthetic_filtered_heatmap.png", width = 5.5, height = 9, units = "in", dpi = 300)
ggsave( "HA_18SCOI_combinedIODPsynthetic_filtered_heatmap.png", width = 5.5, height = 9, units = "in", dpi = 300)
#################################### pygosclis antarcticus 18S v. COI ##################################
pa_reads <- all_assignment %>%
  filter(read_name %in% c("KP875235", "EU525471")) %>%
  mutate(read_name = case_when(
    read_name == "KP875235" ~ "18S - KP875235",
    read_name == "EU525471"~ "COI - EU525471",
    TRUE ~ read_name # In case of any other values, retain them unchanged
  ))

pa_reads %>%
  #filter(read_name == "18S - KP875235") %>%
  #filter(read_name == "COI - EU525471") %>% 
ggplot(aes(assignment_level, file_name, fill = read_count)) + 
  geom_tile()+
  scale_fill_viridis_c(option = "plasma", direction = -1)+
  scale_x_discrete(expand = c(0, 0)) +  # remove x padding around tiles
  facet_grid(~read_name, scales = "free", space = "free") + # enables varying widths
  mytheme_bigheatmap_facetgrid+ 
  labs(
    x = "MALT Taxonomic Assignment", 
    y = "Damage and Fragment Input Run", 
    #title = "Harpagifer antarcticus - 18S vs. COI taxonomic assignment sytnthic datasets", #title for synthetic_only data
    title = "Pygoscelis antarcticus - 18S vs. COI taxonomic assignment combined IODP + synthetic datasets", #title for merged_combo data
    fill = "Total count\nfragment assignment"
  )
ggsave("PA_18SCOI_synthetic_fullheatmap.png", width = 12, height = 7, units = "in", dpi = 300)
ggsave("PA_18SCOI_combinedIODPsynthetic_fullheatmap.png", width = 12, height = 7, units = "in", dpi = 300)
##### same filtered thing as nefore
# Filter to include only assignment levels present in all file_names per read_name
filtered_pa_reads <- pa_reads %>%
  group_by(read_name) %>%
  mutate(total_files = n_distinct(file_name)) %>%
  group_by(read_name, assignment_level) %>%
  mutate(n_files = n_distinct(file_name)) %>%
  ungroup() %>%
  filter(n_files == total_files)

# Plot
filtered_pa_reads %>%
  ggplot(aes(assignment_level, file_name, fill = read_count)) + 
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", direction = -1)+ #more colourblind accessible
  # scale_fill_gradient(low = "purple", high = "yellow")+
  #scale_fill_gradient(low = "#440154", high = "#FDE725")+
  facet_wrap(~read_name, scales = "free") + 
  theme_minimal() +
  mytheme_test +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)) + 
  labs(
    x = "MALT Taxonomic Assignment", 
    y = "Damage and Fragment Input Run", 
    #title = "Pygoscelis antarcticus - 18S vs. COI taxonomic assignment",
    fill = "Total count\nfragment assignment"
  )
ggsave( "PA_18SCOI_synthetic_filtered_heatmap.png", width = 5.5, height = 9, units = "in", dpi = 300)
ggsave( "PA_18SCOI_combinedIODPsynthetic_filtered_heatmap.png", width = 5.5, height = 9, units = "in", dpi = 300)
    
############################################# testing species_only assignment across all runs#############

library(reshape2)
#################### look at species only assignments ############
# 
# species_assigment <- merged_data %>% 
#   filter(str_detect(assignment_level, "_")) %>%
#   filter(!assignment_level %in% c("cellular_organisms", "No_hits", "Not_assigned")) %>%
#   group_by(assignment_level, read_name, file_name) %>%
#   mutate(read_count = n())
# View(species_assigment)
# 
# ### for spiked datasets 
# species_assigment <- merged_combo_data %>% 
#   filter(str_detect(assignment_level, "_")) %>%
#   filter(!assignment_level %in% c("cellular_organisms", "No_hits", "Not_assigned")) %>%
#   group_by(assignment_level, read_name, file_name) %>%
#   mutate(read_count = n())
# View(species_assigment)
# 
# 
# # Filter to just the observed data first (non-zero read_count)
# nonzero_species <- species_assigment %>%
#   filter(read_count > 0)
# 
# # Get total file count per read_name
# total_files_df <- nonzero_species %>%
#   group_by(read_name) %>%
#   summarise(total_files = n_distinct(file_name), .groups = "drop")
# 
# # Count how many files each assignment_level appears in (with non-zero reads)
# assignments_per_read <- nonzero_species %>%
#   group_by(read_name, assignment_level) %>%
#   summarise(n_files = n_distinct(file_name), .groups = "drop")
# 
# # Filter to species present (with counts) in all file_names
# assignments_present_in_all <- assignments_per_read %>%
#   left_join(total_files_df, by = "read_name") %>%
#   filter(n_files == total_files)
# 
# # Filter the full dataset (including zeroes if needed for tile structure)
# filtered_species_reads <- species_assigment %>%
#   semi_join(assignments_present_in_all, by = c("read_name", "assignment_level"))
# 
# 
# # Plot
# filtered_species_reads %>%
#            filter(read_count > 0, 
#                   !assignment_level %in% c("Harpagifer_bispinis", "Hesseltinella_vesiculosa", "Trematomus_loennbergii")) %>% #species with zero counts breaking the plot
#   ggplot(aes(assignment_level, file_name, fill = read_count)) + 
#   geom_tile(color = "white") +
#   scale_fill_viridis_c(option = "plasma", direction = -1) +
#   # coord_fixed(ratio = 1) +
#   theme_minimal() +
#   mytheme_test_nofacet + 
#   labs(
#     x = "MALT Taxonomic Assignment", 
#     y = "Damage and Fragment Input Run", 
#     title = "Species-level taxonomic assignment across synthetic datasets"
#   )
# 
# 
# 
# plot_data <- filtered_species_reads %>%
#   select(read_name, assignment_level, file_name, read_count) %>%
#   filter(read_name == unique(read_name)[1])
# 
# ggplot(plot_data, aes(assignment_level, file_name, fill = read_count)) +
#   geom_tile(color = "white") +
#   scale_fill_viridis_c(option = "plasma", direction = -1) +
#   theme_minimal() +
#   mytheme_test_nofacet +
#   labs(
#     x = "MALT Taxonomic Assignment", 
#     y = "Damage and Fragment Input Run", 
#     title = "Species-level taxonomic assignment (preview)"
#   )

################## filtering and renaming accession numbers ############### 
accession <- c("EU525299", "EU525303","OK493624",
"EF609319","MW829388","ON000293",
"FJ582589","FJ582593","MG739962",
"MG740185","EU525360","ON891147",
"AF518190","OL339430","ON891170",
"ON891164","MK262211","MK261839",
"EU525471","KP875235","OK493720",
"OK493634"
)

all_assignment_syn  <- all_assignment %>%
  filter(read_name %in% accession) %>%
  mutate(read_name = case_when(
    read_name == "AF518190" ~ "18S - H. antarcticus",
    read_name == "ON891147" ~ "COI - H. antarcticus",
    read_name == "EU525299" ~ "Aptenodytes forsteri", 
    read_name == "EU525303" ~ "Aptenodytes patagonicus",
    read_name =="OK493624" ~ "Artedidraco lonnbergi",
    read_name =="EF609319" ~ "Champsocephalus gunnari",
    read_name =="MW829388" ~ "Dissostichus eleginoides",
    read_name =="ON000293" ~ "Dissostichus mawsoni",
    read_name =="FJ582589" ~ "Eudyptes chrysocome",
    read_name =="FJ582593" ~ "Eudyptes chrysolophus",
    read_name =="MG739962" ~ "Eudyptes filholi",
    read_name =="MG740185" ~ "Eudyptes schlegeli",
    read_name =="EU525360" ~ "Eudyptula minor",
    read_name =="OL339430" ~ "Harpagifer bispinis",
    read_name =="ON891170" ~ "Harpagifer geogianus",
    read_name == "ON891164" ~ "Harpagifer kerguelensis",
    read_name =="MK262211" ~ "Megadyptes antipodes",
    read_name == "MK261839" ~ "Pygoscelis adeliae",
    read_name =="EU525471" ~ "COI - P. antarcticus",
    read_name =="KP875235" ~ "18S - P. antarcticus", 
    read_name =="OK493720" ~ "Trematomus loennbergii",
    read_name =="OK493634" ~ "Trematomus scotti",
    TRUE ~ read_name # In case of any other values, retain them unchanged
  ))

all_assignment_combo_syn  <- all_assignment %>%
  filter(read_name %in% accession) %>%
  mutate(read_name = case_when(
    read_name == "AF518190" ~ "18S - H. antarcticus",
    read_name == "ON891147" ~ "COI - H. antarcticus",
    read_name == "EU525299" ~ "Aptenodytes forsteri", 
    read_name == "EU525303" ~ "Aptenodytes patagonicus",
    read_name =="OK493624" ~ "Artedidraco lonnbergi",
    read_name =="EF609319" ~ "Champsocephalus gunnari",
    read_name =="MW829388" ~ "Dissostichus eleginoides",
    read_name =="ON000293" ~ "Dissostichus mawsoni",
    read_name =="FJ582589" ~ "Eudyptes chrysocome",
    read_name =="FJ582593" ~ "Eudyptes chrysolophus",
    read_name =="MG739962" ~ "Eudyptes filholi",
    read_name =="MG740185" ~ "Eudyptes schlegeli",
    read_name =="EU525360" ~ "Eudyptula minor",
    read_name =="OL339430" ~ "Harpagifer bispinis",
    read_name =="ON891170" ~ "Harpagifer geogianus",
    read_name == "ON891164" ~ "Harpagifer kerguelensis",
    read_name =="MK262211" ~ "Megadyptes antipodes",
    read_name == "MK261839" ~ "Pygoscelis adeliae",
    read_name =="EU525471" ~ "COI - P. antarcticus",
    read_name =="KP875235" ~ "18S - P. antarcticus", 
    read_name =="OK493720" ~ "Trematomus loennbergii",
    read_name =="OK493634" ~ "Trematomus scotti",
    TRUE ~ read_name # In case of any other values, retain them unchanged
  ))

################## try to create sensitivity (relative abundance - taxonomic assignment) read counts #####################
#### for synthetic only visualization normalized read counts ###### 
all_assignment_syn <- all_assignment_syn %>% 
  mutate(file_number = as.numeric(str_extract(file_name, "\\d+")))

all_assignment_normal <- all_assignment_syn %>% 
  mutate(nor_read_count = read_count/file_number) 

ggplot(all_assignment_normal, aes(x = file_name, y = read_name, fill = nor_read_count)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  theme_minimal() +
  mytheme_test_nofacet +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)) + 
  labs(
   #title = "Heatmap of Normalised Read Counts across Input Accessions",
    x = "Dataset",
    y = "Input Species", 
    fill = "Proportion fragments\nassigned correctly"
  )

ggsave("normal_synthetic_fullheatmap.png", width = 12, height = 7, units = "in", dpi = 300)
########## combined visualization normalized read counts #########
all_assignment_combo_test <- all_assignment_combo_syn %>% 
  mutate(file_number = str_extract(file_name, "\\d+"))

all_assignment_normal_combo <- all_assignment_combo_test %>% 
  mutate(nor_read_count = read_count/as.numeric(file_number)) 

ggplot(all_assignment_normal_combo, aes(x = file_name, y = read_name, fill = nor_read_count)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  theme_minimal() +
  mytheme_test_nofacet +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)) + 
  labs(
    #title = "Heatmap of Normalised Read Counts across Input Accessions",
    x = "Dataset",
    y = "Input Species", 
    fill = "Proportion fragments\nassigned correctly"
  )

ggsave("normal_combo_fullheatmap.png", width = 12, height = 7, units = "in", dpi = 300)


