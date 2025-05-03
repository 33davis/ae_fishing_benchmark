# Load necessary library
library(seqinr)

# Function to read a FASTA file and return a list of sequences with headers
read_fasta <- function(file_path) {
  fasta_data <- read.fasta(file = file_path, seqtype = "DNA", as.string = TRUE)
  return(fasta_data)
}

# Function to write selected sequences to a new FASTA file
write_fasta <- function(sequences, output_path) {
  write.fasta(sequences, names = names(sequences), file.out = output_path, as.string = TRUE)
}

# Function to number sequences and select a random subset
select_random_sequences <- function(fasta_data, n) {
  total_sequences <- length(fasta_data)
  if (n > total_sequences) {
    stop("Number of sequences to select exceeds the total number of sequences available.")
  }
  
  selected_indices <- sample(seq_len(total_sequences), n, replace = FALSE)
  selected_sequences <- fasta_data[selected_indices]
  non_selected_sequences <- fasta_data[-selected_indices]
  
  return(list(selected = selected_sequences, non_selected = non_selected_sequences))
}

# Main script
input_dir <- "./test_setseed_lowersd/"  # Directory containing input FASTA files
output_dir <- "./test_setseed_lowersd/combo_files/"  # Directory to save output files
percentages <- c(0.03, 0.07, 0.18)  # Vector of percentages to select

# List all FASTA files in the input directory
fasta_files <- list.files(input_dir, pattern = "\\.fasta$", full.names = TRUE)

# Check if output directory exists; if not, create it
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop through each FASTA file
for (fasta_file in fasta_files) {
  # Get the base name of the FASTA file (without directory path)
  base_name <- basename(fasta_file)
  file_name <- tools::file_path_sans_ext(base_name)
  
  # Read the FASTA file
  fasta_data <- read_fasta(fasta_file)
  
  # Loop through each percentage
  for (percentage in percentages) {
    # Calculate the number of random sequences to select
    total_sequences <- length(fasta_data)
    n <- round(percentage * total_sequences)
    
    # Select random sequences and get non-selected sequences
    sequences <- select_random_sequences(fasta_data, n)
    selected_sequences <- sequences$selected
    non_selected_sequences <- sequences$non_selected
    
    # Define output file names
    output_selected_file <- file.path(output_dir, paste0(file_name, "_selected_", percentage * 100, ".fasta"))
    output_non_selected_file <- file.path(output_dir, paste0(file_name, "_non_selected_", percentage * 100, ".fasta"))
    
    # Write the selected sequences to the output FASTA file
    write_fasta(selected_sequences, output_selected_file)
    
    # Write the non-selected sequences to another output FASTA file
    write_fasta(non_selected_sequences, output_non_selected_file)
    
    cat("Random selection of", percentage * 100, "% sequences from", base_name, "written to", output_selected_file, "\n")
    cat("Non-selected sequences from", base_name, "written to", output_non_selected_file, "\n")
  }
}


