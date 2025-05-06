# Read the FASTA file
fasta_file <- "ITS_58S_ITS.fasta"  # Replace with the path to your FASTA file
fasta_data <- readLines(fasta_file)

# Create a directory to save individual sequence files
output_dir <- "ITS_58S_sequences"
dir.create(output_dir, showWarnings = FALSE)

# Function to generate a new filename
generate_filename <- function(seq_id) {
  # Modify this function to generate filenames as per your requirement
  seq_id_cleaned <- gsub("[^A-Za-z0-9 ]", "", seq_id)  # Remove non-alphanumeric characters
  first_space_index <- regexpr(" ", seq_id_cleaned)  # Find the index of the first space
  if (first_space_index == -1) {
    return(paste0(seq_id_cleaned, ".fasta"))
  } else {
    return(paste0(substr(seq_id_cleaned, 1, first_space_index - 1), ".fasta"))
  }
}

# Function to write a sequence to a FASTA file
write_fasta <- function(seq_id, sequence, output_dir) {
  output_file <- file.path(output_dir, paste0(gsub("[^A-Za-z0-9]", "_", seq_id), ".fasta"))
  writeLines(c(paste0(">", seq_id), sequence), output_file)
}

# Split sequences into separate files
current_seq_id <- NULL
current_seq <- NULL
for (line in fasta_data) {
  if (substr(line, 1, 1) == ">") {
    # Start of a new sequence
    if (!is.null(current_seq_id)) {
      write_fasta(current_seq_id, current_seq, output_dir)
    }
    current_seq_id <- gsub("^>", "", line)
    current_seq <- character()
  } else {
    # Add sequence lines
    current_seq <- c(current_seq, line)
  }
}

# Write the last sequence to a file
if (!is.null(current_seq_id)) {
  write_fasta(current_seq_id, current_seq, output_dir)
}

cat("Splitting complete. Individual sequence files saved in", output_dir, "\n")
