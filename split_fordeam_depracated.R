
### this code will split fasta files into percentages (blind by species) to percentages
#This function will split each sequence in the input file according to the specified percentages and write the fragments to two output files. 
#Adjust the input_file, output_dir, percentage_1, and percentage_2 variables as needed for your specific case.
split_random_percentage <- function(input_file, output_dir, percentage_1, percentage_2) {
  # Check if percentages sum to 1
  if (percentage_1 + percentage_2 != 1) {
    stop("Error: The sum of percentage_1 and percentage_2 must equal 1.")
  }
  
  # Read FASTA file
  fasta_data <- readLines(input_file)
  
  # Process each sequence
  process_sequence <- function(sequence) {
    current_fragment <- NULL
    fragments <- character()
    
    for (line in sequence) {
      if (substr(line, 1, 1) == ">") {
        if (!is.null(current_fragment)) {
          fragments <- c(fragments, current_fragment)
        }
        current_fragment <- ""
      } else {
        current_fragment <- paste0(current_fragment, line)
      }
    }
    
    if (!is.null(current_fragment)) {
      fragments <- c(fragments, current_fragment)
    }
    
    n_fragments <- length(fragments)
    n_frag_1 <- floor(percentage_1 * n_fragments)
    selected_indices_1 <- sample(length(fragments), size = n_frag_1)
    fragments_1 <- fragments[selected_indices_1]
    fragments_2 <- fragments[-selected_indices_1]
    
    list(fragments_1 = fragments_1, fragments_2 = fragments_2)
  }
  
  # Split sequences using lapply
  split_sequences <- lapply(split(fasta_data, cumsum(substr(fasta_data, 1, 1) == ">")), process_sequence)
  
  # Get input file base name
  input_basename <- tools::file_path_sans_ext(basename(input_file))
  
  # Create output file names with the percentages and input filename
  output_file1 <- paste0(output_dir, "/", input_basename, "_", percentage_1*100, "%_1.fasta")
  output_file2 <- paste0(output_dir, "/", input_basename, "_", percentage_2*100, "%_2.fasta")
  
  # Open the output files for writing
  con1 <- file(output_file1, "w")
  con2 <- file(output_file2, "w")
  
  # Write fragments to the output files
  for (seq in split_sequences) {
    for (i in seq$fragments_1) {
      cat(i, "\n", file = con1)
    }
    for (i in seq$fragments_2) {
      cat(i, "\n", file = con2)
    }
  }
  
  # Close the output files
  close(con1)
  close(con2)
}

# Example usage
input_file <- "~/Bioinformatics_Documentation/Dataset_Setup/datafiles/U1538/test_split/output_1000geomedian.fasta"
output_dir <- "~/Bioinformatics_Documentation/Dataset_Setup/datafiles/U1538/split_4_deamSim/"
percentage_1 <- 0.97  # 67% of fragments go to output1
percentage_2 <- 0.03  # 33% of fragments go to output2

split_random_percentage(input_file, output_dir, percentage_1, percentage_2)
