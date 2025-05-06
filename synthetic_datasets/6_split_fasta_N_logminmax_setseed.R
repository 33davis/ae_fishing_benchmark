nsplit_fasta_files_lognormal <- function(input_dir, output_file, meanlog, sdlog, min_frag_length, max_frag_length, n, seed = NULL) {
  # Set seed for reproducibility
  if (!is.null(seed)) set.seed(seed)
  
  # Get list of FASTA files in the input directory
  fasta_files <- list.files(input_dir, pattern = "\\.fasta$", full.names = TRUE)
  
  # Process each n value
  for (current_n in n) {
    # Create dynamic output file name based on the current n value
    output_file_n <- gsub("n[0-9]+", paste0("n", current_n), output_file)
    
    # Open the output file for writing
    con <- file(output_file_n, "w")
    
    # Process each FASTA file
    for (file in fasta_files) {
      # Read FASTA file
      fasta_data <- readLines(file)
      
      # Process each sequence
      current_sequence <- NULL
      current_header <- NULL
      fragments_count <- 0
      
      while (fragments_count < current_n) {
        for (line in fasta_data) {
          if (substr(line, 1, 1) == ">") {
            # Save previous sequence
            if (!is.null(current_sequence) && nchar(current_sequence) > 0) {
              # Split the sequence into fragments
              nchar_remaining <- nchar(current_sequence)
              while (nchar_remaining > 0 && fragments_count < current_n) {
                fragment_length <- max(min(rlnorm(1, meanlog = meanlog, sdlog = sdlog), max_frag_length), min_frag_length)
                frag_len <- min(fragment_length, nchar_remaining)
                fragment <- substring(current_sequence, 1, frag_len)
                current_sequence <- substring(current_sequence, frag_len + 1)
                nchar_remaining <- nchar(current_sequence)
                fragments_count <- fragments_count + 1
                
                # Write each fragment to the output file
                cat(current_header, "\n", fragment, "\n", file = con, sep = "")
              }
            }
            
            # Start a new sequence
            current_header <- line
            current_sequence <- ""
          } else {
            # Concatenate sequence lines
            current_sequence <- paste0(current_sequence, line)
          }
        }
        
        # Save the last sequence
        if (!is.null(current_sequence) && nchar(current_sequence) > 0) {
          # Split the last sequence into fragments
          nchar_remaining <- nchar(current_sequence)
          while (nchar_remaining > 0 && fragments_count < current_n) {
            fragment_length <- max(min(rlnorm(1, meanlog = meanlog, sdlog = sdlog), max_frag_length), min_frag_length)
            frag_len <- min(fragment_length, nchar_remaining)
            fragment <- substring(current_sequence, 1, frag_len)
            current_sequence <- substring(current_sequence, frag_len + 1)
            nchar_remaining <- nchar(current_sequence)
            fragments_count <- fragments_count + 1
            
            # Write each fragment to the output file
            cat(current_header, "\n", fragment, "\n", file = con, sep = "")
          }
        }
        
        if (fragments_count >= current_n) {
          break  # Exit while loop if n fragments have been generated
        }
      }
    }
    
    # Close the output file
    close(con)
  }
}

# Example usage
setwd('~/Bioinformatics_Documentation/Dataset_Setup/datafiles/U1538/')
input_directory <- "U1538"
output_file <- "~/Bioinformatics_Documentation/Dataset_Setup/datafiles/synethic/test_setseed_lowersd//n500_iterative.fasta"
# meanlog <-  4.0432176
# sdlog <- 0.41833069
#meanlog <- 4.02
meanlog <-  3.99999
#sdlog <- 0.31033069 for cutoff run 150 
sdlong <- 0.19 #0.21033069 #for cutoff run 100 bp 
min_frag_length <- 25
max_frag_length <- 100
n <- c(25, 50, 100, 250, 500, 1000, 2500, 5000)  # Vector of n values

nsplit_fasta_files_lognormal(input_directory, output_file, meanlog, sdlog, min_frag_length, max_frag_length, n)
