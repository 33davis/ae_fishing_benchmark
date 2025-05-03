library(ggplot2)


mytheme_nofacet <-   theme(legend.position = "bottom", 
                           axis.text.x = element_text(size = 12), 
                           legend.text = element_text(size = 12), 
                           legend.title = element_text(size = 14), 
                           axis.text.y = element_text(size = 12), 
                           axis.title = element_text(size =12), 
                           legend.key.width = unit(0.5, "cm"),
                           legend.spacing.x = unit(0.2, 'cm'))

# Set log-normal parameters
meanlog <- 3.99999  # Geometric mean (log-transformed)
sdlog <- 0.19  # Standard deviation (log-transformed)

min_frag_length <- 25
max_frag_length <- 100

# Generate a sample of fragment lengths
set.seed(42)  # For reproducibility
sample_lengths <- rlnorm(1000, meanlog = meanlog, sdlog = sdlog)

# Apply fragment length constraints
sample_lengths <- pmin(pmax(sample_lengths, min_frag_length), max_frag_length)

# Create the histogram
hist(sample_lengths, breaks = 30, main = "Fragment Length Distribution", 
     xlab = "Fragment Length", col = "lightblue")

# Add a legend to the plot, showing the input meanlog and sdlog values
legend("topright", 
       legend = c(paste("meanlog =", round(meanlog, 5)), 
                  paste("sdlog =", round(sdlog, 5))),
      # col = c("red", "blue"), 
      # lty = 2, 
      # lwd = 2,
       cex = 1,  # Scale the legend text
       bty = "n")  # No border for the legend


ggplot(data.frame(sample_lengths), aes(x = sample_lengths)) + 
  geom_histogram(binwidth = 2, fill = "lightblue", color = "black", aes(y = ..density..)) + 
  labs(
    title = "Fragment Length Distribution",
    x = "Fragment Length",
    y = "Density"
  ) + 
  annotate("label", 
           x = 85, 
           y = 0.04, 
           label = "U1538 geometric mean = 3.9999\nstandard deviation = 0.19", 
           hjust = 0.5) + 
  theme_minimal() +
  mytheme_nofacet

  