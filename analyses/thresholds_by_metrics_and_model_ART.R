# =========================================
# Title: Extract thresholds from ART results by model type
# Input: tidied ART results table (art_all)
# Output: CSV of thresholds by metric and model
# =========================================

# Load packages
library(dplyr)
library(readr)

# ------------------------------
# Helper function: Threshold extraction
# ------------------------------
find_threshold_art <- function(art_tbl, metric_name, model_type = "pooled", ref_complexity = "5000") {
  
  art_tbl %>%
    filter(metric == metric_name, model == model_type) %>%  # filter by metric and model
    mutate(contrast = as.character(contrast)) %>%           # ensure character
    filter(grepl(ref_complexity, contrast)) %>%            # keep contrasts involving reference
    filter(p.value >= 0.05) %>%                            # non-significant comparisons
    rowwise() %>%
    mutate(candidate = strsplit(contrast, " - ")[[1]][
      which(strsplit(contrast, " - ")[[1]] != paste0("complexity", ref_complexity))
    ]) %>%
    ungroup() %>%
    mutate(candidate = as.numeric(gsub("complexity", "", candidate))) %>%
    summarise(threshold = min(candidate, na.rm = TRUE)) %>%
    pull(threshold)
}

# ------------------------------
# Load ART results table
# ------------------------------
# Example: art_all.csv exported previously
art_all <- read_csv("art_results.csv")

# ------------------------------
# Extract thresholds for each metric & model
# ------------------------------
metrics <- c("precision", "sensitivity", "F1")
models  <- c("pooled", "paired")

thresholds_art <- expand.grid(metric = metrics, model = models) %>%
  rowwise() %>%
  mutate(threshold = find_threshold_art(art_all, metric, model)) %>%
  ungroup() %>%
  arrange(metric, model)

# ------------------------------
# Save thresholds
# ------------------------------
write_csv(thresholds_art, "Thresholds_by_metric_and_model_ART.csv")

# Print to console
print(thresholds_art)
