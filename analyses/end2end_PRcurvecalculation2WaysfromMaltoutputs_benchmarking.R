# ================================
# End-to-End Benchmarking PR curve calculation for MALT outputs
# ================================
library(tidyverse)
library(PRROC)     # For PR curves
library(viridis)   # For plasma colormap

# ----------------
# 1. Import Data
# ----------------
# Assume you ran: rma2info -i sampleX.rma6 -o sampleX.tsv -r2c Taxonomy
# Each TSV has columns like: read_id, taxon_id, confidence

# Load all sample TSVs
file_list <- list.files("results/malt_exports/", pattern = "*.tsv", full.names = TRUE)

malt_df <- file_list %>%
  map_dfr(~ read_tsv(.x) %>%
            mutate(sample_id = basename(.x) %>% str_remove(".tsv"))
  )

# Example structure after binding:
# read_id | taxon_id | confidence | sample_id

# -------------------------
# 2. Add Ground Truth Labels
# -------------------------
# Suppose you have a table of synthetic truth per read
# e.g. truth_df: read_id | true_taxon
truth_df <- read_tsv("data/ground_truth_reads.tsv")

# Join into malt_df
joined_df <- malt_df %>%
  left_join(truth_df, by = "read_id") %>%
  rename(predicted_taxon = taxon_id)

# -------------------------------
# 3. READ-LEVEL Precision/Recall
# -------------------------------
# Define TP if predicted_taxon == true_taxon
# Otherwise FP if predicted, FN if not predicted

read_eval <- joined_df %>%
  mutate(correct = predicted_taxon == true_taxon) 

# Create a PR curve object using confidence as score
fg <- read_eval %>% filter(correct) %>% pull(confidence)  # scores of correct predictions
bg <- read_eval %>% filter(!correct) %>% pull(confidence) # scores of incorrect predictions

pr_read <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)

# -------------------------------
# 4. SAMPLE-LEVEL Precision/Recall
# -------------------------------
# Suppose you also have ground truth taxa per sample at each rank
truth_taxa <- read_tsv("data/ground_truth_taxa.tsv")
# schema: sample_id | rank | true_taxa (comma separated)

# Collapse predictions by sample + rank
pred_taxa <- joined_df %>%
  group_by(sample_id, predicted_taxon) %>%
  summarise(max_conf = max(confidence), .groups = "drop") %>%
  mutate(rank = "Genus")  # or parse rank if available

# Convert to presence/absence set at chosen confidence cutoff
cutoff <- 0.9
pred_sets <- pred_taxa %>%
  filter(max_conf >= cutoff) %>%
  group_by(sample_id, rank) %>%
  summarise(pred_taxa = list(unique(predicted_taxon)), .groups = "drop")

# Join with truth
sample_eval <- truth_taxa %>%
  left_join(pred_sets, by = c("sample_id", "rank"))

# Compute precision/recall
sample_eval <- sample_eval %>%
  rowwise() %>%
  mutate(
    tp = length(intersect(true_taxa, pred_taxa)),
    fp = length(setdiff(pred_taxa, true_taxa)),
    fn = length(setdiff(true_taxa, pred_taxa)),
    precision = ifelse(tp + fp == 0, NA, tp / (tp + fp)),
    recall    = ifelse(tp + fn == 0, NA, tp / (tp + fn))
  )

# -------------------
# 5. Plotting PR Curves
# -------------------

# READ-level PR curve
autoplot_pr <- function(pr_obj, title) {
  pr_df <- data.frame(recall = pr_obj$curve[,1],
                      precision = pr_obj$curve[,2],
                      threshold = pr_obj$curve[,3])
  
  ggplot(pr_df, aes(x = recall, y = precision, color = threshold)) +
    geom_line(size = 1.2) +
    scale_color_viridis(option = "C") +
    theme_minimal(base_size = 14) +
    labs(title = title, x = "Recall", y = "Precision", color = "Threshold")
}

p1 <- autoplot_pr(pr_read, "Read-level PR Curve")

# SAMPLE-level PR curve (need to vary cutoff)
sample_pr <- map_dfr(seq(0, 1, by = 0.1), function(cut) {
  pred_sets <- pred_taxa %>%
    filter(max_conf >= cut) %>%
    group_by(sample_id, rank) %>%
    summarise(pred_taxa = list(unique(predicted_taxon)), .groups = "drop")
  
  eval_df <- truth_taxa %>%
    left_join(pred_sets, by = c("sample_id", "rank")) %>%
    rowwise() %>%
    mutate(
      tp = length(intersect(true_taxa, pred_taxa)),
      fp = length(setdiff(pred_taxa, true_taxa)),
      fn = length(setdiff(true_taxa, pred_taxa)),
      precision = ifelse(tp + fp == 0, NA, tp / (tp + fp)),
      recall    = ifelse(tp + fn == 0, NA, tp / (tp + fn)),
      cutoff = cut
    )
  eval_df
})

p2 <- ggplot(sample_pr, aes(x = recall, y = precision, color = cutoff, group = cutoff)) +
  geom_line(size = 1.2) +
  scale_color_viridis(option = "C") +
  theme_minimal(base_size = 14) +
  labs(title = "Sample-level PR Curve", x = "Recall", y = "Precision", color = "Cutoff")

# ----------------
# 6. Save Plots
# ----------------
ggsave("read_level_pr.png", p1, width = 6, height = 5, dpi = 300)
ggsave("sample_level_pr.png", p2, width = 6, height = 5, dpi = 300)

