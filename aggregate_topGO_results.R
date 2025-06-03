# Clear environment
rm(list = ls())

# Load packages
library(tidyverse)

# Define input directory
results_dir <- "07_GeneEnrichment"

# List all topGO results
go_files <- list.files(results_dir, pattern = "^topGO_results_.*\\.csv$", full.names = TRUE)

# Helper to parse metadata from filename
parse_filename <- function(fname) {
  parts <- stringr::str_match(
    basename(fname),
    "^topGO_results_(BP|MF|CC)_(.+?)_ref_(.+?)_(oe|ue|deg)\\.csv$"
  )
  tibble(
    ontology = parts[, 2],
    strain_condition = parts[, 3],
    ref_condition = parts[, 4],
    deg_type = parts[, 5]
  )
}

# Read and bind all topGO result files safely
go_results <- map_dfr(go_files, function(file) {
  df <- read_csv(file, show_col_types = FALSE)
  
  # Coerce problematic columns to correct types if needed
  df <- df %>%
    mutate(
      Fisher = as.numeric(Fisher),
      Fisher_raw = as.numeric(Fisher_raw),
      Fisher_adj = as.numeric(Fisher_adj),
      num_DEGs = as.integer(num_DEGs)
    )
  
  meta <- parse_filename(file)
  bind_cols(meta, df)
})

# Save full aggregated dataset
write_csv(go_results, file.path(results_dir, "aggregated_topGO_results.csv"))
message("✅ Aggregation complete: ", nrow(go_results), " rows.")

