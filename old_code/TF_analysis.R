# Clear environment (optional—consider if you really need to remove everything)
rm(list = ls())

# Load needed packages
library(tidyverse)
library(pheatmap)

# Define and validate paths
raw_data_dir <- "00_RawData"
tidy_data_dir <- "01_TidyData"
qc_dir <- "02_QualityCheck"
dds_dir <- "03_DESeqObjects"
deg_dir <- "04_DEGanalysis"
tidy_FC_dir <- "05_TidyFC"
out_dir <- "06_TFanalysis"

# Select export type
# export_file <- "png"
export_file <- "pdf"

# Define the plotting order
# strain_condition_order <- c("CBS_464_glycerol_5d",
#                             "CBS_464_glucose_5d",
#                             "CBS_464_mannose_5d",
#                             "CBS_464_cellulose_5d",
#                             "CBS_464_guargum_5d",
#                             "CBS_464_xylan_5d",
#                             "bmtr3_35.2_glycerol_5d",
#                             "bmtr3_35.2_glucose_5d",
#                             "bmtr3_35.2_mannose_5d",
#                             "bmtr3_35.2_cellulose_5d",
#                             "bmtr3_35.2_guargum_5d",
#                             "bmtr3_35.2_xylan_5d")
strain_condition_order <- c("CBS_464_glycerol_5d","bmtr3_35.2_glycerol_5d",
                            "CBS_464_glucose_5d","bmtr3_35.2_glucose_5d",
                            "CBS_464_mannose_5d","bmtr3_35.2_mannose_5d",
                            "CBS_464_cellulose_5d","bmtr3_35.2_cellulose_5d",
                            "CBS_464_guargum_5d","bmtr3_35.2_guargum_5d",
                            "CBS_464_xylan_5d","bmtr3_35.2_xylan_5d" )


# Ensure input directories exist
required_dirs <- c(raw_data_dir, tidy_data_dir, qc_dir, dds_dir, deg_dir,tidy_FC_dir)
lapply(required_dirs, function(d) if (!dir.exists(d)) stop("Directory not found: ", d))

# Create output directory if missing
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load expression matrix
data_mat <- read.csv(file.path(tidy_FC_dir, "tidy_all_strains_vst_avg.csv"), row.names = 1) %>% # change to other transformation if you want
  select(all_of(strain_condition_order))

data_labels <- read.csv(file.path(tidy_FC_dir, "tidy_all_strains_labels.csv"), row.names = 1)
ref_conditions <- unique(data_labels$condition)

# Load the gene list
sugarGenes <- read.csv(file.path(raw_data_dir,"Dsqul_masterList_genes_v5.csv")) %>%
  mutate(protein_id = as.character(protein_id))

TFGenes <- sugarGenes %>%
  filter(predicted_function == "TranscriptionFactors",
         note != "Remove_Corrected") %>%
  mutate(label = paste(protein_id, family, sep = " - "))

family_order <- sort(unique(TFGenes$family))

# Create the hm matrix
hm_data_mat <- data_mat[TFGenes$protein_id, , drop = FALSE]
hm_row_labels <- TFGenes %>%
  select (protein_id, label) %>%
  column_to_rownames("protein_id") 

# Store enrichment results
enrichment_results <- data.frame()

# condition <- ref_conditions[5] # for debugging purposes

# Loop over conditions
for (condition in ref_conditions) {
  message("Generating heatmap for ref condition ", condition)
  
  DEG_mat <- read.csv(file.path(tidy_FC_dir, paste0("tidy_all_strains_DEGmat_ref_", condition, ".csv")), row.names = 1)
  DEG_mat[is.na(DEG_mat)] <- " "
  hm_DEG_mat <- DEG_mat[rownames(hm_row_labels), colnames(hm_data_mat), drop = FALSE]
  
  # Heatmap output
  hm_file_name <- file.path(out_dir, paste0("heatmap_TF_ref_", condition, ".", export_file))
  
  pheatmap(hm_data_mat,
           # scale = "row",
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           clustering_distance_rows = "euclidean",
           display_numbers = hm_DEG_mat,
           show_rownames = TRUE,
           labels_row = hm_row_labels[rownames(hm_data_mat),],
           color = colorRampPalette(c("navy", "white", "red"))(50),
           border_color = NA,
           fontsize_row = 6,
           fontsize_col = 6,
           height = nrow(hm_data_mat)/7,
           gaps_col = seq(from = 2, to = ncol(hm_data_mat)-2, by = 2),
           filename = hm_file_name)
}
