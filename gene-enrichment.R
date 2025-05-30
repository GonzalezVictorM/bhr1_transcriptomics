# Clear environment (optional—consider if you really need to remove everything)
rm(list = ls())

# Load needed packages
library(tidyverse)

# Define and validate paths
raw_data_dir <- "00_RawData"
func_anot_dir <- file.path(raw_data_dir,"functional_annotation")
tidy_data_dir <- "01_TidyData"
qc_dir <- "02_QualityCheck"
dds_dir <- "03_DESeqObjects"
deg_dir <- "04_DEGanalysis"
tidy_FC_dir <- "05_TidyFC"
out_dir <- "07_GeneEnrichment"


# Define your mutant strain
mut_strain <- "bmtr3_35-2"

# Ensure input directories exist
required_dirs <- c(raw_data_dir, func_anot_dir, tidy_data_dir, qc_dir, dds_dir, deg_dir, tidy_FC_dir)
lapply(required_dirs, function(d) if (!dir.exists(d)) stop("Directory not found: ", d))

# Create output directory if missing
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load expression matrix
data_mat <- read.csv(file.path(tidy_FC_dir, "tidy_all_strains_vst_avg.csv"), row.names = 1) %>% # change to other transformation if you want
  select( # sadly I have not figured out how to automate this part
    CBS_464_glycerol_5d,
    CBS_464_glucose_5d,
    CBS_464_mannose_5d,
    CBS_464_cellulose_5d,
    CBS_464_guargum_5d,
    CBS_464_xylan_5d,
    bmtr3_35.2_glycerol_5d,
    bmtr3_35.2_glucose_5d,
    bmtr3_35.2_mannose_5d,
    bmtr3_35.2_cellulose_5d,
    bmtr3_35.2_guargum_5d,
    bmtr3_35.2_xylan_5d)

data_labels <- read.csv(file.path(tidy_FC_dir, "tidy_all_strains_labels.csv"), row.names = 1)
ref_conditions <- unique(data_labels$condition)

# Load the gene list with gene ontology terms
go_annot <- read_tsv(file.path(func_anot_dir, "goinfo_FilteredModels2.tab"), col_names = TRUE) %>%
  dplyr::rename(
    protein_id = proteinId,
    go_id = goAcc,
    go_term = goName,
    go_type = gotermType
  )

# Loop over conditions
for (condition in ref_conditions) {
  message("Generating Gene enrichment for ref condition ", condition)
  
  DEG_mat <- read.csv(file.path(tidy_FC_dir, paste0("tidy_all_strains_DEGmat_ref_", condition, ".csv")), row.names = 1)
  DEG_mat[is.na(DEG_mat)] <- " "
  
 # calculate gene enrihment here
}