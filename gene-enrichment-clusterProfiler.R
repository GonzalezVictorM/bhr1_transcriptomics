# Clear environment
rm(list = ls())

# Load required packages
library(clusterProfiler)
library(tidyverse)

# Load organism-specific annotation package if available (e.g., org.Sc.sgd.db for yeast)
# library(org.Sc.sgd.db)  # Example
# If not available, provide TERM2GENE manually

# Define paths
data_dir <- "05_TidyFC"
func_anot_dir <- "00_RawData/functional_annotation"
out_dir <- "07_GeneEnrichment"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load GO annotation file (custom mapping)
go_annotations <- read_tsv(file.path(func_anot_dir, "goinfo_FilteredModels2.tab"), col_names = TRUE) %>%
  rename(protein_id = proteinId, go_id = goAcc, go_term = goName, go_type = gotermType)

# Prepare TERM2GENE and TERM2NAME data.frames
TERM2GENE <- go_annotations %>% select(go_id, protein_id)
TERM2NAME <- go_annotations %>% select(go_id, go_term) %>% distinct()

# Load expression metadata
data_labels <- read_csv(file.path(data_dir, "tidy_all_strains_labels.csv"), show_col_types = FALSE)
ref_conditions <- unique(data_labels$condition)
strain_conditions <- gsub("-",".",unique(data_labels$strain_condition))


ref_condition <- ref_conditions[1]
strain_condition <- strain_conditions[2]
deg_type <- "oe"

# Loop through reference and strain conditions
for (ref_condition in ref_conditions) {
  message("Ref condition: ", ref_condition)
  
  deg_mat_path <- file.path(data_dir, paste0("tidy_all_strains_DEGmat_ref_", ref_condition, ".csv"))
  if (!file.exists(deg_mat_path)) {
    warning("Missing DEG matrix for ref: ", ref_condition)
    next
  }
  
  deg_mat <- read.csv(deg_mat_path, row.names = 1)
  deg_mat[is.na(deg_mat)] <- " "
  
  for (strain_condition in strain_conditions) {
    message("  → Strain condition: ", strain_condition)
    
    # Skip if strain condition is same as reference
    if (stringr::str_ends(strain_condition, ref_condition)) {
      message("Skipping self-comparison: ", strain_condition, " vs ref ", ref_condition)
      next
    }
    
    if (!strain_condition %in% colnames(deg_mat)) {
      warning("Missing column in DEG matrix: ", strain_condition)
      next
    }
    
    # Gene universe and DEG types
    all_genes <- rownames(deg_mat)
    gene_status <- deg_mat[[strain_condition]]
    oe_genes <- rownames(deg_mat)[gene_status == "+"]
    ue_genes <- rownames(deg_mat)[gene_status == "-"]
    all_deg <- union(oe_genes, ue_genes)
    
    deg_types <- list(oe = oe_genes, ue = ue_genes, deg = all_deg)
    
    for (deg_type in names(deg_types)) {
      genes_of_interest <- deg_types[[deg_type]]
      
      if (length(genes_of_interest) < 10) {
        message("    Skipping ", deg_type, ": too few DEGs")
        next
      }
      
      for (ontology in c("BP", "MF", "CC")) {
        message("    → Ontology: ", ontology, " | DEG: ", deg_type)
        
        ego <- enrichGO(
          gene          = genes_of_interest,
          universe      = all_genes,
          OrgDb         = NULL,
          keyType       = "GENEID",
          ont           = ontology,
          pAdjustMethod = "BH",
          pvalueCutoff  = 1,
          qvalueCutoff  = 1,
          readable      = FALSE,
          TERM2GENE     = TERM2GENE,
          TERM2NAME     = TERM2NAME
        )
        
        result_df <- as.data.frame(ego)
        
        # Save to CSV
        out_file <- file.path(out_dir, paste0("clusterProfiler_results_", ontology, "_", strain_condition, "_ref_", ref_condition, "_", deg_type, ".csv"))
        write.csv(result_df, out_file, row.names = FALSE)
      }
    }
  }
}
