# Clear environment
rm(list = ls())

# Load packages
library(topGO)
library(tidyverse)

# Define paths
raw_data_dir <- "00_RawData"
func_anot_dir <- file.path(raw_data_dir, "functional_annotation")
tidy_FC_dir <- "05_TidyFC"
out_dir <- "07_GeneEnrichment"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load GO annotations
go_annotations_raw <- read_tsv(file.path(func_anot_dir, "goinfo_FilteredModels2.tab"), col_names = TRUE) %>%
  dplyr::rename(
    protein_id = proteinId,
    go_id = goAcc,
    go_term = goName,
    go_type = gotermType
  )

# Build gene2GO list: named list of gene_id -> vector of GO terms
gene2go <- go_annotations_raw %>%
  group_by(protein_id) %>%
  summarise(go_terms = list(unique(go_id))) %>%
  deframe()

# Load expression metadata
data_labels <- read.csv(file.path(tidy_FC_dir, "tidy_all_strains_labels.csv"), row.names = 1)
ref_conditions <- unique(data_labels$condition)
strain_conditions <- gsub("-",".",unique(data_labels$strain_condition))

# Start logging
log_file <- file.path(out_dir, "topGO_log.txt")
sink(log_file, append = TRUE)

 ref_condition <- strain_conditions[2]
 strain_condition <- strain_conditions[7]
 deg_type<- "oe"
 ontology <- "BP"

# Loop over each ref and strain condition
for (ref_condition in strain_conditions) {
  message("Ref condition: ", ref_condition)
  
  deg_mat_path <- file.path(tidy_FC_dir, paste0("tidy_all_strains_DEGmat_ref_", ref_condition, ".csv"))
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
      warning("Column missing in DEG matrix: ", strain_condition)
      next
    }
    
    # Get gene universe (only those with GO annotations)
    gene_universe <- intersect(rownames(deg_mat), names(gene2go))
    
    # Get DEG subsets
    oe_genes <- rownames(deg_mat)[deg_mat[[strain_condition]] == "+"]
    ue_genes <- rownames(deg_mat)[deg_mat[[strain_condition]] == "-"]
    all_deg_genes <- union(oe_genes, ue_genes)
    
    deg_types <- list(
      oe = oe_genes,
      ue = ue_genes,
      deg = all_deg_genes
    )
    
    for (deg_type in names(deg_types)) {
      genes_of_interest <- deg_types[[deg_type]]
      
      gene_list <- factor(as.integer(gene_universe %in% genes_of_interest))
      names(gene_list) <- gene_universe
      
      if (sum(as.numeric(as.character(gene_list))) < 10) {
        message("    Skipping ", deg_type, ": not enough DEGs")
        next
      }
      
      for (ontology in c("BP", "MF", "CC")) {
        message("    → Ontology: ", ontology, " | DEG: ", deg_type)
        
        GOdata <- new("topGOdata",
                      ontology = ontology,
                      allGenes = gene_list,
                      annot = annFUN.gene2GO,
                      gene2GO = gene2go,
                      nodeSize = 10)
        
        result_fisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
        
        # Extract all GO terms tested, not just top N
        all_go_terms <- usedGO(GOdata)
        
        # Build full result table
        all_go_table <- GenTable(
          GOdata,
          Fisher = result_fisher,
          topNodes = length(all_go_terms)  # use all GO terms tested
        )
        
        # Add raw p-values
        raw_pvals <- score(result_fisher)[all_go_table$GO.ID]  # Ensure alignment
        all_go_table$Fisher_raw <- raw_pvals
        
        # Add adjusted p-values (FDR)
        all_go_table$Fisher_adj <- p.adjust(all_go_table$Fisher_raw, method = "BH")
        
        # Add DEG count (only genes that have GO terms)
        all_go_table$num_DEGs <- sum(as.numeric(as.character(gene_list)))
        
        # Save full table
        out_path <- file.path(
          out_dir,
          paste0("topGO_results_", ontology, "_", strain_condition, "_ref_", ref_condition, "_", deg_type, ".csv")
        )
        write.csv(all_go_table, out_path, row.names = FALSE)
        
      }
    }
  }
}

# Stop logging
sink()
