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
out_dir <- "06_STanalysis"

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
strain_condition_order <- c(#"CBS_464_glycerol_5d","bmtr3_35.2_glycerol_5d",
                            "CBS_464_glucose_5d","bmtr3_35.2_glucose_5d",
                            "CBS_464_mannose_5d","bmtr3_35.2_mannose_5d",
                            "CBS_464_cellulose_5d","bmtr3_35.2_cellulose_5d",
                            "CBS_464_guargum_5d","bmtr3_35.2_guargum_5d",
                            "CBS_464_xylan_5d","bmtr3_35.2_xylan_5d" )

# Ensure input directories exist
required_dirs <- c(raw_data_dir, tidy_data_dir, qc_dir, dds_dir, deg_dir, tidy_FC_dir)
lapply(required_dirs, function(d) if (!dir.exists(d)) stop("Directory not found: ", d))

# Create output directory if missing
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load expression matrix
data_mat <- read.csv(file.path(tidy_FC_dir, "tidy_all_strains_vst_avg.csv"), row.names = 1) %>% # change to other transformation if you want
  select(all_of(strain_condition_order))

data_labels <- read.csv(file.path(tidy_FC_dir, "tidy_all_strains_labels.csv"), row.names = 1)
ref_conditions <- unique(data_labels$condition)
ref_strains <- gsub("-",".",unique(data_labels$strain))
deg_types <- c("oe", "ue")

# Load the gene list
sugarGenes <- read.csv(file.path(raw_data_dir,"Dsqul_masterList_genes_v5.csv")) %>%
  mutate(protein_id = as.character(protein_id))

STGenes <- sugarGenes %>%
  filter(predicted_function == "SugarTransporter") %>%
  dplyr::rename(specificity = substrate) %>%
  mutate(label = paste(protein_id, specificity, sep = " - "))

specificity_order <- c( "cellodextrin/lactose","D-galacturonic, quinic acid",
                      "glucose, pentose","hexose","inositol, fructose",
                      "pentose, glycerol, arabitol","xylose","unknown group")

# Create the hm matrix
hm_data_mat <- data_mat[STGenes$protein_id, , drop = FALSE]
hm_row_labels <- STGenes %>%
  select (protein_id, label) %>%
  column_to_rownames("protein_id") 

# Store enrichment results
enrichment_results <- data.frame()

condition <- ref_conditions[5] # for debugging purposes

# Loop over conditions
for (strain in ref_strains) {
  message("Generating heatmap for ref strain ", strain)
  hm_DEG_mat <- data.frame(protein_id = rownames(hm_row_labels))
  
  for (condition in ref_conditions) {
    DEG_mat <- read.csv(file.path(tidy_FC_dir, paste0("tidy_all_strains_DEGmat_ref_",strain, "_", condition, ".csv")), row.names = 1) %>%
      rownames_to_column("protein_id") %>%
      select(protein_id, grep(condition,colnames(hm_data_mat), value = T))
    DEG_mat[is.na(DEG_mat)] <- " "
    hm_DEG_mat <- left_join(hm_DEG_mat, DEG_mat)
  
    # === Hypergeometric test (hgt)===
    hgt_total_data <- DEG_mat %>%
      pivot_longer(-protein_id, names_to = "strain_condition", values_to = "Expression") %>%
      # Map Expression → symbol
      mutate(
        ref_strain = strain,
        ref_sub = condition,
        Expression = case_when(
          Expression == "+" ~ "Overexpressed",
          Expression == "-" ~ "Underexpressed",
          TRUE                           ~ "Lowly expressed"
        )) %>%
      group_by(ref_strain, ref_sub, strain_condition) %>%
      summarise(total_genes = n(),
                total_oeg = sum(Expression == "Overexpressed"),
                total_ueg = sum(Expression == "Underexpressed")) %>%
      ungroup()
    
    hgt_set_data <- hm_DEG_mat %>%
      select(protein_id, grep(condition,colnames(hm_data_mat), value = T)) %>%
      pivot_longer(-protein_id, names_to = "strain_condition", values_to = "Expression") %>%
      left_join(STGenes[c("protein_id","specificity")]) %>%
      # Map Expression → symbol
      mutate(
        ref_strain = strain,
        ref_sub = condition,
        Expression = case_when(
          Expression == "+" ~ "Overexpressed",
          Expression == "-" ~ "Underexpressed",
          TRUE                           ~ "Lowly expressed"
        )) %>%
      group_by(ref_strain, ref_sub, strain_condition, specificity) %>%
      summarise(set_genes = n(),
                set_oeg = sum(Expression == "Overexpressed"),
                set_ueg = sum(Expression == "Underexpressed"))%>%
      ungroup()
    
    hgt_data <- full_join(hgt_total_data,hgt_set_data)
    
    enrichment_results <- rbind(enrichment_results,hgt_data)
  }
  
  # Heatmap output
  hm_file_name <- file.path(out_dir, paste0("heatmap_ST_ref_", strain, ".", export_file))
  
  hm_DEG_mat <- hm_DEG_mat %>%
    column_to_rownames("protein_id")
  
  pheatmap(hm_data_mat,
           # scale = "row",
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           clustering_distance_rows = "euclidean",
           display_numbers = hm_DEG_mat[rownames(hm_data_mat), colnames(hm_data_mat)],
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

enrichment_data <- enrichment_results %>%
  rowwise() %>%
  mutate(
    pval_oe = phyper(set_oeg - 1, total_oeg, total_genes - total_oeg, set_genes, lower.tail = FALSE),
    pval_ue = phyper(set_ueg - 1, total_ueg, total_genes - total_ueg, set_genes, lower.tail = FALSE)
  ) %>%
  ungroup() %>%
  # Adjust p-values for multiple testing (FDR)
  mutate(
    padj_oe = p.adjust(pval_oe, method = "fdr"),
    padj_ue = p.adjust(pval_ue, method = "fdr"),
    # Significance indicator for overexpression
    sig_oe = case_when(
      padj_oe < 0.05 ~ "Significant",
      TRUE ~ "Non-significant"),
    sig_ue = case_when(
      padj_ue < 0.05 ~ "Significant",
      TRUE ~ "Non-significant")
  ) %>%
  # Order strain_condition by data_mat columns for consistency
  mutate(strain_condition = factor(strain_condition, levels = colnames(data_mat)),
         specificity = factor(specificity, levels = rev(specificity_order)))

#substrate <- ref_conditions[5] # for debugging
#deg_type <- deg_types[1] # for debugging

# Loop through the reference conditions and the deg type
for (strain in ref_strains) {
  set_data <- enrichment_data %>%
    filter(ref_strain == strain, !grepl(strain, strain_condition))
  
  for (deg_type in deg_types) {
    plot_data <- set_data %>%
      mutate(Significant = case_when(
        deg_type == "oe" ~ set_oeg,
        TRUE ~ set_ueg),
        sig_Fisher_adj = case_when(
          deg_type == "oe" ~ sig_oe,
          TRUE ~ sig_ue),
        Fisher_adj = case_when(
          deg_type == "oe" ~ padj_oe,
          TRUE ~ padj_ue),
        logFDR = -log10(Fisher_adj))
    
    p <- plot_data %>%
      ggplot(aes(y = specificity, x = strain_condition, size = Significant, color = logFDR)) +
      geom_point(aes(shape = sig_Fisher_adj), alpha = 0.8) +
      # Custom color scale for -log10 adjusted p-values
      scale_color_viridis_c(
        option = "D",  # Viridis D palette (magma)
        name = "-log10 FDR",
        limits = c(0, max(plot_data$logFDR, na.rm = TRUE)),
        breaks = seq(0, max(plot_data$logFDR, na.rm = TRUE), 
                     by = ceiling(max(plot_data$logFDR, na.rm = TRUE) / 5))
      ) +
      # Custom size scale for overexpressed gene counts
      scale_size_continuous(
        name = paste0(deg_type," genes"),
        range = c(2, 18),
        breaks = c(0, 2, 4),
        labels = c("0", "2", "4+")
      ) +
      # Shape for significance
      scale_shape_manual(
        name = "Significance (FDR < 0.05)",
        values = c("Significant" = 16, "Non-significant" = 1)
      ) +
      # Improve axis labels and theme
      labs(
        title = paste0("Gene Set Enrichment for sugar transporter ", deg_type, " genes (Ref: ",strain,")"),
        y = "Predicted pathway",
        x = "Strain and Condition"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        legend.box = "vertical",
        legend.margin = margin(t = 0, r = 10, b = 0, l = 10)
      )
    
    # Save plot
    ggsave(
      filename = file.path(out_dir, paste0("enrichment_plot_ref_",strain,"_",deg_type,".", export_file)),
      plot = p,
      width = 20,
      height = 16,
      device = export_file
    )
  }
}

