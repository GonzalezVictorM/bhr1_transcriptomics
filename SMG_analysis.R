# Clear environment (optional—consider if you really need to remove everything)
rm(list = ls())

# Load needed packages
library(tidyverse)

# Define and validate paths
raw_data_dir <- "00_RawData"
tidy_data_dir <- "01_TidyData"
qc_dir <- "02_QualityCheck"
dds_dir <- "03_DESeqObjects"
deg_dir <- "04_DEGanalysis"
tidy_FC_dir <- "05_TidyFC"
out_dir <- "06_SMGanalysis"

# Define your mutant strain
mut_strain <- "bmtr3_35-2"

# Ensure input directories exist
required_dirs <- c(raw_data_dir, tidy_data_dir, qc_dir, dds_dir, deg_dir)
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

# Load the gene list
sugarGenes <- read.csv(file.path(raw_data_dir,"Dsqul_masterList_genes_v5.csv")) %>%
  mutate(protein_id = as.character(protein_id))

SMGenes <- sugarGenes %>%
  filter(predicted_function == "SugarMetabolicGenes") %>%
  mutate(pathway = gsub("tricarboxylic acid cycle","TCA", pathway),
         pathway = gsub("pentose catabolic","PCP", pathway))

SMGenes_short <- SMGenes %>%
  ungroup() %>%
  select(protein_id, pathway) %>%
  group_by(protein_id) %>%
  summarise(pathway = paste(pathway, collapse = "; ")) %>%
  mutate(label = paste(protein_id, pathway, sep = " - "))

pathway_order <- sort(unique(SMGenes$pathway))

# Create the hm matrix
hm_data_mat <- data_mat[SMGenes_short$protein_id, , drop = FALSE]
hm_row_labels <- SMGenes_short %>%
  select (protein_id, label) %>%
  column_to_rownames("protein_id")

# Store enrichment results
enrichment_results <- data.frame()

condition <- ref_conditions[5] # for debugging purposes

# Loop over conditions
for (condition in ref_conditions) {
  message("Generating heatmap for ref condition ", condition)
  
  DEG_mat <- read.csv(file.path(tidy_FC_dir, paste0("tidy_all_strains_DEGmat_ref_", condition, ".csv")), row.names = 1)
  DEG_mat[is.na(DEG_mat)] <- " "
  hm_DEG_mat <- DEG_mat[rownames(hm_row_labels), colnames(hm_data_mat), drop = FALSE]
  
  # Heatmap output
  hm_file_name <- file.path(out_dir, paste0("heatmap_SMG_ref_", condition, ".pdf"))
  
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
           filename = hm_file_name)
  
  # === Hypergeometric test (hgt)===
  hgt_total_data <- DEG_mat %>%
    rownames_to_column("protein_id") %>%
    pivot_longer(-protein_id, names_to = "strain_condition", values_to = "Expression") %>%
    # Map Expression → symbol
    mutate(
      ref_sub = condition,
      Expression = case_when(
        Expression == "+" ~ "Overexpressed",
        Expression == "-" ~ "Underexpressed",
        TRUE                           ~ "Lowly expressed"
      )) %>%
    group_by(ref_sub, strain_condition) %>%
    summarise(total_genes = n(),
              total_oeg = sum(Expression == "Overexpressed"),
              total_ueg = sum(Expression == "Underexpressed")) %>%
    ungroup()
  
  hgt_set_data <- hm_DEG_mat %>%
    rownames_to_column("protein_id") %>%
    pivot_longer(-protein_id, names_to = "strain_condition", values_to = "Expression") %>%
    # right_join(SMGenes[c("protein_id","pathway")], relationship = "many-to-many") %>%
    right_join(unique(SMGenes[c("protein_id","pathway")]), relationship = "many-to-many") %>%
    # Map Expression → symbol
    mutate(
      ref_sub = condition,
      Expression = case_when(
        Expression == "+" ~ "Overexpressed",
        Expression == "-" ~ "Underexpressed",
        TRUE                           ~ "Lowly expressed"
      )) %>%
    group_by(ref_sub, strain_condition, pathway) %>%
    summarise(set_genes = n(),
              set_oeg = sum(Expression == "Overexpressed"),
              set_ueg = sum(Expression == "Underexpressed"))%>%
    ungroup()
  
  hgt_data <- full_join(hgt_total_data,hgt_set_data)
  
  enrichment_results <- rbind(enrichment_results, hgt_data)
  
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
  mutate(strain_condition = factor(strain_condition, levels = rev(colnames(data_mat)))) %>%
  # Order gene_set by total overexpressed genes for better visualization
  arrange(desc(set_oeg)) %>%
  mutate(pathway = factor(pathway, levels = pathway_order))

for (substrate in ref_conditions) {
  p_data <- enrichment_data %>%
    filter(ref_sub == substrate) 
  p <-  p_data %>%
    ggplot(aes(x = pathway, y = strain_condition, size = -log10(padj_oe), color = set_oeg)) +
    geom_point(aes(shape = sig_oe), alpha = 0.8) +
    # Custom color scale for overexpressed gene counts
    scale_color_gradientn(
      colors = c("blue", "white", "red"),
      name = "Overexpressed Genes",
      limits = c(0, max(p_data$set_oeg, na.rm = TRUE)),
      breaks = seq(0, max(p_data$set_oeg, na.rm = TRUE), by = ceiling(max(p_data$set_oeg, na.rm = TRUE)/5))
    ) +
    # Custom size scale for -log10 adjusted p-values
    scale_size_continuous(
      name = "-Log10 Adjusted P-value (Overexpression)",
      range = c(2, 10),
      breaks = c(0, 1, 2, 3, 4),
      labels = c("0", "1", "2", "3", "4+")
    ) +
    # Shape for significance
    scale_shape_manual(
      name = "Significance (FDR < 0.05)",
      values = c("Significant" = 16, "Non-significant" = 1)
    ) +
    # Improve axis labels and theme
    labs(
      title = paste0("Gene Set Enrichment for Overexpressed CAZy Genes (Ref: ",substrate,")"),
      x = "ST predicted specificity",
      y = "Strain and Condition"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      legend.box = "vertical",
      legend.margin = margin(t = 0, r = 10, b = 0, l = 10)
    ) +
    # Ensure y-axis labels are clear
    scale_y_discrete(labels = function(x) gsub("\\.", "-", x))
  
  # Save plot
  ggsave(
    filename = file.path(out_dir, paste0("enrichment_plot_ref_",substrate,".pdf")),
    plot = p,
    width = 10,
    height = 8,
    device = "pdf"
  )
}
