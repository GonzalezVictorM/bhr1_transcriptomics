# Clear environment
rm(list = ls())

# Load packages
library(tidyverse)

# Define input/output paths
enrich_dir <- "07_GeneEnrichment"
input_file <- file.path(enrich_dir, "aggregated_topGO_results.csv")
output_dir <- "08_GOplots"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Define the plotting order
strain_condition_order <- c("CBS_464_glycerol_5d",
                            "CBS_464_glucose_5d",
                            "CBS_464_mannose_5d",
                            "CBS_464_cellulose_5d",
                            "CBS_464_guargum_5d",
                            "CBS_464_xylan_5d",
                            "bmtr3_35.2_glycerol_5d",
                            "bmtr3_35.2_glucose_5d",
                            "bmtr3_35.2_mannose_5d",
                            "bmtr3_35.2_cellulose_5d",
                            "bmtr3_35.2_guargum_5d",
                            "bmtr3_35.2_xylan_5d")
strain_condition_order <- c("CBS_464_glycerol_5d","bmtr3_35.2_glycerol_5d",
                            "CBS_464_glucose_5d","bmtr3_35.2_glucose_5d",
                            "CBS_464_mannose_5d","bmtr3_35.2_mannose_5d",
                            "CBS_464_cellulose_5d","bmtr3_35.2_cellulose_5d",
                            "CBS_464_guargum_5d","bmtr3_35.2_guargum_5d",
                            "CBS_464_xylan_5d","bmtr3_35.2_xylan_5d" )


# Load aggregated GO enrichment results
go_data <- read.csv(input_file) %>%
  mutate(
    strain_condition = factor(strain_condition, levels = strain_condition_order),
    # Significance indicator for overexpression
    sig_Fisher_raw = case_when(
      Fisher_raw < 0.05 ~ "Significant",
      TRUE ~ "Non-significant"),
    sig_Fisher_adj = case_when(
      Fisher_adj < 0.05 ~ "Significant",
      TRUE ~ "Non-significant"),
    logFisher = -log10(Fisher_raw),
    logFDR = -log10(Fisher_adj),
    GO_label = paste(ontology,GO.ID,Term, sep = " - ")
  )

# Parameters to customize
ref_conditions <- unique(go_data$ref_condition)
deg_types <- c("oe", "ue", "deg")
p_val_cutoff <- 0.05   # cutoff for significance

# Loop through the reference conditions and
for (ref_condition in ref_conditions) {
  go_ref_data <- go_data %>%
    filter(ref_condition == !!ref_condition)
  
  for (deg_type in deg_types) {
    go_set_data <- go_ref_data %>%
      filter(deg_type == !!deg_type)
    
    plot_GOs <- go_set_data %>% 
      filter(Fisher_adj <= p_val_cutoff) %>%
      arrange(desc(logFDR)) %>%
      .$GO.ID %>% unique
    plot_data <- go_set_data %>% 
      filter(GO.ID %in% plot_GOs)
    
    # Check if enough data to plot
    if (nrow(plot_data) == 0) {
      message("No enriched GO terms found for the specified conditions.")
    } else {
      p <-  plot_data %>%
        ggplot(aes(y = reorder(GO_label, logFDR), x = strain_condition, size = Significant, color = logFDR)) +
        geom_point(aes(shape = sig_Fisher_adj), alpha = 0.8) +
        # Custom color scale for -log10 FDR
        scale_color_viridis_c(
          option = "D",  # Viridis D palette (magma)
          name = "-log10 FDR",
          limits = c(0, max(plot_data$logFDR, na.rm = TRUE)),
          breaks = seq(0, max(plot_data$logFDR, na.rm = TRUE), 
                       by = ceiling(max(plot_data$logFDR, na.rm = TRUE) / 5))
        ) +
        # Custom size scale for 
        scale_size_continuous(
          name = paste0(deg_type," genes"),
          range = c(2, 20),
          breaks = c(0, 5, 10, 15, 20),
          labels = c("0", "5", "10", "15", "20+")
        ) +
        # Shape for significance
        scale_shape_manual(
          name = "Significance (FDR < 0.05)",
          values = c("Significant" = 16, "Non-significant" = 1)
        ) +
        # Improve axis labels and theme
        labs(
          title = paste0("GO Enrichment for ", deg_type, " Genes (Ref: ",ref_condition,")"),
          y = "GO Term",
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
        filename = file.path(output_dir, paste0("topGO_enrichment_plot_ref_",ref_condition,"_",deg_type,".pdf")),
        plot = p,
        width = 20,
        height = 10,
        device = "pdf")
    }
  }
}
