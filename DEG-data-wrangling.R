# Clear environment (optional—consider if you really need to remove everything)
rm(list = ls())

# Load needed packages
library(tidyverse)

# Define and validate paths
tidy_data_dir <- "01_TidyData"
qc_dir <- "02_QualityCheck"
dds_dir <- "03_DESeqObjects"
deg_dir <- "04_DEGanalysis"
out_dir <- "05_TidyFC"

# Define your mutant strain
mut_strain <- "bmtr3_35-2"

# Ensure input directory exists
if (!dir.exists(tidy_data_dir)) {
  stop("Directory not found: ", tidy_data_dir)
}
if (!dir.exists(qc_dir)) {
  stop("Directory not found: ", qc_dir)
}
if (!dir.exists(dds_dir)) {
  stop("Directory not found: ", dds_dir)
}
if (!dir.exists(deg_dir)) {
  stop("Directory not found: ", deg_dir)
}

# Create output directory if missing
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Load data
all_labels <- read.csv(file.path(tidy_data_dir, "all_strains_labels.csv")) %>%
  mutate(condition = paste(cSource, time, sep = "_"),
         strain_condition = paste(strain, condition, sep = "_"))

all_counts <- read.csv(file.path(tidy_data_dir, "all_strains_counts.csv")) %>%
  mutate(protein_id = as.character(protein_id))

all_fpkms <- read.csv(file.path(tidy_data_dir, "all_strains_fpkms.csv")) %>%
  mutate(protein_id = as.character(protein_id))
  

all_vst <- read.csv(file.path(dds_dir, "vst_CBS_464_all_strains.csv")) %>%
  dplyr::rename(protein_id = X) %>%
  mutate(protein_id = as.character(protein_id))

## (Addition) filter only the labels that you are interested in
strains <- gsub("-",".",unique(all_labels$strain))

conditions_in_mutant <- all_labels %>% 
  filter(strain == mut_strain) %>% 
  .$condition %>% 
  as.character() %>% unique

all_labels <- all_labels %>%
  filter(condition %in% conditions_in_mutant)

all_counts <- all_counts %>% 
  select(protein_id, all_labels$sample)

all_fpkms <- all_fpkms %>% 
  select(protein_id, all_labels$sample)

# Load outliers (from previous script)
outlier_file <- file.path(qc_dir, "replicates_low_corr_to_all_peers.csv")

# If outliers exist, read and remove from data
if (file.exists(outlier_file)) {
  outliers <- read.csv(outlier_file)
  outlier_ids <- outliers$sample1
  message("Outliers detected and will be removed:")
  print(outlier_ids)
  
  all_labels <- all_labels %>%
    filter(!sample %in% outlier_ids)
  
  all_counts <- all_counts %>%
    select(-all_of(outlier_ids))
  
  all_fpkms <- all_fpkms %>%
    select(-all_of(outlier_ids))
  
} else {
  message("No outlier file found — using full dataset.")
}

if (!all(colnames(all_counts) == colnames(all_vst))) {
  stop("Column names of counts and VST don’t match!")
}

# Save the files
write.csv(all_counts, file = file.path(out_dir,"tidy_all_strains_counts.csv"))
write.csv(all_vst, file = file.path(out_dir,"tidy_all_strains_vst.csv"))
write.csv(all_fpkms, file = file.path(out_dir,"tidy_all_strains_fpkms.csv"))
write.csv(all_labels, file = file.path(out_dir,"tidy_all_strains_labels.csv"))

# Create the averages matrix ######
all_counts_avg <- all_counts %>%
  pivot_longer(cols = -protein_id, names_to = "sample", values_to = "raw_counts") %>%
  left_join(all_labels) %>%
  group_by(protein_id, strain_condition) %>%
  dplyr::summarise(mean_raw_counts = mean(raw_counts)) %>%
  ungroup() %>%
  mutate(strain_condition = gsub("-",".",strain_condition)) %>%
  pivot_wider(names_from = strain_condition, values_from = mean_raw_counts) %>%
  column_to_rownames("protein_id")

all_vst_avg <- all_vst %>%
  pivot_longer(cols = -protein_id, names_to = "sample", values_to = "vst_counts") %>%
  left_join(all_labels) %>%
  group_by(protein_id,strain_condition) %>%
  dplyr::summarise(mean_vst_counts = mean(vst_counts)) %>%
  ungroup() %>%
  mutate(strain_condition = gsub("-",".",strain_condition)) %>%
  pivot_wider(names_from = strain_condition, values_from = mean_vst_counts) %>%
  column_to_rownames("protein_id")

all_fpkms_avg <- all_fpkms %>%
  pivot_longer(cols = -protein_id, names_to = "sample", values_to = "fpkm_counts") %>%
  left_join(all_labels) %>%
  group_by(protein_id,strain_condition) %>%
  dplyr::summarise(mean_fpkm_counts = mean(fpkm_counts)) %>%
  ungroup() %>%
  mutate(strain_condition = gsub("-",".",strain_condition)) %>%
  pivot_wider(names_from = strain_condition, values_from = mean_fpkm_counts) %>%
  column_to_rownames("protein_id")

# Save the files
write.csv(all_counts_avg, file = file.path(out_dir,"tidy_all_strains_counts_avg.csv"))
write.csv(all_vst_avg, file = file.path(out_dir,"tidy_all_strains_vst_avg.csv"))
write.csv(all_fpkms_avg, file = file.path(out_dir,"tidy_all_strains_fpkms_avg.csv"))

# Load merged “all_strains” DEGs
all_DEG_data0 <- read.csv(file.path(deg_dir, "merged_all_strains_DEGs.csv")) 

all_DEG_data <- all_DEG_data0 %>%
  dplyr::rename(protein_id = ProteinID) %>%
  separate(
    ref_cond,
    into = c("ref_strain", "ref_condition"),
    sep  = "_(?=[^_]+_[^_]+$)",  # split at the underscore before the last two tokens
    remove = FALSE               # keep original column if you want
  ) %>%
  separate(
    condition,
    into = c("sample_strain", "sample_condition"),
    sep  = "_(?=[^_]+_[^_]+$)",  # split at the underscore before the last two tokens
    remove = FALSE               # keep original column if you want
  )

# sanity‐check: did renaming and splitting work?
if (!"protein_id" %in% names(all_DEG_data)) {
  stop("protein_id column missing after rename!")
}
if (any(is.na(all_DEG_data$ref_strain_file) |
        is.na(all_DEG_data$ref_condition_file) |
        is.na(all_DEG_data$sample_strain) |
        is.na(all_DEG_data$sample_condition))) {
  warning("Some rows failed the filename‐split into strain/condition.")
}

# Filter to the self‐vs‐self contrast for the chosen reference
comparison = conditions_in_mutant[1] # for debugging
for (comparison in conditions_in_mutant) {
  hm_DEG_mat <- all_DEG_data %>%
    filter(ref_condition == comparison,
           ref_strain == sample_strain) %>%
    # Map Expression → symbol
    mutate(DEG = case_when(
      Expression == "Overexpressed"  ~ "+",
      Expression == "Underexpressed" ~ "-",
      TRUE                           ~ " "
    )) %>%
    # Pivot to wide matrix form
    select(protein_id,condition, DEG) %>%
    pivot_wider(names_from = condition,
                values_from = DEG) %>%
    column_to_rownames("protein_id") %>%
    # add the reference columns
    mutate(!!!set_names(c(" "," "), paste0(strains,"_",comparison))) %>%
    as.matrix()
  
  write.csv(hm_DEG_mat, file = file.path(out_dir,paste0("tidy_all_strains_DEGmat_ref_",comparison,".csv")))
}