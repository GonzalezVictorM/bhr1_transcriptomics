# Clear environment (optional—consider if you really need to remove everything)
rm(list = ls())

# Load needed packages
library(tidyverse)

# Define and validate paths 
deg_dir <- "04_DEGanalysis"
out_dir <- deg_dir   # Change if you want merged files elsewhere

# Ensure input directory exists
if (!dir.exists(deg_dir)) {
  stop("Directory not found: ", deg_dir)
}

# Create output directory if missing
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Discover DEG files 
# Pattern assumes filenames like
#   dds_<all|single>_strains_ref_<strain>_<cond>_DEG_p<padj>_fc<fc>.csv
deg_files <- list.files(
  deg_dir,
  pattern = "^dds_(all_strains|single_strain)_ref_.+_DEG_.*\\.csv$",
  full.names = TRUE
)

if (length(deg_files) == 0) {
  stop("No DEG CSVs found in ", deg_dir)
}

# Parse metadata from filenames
meta <- tibble(
  path = deg_files,
  file = basename(deg_files)
) %>%
  extract(
    file,
    into = c("set", "strain_ref", "condition"),
    regex = "^dds_(all_strains|single_strain)_ref_(.+)_(.+_.+)_DEG_.*\\.csv$",
    remove = FALSE
  )

# Merge “all_strains” contrasts
all_paths <- meta %>% filter(set == "all_strains") %>% pull(path)

if (length(all_paths) > 0) {
  merged_all <- all_paths %>%
    set_names() %>%       # so .id is filename
    map_dfr(
      ~ read_csv(.x, show_col_types = FALSE, progress = FALSE) %>% 
        mutate(source = basename(.x)),
      .id = "source_file"
    ) %>%
    mutate(
      # basename without the “dds_all_strains_ref_…_DEG_…” wrapper
      contrast = source %>%
        # 1) remove the prefix up to and including the last underscore before the contrast
        str_remove("^dds_all_strains_ref_") %>%
        # 2) remove the suffix starting at _DEG_
        str_remove("_DEG_.*$"),
      set = "all_strains"
    ) %>%
    # split the cSource column into two sides, without clobbering file_condition
    separate(
      cSource,
      into = c("condition", "ref_cond"),
      sep = "_vs_",
      remove = FALSE
    ) %>%
    select(-source_file, -cSource, -source)
  
  write_csv(merged_all, file.path(out_dir, "merged_all_strains_DEGs.csv"))
  message("✔️ merged_all_strains_DEGs.csv written (", nrow(merged_all), " rows).")
} else {
  message("ℹ️  No all_strains files to merge.")
}

# Merge “single_strain” by strain_ref
single_meta <- meta %>% filter(set == "single_strain")

if (nrow(single_meta) > 0) {
  single_meta %>%
    split(.$strain_ref) %>%
    iwalk(function(df, strain) {
      paths <- df$path
      
      merged <- paths %>%
        set_names() %>%
        map_dfr(
          ~ read_csv(.x, show_col_types = FALSE, progress = FALSE) %>% 
            mutate(source = basename(.x)),
          .id = "source_file"
        ) %>%
        mutate(
          # basename without the “dds_all_strains_ref_…_DEG_…” wrapper
          contrast = source %>%
            # 1) remove the prefix up to and including the last underscore before the contrast
            str_remove("^dds_single_strain_ref_") %>%
            # 2) remove the suffix starting at _DEG_
            str_remove("_DEG_.*$"),
          set = "single_strain",
          strain_ref = strain,
          cSource = gsub("condition_","",cSource)
        ) %>%
        separate(
          cSource,
          into = c("condition", "ref_cond"),
          sep = "_vs_",
          remove = FALSE
        ) %>%
        mutate(condition = paste0(strain,"_",condition),
               ref_cond = paste0(strain,"_",ref_cond)) %>%
        select(-source_file, -cSource, -source)
      
      out_file <- file.path(out_dir, paste0("merged_single_strain_", strain, "_DEGs.csv"))
      write_csv(merged, out_file)
      message("✔️ ", basename(out_file), " written (", nrow(merged), " rows).")
    })
} else {
  message("ℹ️  No single_strain files to merge.")
}

# Done
message("✅ All merges complete.")
