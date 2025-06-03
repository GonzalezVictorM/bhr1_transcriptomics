# Package installation and version logging script

# List of packages
cran_packages <- c("tidyverse", "Hmisc", "pheatmap")
bioc_packages <- c("DESeq2", "vsn", "apeglm", "topGO", "clusterProfiler")

# Ensure BiocManager is available
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install CRAN packages
install.packages(cran_packages)

# Install Bioconductor packages
BiocManager::install(bioc_packages)

# Create log
log_df <- data.frame(
  Date = Sys.Date(),
  Package = c(cran_packages, bioc_packages),
  Version = sapply(c(cran_packages, bioc_packages), function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) as.character(packageVersion(pkg)) else "Not Installed"
  }),
  stringsAsFactors = FALSE
)

# Write log to file
log_file <- paste0("package_install_log_", Sys.Date(), ".txt")
write.table(log_df, file = log_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Installation complete. Log written to", log_file, "\n")