# bhr1_transcriptomics
Analysis of the transcriptomic data from bhr1 TF

## data-wrangling.R
Reorganizes counts and sample master data (MD) keeping all conditions present in either strains.

## data-exploration.R
Takes the tidy data and runs QC to identify outliers and printing a PCA of the cleaned up data. A section marked with (Added) can be run to keep only conditions present in both strains. Remember to double check the reference conditions used in the code.

## dds-pipeline.R
This is a robust code that will run all single strain and all strain comparisons. It is a bit of an overkill but it will save a lot of time later if you come up with new conditions to test. The "Define the strains and conditions "[..]" can be modified to select specific conditions and not run the whole set. Once again, the section marked with (Added) can be run to keep only conditions present in both strains.

## fold-change.R and merge-fold-change.R
These to codes run the shrinking, extraction of the fold change (FC) and false discovery rate (FDR), and classification of genes as DEG based on FC and FDR cut-offs.It takes quite long to run but again, it is worth it in the end.

## DEG-data-wrangling.R
It takes the counts, fpkms, vst counts, and labels and removes outliers and conditions not present in both strains, as well as calculates the avg of the biological replicates. Finally, it generates the DEG matrixes with + for overexpressed and - for underexpressed according to the classification done in fold-change.R

## CAZy_analysis.R
Using the previously tidied data, it makes heat maps and gene set enrichment analysis depending on the type of enzymes (eg. cellulases, hemicellulases).

## ST_analysis
A modified version of the CAZy_analysis code to perform the same tasks on sugar transporters.